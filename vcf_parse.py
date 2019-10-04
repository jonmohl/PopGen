#!/usr/bin/python
#Created by Jon Mohl, Nov 27, 2018
#20190328: Added chromosomal output to include both "A & B". JEM
#20190404: Added chromosomal sort function. JEM
#

#Popmap format:   sampleid<tab>pop<tab>sex
#
#Example command: python ~/scripts/popgenome/vcf_parse_v1.3.py -i jm_test/PROJECT.sorted.bam.vcf -o jm_test_duck -p popmap2 -f ../references/bait-redo.fasta -mwpc 0.1 -af 0.05 -ac 2
#
#

import os, re, sys, random, operator
import argparse, random

def nc_nucleotide(nuc):
   if 'A' in nuc and 'G' in nuc:
      return 'R'
   elif 'C' in nuc and 'T' in nuc:
      return 'Y'
   elif 'G' in nuc and 'C' in nuc:
      return 'S'
   elif 'A' in nuc and 'T' in nuc:
      return 'W'
   elif 'T' in nuc and 'G' in nuc:
      return 'K'
   elif 'A' in nuc and 'C' in nuc:
      return 'M'
   else:
      return 'N'

def create_bait(bait, sn):
   bait_d = [bait[0],len(bait[3]),list(bait[3]),[],[],[],[],0,bait[1]]
   tmp = ' '.join(bait_d[2])
   for n in range(0,len(sn)):
      bait_d[3].append(tmp.split(' '))
      bait_d[4].append(tmp.split(' '))
      bait_d[5].append(tmp.split(' '))
   return bait_d

def pop_checks(vl, sl, pops, tc):
   pc = {}
   for x in range(0,len(vl)):
      if int(vl[x].split(':')[1]) >= min_seq and sl[x] in pops.keys(): 
         if pops[sl[x]] in pc.keys():
            pc[pops[sl[x]]] = pc[pops[sl[x]]] + 1
         else:
            pc[pops[sl[x]]] = 1
   i=0
   for k in pc.keys():
      p = float(pc[k])/tc[k]
      if p <= mwpc:
         return False
      else:
         i+=1
   if i == 0:
      return False
   elif float(i)/len(pc.keys()) <= mpc:
      return False
   return True

#Defaults
min_seq = 5
un_cut = 0.75
mwpc = 0.75
mpc = 0.75
cutoff = 0.05
ign_list = ['BLASHAHASALHFLKAHSDLKHJASLKHDLAKHSDLKHADLKJH']
all_count = 2

parser = argparse.ArgumentParser(description='Filtering of VCF file.')
parser.add_argument('-i', '--infile', metavar="FILENAME", dest="fn_in", required=True, help='Required: input file is necessary to run vcf_parse.')
parser.add_argument('-o', '--outdir', metavar="DIRNAME", dest="direc", required=True, help='Required: output directory in which all output is stored.')
parser.add_argument('-p', '--popmap', metavar="FILENAME", dest="pop_file",required=True, help='Required: population mapping file containing sample name from VCF file tab separated from population information.')
parser.add_argument('-f', '--fasta_ref', metavar="FILENAME", dest="fasta_file", required=True, help='Required: fasta file containing the reference used to generate the VCF file.')
parser.add_argument('-af', '--alt_allele_freq', metavar="Float", dest="cutoff", type=float, help='The minimal alternative allele frequency to inlcude the alternative allele. Default=0.05')
parser.add_argument('-ac', '--alt_allele_count', metavar="Float", dest="all_count", type=float, help='The minimal alternative allele sequence depth to include the alternative allele. Default=2')  
parser.add_argument('-ms', '--min_seq_depth', metavar="INT", dest="min_seq", type=float, help='The minimal sequence depth to maintain the GSV. Default=5')
parser.add_argument('-sc', '--seq_cutoff', metavar="FLOAT", dest="un_cut", type=float, help='The minimal amount of unknowns that can be within a sequence before the sample sequence is dropped. Default=0.75')
parser.add_argument('-mwpc', '--min_within_pop_coverage', metavar="FLOAT", dest="mwpc", type=float, help='The minimal percentage of samples within a population before the site is dropped. Default=0.75')
parser.add_argument('-mpc',  '--min_pop_coverage', metavar="FLOAT", dest="mpc", type=float, help='The minimal percentage of populations before the site is dropped. Default=0.75')
parser.add_argument('-snp', '--snp_selection', metavar="STR", dest="snp_select", help='The number of SNPs returned per locus. Can either be first, random, or all.  Default=all', )
parser.add_argument('-ic', '--ignore_chromosome', metavar="LIST", dest="ign_list", nargs='+', help='A list of chromosomes to ignore while processing the VCF file.  Example -ic chr1 chr2')
parser.add_argument('-cf', '--chromosome_file', metavar="FILENAME", dest="chromo_file", help='A file containing the order of the chromosomes for the chromosomal sort function')

args = parser.parse_args()

#Create variables
fn_in = args.fn_in
direc = args.direc
pop_file = args.pop_file
fasta_file = args.fasta_file
if args.min_seq:
   min_seq = args.min_seq
if args.un_cut:
   un_cut = args.un_cut
if args.mwpc:
   mwpc = args.mwpc
if args.mpc:
   mpc = args.mpc
if args.cutoff:
   cutoff = args.cutoff
if args.ign_list:
   ign_list = args.ign_list
if args.all_count:
   all_count = args.all_count
if args.snp_select:
   if args.snp_select in ['all','first','random']:
      snp_select = args.snp_select
   else:
      sys.exit('Error: SNP Selection is not valid.  Please use either: all, first, or random')
      
else:
   snp_select = 'all'

print('Selecting'+ snp_select + 'snps')

chromo_order = []
if args.chromo_file:
   f = open(args.chromo_file,'r')
   lines = f.readlines()
   f.close()
   for l in lines:
      chromo_order.append(l.strip())

if '/' in fn_in:
   file_root = fn_in.split('/')[-1][0:-4]
else:
   file_root = fn_in[0:-4]       #Creates root output file name

os.mkdir(direc)

bait_dict = {}
chromo_list = {}
bait2chromo = {}
ff = open(fasta_file,'r')
lines = ff.readlines()
ff.close()
boo = 0
seq = ''
ch = 'AHCALHSLKHJALKJHGLKDHGLSKHDFLKHJSDLKFHSLDKFH'
chromo = '-1'

#Create fasta sequences for baits that could be used
ref_seq = {}
for line in lines:
   if line[0] == '>':
      if boo == 1 and not any(s in ch for s in ign_list) and not any(s in chromo for s in ign_list):
         bait2chromo[ch] = chromo
         ref_seq[ch] = [chromo, start, ch, seq]
         ch = line[1:].split()[0]
         if 'CHROMO=' in line:
            chromo = line.strip().split('CHROMO=')[1].split(' ')[0]
         else:
            chromo = 'UN'
         if 'START=' in line:
            start = int(line.strip().split('START=')[1].split(' ')[0])
         else:
            start = 0
         seq = ''
      elif boo == 0 and not any(s in ch for s in ign_list) and not any(s in chromo for s in ign_list):
         bait2chromo[ch] = chromo
         ch = line[1:].split()[0]
         if 'CHROMO=' in line:
            chromo = line.strip().split('CHROMO=')[1].split(' ')[0]
         else:
            chromo = 'UN'
         if 'START=' in line:
            start = int(line.strip().split('START=')[1].split(' ')[0])
         else:
            start = 0
         seq = ''
         boo = 1
      else:
         bait2chromo[ch] = chromo
         seq = ''
         ch = line[1:].split()[0]
         if 'CHROMO=' in line:
            chromo = line.strip().split('CHROMO=')[1].split(' ')[0]
         else:
            chromo = 'UN'
         if 'START=' in line:
            start = int(line.strip().split('START=')[1].split(' ')[0])
         else:
            start = 0
   else:
      seq = seq + line.strip()
bait2chromo[ch] = chromo
if not any(s in ch for s in ign_list) and not any(s in chromo for s in ign_list):
   ref_seq[ch] = [chromo, start, ch, seq]

#Read in the population file
file = open(pop_file,'r')
lines = file.readlines()
file.close()
popmap = {}
pop_tc = {}
sex = {}
for line in lines:
   sl = line.strip().split()
   popmap[sl[0]] = sl[1]
   if len(sl) == 3:
      sex[sl[0]] = sl[2]
   else:
      sex[sl[0]] = 'U'
   if sl[1] in pop_tc.keys():
      pop_tc[sl[1]] += 1
   else:
      pop_tc[sl[1]] = 1


#Open vcf output files for parsing
vcf_out = open('%s/%s_good.vcf'%(direc,file_root),'w')
vcf_bad_out = open('%s/%s_bad.vcf'%(direc,file_root),'w')

snp_list = []               #Create a list of snps to process
indel_list = []             #Create a list of indels to process

#Sort the different lines of the VCF file
with open(fn_in) as f:
  for line in f:
   if not line.startswith('#'):
      temp = []
      sl = line.split('\t')
      if 'INDEL' not in sl[7] and not any(s in sl[0] for s in ign_list) and not any(s in bait2chromo[sl[0]] for s in ign_list):
         bait = sl[0]
         loc = int(sl[1]) -1
         adl = sl[7].split(';')[1].split(',')[1:]
         mtalc = -1
         for a in adl:
            if int(a) > mtalc:
               mtalc = int(a) 
         if '<*>' != sl[4] and len(sl[4].split(',')) < 3 and mtalc >= all_count and  pop_checks(sl[9:], sample_names, popmap, pop_tc):
            vcf_out.write(line)
            snp_list.append(line)
         else:
            vcf_bad_out.write(line)
      elif 'INDEL' in sl[7] and len(sl[4].split(',')) == 1 and not any(s in sl[0] for s in ign_list) and not any(s in bait2chromo[sl[0]] for s in ign_list) and pop_checks(sl[9:], sample_names, popmap, pop_tc):
         indel_list.append(line)
   elif line.startswith('#C'):
      vcf_out.write(line)
      vcf_bad_out.write(line)
      sample_names = [n.split('/')[-1] for n in line.strip().split('\t')[9:]]
   elif line.startswith('##contig='):
      vcf_out.write(line)
      vcf_bad_out.write(line)
   else:
      vcf_out.write(line)
      vcf_bad_out.write(line)

#Close vcf files
vcf_out.close()
vcf_bad_out.close()

cs = 0
for x in popmap.keys():
   if x in sample_names:
      cs += 1

print('Number of samples: %s/%s used within VCF file.'%(cs,len(sample_names)))

#Process the INDELs
good_indel_count = 0
print('Processing the INDELs...%s INDELs'%len(indel_list))
for line in indel_list:
   sl = line.split('\t')
   if sl[0] not in bait_dict.keys():
      bait_dict[sl[0]] = create_bait(ref_seq[sl[0]], sample_names)
      if  bait_dict[sl[0]][0] not in chromo_list.keys():
         chromo_list[bait_dict[sl[0]][0]] = [sl[0]]
      else:
         chromo_list[bait_dict[sl[0]][0]].append(sl[0])
   if bait_dict[sl[0]][7] == 0:
      bait_dict[sl[0]][7] = 1
   for x in range(0,len(sample_names)):
      if sample_names[x] in popmap:
         il = sl[9+x].split(':')
         cl = il[2].split(',')
         if int(cl[0]) > min_seq and int(cl[1]) > all_count and (float(cl[1])/float(il[1])) >= (1 - cutoff):
            for n in range((int(sl[1])+len(sl[4].split(',')[0])-1),(int(sl[1])+len(sl[3])-1)):
               tmp = bait_dict[sl[0]][3][x]
               tmp[n] = '-'
               bait_dict[sl[0]][3][x] = tmp
               tmp = bait_dict[sl[0]][4][x]
               tmp[n] = '-'
               bait_dict[sl[0]][4][x] = tmp
               tmp = bait_dict[sl[0]][5][x]
               tmp[n] = '-'
               bait_dict[sl[0]][5][x] = tmp
               if n not in bait_dict[sl[0]][6]:
                  bait_dict[sl[0]][6].append(n)
         elif int(cl[0]) > min_seq  and int(cl[1]) > all_count and (float(cl[1])/float(il[1])) >= cutoff:
            r = random.uniform(0, 1)
            tmp = bait_dict[sl[0]][3][x]
            for n in range((int(sl[1])+len(sl[4].split(',')[0])-1),(int(sl[1])+len(sl[3])-1)):
               tmp = bait_dict[sl[0]][3][x]
               tmp[n] = tmp[n].lower()
               bait_dict[sl[0]][3][x] = tmp
               if r > 0.5:
                  tmp = bait_dict[sl[0]][4][x]
                  tmp[n] = '-'
                  bait_dict[sl[0]][4][x] = tmp
               else:
                  tmp = bait_dict[sl[0]][5][x]
                  tmp[n] = '-'
                  bait_dict[sl[0]][5][x] = tmp
               if n not in bait_dict[sl[0]][6]:
                  bait_dict[sl[0]][6].append(n)

#Processing the SNPs
good_snp_count = 0
print('Processing the SNPs... %s SNPs'%len(snp_list))
for line in snp_list[0:]:
   sl = line.split('\t')
   if sl[0] not in bait_dict.keys():
      bait_dict[sl[0]] = create_bait(ref_seq[sl[0]], sample_names)
      if  bait_dict[sl[0]][0] not in chromo_list.keys():
         chromo_list[bait_dict[sl[0]][0]] = [sl[0]]
      else:
         chromo_list[bait_dict[sl[0]][0]].append(sl[0])
   #if bait_dict[sl[0]][7] == 0:
   #   bait_dict[sl[0]][7] = 1
   n = int(sl[1]) - 1
   for x in range(0,len(sample_names)):
      if sample_names[x] in popmap.keys():
         il = sl[9+x].split(':')
         cl = il[2].split(',')
         sam_dp = int(il[1])
         if int(cl[0]) > min_seq and int(cl[1]) > all_count and (float(cl[1])/float(il[1])) >= (1-cutoff):
            if '-' not in bait_dict[sl[0]][4][x][n] and '-' not in bait_dict[sl[0]][5][x][n]:
               tmp = bait_dict[sl[0]][3][x]
               tmp[n] = sl[4].split(',')[0]
               bait_dict[sl[0]][3][x] = tmp
               tmp = bait_dict[sl[0]][4][x]
               tmp[n] = sl[4].split(',')[0]
               bait_dict[sl[0]][4][x] = tmp
               tmp = bait_dict[sl[0]][5][x]
               tmp[n] = sl[4].split(',')[0]
               bait_dict[sl[0]][5][x] = tmp
            elif '-' in bait_dict[sl[0]][5][x][n]:
               tmp = bait_dict[sl[0]][4][x]
               tmp[n] = sl[4].split(',')[0]
               bait_dict[sl[0]][4][x] = tmp
            elif '-' in bait_dict[sl[0]][4][x][n]:
               tmp = bait_dict[sl[0]][5][x]
               tmp[n] = sl[4].split(',')[0]
               bait_dict[sl[0]][5][x] = tmp
            if n not in bait_dict[sl[0]][6]:
               bait_dict[sl[0]][6].append(n)
               bait_dict[sl[0]][7] = 1
         elif int(cl[0]) > min_seq and int(cl[1]) > all_count and (float(cl[1])/float(il[1])) >= cutoff:
            if '-' not in bait_dict[sl[0]][4][x][n] and '-' not in bait_dict[sl[0]][5][x][n]:
               tmp =bait_dict[sl[0]][3][x]
               tmp[n] = nc_nucleotide([sl[3],sl[4].split(',')[0]])
               bait_dict[sl[0]][3][x] = tmp
               r = random.uniform(0, 1)
               if r > 0.5:
                  tmp = bait_dict[sl[0]][4][x]
                  tmp[n] = sl[4].split(',')[0]
                  bait_dict[sl[0]][4][x] = tmp
               else:
                  tmp = bait_dict[sl[0]][5][x]
                  tmp[n] = sl[4].split(',')[0]
                  bait_dict[sl[0]][5][x] = tmp
            elif '-' in bait_dict[sl[0]][5][x][n]:
               tmp = bait_dict[sl[0]][4][x]
               tmp[n] = sl[4].split(',')[0]
               bait_dict[sl[0]][4][x] = tmp
            elif '-' in bait_dict[sl[0]][4][x][n]:
               tmp = bait_dict[sl[0]][5][x]
               tmp[n] = sl[4].split(',')[0]
               bait_dict[sl[0]][5][x] = tmp
            if n not in bait_dict[sl[0]][6]:
               bait_dict[sl[0]][6].append(n)
               bait_dict[sl[0]][7] = 1
         elif il[1] == '0' and '-' not in bait_dict[sl[0]][4][x][n] and '-' not in bait_dict[sl[0]][5][x][n]:
            tmp = bait_dict[sl[0]][3][x]
            tmp[n] = '?'
            bait_dict[sl[0]][3][x] = tmp
            tmp = bait_dict[sl[0]][4][x]
            tmp[n] = '?'
            bait_dict[sl[0]][4][x] = tmp
            tmp = bait_dict[sl[0]][5][x]
            tmp[n] = '?'
            bait_dict[sl[0]][5][x] = tmp
         elif int(cl[0]) > min_seq and int(cl[1]) > all_count and (float(cl[1])/float(il[1])) >= cutoff and (bait_dict[sl[0]][4][x][n] == '-' or bait_dict[sl[0]][5][x][n] == '-'):
            print(line)

print('Finished processing input file...\nStarting to generate output files...')

print('Before bait sort')
if len(chromo_list) != 1 and not chromo_list.keys()[0] == 'UN':
   for chr in chromo_list:
      b = chromo_list[chr]
      #b_sorted = sorted(b, key=lambda b : b[7])
      b_sorted = sorted(b)
      chromo_list[chr] = b_sorted

print('Before chromosome sort')
if len(chromo_list.keys()) != 1:
   if len(chromo_order) > 0:
      sor_chromo_list = []
      for c in chromo_order:
         if c in chromo_list.keys():
            sor_chromo_list.append(c)
   else:
      sor_chromo_list = sorted(chromo_list)
else:
   sor_chromo_list = chromo_list.keys()


os.mkdir('%s/bait_fasta'%direc)
os.mkdir('%s/bait_single_fasta'%direc)

no_bait = []
##Create fasta files
for bait in bait_dict.keys():
   if bait_dict[bait][7] == 1:
      fn = '_'.join(bait.split(':'))
      fa_bait = open('%s/bait_fasta/%s.fa'%(direc,fn),'w')
      fa_sin_bait = open('%s/bait_single_fasta/%s.fa'%(direc,fn),'w')
      for y in range(0,len(sample_names)):
         if sample_names[y] in popmap:
            if bait_dict[bait][4][y].count('?') == 0:
#            if bait_dict[bait][4][y].count('?') <= len(bait_dict[bait][4][y])*un_cut:
               fa_bait.write('>%sa\n'%sample_names[y])
               fa_bait.write(''.join(bait_dict[bait][4][y]).replace('.',''))
               fa_bait.write('\n>%sb\n'%sample_names[y])
               fa_bait.write(''.join(bait_dict[bait][5][y]).replace('.',''))
               fa_bait.write('\n')
               fa_sin_bait.write('>%s\n%s\n'%(sample_names[y],''.join(bait_dict[bait][3][y]).replace('.','')))
      fa_bait.close()
      fa_sin_bait.close()
   else:
      no_bait.append(bait)

print('Number of baits without reported information: %s'%len(no_bait))


os.mkdir('%s/chromosome_fasta'%direc)
os.mkdir('%s/chromosome_fasta_ab'%direc)
for chromo in sor_chromo_list:
   seq = []
   seqa = []
   seqb = []
   sam_num = []
   fo = open('%s/chromosome_fasta/%s.fa'%(direc,chromo),'w')
   fo2 = open('%s/chromosome_fasta_ab/%s_ab.fa'%(direc,chromo),'w')
   for x in range(0,len(sample_names)):
      if sample_names[x] in popmap.keys():
         if '/' in sample_names[x]:
            seq.append(['>%s\n'%sample_names[x].split('/')[1]])
            seqa.append(['>%sa\n'%sample_names[x].split('/')[1]])
            seqb.append(['>%sb\n'%sample_names[x].split('/')[1]])
         else:
            seq.append(['>%s\n'%sample_names[x]])
            seqa.append(['>%sa\n'%sample_names[x]])
            seqb.append(['>%sb\n'%sample_names[x]])
         sam_num.append(x)
   for cb in chromo_list[chromo]:
      if bait_dict[cb][7] == 1:
         for x in sam_num:
            if sample_names[x] in popmap.keys():
               if bait_dict[bait][4][y].count('?') == 0:
                  seq[sam_num.index(x)].append(''.join(bait_dict[cb][3][x]).replace('.',''))
                  seqa[sam_num.index(x)].append(''.join(bait_dict[cb][4][x]).replace('.',''))
                  seqb[sam_num.index(x)].append(''.join(bait_dict[cb][5][x]).replace('.',''))
               else:
                  seq[sam_num.index(x)].append('?'*len(bait_dict[cb][3][x]))
                  seqa[sam_num.index(x)].append('?'*len(bait_dict[cb][3][x]))
                  seqb[sam_num.index(x)].append('?'*len(bait_dict[cb][3][x]))
   for x in sam_num:
      if sample_names[x] in popmap.keys():
         fo.write(''.join(seq[sam_num.index(x)]))
         fo.write('\n')
         fo2.write(''.join(seqa[sam_num.index(x)]))
         fo2.write('\n')
         fo2.write(''.join(seqb[sam_num.index(x)]))
         fo2.write('\n')
   fo.close()

#Create ped and map files
ped = []
mapp = []

for x in range(0,len(sample_names)):
   if sample_names[x] in popmap.keys():
      ped_temp = ''
      if sex[sample_names[x]] == 'M':
         s = '1'
      elif sex[sample_names[x]] == 'F':
         s = '2'
      else:
         s = 'U'
      ped_temp = '%s\t%s\t0\t0\t%s\t0'%(popmap[sample_names[x]],sample_names[x],s)
      ped.append(ped_temp)

for chromo in sor_chromo_list:
   for cb in chromo_list[chromo]:
      if snp_select == 'all':
         for i in bait_dict[cb][6]:
            loc = bait_dict[cb][7] + i
            ap = ref_seq[cb][1] + loc
            mapp.append('%s\t%s_%s\t0\t%s'%(chromo,cb,loc,ap))
            for x in range(0,len(sample_names)):
               if sample_names[x] in popmap.keys():
                  if bait_dict[cb][4][x][i] != '?':
                     ped[sam_num.index(x)] = ped[sam_num.index(x)] + '\t' + bait_dict[cb][4][x][i] + '\t'
                  else:
                     ped[sam_num.index(x)] = ped[sam_num.index(x)] + '\t0\t'
                  if bait_dict[cb][5][x][i] != '?':
                     ped[sam_num.index(x)] = ped[sam_num.index(x)] + bait_dict[cb][5][x][i]
                  else:
                     ped[sam_num.index(x)] = ped[sam_num.index(x)] + '0'
      elif snp_select == 'first':
         boo = 0
         j = 0
         if len(bait_dict[cb][6]) != 0:
            while boo == 0 and j < len(bait_dict[cb][6]):
               count = 0
               m = []
               i = bait_dict[cb][6][j]
               for x in range(0,len(sample_names)):
                  if sample_names[x] in popmap.keys():
                     if (bait_dict[cb][4][x][i] != bait_dict[cb][2][i] or bait_dict[cb][5][x][i] != bait_dict[cb][2][i]) and bait_dict[cb][4][x][i] != '?':
                        count +=1
                        if bait_dict[cb][4][x][i] != bait_dict[cb][2][i] and bait_dict[cb][4][x][i] not in m and bait_dict[cb][4][x][i] != '?':
                           m.append(bait_dict[cb][4][x][i])
                        if bait_dict[cb][5][x][i] != bait_dict[cb][2][i] and bait_dict[cb][5][x][i] not in m and bait_dict[cb][4][x][i] != '?':
                           m.append(bait_dict[cb][5][x][i])
               af = count/float(len(sample_names))
               if af >= cutoff and len(m) == 1:
                  print(j)
                  boo = 1
               else:
                  j += 1
            if boo == 1:
               lloc = bait_dict[cb][6][j]
               loc = bait_dict[cb][7] + lloc
               ap = ref_seq[cb][1] + loc
               mapp.append('%s\t%s_%s\t0\t%s'%(chromo,cb,loc,ap))
               for x in range(0,len(sample_names)):
                  if sample_names[x] in popmap.keys():
                     if bait_dict[cb][4][x][lloc] != '?':
                        ped[sam_num.index(x)] = ped[sam_num.index(x)] + '\t' + bait_dict[cb][4][x][lloc] + '\t'
                     else:
                        ped[sam_num.index(x)] = ped[sam_num.index(x)] + '\t0\t'
                     if bait_dict[cb][5][x][lloc] != '?':
                        ped[sam_num.index(x)] = ped[sam_num.index(x)] + bait_dict[cb][5][x][lloc]
                     else:
                        ped[sam_num.index(x)] = ped[sam_num.index(x)] + '0'
      elif snp_select == 'random':
         j = 0
         l = []
         if len(bait_dict[cb][6]) != 0:
            while j < len(bait_dict[cb][6]):
               count = 0
               m = []
               i = bait_dict[cb][6][j]
               for x in range(0,len(sample_names)):
                  if sample_names[x] in popmap.keys():
                     if (bait_dict[cb][4][x][i] != bait_dict[cb][2][i] or bait_dict[cb][5][x][i] != bait_dict[cb][2][i]) and bait_dict[cb][4][x][i] != '?':
                        count +=1
                        if bait_dict[cb][4][x][i] != bait_dict[cb][2][i] and bait_dict[cb][4][x][i] not in m and bait_dict[cb][4][x][i] != '?':
                           m.append(bait_dict[cb][4][x][i])
                        if bait_dict[cb][5][x][i] != bait_dict[cb][2][i] and bait_dict[cb][5][x][i] not in m and bait_dict[cb][4][x][i] != '?':
                           m.append(bait_dict[cb][5][x][i])
               af = count/float(len(sample_names))
               print(cutoff,af,m)
               if af >= cutoff and len(m) == 1:
                  l.append(bait_dict[cb][6][j])
               j += 1
         #Radnom select snp
         if len(l) > 0:
            lloc = random.choice(l)
            loc = bait_dict[cb][7] + lloc
            ap = ref_seq[cb][1] + loc
            mapp.append('%s\t%s_%s\t0\t%s'%(chromo,cb,loc,ap))
            for x in range(0,len(sample_names)):
               if sample_names[x] in popmap.keys():
                  if bait_dict[cb][4][x][lloc] != '?':
                     ped[sam_num.index(x)] = ped[sam_num.index(x)] + '\t' + bait_dict[cb][4][x][lloc] + '\t'
                  else:
                     ped[sam_num.index(x)] = ped[sam_num.index(x)] + '\t0\t'
                  if bait_dict[cb][5][x][lloc] != '?':
                     ped[sam_num.index(x)] = ped[sam_num.index(x)] + bait_dict[cb][5][x][lloc]
                  else:
                     ped[sam_num.index(x)] = ped[sam_num.index(x)] + '0'

fo = open('%s/%s.ped'%(direc,file_root),'w')
fo.write('\n'.join(ped))
fo.close()

fo = open('%s/%s.map'%(direc,file_root),'w')
fo.write('\n'.join(mapp))
fo.close()

temp = ''
fo = open('%s/%s_info.txt'%(direc,file_root),'w')
fo.write('Command: python %s\n'%(' '.join(sys.argv)))
fo.write('Output directory: %s\n'%direc)
fo.write('Number of samples: %s/%s used in VCF file\n'%(len(popmap.keys()),len(sample_names)))
fo.write('Number of variants in original file: %s\n'%(len(indel_list) + len(snp_list)))
fo.write('Number of INDELs: %s\n'%len(indel_list))
fo.write('Number of SNPs: %s\n'%len(snp_list))
fo.close()

