#!/bin/python
#Created to process multiple sequencing files including a mapDamage check.
#UpdatedNovember 2, 2018
#Jon Mohl

import os,subprocess
import sys
import argparse

PEAR='/usr/local/pear-0.9.11/bin/pear'
TRIM='java -jar /usr/local/Trimmomatic-0.38/trimmomatic-0.38.jar'
BWA='/usr/local/bwa-0.7.15/bwa'
MAP='/usr/local/mapDamage/mapdamage'
SAM='/usr/local/samtools-1.6/samtools'
BCF='/usr/local/bcftools/bcftools'

#Create arg parser to allow for input of parameters
parser = argparse.ArgumentParser(description='Processing Bait Sequencing Files')
parser.add_argument('-i', '--input_directory', metavar="DIRNAME", dest='raw', required=True, help='Required: Input directory with raw fastq files.')
parser.add_argument('-b', '--bait_reference', metavar="FILENAME", dest='bait_ref', required=True, help='Required: Fasta sequence containing bait sequences.')
parser.add_argument('-o', '--output_dir', metavar="DIRNAME", dest='outdir', required=True, help='Required: Name of output directory to create.')
parser.add_argument('-a', '--ancient_list', metavar="FILENAME", dest='anc_list', help='This file contains the names of the files that may need the MapDamage correction')
parser.add_argument('-t', '--threads', metavar="INT", dest='threads', type=int, default=1)
parser.add_argument('-s', '--start_at', metavar="STRING", dest='start_at', type=str, default='1', help='This option allows the user to start at a particular portion of the analysis pipeline. (1) Full pipeline starting at quality trimming, (2) BWA alignment, (3) Mapdamage, (4) Samtools rmdup, (5) Samtools index and sort, (6) Samtools mpileup. Example: \'-s\' 3 would start the process at the Mapdamage step.') 
parser.add_argument('-PE', '--PairedEnd', dest='PE', action='store_true',help='Is used for PE fastq files')
args = parser.parse_args()

raw = args.raw
bait_ref = args.bait_ref
outdir = args.outdir
threads = args.threads
ancient = args.anc_list
if args.start_at:
   start_at = int(args.start_at)
else:
   start_at = 1
if args.PE:
   PE = True
else:
   PE = False


#Create output directories
if start_at == 1:
   #Check to make sure directory is not there
   if os.path.exists(outdir):
      print("Error: Output directory already exists.  Please choose different directory or remove old directory.")
      quit()
   os.mkdir(outdir)
   os.mkdir('%s/trimmed'%outdir)

mapdamage_log = open('%s/bait_run_log.txt'%outdir,'a')

print('Starting to process sequences.')

file_names = os.listdir(raw)
acc = []

if PE:
   os.mkdir('%s/merged'%outdir)
   os.mkdir('%s/merged/tmp'%outdir)
   for fn in file_names:
      if fn[0] != '.':
         if ('_R1' in fn) and fn[-3:] != '.gz':
            acc.append(fn.split('_R1')[0])
            fn2 = fn.split('_R1')[0]+'_R2.fastq'
            comm = '/usr/local/pear-0.9.11/bin/pear -f %s/%s -r %s/%s -p 0.0001 -v 10 -m 0 -n 0 -q 25 -i -j %s -o merged/%s_merged.fastq'%(raw,fn,raw,fn2,threads,outdir,fn.split('_R1')[0])
         elif '_R1' in fn:
            acc.append(fn.split('_R1')[0])
            fn2 = fn.replace('_R1','_R2')
            if start_at == 1:
               comm = 'cp %s/%s %s/merged/tmp/'%(raw,fn,outdir)
               process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
               process.wait()
               comm = 'cp %s/%s %s/merged/tmp/'%(raw,fn2,outdir)
               process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
               process.wait()
               comm = 'gunzip %s/merged/tmp/*.gz'%(outdir)
               process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
               process.wait()
               comm = '/usr/local/pear-0.9.11/bin/pear -f %s/merged/tmp/%s -r %s/merged/tmp/%s -p 0.0001 -v 10 -m 0 -n 0 -q 25 -i -j %s -o %s/merged/%s_merged.fq'%(outdir,fn[0:-3],outdir,fn2[0:-3],threads,outdir,fn.split('_R1')[0])
         print(comm)
         if '_R1' in fn and start_at == 1:
            process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            process.wait()
            comm = 'gzip %s/merged/*.fq'%outdir
            print(comm)
            process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            process.wait()
            mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
            mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))
            comm = '%s SE -threads %s -phred33 %s/merged/%s_merged.fq.gz %s/trimmed/%s.trimmed.fastq.gz SLIDINGWINDOW:5:30 MINLEN:30'%(TRIM,threads,outdir,acc[-1],outdir,acc[-1])
            print(comm)
            process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            process.wait()
            mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
            mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))
   if start_at == 1:
      comm = 'gzip %s/merged/tmp/*'
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()

else:
   for fn in file_names:
      if 'fastq' in fn.split('.'):
         acc.append(fn.split('.')[0])
         comm = '%s SE -threads %s -phred33 %s/%s %s/trimmed/%s.trimmed.fastq.gz SLIDINGWINDOW:5:30 MINLEN:30'%(TRIM,threads,raw,fn,outdir,fn.split('.')[0])
      else:
         new_acc = fn[0:-3]
         acc.append(new_acc)
         comm = '%s SE -threads %s -phred33 %s/%s %s/trimmed/%s.trimmed.fastq.gz SLIDINGWINDOW:5:30 MINLEN:30'%(TRIM,threads,raw,fn,outdir,new_acc)
      if start_at == 1:
         process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
         process.wait()
         mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
         mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))


print('There are %s samples to process.'%len(acc))

mpileup_comm = '%s mpileup -Ou -Q 30 -q 25 -a AD,DP,SP,INFO/AD -f %s '%(BCF,bait_ref)

ancient_list = []

try:
   file_ancient = open(ancient,'r')
   lines = file_ancient.readlines()
   file_ancient.close()
   for line in lines:
      ancient_list.append(line.strip())
except:
   print('No ancient samples detected.')
   mapdamage_log.write('No ancient samples detected.\n')

i = 1
for a in acc:
   print('Processing %s of %s files.'%(i,len(acc)))
   mapdamage_log.write('Processing %s of %s files.'%(i,len(acc)))

   if start_at <= 2:
      if not os.path.exists('%s/aligned'%outdir):
         os.mkdir('%s/aligned'%outdir)

      comm = '%s mem -t %s %s -T 30 -I 1024 -w O %s/trimmed/%s.trimmed.fastq.gz > %s/aligned/%s.bwa.sam'%(BWA,threads,bait_ref,outdir,a,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      mapdamage_log.write('bwa command: %s\n'%comm)
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))


   if a in ancient_list and start_at <= 3:
      if not os.path.exists('%s/aligned_prescale'%outdir):
         os.mkdir('%s/aligned_prescale'%outdir)
      if not os.path.exists('%s/mapdamage'%outdir):
         os.mkdir('%s/mapdamage'%outdir)

      comm = 'mv %s/aligned/%s.bwa.sam %s/aligned_prescale/%s.bwa.sam'%(outdir,a,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      mapdamage_log.write('Move bam file for mapDamage command: %s\n'%comm)
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))


      #comm = '%s -i %s/aligned_prescale/%s.bwa.sam -r %s --rescale -d %s/mapdamage/%s --merge-reference-sequences'%(MAP,outdir,a,bait_ref,outdir,a)
      comm = '%s -i %s/aligned_prescale/%s.bwa.sam --rescale -r %s -d %s/mapdamage/%s --merge-reference-sequences'%(MAP,outdir,a,bait_ref,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      pso = [x.decode('utf-8') for x in iter(process.stdout.readlines())]
      pse = [x.decode('utf-8') for x in iter(process.stderr.readlines())]
      mapdamage_log.write('mapDamage command: %s\n'%comm)
      mapdamage_log.write(''.join(pso))
      mapdamage_log.write(''.join(pse))



      map_boo = 1
      for l in pso:
         if 'Warning: DNA damage levels are too low, the Bayesian computation should not be performed' in l:
            map_boo = 0

      if map_boo == 0:
         mapdamage_log.write('%s was not rescaled.\n'%a)
         comm = 'mv %s/aligned_prescale/%s.bwa.sam %s/aligned/%s.bwa.sam'%(outdir,a,outdir,a)
         process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
         process.wait()
         mapdamage_log.write('Return corrected bam file to command: %s\n'%comm)
         print('No adjustment made on sample')
         mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
         mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))
      else:
         mapdamage_log.write('%s was rescaled.\n'%a)
         comm = '%s view -h -o %s/aligned/%s.bwa.sam %s/mapdamage/%s/%s.bwa.rescaled.bam '%(SAM,outdir,a,outdir,a,a)
         process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
         process.wait()
         print('MapDamage adjustment made on sample')
         mapdamage_log.write('samtools view command: %s\n'%comm)
         mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
         mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))

      print('After mapdamage')      

   if start_at <= 4:
      if not os.path.exists('%s/dedup'%outdir):
         os.mkdir('%s/dedup'%outdir)
      comm = '%s rmdup --output-fmt BAM -sS %s/aligned/%s.bwa.sam %s/dedup/%s.dedup.bam'%(SAM,outdir,a,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      mapdamage_log.write('samtools rmdup command: %s\n'%comm)
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))

      comm = 'gzip %s/aligned/*.sam'
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()

   if start_at <= 5:
      if not os.path.exists('%s/sorted'%outdir):
         os.mkdir('%s/sorted'%outdir)
      comm = '%s sort -o %s/sorted/%s.sorted.bam %s/dedup/%s.dedup.bam'%(SAM,outdir,a,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      mapdamage_log.write('samtools sort command: %s\n'%comm)
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))

      comm = '%s index %s/sorted/%s.sorted.bam'%(SAM,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      mapdamage_log.write('samtools index command: %s\n'%comm)
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))


   mpileup_comm += ' %s/sorted/%s.sorted.bam'%(outdir,a)

   i+= 1


if start_at <= 6:
   print('Starting samtools mpileup command...')

   mpileup_comm += '| %s call --threads %s -a GQ,GP -m -o %s/PROJECT.sorted.bam.vcf'%(BCF,threads,outdir)

   mapdamage_log.write('mpileup command: %s\n'%mpileup_comm)
   process = subprocess.Popen(mpileup_comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
   process.wait()
   mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
   mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))

mapdamage_log.write('Completed.')
mapdamage_log.close()
print('Completed')
