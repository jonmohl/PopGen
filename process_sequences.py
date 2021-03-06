#!/bin/python
#Created to process multiple sequencing files including a mapDamage check.
#UpdatedNovember 2, 2018
#Jon Mohl

import os,subprocess
import sys
import argparse

TRIM='java -jar /usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'
BWA='/usr/local/bwa-0.7.15/bwa'
MAP='/usr/bin/mapDamage'
SAM='/usr/local/samtools-1.6/samtools'
BCF='/usr/local/bin/bcftools'

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
start_at = int(args.start_at)
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
for fn in file_names:
   if fn[0] != '.':
      if PE and ('_R1' in fn):
         acc.append(fn.split('_R1')[0])
         if start_at == 1:
            fn2 = fn.split('_R1')[0]+'_R2.fastq.gz'
            comm = '%s PE -threads %s -phred33 %s/%s %s/%s %s/trimmed/%s_R1.tp.fastq.gz %s/trimmed/%s_R1.tu.fastq.gz %s/trimmed/%s_R2.tp.fastq.gz %s/trimmed/%s_R2.tu.fastq.gz SLIDINGWINDOW:5:30 MINLEN:30'%(TRIM,threads,raw,fn,raw,fn2,outdir,fn.split('_R1')[0],outdir,fn.split('_R1')[0],outdir,fn.split('_R1')[0],outdir,fn.split('_R1')[0])
            process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            process.wait()
            mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
            mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))
            #mapdamage_log.write(''.join(process.stdout.readlines()))
            #mapdamage_log.write(''.join(process.stderr.readlines()))
      elif not PE:
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
            #mapdamage_log.write(''.join(process.stdout.readlines()))
            #mapdamage_log.write(''.join(process.stderr.readlines()))

print('There are %s samples to process.'%len(acc))

mpileup_comm = '%s mpileup -a -g -Q 30 -q 25 -t AD,DP,INFO/AD -f %s '%(SAM,bait_ref)

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
      if PE:
         comm = '%s mem -t %s %s -T 30 -I 1024 -w O %s/trimmed/%s_R1.tp.fastq.gz %s/trimmed/%s_R2.tp.fastq.gz> %s/aligned/%s.bwa.sam'%(BWA,threads,bait_ref,outdir,a,outdir,a,outdir,a)
      else:
         comm = '%s mem -t %s %s -T 30 -I 1024 -w O %s/trimmed/%s.trimmed.fastq.gz > %s/aligned/%s.bwa.sam'%(BWA,threads,bait_ref,outdir,a,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      mapdamage_log.write('bwa command: %s\n'%comm)
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))
      #mapdamage_log.write(''.join(process.stdout.readlines()))
      #mapdamage_log.write(''.join(process.stderr.readlines()))

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
      #mapdamage_log.write(''.join(process.stdout.readlines()))
      #mapdamage_log.write(''.join(process.stderr.readlines()))

      #comm = '%s -i %s/aligned_prescale/%s.bwa.sam -r %s --rescale -d %s/mapdamage/%s --merge-reference-sequences'%(MAP,outdir,a,bait_ref,outdir,a)
      comm = '%s -i %s/aligned_prescale/%s.bwa.sam --rescale -r %s -d %s/mapdamage/%s --merge-reference-sequences'%(MAP,outdir,a,bait_ref,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      pso = [x.decode('utf-8') for x in iter(process.stdout.readlines())]
      pse = [x.decode('utf-8') for x in iter(process.stderr.readlines())]
      #pso = process.stdout.readlines()
      #pse = process.stderr.readlines()
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
         #mapdamage_log.write(''.join(process.stdout.readlines()))
         #mapdamage_log.write(''.join(process.stderr.readlines()))
      else:
         mapdamage_log.write('%s was rescaled.\n'%a)
         comm = '%s view -h -o %s/aligned/%s.bwa.sam %s/mapdamage/%s/%s.bwa.rescaled.bam '%(SAM,outdir,a,outdir,a,a)
         process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
         process.wait()
         print('MapDamage adjustment made on sample')
         mapdamage_log.write('samtools view command: %s\n'%comm)
         mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
         mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))
         #mapdamage_log.write(''.join(process.stdout.readlines()))
         #mapdamage_log.write(''.join(process.stderr.readlines()))
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
      #mapdamage_log.write(''.join(process.stdout.readlines()))
      #mapdamage_log.write(''.join(process.stderr.readlines()))

   if start_at <= 5:
      if not os.path.exists('%s/sorted'%outdir):
         os.mkdir('%s/sorted'%outdir)
      comm = '%s sort -o %s/sorted/%s.sorted.bam %s/dedup/%s.dedup.bam'%(SAM,outdir,a,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      mapdamage_log.write('samtools sort command: %s\n'%comm)
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))
      #mapdamage_log.write(''.join(process.stdout.readlines()))
      #mapdamage_log.write(''.join(process.stderr.readlines()))

      comm = '%s index %s/sorted/%s.sorted.bam'%(SAM,outdir,a)
      process = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      process.wait()
      mapdamage_log.write('samtools index command: %s\n'%comm)
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
      mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))
      #mapdamage_log.write(''.join(process.stdout.readlines()))
      #mapdamage_log.write(''.join(process.stderr.readlines()))

   mpileup_comm += ' %s/sorted/%s.sorted.bam'%(outdir,a)

   i+= 1


if start_at <= 6:
   print('Starting samtools mpileup command...')

   mpileup_comm += '| %s call --threads %s  -m -O v -o %s/PROJECT.sorted.bam.vcf'%(BCF,threads,outdir)

   mapdamage_log.write('mpileup command: %s\n'%mpileup_comm)
   process = subprocess.Popen(mpileup_comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
   process.wait()
   mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stdout.readlines())]))
   mapdamage_log.write('\n'.join([x.decode('utf-8').rstrip('\n') for x in iter(process.stderr.readlines())]))
   #mapdamage_log.write(''.join(process.stdout.readlines()))
   #mapdamage_log.write(''.join(process.stderr.readlines()))

mapdamage_log.write('Completed.')
mapdamage_log.close()
print('Completed')
