#!usr/bin/env python
#Created by Armando Lerma on December 2, 2015 to filter sequences from a E. vaginatum  fasta file using a list of target Locus Ids.
import re, sys

fastafile = sys.argv[1]  #complete full fasta file 
idsfile = sys.argv[2]    #ids file of target SNPs
seqfile = sys.argv[3]    #output fasta file

f2 = open(idsfile,'r')
f1 = open(fastafile,'r')
f3 = open(seqfile,'w')

#Store Ids from ids file in a dictionary
AI_DICT = {}
for line in f2:
    AI_DICT[line[:-1]] = 1
f2.close()
print(len(AI_DICT.keys()))


lines = f1.readlines()
f1.close()
i=0
while i < len(lines):
   locusid = lines[i].split(' ')[0][1:]
   if locusid in AI_DICT.keys():
      f3.write(lines[i])
      f3.write(lines[i+1])
   i+=2 

f3.close()
