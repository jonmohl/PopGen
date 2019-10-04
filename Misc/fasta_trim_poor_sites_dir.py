#Usage:  python fasta_trim_poor_sites_dir.py fastadir 0.7

import os, sys

dir = sys.argv[1]
cutoff = float(sys.argv[2])

files = os.listdir(dir)

for file in files:
   if '.fasta' in file or '.fa' in file:
      f = open('%s/%s'%(dir,file), 'r')
      lines = f.readlines()
      f.close()

      i = 0
      seq_list = []
      while i < len(lines):
          seq_list.append([lines[i],list(lines[i+1].strip())])
          i+=2

      toprint = []
      for i in range(0,len(seq_list[0][1])):
         count = 0
         for seq in seq_list:
            if seq[1][i] != 'N':
               count+=1
         if float(count)/len(seq_list) >= cutoff:
            toprint.append(1)
         else:
            toprint.append(0)
            
      new_seq = []
      for seq in seq_list:
         new_seq.append([seq[0],[]])

      for i in range(0,len(seq_list[0][1])):
         if toprint[i] == 1:
            for j in range(0,len(seq_list)):
               new_seq[j][1].append(seq_list[j][1][i])


      fo = open('%s/%s.trimmed.fasta'%(dir,file[:-6]),'w')
      for seq in new_seq:
         fo.write('%s%s\n'%(seq[0],''.join(seq[1])))
      fo.close()
