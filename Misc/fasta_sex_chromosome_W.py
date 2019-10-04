import os, sys

directory = sys.argv[1]
popmap = sys.argv[2]

f = open(popmap,'r')
lines = f.readlines()
f.close()
pop = {}
for line in lines:
   sl = line.strip().split('\t')
   pop[sl[0]] = [sl[1],sl[2]]

files = os.listdir(directory)
for f in files:
   if '.fa' in f or '.fasta' in f:
      fi = open('%s/%s'%(directory,f), 'r')
      lines = fi.readlines()
      fi.close()

      i = 0

      fo = open('%s/%s_W.fa'%(directory,'.'.join(f.split('.')[0:-1])),'w')
      while i < len(lines):
         if i<len(lines) and pop[lines[i].strip()[1:-1]][1] == 'F' and lines[i].strip()[-1]  == 'a':
            new_header = lines[i]
            fo.write(new_header)
            fo.write(lines[i+1])
         i += 2

         
      fo.close()
