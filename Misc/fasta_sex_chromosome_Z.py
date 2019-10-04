import os, sys

directory = sys.argv[1]
ic = sys.argv[2]

files = os.listdir(directory)
for f in files:
   if '.fa' in f or '.fasta' in f:
      fi = open('%s/%s'%(directory,f), 'r')
      lines = fi.readlines()
      fi.close()

      i = 0

      fo = open('%s/%s_sc.fa'%(directory,'.'.join(f.split('.')[0:-1])),'w')

      while i < len(lines):
         new_header = lines[i]
         fo.write(new_header)
         fo.write(lines[i+1])
         i += 2
         if i<len(lines) and lines[i-2][0:-1] == new_header[0:-1]:
            new_header = lines[i]
            fo.write(new_header)
            for n in range(0,len(lines[i+1].strip())):
               fo.write(ic)
            fo.write('\n')
            i += 2
         
      fo.close()
