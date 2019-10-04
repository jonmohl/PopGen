import os, sys

file = sys.argv[1]

if '.fasta' in file:
   f = open(%s'%(file), 'r')
   lines = f.readlines()
   f.close()

   i = 0

   fo = open('%s.headers.txt'%(di,file[0:-6]),'w')

   while i < len(lines):
      new_header = lines[i].strip().split()[0] + 'a\n'
      fo.write(new_header)
      fo.write(lines[i+1])
      i += 2
      if i<len(lines) and lines[i-2].split()[0] in lines[i]:
         new_header = lines[i].split()[0] + 'b\n'
         fo.write(new_header)
         fo.write(lines[i+1])
         i += 2
         
   fo.close()
