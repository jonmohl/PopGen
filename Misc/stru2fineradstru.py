#Created by Jon Mohl, 1/17/19
#Usage:  python stru2fineradstru.py test.stru test.txt
#test.stru --> Input structure file with only population as extra column
#test.txt  --> fineRADstructure as output file
#header    --> If header is present put header as next argument


import sys

if len(sys.argv) == 4 and sys.argv[3] == 'header':
   header = 1
else:
   header = 0


f = open(sys.argv[1],'r')
lines = f.readlines()
f.close()

i2d = {'1':'A', '4':'C', '3':'G', '2':'T'}

toOut = []

if header == 1:
   j=1
else:
   j = 0

while j < len(lines):
   sla = lines[j].rstrip().split('\t')
   slb = lines[j+1].rstrip().split('\t')
   temp = [sla[0]]
   i = 2
   while i < len(sla):
      if sla[i] == '-9':
         temp.append('')
      else:
         temp.append(i2d[sla[i]]+'/'+i2d[slb[i]])
      i+=1
   j+=2
   toOut.append(temp)


if header == 1:
   loci_acc = lines[0].rstrip().split('\t')[2:]


trans = []
for i in range(0,len(toOut[0])):
   temp = []
   for j in range(0,len(toOut)):
      temp.append(toOut[j][i])
   trans.append(temp)


fo = open(sys.argv[2],'w')


fo.write('Chr\t' + '\t'.join(trans[0]) + '\n')

i=1
for t in trans[1:]:
   if header == 1:
      fo.write(loci_acc[i-1] + '\t' + '\t'.join(t) + '\n')
   else:
      fo.write(str(i) + '\t' + '\t'.join(t) + '\n')
   i+=1


fo.close()
