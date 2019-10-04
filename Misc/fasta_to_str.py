#Usage:  python fasta_to_str.py file.fasta popmap

import os, sys

fasta_name = sys.argv[1]
popmap_name = sys.argv[2]

file = open(fasta_name,'r')
fasta_in = file.readlines()
file.close()

file = open(popmap_name,'r')
lines = file.readlines()
file.close()

pop_list = []
popmap = {}
for line in lines:
   sl = line.strip().split('\t')
   popmap[sl[0]] = sl[1]
   if sl[1] not in pop_list:
      pop_list.append(sl[1])


stru = []

to_str = []
#Create fasta list
i = 0
fasta_list = []
while i < len(fasta_in):
   if fasta_in[i][0:-3] == fasta_in[i+2][0:-3]:
      fasta_list.append([fasta_in[i].strip(),list(fasta_in[i+1].strip()),list(fasta_in[i+3].strip())])
      i+=4
   else:
      fasta_list.append([fasta_in[i].strip(),list(fasta_in[i+1].strip())])
      i+=2

#Determine sites to keep

for i in range(0,len(fasta_list[1][1])):
   nuc_list = []
   for fasta in fasta_list:
      if fasta[1][i] not in nuc_list:
         nuc_list.append(fasta[1][i])
      if len(fasta) > 2 and fasta[2][i] not in nuc_list:
         nuc_list.append(fasta[1][i])
   if len(nuc_list) > 1:
      to_str.append(1)
   else:
      to_str.append(0)

str_temp = '\t'

temp = fasta_name.replace(',','_')
for x in range(0,len(to_str)):
   if to_str[x] == 1:
      str_temp = str_temp + '\t%s_%s'%(temp,x)
stru.append(str_temp)

for fasta in fasta_list:
   str_temp = '%s\t%s'%(fasta[0][1:-1],pop_list.index(popmap[fasta[0][1:-1]]))
   str_temp2 = '%s\t%s'%(fasta[0][1:-1],pop_list.index(popmap[fasta[0][1:-1]]))
   for i in range(0,len(to_str)):
      if to_str[i] == 1:
         if fasta[1][i] == '?':
            str_temp = str_temp + '\t-9'
            str_temp2 = str_temp2 + '\t-9'
         else:
            str_temp = str_temp + '\t%s'%fasta[1][i]
            str_temp2 = str_temp2 + '\t%s'%fasta[2][i]
   stru.append(str_temp)
   stru.append(str_temp2)


fo = open('%s.stru'%(fasta_name),'w')
fo.write('\n'.join(stru))
fo.close()

