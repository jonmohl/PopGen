#Usage:  python fasta_to_ped.py file.fasta popmap 0.5

import os, sys

fasta_name = sys.argv[1]
popmap_name = sys.argv[2]
maf = float(sys.argv[3])

file = open(fasta_name,'r')
fasta_in = file.readlines()
file.close()

file = open(popmap_name,'r')
lines = file.readlines()
file.close()

popmap = {}
for line in lines:
   sl = line.strip().split('\t')
   popmap[sl[0]] = sl[1]

ped = []
mapp = []
to_ped = []
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
   nuc_list = {}
   for fasta in fasta_list:
      if fasta[1][i] not in nuc_list.keys() and fasta[1][i] != '?' and fasta[1][i] != 'N'and fasta[1][i] != '-':
         nuc_list[fasta[1][i]] = 1
      elif fasta[1][i] != '?' and fasta[1][i] != 'N'and fasta[1][i] != '-':
         nuc_list[fasta[1][i]] += 1
      
      if fasta[2][i] not in nuc_list.keys() and fasta[2][i] != '?' and fasta[2][i] != 'N'and fasta[2][i] != '-':
         nuc_list[fasta[2][i]] = 1
      elif fasta[2][i] != '?' and fasta[2][i] != 'N'and fasta[2][i] != '-':
         nuc_list[fasta[2][i]] += 1
   if len(nuc_list) == 2:
      boo = 0
      tco = 0
      for k in nuc_list.keys():
         tco += nuc_list[k]
      for k in nuc_list.keys():
         if (float(nuc_list[k])/tco) > maf:
            boo += 1
      if boo == 2:
         to_ped.append(1)
      else:
         to_ped.append(0)
   else:
      to_ped.append(0)


for x in range(0,len(fasta_list)):
   ped_temp = ''
   ped_temp = '%s\t%s\t0\t0\tU\t0'%(popmap[fasta_list[x][0][1:-1]],fasta_list[x][0][1:-1])
   ped.append(ped_temp)

for i in range(0,len(fasta_list[0][1])):
   if to_ped[i] == 1:
      mapp.append('%s\t%s_%s\t0\t%s'%(fasta_name,fasta_name,i,i))
      for x in range(0,len(fasta_list)):
         if fasta_list[x][1][i] == '?' or fasta_list[x][1][i] == 'N' or fasta_list[x][1][i] == '-'or fasta_list[x][2][i] == '?' or fasta_list[x][2][i] == 'N'or fasta_list[x][2][i] == '-':
            ped[x] = ped[x] + '\t0\t0'
         else:
            ped[x] = ped[x] + '\t' + fasta_list[x][1][i] + '\t'  + fasta_list[x][2][i]


fo = open('%s.ped'%(fasta_name),'w')
fo.write('\n'.join(ped))
fo.close()

fo = open('%s.map'%(fasta_name),'w')
fo.write('\n'.join(mapp))
fo.close()
