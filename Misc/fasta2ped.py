import os,sys,random

dir_in = sys.argv[1]
snp_choice = sys.argv[2]
if snp_choice not in ['all','random','first']:
   sys.exit('Not a valid snp selection')
popmap_name = sys.argv[3]
maf = float(sys.argv[4])

popmap = {}
sample_names = []
sex = {}
f = open(popmap_name,'r')
lines = f.readlines()
f.close()
for line in lines:
   sl = line.strip().split('\t')
   popmap[sl[0]] = sl[1]
   sample_names.append(sl[0])
   sex[sl[0]] = sl[2]

fasta_files = os.listdir(dir_in)

#Create ped and map files
ped = []
mapp = []

for x in range(0,len(sample_names)):
   if sample_names[x] in popmap.keys():
      ped_temp = ''
      if sex[sample_names[x]] == 'M':
         s = '1'
      elif sex[sample_names[x]] == 'F':
         s = '2'
      else:
         s = 'U'
      ped_temp = '%s\t%s\t0\t0\t%s\t0'%(popmap[sample_names[x]],sample_names[x],s)
      ped.append(ped_temp)

fasta_files.sort()


for ff in fasta_files:
   if '.fasta' in ff or '.fa' in ff:
      f = open('%s/%s'%(dir_in,ff), 'r')
      lines = f.readlines()
      f.close()
      seq_dic = {}
      v = 0
      while v < len(lines):
        name = lines[v].strip()[1:-1]
        seq_dic[name] = [lines[v+1].strip().replace('N','?').replace('-','?'),lines[v+3].strip().replace('N','?').replace('-','?')]
        v+=4
      i = 0
      snp = []

      for i in range(0,len(seq_dic[seq_dic.keys()[0]][0])):
         m = []
         mc = []
         for x in range(0,len(sample_names)):
            if sample_names[x] in seq_dic.keys(): 
               if seq_dic[sample_names[x]][0][i] not in m and seq_dic[sample_names[x]][0][i] != '?':
                  m.append(seq_dic[sample_names[x]][0][i])
                  mc.append(1)
               elif seq_dic[sample_names[x]][0][i] != '?':
                  mc[m.index(seq_dic[sample_names[x]][0][i])] += 1
               if seq_dic[sample_names[x]][1][i] not in m and seq_dic[sample_names[x]][1][i] != '?':
                  m.append(seq_dic[sample_names[x]][0][i])
                  mc.append(1)
               elif seq_dic[sample_names[x]][1][i] != '?':
                  mc[m.index(seq_dic[sample_names[x]][1][i])] += 1
         if len(m) == 2 and float(min(mc))/sum(mc) > maf:
            snp.append(i)

      if len(snp) > 0:
         if snp_choice == 'first':
            snp = [snp[0]]
         elif snp_choice == 'random':
            snp = [random.choice(snp)]

         for s in snp:
            ffs = ff.split('.')[0].split('_')
            ap = 1
            if len(ffs) == 2:
               try:
                  ap = int(ffs[1]) + s
               except:
                  ap = 1 
            mapp.append('%s\t%s_%s\t0\t%s'%(ffs[0],ff.split('.')[0],i,ap))
            for x in range(0,len(sample_names)):
               if sample_names[x] in seq_dic.keys():
                  if seq_dic[sample_names[x]][0][s] == '?' or seq_dic[sample_names[x]][1][s] == '?':
                     ped[x] = ped[x] + '\t0\t0'
                  else:
                     ped[x] = ped[x] + '\t' + seq_dic[sample_names[x]][0][s] + '\t'+ seq_dic[sample_names[x]][1][s]
               else:
                     ped[x] = ped[x] + '\t0\t0'

#Write results to files
fo = open('%s.ped'%(dir_in),'w')
fo.write('\n'.join(ped))
fo.close()

fo = open('%s.map'%(dir_in),'w')
fo.write('\n'.join(mapp))
fo.close()

