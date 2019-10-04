import sys

fin = open(sys.argv[1],'r')
lines = fin.readlines()
fin.close()

num_samples = int(sys.argv[2])
min_cutoff = int(sys.argv[3])


mat = []

for line in lines[1:]:
   count = 0
   sl = line.split('\t')
   temp = [sl[2]]
   for n in range(0,num_samples):
      temp.append('0')
   sam_list = sl[8].split(',')
   for sam in sam_list:
      n = int(sam.split('_')[0])
      temp[n+1] = '1'
      count += 1
   if count >= min_cutoff:
      mat.append(temp)

fo = open('test2.csv','w')
for m in mat:
   fo.write(','.join(m) + '\n')

fo.close()
   
  
