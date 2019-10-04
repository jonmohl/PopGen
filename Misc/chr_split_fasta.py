import sys, collections

fn_in = sys.argv[1]
ln_in = sys.argv[2]

file = open(ln_in, 'r')
lines = file.readlines()
file.close()

chr_list = []
chr_st2 = {}

#First line
sl = lines[0].split('\t')
pos = int(sl[0]) - 1
if 'apl' in sl[1]:
   chr = sl[1].split(':')[0]
elif 'mtDNA' in sl[1]:
   chr = 'mtDNA'
else:
   print line

st2 = [pos,0]
chr_list.append(chr)
chr_run = chr
chr_st2[chr] = []
c=0
for line in lines[1:]:
   sl = line.split('\t')
   pos = int(sl[0]) - 1
   if 'apl' in sl[1]:
      chr = sl[1].split(':')[0]
   elif 'mtDNA' in sl[1]:
      chr = 'mtDNA'
   else:
      print line

   if chr not in chr_list:
      chr_list.append(chr)
      chr_st2[chr] = []

   if chr != chr_run:
      chr_st2[chr_run].append(st2)
      st2 = [pos,0]
      chr_run = chr
   else:
      st2[1] = pos

chr_st2[chr_run].append(st2)


file = open(fn_in,'r')
lines = file.readlines()
file.close()

file_root = '.'.join(fn_in.split('.')[0:-1])

for chromo in chr_list:
   file = open('chromosome_split/%s_%s.fasta'%(file_root,chromo),'w')
   for x in xrange(0,len(lines),2):
      file.write(lines[x])
      for y in range(0,len(chr_st2[chromo])):

         try:
            file.write('%s'%lines[x+1][chr_st2[chromo][y][0]:chr_st2[chromo][y][1]])
         except:
            print x
            print lines[x]
      file.write('\n')
   file.close()

