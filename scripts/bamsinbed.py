import sys, numpy, pysam, pickle, optparse
p = optparse.OptionParser()
p.add_option("-t", "--tab", action="store_true", dest="t", help="Output a tab deliminated array.Default is pickles")
options, args = p.parse_args()
bed = open(args[0],'rU')
w = open(args[-1],'w')
bedheight = 0
bedwidth = 0
dataformat='p'
if options.t==True:
    dataformat='t'
for line in bed:
    if bedwidth == 0:
	line=line.split('\t')
	bedwidth=int(line[2])-int(line[1])+1
    bedheight+=1
a = numpy.zeros((bedheight,bedwidth))
filecounter = 1
while args[filecounter] != args[-1]:
    bam = pysam.Samfile(args[filecounter],'rb')
    bed.seek(0)
    for num, line in enumerate(bed):
	line = line.split('\t')
	for entry in bam.fetch(line[0], int(line[1]), int(line[2])):
	    pos = entry.pos
	    a[num,pos-int(line[1])]+=1
    filecounter += 1
if dataformat=='p':
    pickle.dump(a,w)
if dataformat=='t':
    for row in a:
	for element in row:
	    w.write(str(element)+'\t')
	w.write('\n')
