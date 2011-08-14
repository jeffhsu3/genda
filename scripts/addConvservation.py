'''
Appends the average conservation score of a region to each entry in a bed file
'''
import sys, tabix, optparse
p = optparse.OptionParser()
p.add_option('-p', '--padded',action = 'store', dest = 'pad', help = 'Program will remove this amount of padding when checking for scores, but will remain in the output bed file', default = 0)
options, args = p.parse_args()
# PWM bed
bed = open(args[0], 'rU')
# PhastCons
try:
    fc = tabix.Tabix(args[2])
except:
    fc = tabix.Tabix('fastcons44.bed.gz')
#Output
w = open(args[1], 'w')

for line in bed:
    scores = []
    l1 = line.strip().split('\t')
    chrom = l1[0]
    #Positive strand
    if int(l1[1]) > 0:
	start = int(l1[1])
	end = int(l1[2])
    #Negative strand
    else:
	diff = int(l1[1]) - int(l1[2])
	diff *= -1
	start = int(l1[1])*-1
	end = start+diff
    start += int(options.pad)
    end -= int(options.pad)
    try:
	#Finds all matching results
	res = fc.query(chrom, start, end)
	#Avereges the results
	total = 0
	howmany = 0
	for num,entry in enumerate(res):
	    howmany = num+1
	    total += int(entry[-1])
	if howmany == 0:
	    avg =0
	else:
	   avg = total/howmany
    #If a certain chromosome or haplotype is not found
    except tabix.error:
	avg = 0
    w.write('\t'.join(l1)+'\t'+str(avg)+'\n')
