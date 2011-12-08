'''
Appends the average conservation score of a region to each entry in a bed file

usage:
    python conservationBed.py [OPTIONS] bedfile outputfile <conservationfile> 

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
    fc = tabix.Tabix('/gen_local/hsuj/ref/PWM/fastcons44.bed.gz')
#Output
w = open(args[1], 'w')

def main():
    for line in bed:
        scores = []
        l1 = line.strip().split('\t')
        chrom = l1[0]
        start = int(l1[1])
        end = int(l1[2])
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

if __name__ == '__main__':
    main()
