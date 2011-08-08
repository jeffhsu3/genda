"""
    python walker.py [OPTIONS] <bed/vcf/gff/tab> <bamfile1> [bamfile2, ...]

    Right now takes a GTF produced by cufflinks (cuffcompare) output from
    cuffcompare and outputs a sql database of their transcripts and their exon
    positions 

"""
import optparse, sys, csv, os, sqlite3
from pySeq.common_funcs import regionParse

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-N", "--novel", action="store_true", dest="novel",\
	    help="Only check at novel loci")
    p.add_option("-i", "--identifier", dest="i",\
	    help="Identifier to use")
    p.add_option("-o", "--outupt", dest="o",\
	    help="Output database")
    options, args = p.parse_args()

    fh = csv.reader(open(args[0], "rU"), delimiter="\t")
    # Things to look for in the detail column 
    ID = "transcript_id"
    class_code = "class_code"

    padding = 2

    debug = 0 
    seen = {}
    for line in fh:
	#parsing the data column
	start=line[8].find(ID)
	start+=len(ID)+padding
	length = line[8][start:].find(";") 
	ident = line[8][start:start+length].strip('"')
	start=line[8].find(class_code)
	start+=len(class_code) + padding
	length = line[8][start:].find(";")
	tclass = line[8][start:start+length].strip('"')

	if ident in seen:
	    # :TODO The exons probably need to be sorted
	    seen[ident][3] += ","+str(line[3])
	    seen[ident][4] += ","+str(line[4])
	else:
	    seen[ident] = [ident,tclass, line[0], str(line[3]), str(line[4])]
	
	if options.D:
	    debug += 1
	    if debug > 40:
		break

    transcript = map(tuple, seen.values())

    if options.o:
	conn = sqlite3.connect(options.o) 
    else:
	conn = sqlite3.connect(os.join(os.getcwd(), args.split('.')[0]))

    c = conn.cursor()
    c.execute('''create table transcripts\
	    (tcons text,class text, region text, start int, end int)''')
    c.executemany("insert into transcripts(tcons,class, region, start, end) values \
	    (?,?,?,?,?)",transcript)
    conn.commit()
    c.close()




if __name__=='__main__':
    main()

