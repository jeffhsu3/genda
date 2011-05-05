"""
    python walker.py [OPTIONS] tconsdatabase cufflinks.diff <bamfile1> [bamfile2, ...]
    Goes through a cuffdiff diff file and gets counts at the transcripts
"""

import optparse, sys, csv, os, sqlite3, pysam, re
import numpy as np
from pySeq.pysam_callbacks.geneCounter import geneCounter


def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-d", "--database", dest="d",\
	    help="database")
    p.add_option("-T", "--transcript", dest="t",\
	    help="Transcript_id")
    p.add_option("-C", "--class", dest="C",\
	    help="Look only at a particular Cufflinks class of transcript")

    options, args = p.parse_args()
    fh = csv.reader(open(args[1], "rU"), delimiter="\t")
    debug=0
    conn = sqlite3.connect(args[0]) 
    # Don't return unicode
    conn.text_factory = str
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    bamfile=pysam.Samfile(args[2], "rb")	

    # Print Header
    header=fh.next()
    header.extend(["length", "counts"])
    print("\t".join(header))
    unanno = re.compile('^RP[0-9]')

    for line in fh:
	#if entry['class'] == options.C or unanno.match(line[4]):
	if line[11] == "yes" and (unanno.match(line[1]) or line[1]=="-"):
	    t = (line[0],)
	    c.execute('select * from transcripts where tcons=?', t) 
	    entry = c.fetchone()
	    transcript = geneCounter(line[0]) 
	    if unanno.match(line[1]) or (line[1]=="-" and entry['class']=='u'):
		try:
		    region, start, end = entry["region"], entry["start"].split(','),\
			    entry["end"].split(',')
		    start = np.array(map(int, start))	
		    end = np.array(map(int, end))
		    for i, j in zip(start, end):
			bamfile.fetch(region,i,j, callback=transcript)
		    line.append(str(sum(end-start)))
		    line.append(str(transcript.mCounts))
		    print("\t".join(line))
		except AttributeError:
		    # For single Exon Entries
		    region, start, end = entry["region"], int(entry["start"]), \
			    int(entry["end"]) 
		    bamfile.fetch(region,start,end, callback=transcript)
		    print("\t".join([transcript.gene, str(transcript.mCounts),\
			    str(end-start), line[4]]))


		if options.D:
		    debug += 1
		    if debug > 40:
			break
	else: pass
    



if __name__=='__main__':
    main()

