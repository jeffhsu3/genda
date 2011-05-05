"""
    Functions for parsing GTF files
"""

def parseID(GTFline, ID):
	padding = 2
	start=GTFline[8].find(ID)
	start+=len(ID)+padding
	length = GTFline[8][start:].find(";") 
	return(GTFline[8][start:start+length].strip('"'))
