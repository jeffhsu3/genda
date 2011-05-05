def mirFamDict( filename ):
    """
    Create a dictionary from the miRBase miRNA family database in memory
    """
    fam_dict = {}
    with open(filename, "rU") as filehandle:
	for line in filehandle:
	    line = line.strip('\n').split('  ')
	    if line[0] == 'ID':
		current = line[1].strip(' ')
	    else: pass
	    if line[0] == 'MI':
		fam_dict[line[2]]=current
	    else: pass
        return(fam_dict)

