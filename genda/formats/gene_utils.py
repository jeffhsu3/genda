"""
"""
import json
import requests
import gzip


def get_symbol_ensembl(ensembl):
    """Get symbol from ensembl ID from 
    ensembl REST API   
    """

    server = "http://rest.ensembl.org"  
    ext = "/lookup/id/{0}?expand=1"  
    r = requests.get(server + ext.format(str(ensembl)), 
                    headers={ "Content-Type" : "application/json"})
    if not r.ok:   
        return('NA')
    else:
        decoded = r.json()
        try:
            return(decoded['display_name']) 
        except KeyError:
            return('NA')


def grab_gene_location(gid, ann_file = None, cis_buffer=0):
    """ Grab the gene boundries as a tuple buffered by 
    cis_buffer

    returns (chr, pos0, pos1)
    """

    if ann_file:
        with gzip.open(ann_file, mode='rt') as annot:
            for line in annot:
                if gid in line:
                    line = line.split("\t")
                    return(line[0], int(line[1]) - cis_buffer, 
                            int(line[2]) + cis_buffer, str(line[4]))
    else:
        raise NotImplementedError 


def get_transcript_ids(gene, ref37=False):
    """ Get transcript ids from ensembl 
    REST API from ensemble gene ID
    """
    if ref37:
        server = "http://grch37.rest.ensembl.org/"
    else:
        server = "http://rest.ensembl.org"
    ext = "/lookup/id/{0}?species=homo_sapiens;expand=1"
    r = requests.get(server+ext.format(gene), headers={ "Content-Type" :
        "application/json"})
    out_transcripts = []
    decoded = r.json()
    for i in decoded['Transcript']:
        out_transcripts.append(i['id'])
    return(out_transcripts)
