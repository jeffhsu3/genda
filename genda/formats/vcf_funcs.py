import numpy as np

def skip_headers(VCF):
    """ Skips header of a vcf file.  Returns all the meta information to a
    string, the sample identifiers, and the number of meta_lines.  
    """
    meta = ''
    in_header = True
    line = fh.readline()
    l_c = 0
    while in_header:
        if line[0:2] == "##":
            meta += line
            line = fh.readline()
            l_c += 1
        elif line[0] == "#":
            s_l = line.rstrip("\n")
            samples = line.rstrip("\n").split("\t")[9:]
            in_header = False
    return (meta, samples, l_c)


def parse_genotype(line):
    """ Returns the genotype for the VCF file for each line
    """
    def parse_vcf_g(entry):
        g = entry.split(":")[0]
        if g == "0/0" or g == "0|0":
            return 0
        elif g == "1/1" or g == "1|1":
            return 2
        elif g == "0/1" or g == "0|1":
            return 1
        else:
            return np.NaN
    genotypes = [parse_vcf_g(i) for i in line.rstrip("\n").split("\t")[9:]]
    return genotypes
