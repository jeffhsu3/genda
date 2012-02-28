""" Functions for working with VCF files.  Each function either takes a VCF or
a line from a VCF file
"""

import sys
import numpy as np
import pandas as pd

# Definitions taken from www.1000genomes.org/node/101

class VCF(object):
    """ Class containing genotyping a various metrics.
    """


    def __init__(self, vcf_file, chunksize = None, allelic_ratio=False):
        """ Initialize with a VCF file with genotyping information
        """

        fh = open(vcf_file)
        meta, samples, num_meta = self._skip_headers(fh)
        fh.close()
        self._meta = meta
        self._parse_meta()
        self.samples = samples
        self.info_dict = []
        def parse_geno(entry):
            """ Need to somehow keep phase, maybe use multi_index, but only
            works if haplotype is contigous across whole section

            Or maybe haplotypes should be set in a separate dataframe
            """
            g = entry.split(":")[0]
            if g == "0/0" or g == "0|0":
                return 0
            elif g == "1/1" or g == "1|1":
                return 2
            elif g == "0/1" or g == "0|1":
                return 1
            else:
                return np.NaN


        def allele_ratio(entry):
            try:
                t_i = sum([int(i) for i in\
                                    entry.split(":")[1].split(",")])
                alt_allele = int(entry.split(":")[1].split(",")[1])
                allele_ratio = alt_allele/float(t_i)
                return allele_ratio
            except IndexError:
                return np.NaN
            except ValueError:
                return np.NaN
            except ZeroDivisionError:
                return np.NaN


        def parse_chr(entry):
            try:
                chrom = int(entry.lstrip("chr"))
                return(chrom)
            except ValueError:
                print("VCF not a chrom")

        #convertor = dict((k + 9 , parse_geno) for k in range(len(samples)))
        #convertor = dict((k + 9 , allele_ratio) for k in range(len(samples)))
        convertor = {}
        convertor[0] = parse_chr
        #convertor["INFO"] = self._parse_info
        self.vcf = pd.read_table(vcf_file, sep="\t",
                                 skiprows = xrange(num_meta),
                                 converters = convertor,
                                 chunksize = chunksize,
                                )

        # Index using ID + unique identifier based on position for novel SNPS,
        # INDELS, SVS
        rsID = self.vcf["ID"]
        novel = [str(i) + "_" + str(j) + "_" + str(k) for (i, j, k)\
                 in zip(list(self.vcf["#CHROM"][rsID=="."]),
                        list(self.vcf["POS"][rsID == "."]),
                        list(self.vcf["ALT"][rsID == "."])
                       )
                ]
        self.novel = novel
        rsID[rsID == "."] = np.asarray(novel)
        self.vcf.index = pd.Index(rsID)

        info_dict = self._split_info(self.vcf.INFO)
        def flag_parse(info_dict, key):
            try:
                info_dict[key]
                return(True)
            except KeyError:
                return(False)

        for i in self.info:
            if i[1] == np.bool:
                info_field = pd.Series(np.asarray(
                    [flag_parse(k, i[0]) for k in info_dict],
                                dtype = i[1]), index = self.vcf.index)
                self.vcf[i[0]] = info_field
            elif i[0] == 'Dels':
                # For some reason this is null in some fields
                pass
            else:
                info_field = pd.Series(np.asarray([k[i[0]] for k in info_dict],
                                                dtype = i[1]),
                                    index = self.vcf.index)
                self.vcf[i[0]] = info_field


    def _parse_meta(self):
        """ Parses meta information and returns a list of the INFO and FORMAT
        fields.
        """
        convert = {
            'Integer': np.int,
            'Float': np.float64,
            'Flag': np.bool,
        }
        self.info = []
        self.gformat = []
        meta_lines = self._meta.split("\n")
        meta_lines = [i[2:] for i in meta_lines]
        for line in meta_lines:
            if line[0:4] == 'INFO':
                line = line.lstrip("INFO=<").rstrip(">")
                line = line.split(",")
                try:
                    self.info.append((line[0].split("=")[1],
                                 convert[line[2].split("=")[1]]))
                except KeyError:
                    print('Improperly formatted VCF meta file')
                #info.append(line[line.find('<ID=')+4:line.find(',')])
            elif line[0:6] == 'FORMAT':
                self.gformat.append(line[line.find('<ID=')+4:line.find(',')])
                #:TODO finish this

    def _parse_info(self, entry):
        """ Split the INFO field into a dictionary.  Creates an instance
        of the class varaible ans saves the dictionary to it.
        """
        def parse_field(field):
            try:
                field = field.split("=")
                return((field[0], field[1]))
            except IndexError:
                return((field[0], 0))

        fields = entry.split(";")
        self.info_dict.append(dict(parse_field(k) for k in fields))
        return(int(self.info_dict['AC']))

    def _split_info(self, INFO):
        """ Split the INFO into a dictionary
        """
        def parse_field(field):
            try:
                field = field.split("=")
                return((field[0], field[1]))
            except IndexError:
                return((field[0], 0))
        def make_dict(entry):
            fields = entry.split(";")
            info_dict = dict(parse_field(k) for k in fields)
            return(info_dict)
        out_dict = [make_dict(i) for i in INFO]
        return(out_dict)

    def _parse_ids(self , entry):
        if not entry == ".":
            return(entry)
        else:
            pass


    def HWE(self, df):
        """ Calculates Hardy-Weinberg Equilibrium at each loci in self.ge.  
        Really only detects the blatently wrong calls.  In addition some
        real genotypes can be thrown out, especially loci that are
        currently or have been undergoing selection.
        """


        n = df.shape[1]
        n_p = (df == 0).sum(axis=1)
        n_het = (df == 1).sum(axis=1)
        n_q = (df == 2).sum(axis = 1)
        p = (2 * n_p + n_het)/(2.0 * (n_p + n_het + n_q))
        q = 1 - p
        E_AA = (p ** 2) * n
        E_Aa = 2 * p * q * n
        E_aa = (q ** 2) * n
        chi_sq = (n_p - E_AA)**2/E_AA + (n_het - E_Aa)**2/E_Aa\
                + (n_q-E_aa)**2/E_aa
        from scipy.stats.stats import chisqprob
        return chisqprob(chi_sq, df=1)

    def distance_based(self):
        pass

    def suspicious_loci(self):
        pass

    def ti_tv_ratio(self):
        # 0 = Transition
        # 1 = Transversion
        # Humans expected to be around ~2.1
        table = {'AG' : 0, 'GA' : 0,
                 'CT' : 0, 'TC' : 0,
                 'AC' : 1, 'CA' : 1,
                 'AT' : 1, 'TA' : 1,
                 'GC' : 1, 'CG' : 1,
                 'GT' : 1, 'TG' : 1,
                }
        pass


    def _skip_headers(self, fh):
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

    def _parse_info(self, vcf_line):
        ac = 'k'
        return(ac)

    def _parse_genotypes(self, vcf_line):
        pass


    def _allelic_ratio(self, genotype):
        def no_counts(data):
            try:
                t_i = sum([int(i) for i in\
                                    data.split(":")[1].split(",")])
                alt_allele = int(data.split(":")[1].split(",")[1])
                allele_ratio = alt_allele/float(t_i)
                return allele_ratio
            except IndexError:
                return np.NaN
            except ValueError:
                return np.NaN
            except ZeroDivisionError:
                return np.NaN
        allele_counts = [no_counts(i) for i in line]
        return(allele_counts)

def pos_line_convert(line):
    """ Convert chr:pos into a integer
    """
    line = line.rstrip("\n").split("\t")
    try:
        return int(line[1]) + (int(line[0].lstrip("chr"))*200000000)
    except ValueError:
        print("Malformed VCF file")


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


def assign_novel(fh):
    """ This needs to return UNIQUE values for an index.  Use a hash function?
    """
    ident = []
    novel = []
    counter = 1
    for line in fh:
        line = line.rstrip("\n").split("\t")
        if line[2] == ".":
            novel.append(True)
            rsID = "ns"+line[0]+"_"+line[1]
            if rsID in ident:
                rsID = rsID + "_" + str(counter)
                counter += 1
            else:
                counter = 1
        else:
            rsID = line[2]
            novel.append(False)
        ident.append(rsID)
    return(ident, novel)


def count_alleles(line):
    line = line.rstrip("\n").split("\t")
    line = line[9:]
    def no_counts(data):
        try:
            allele_count = tuple([int(i) for i in\
                                  data.split(":")[1].split(",")])
            return allele_count
        except IndexError:
            return None
        except ValueError:
            return None
    allele_counts = [no_counts(i) for i in line]
    return(allele_counts)


def allelic_ratio(line):
    line = line.rstrip("\n").split("\t")
    line = line[9:]
    def no_counts(data):
        try:
            t_i = sum([int(i) for i in\
                                  data.split(":")[1].split(",")])
            alt_allele = int(data.split(":")[1].split(",")[1])
            allele_ratio = alt_allele/float(t_i)
            return allele_ratio
        except IndexError:
            return np.NaN
        except ValueError:
            return np.NaN
        except ZeroDivisionError:
            return np.NaN
    allele_counts = [no_counts(i) for i in line]
    return(allele_counts)


def is_indel(line):
    """ Returns if a line is an indel or not
    """
    line = line.rstrip("\n").split("\t")
    if len(line[3]) > 1 or len(line[4]) > 1:
        return True
    else:
        return False


def safe_distance(x, threshold):
    if x < threshold:
        return 1
    else:
        return 0


def print_plot(lines, ratios, genotypes, *args):
    """ Prints the plot of bad samples as well as the phred scaled log
    likelihood of the existance of the alternate allele
    """
    for i, line in enumerate(lines):
        line = line.rstrip("\n").split("\t")
        hets = genotypes.columns[genotypes.ix[i] == 1]
        b_s = ratios[hets].ix[i]
        extra = []
        try:
            for eargs in args:
                extra.append(eargs[i])
        except IndexError:
            print("Extra informal arguments must have same length as\
                  lines, ratios and genotypes")
            sys.exit()
        extra.append(line[0])
        extra.append(line[1])
        extra.append(line[2])
        extra.append(line[5])
        extra.append(str(b_s.mean()))
        print("\t".join(extra))

def print_good(lines, ratios, genotypes, *args):
    for i, line in enumerate(lines):
        line = line.rstrip("\n").split("\t")
        hets = genotypes.columns[genotypes.ix[i] == 1]
        b_s = ratios[hets].ix[i]
        try:
            if b_s.mean() > 0.3 and b_s.mean() < 0.7 and len(hets) < 10:
                print("\t".join(line))
            elif line[2] != ".":
                print("\t".join(line))
            else:
                pass
        except TypeError:
            pass


def main():
    import optparse

    fh = open(sys.argv[1], 'rU')
    bad_samples = ['s02-007', 's02-043', 's02-034','s02-069', 's02-077',
                   's02-088', 's03-072','s05-003', 'S02-013', 'S02-063']
    THRESHOLD = 50
    to_realign = open("targets.intervals", "w")

    lines = []
    meta = ""
    for line in fh:
        if line[0:2] == "##":
            meta += line
        elif line[0] == "#":
            s_l = line.rstrip("\n")
            samples = line.rstrip("\n").split("\t")[9:]
            start_of_samples = fh.tell()
        else:
            lines.append(line.rstrip("\n"))
    # lines are the variant lines within a vcf file
    counter = 0
    rs_ID = [assign_novel(line) for line in lines]
    # Right now rs_ID has duplicate values
    indels = [is_indel(line) for line in lines]
    n_loci = len(lines)
    loci = np.array([pos_line_convert(line) for line in lines])
    genotypes = [parse_genotype(line) for line in lines]
    ratio = [allelic_ratio(line) for line in lines]

    # Load genotypes into a pandas.dataframe
    table = pd.DataFrame(genotypes, columns=samples)
    ratios = pd.DataFrame(ratio, columns=samples)

    #print_plot(lines, ratios, table)
    print(meta)
    print(s_l)
    print_good(lines, ratios, table)
    # Checking SNP/INDEL Density
    loci_i = np.roll(loci,1)
    distance = loci - loci_i
    distance[0] = 1000
    potential_miscalls = [safe_distance(i, THRESHOLD) for i in distance]
    # Only matters if the alternate SNPs are called closely in the same samples
    truly_bad = []
    """
    for i in range(len(potential_miscalls)):
        if potential_miscalls[i] == 0:
            pass
        else:
            if i - 1 > 0:
                bad_samples = [] # Maybe make this a boolean vector
                seen = []
                for samp in samples:
                    if all(table[i-1:i+1][samp] > 0):
                        bad_samples.append(samp)
                    else:
                        pass
                if len(bad_samples) > 0:
                    seen.append(rs_ID[i])
                    if rs_ID[i-1] in seen:
                        pass
                    else:
                        pass
                    b_s = ratios[bad_samples].ix[i]
                    if b_s.mean() > 0.6:
                        truly_bad.append(i)
                    else:
                        pass
    """

if __name__ == '__main__':
    main()
