""" Functions for working with VCF files.  Each function either takes a VCF or
a line from a VCF file

Jeff Hsu
"""

import sys
from itertools import cycle

import functools
import numpy as np
import pandas as pd
def parse_chr(entry):
    try:
        chrom = int(entry.lstrip("chr"))
        return(chrom)
    except ValueError:
        print("Chromosomes must be numeric")



def parse_geno(entry,GT=0):
    """ Need to somehow keep phase, maybe use multi_index, but only
    works if haplotype is contigous across whole section.  Maybe
    groups to denote haplotype regions?

    """
    g = entry.split(":")[GT]
    if g == "0/0" or g == "0|0" or g == '0':
        return 0
    elif g == "1/1" or g == "1|1" or g=='1':
        return 2
    elif g == "0/1" or g == "0|1":
        return 1
    else:
        return np.NaN

class VCF(object):
    """ A pandas dataframe wrapper for VCF files.
    """


    def __init__(self, vcf_file, chunksize = None, allelic_ratio=False,
                haplotypes = False, chrom_converter=None):
        """ Initialize with a VCF file with genotyping information
        """
        fh = open(vcf_file)
        meta, samples, num_meta = self._skip_headers(fh)
        fh.close()
        self._meta = meta
        self._parse_meta()
        self.samples = samples
        self.info_dict = []

        sex_chromosomes = ['X', 'Y']



        convertor = {}
        #convertor = dict((k + 9 , parse_geno) for k in range(len(samples)))
        #convertor = dict((k + 9 , allele_ratio) for k in range(len(samples)))
        if chrom_converter:
            convertor[0] = chrom_converter
        #convertor["INFO"] = self._parse_info
        self.vcf = pd.read_table(vcf_file, sep="\t",
                                 skiprows = xrange(num_meta),
                                 converters = convertor,
                                 chunksize = chunksize,
                                )

        pg=functools.partial(parse_geno,GT = self.vcf.ix[0,8].split(":").index("GT"))
        self.vcf.geno = self.vcf.ix[:,9:].applymap(pg)

        #self.vcf.rename(columns = {'#CHROM': 'CHROM'}, inplace=True)

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
        self.vcf.geno.index = self.vcf.index
        info_dict = self._split_info(self.vcf.INFO)

        # Parsing the INFO
        def flag_parse(info_dict, key):
            try:
                info_dict[key]
                return(True)
            except KeyError:
                return(False)

        def empty_field(info_dict, key):
            try:
                return(info_dict[key])
            except KeyError:
                return(np.NaN)

        def string_field(info_dict, key):
            try:
                return(info_dict[key])
            except KeyError:
                return('')

        # Need a better way to do this
        for i in self.info:
            if i[2] == np.bool:
                info_field = pd.Series(np.asarray(
                    [flag_parse(k, i[0]) for k in info_dict],
                                dtype = i[2]), index = self.vcf.index)
                #self.vcf[i[0]] = info_field
            elif i[1] > 1:
                # :TODO This needs to be fixed
                pass
            elif i[2] == np.float64 or i[2] == np.int:
                # :FIXME  float an integer
                info_field = pd.Series(np.asarray(
                    [empty_field(k, i[0]) for k in info_dict],
                                dtype = np.float64), index = self.vcf.index)
                self.vcf[i[0]] = info_field
            elif i[2] == '|S20':
                info_field = pd.Series(np.asarray(
                    [string_field(k, i[0]) for k in info_dict],
                                dtype = i[2]), index = self.vcf.index)
            elif i[0] == 'Dels':
                # For some reason this is null in some fields
                pass
            elif i[0] in ['CIEND', 'CIPOS', 'END', 'HOMLEN', 'HOMSEQ',
                          'SVLEN', 'SVTYPE', 'AA']:
                # :TODO This needs to be fixed 
                pass
            else:
                info_field = pd.Series(np.asarray([k[i[0]] for k in info_dict],
                                                dtype = i[1]),
                                    index = self.vcf.index)
                #self.vcf[i[0]] = info_field
        del self.vcf['INFO']
        del self.vcf['ID']
        # Parse GENO


        def _parse_haplotype(entry):
            g = entry.split(":")[0]
            try:
                g = g.split("|")
                return((int(g[0]), int(g[1])))
            except ValueError:
                return np.NaN
            except KeyError:
                print("Welp")


        def _allele_ratio(entry):
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

        # Create new dataframes for haplotypes and allele ratios
        if 'AD' in self.gformat:
            # Parse out allele ratios automatically
            self.ar = pd.DataFrame({'temp': np.zeros(len(self.vcf.index),
                                                     dtype=np.float)},
                                   index = self.vcf.index)
        if haplotypes:
            def bicycle(iterable, repeat=1):
                for item in iterable:
                    for _ in xrange(repeat):
                        yield item

            haplo_i = cycle([1,2])
            haplo_index = [(i, haplo_i.next())\
                           for i in bicycle(self.samples, repeat=2)]
            haplo_index = pd.MultiIndex.from_tuples(haplo_index, names =\
                                                 ['Individual',
                                                  'Allele'])
            self.haps = pd.DataFrame(np.zeros((len(self.vcf.index),
                                                     2*len(self.samples)),
                                                    dtype = np.float),
                                     columns = haplo_index,
                                     index = self.vcf.index
                                    )
        """
        for ind in self.samples:
            if haplotypes:
                haps = [_parse_haplotype(i) for i in self.vcf[ind]]
                self.haps[(ind,1)] =\
                        np.asarray([i[0] for i in haps], dtype = np.float)
                self.haps[(ind,2)] =\
                        np.asarray([i[0] for i in haps], dtype = np.float)
                del haps

            else:
                geno = np.asarray([_parse_geno(i) for i in self.vcf[ind]],
                                dtype = np.float64)
            if 'AD' in self.gformat:
                a_r = np.asarray([_allele_ratio(i) for i in self.vcf[ind]],
                                 dtype = np.float64)
                self.ar[ind] = a_r
            if haplotypes:
                pass
            else:
                self.vcf[ind] = geno
        """

        if 'AD' in self.gformat:
            del self.ar['temp']

        if haplotypes:
            pass

    def _parse_meta(self):
        """ Parses meta information and returns a list of the INFO and FORMAT
        fields.

        Okay so all VCF files are slightly different.  Why is the number sometimes a letter?  Jesus
        """
        convert = {
            'Integer': np.int,
            'Float': np.float64,
            'Flag': np.bool,
            'String': '|S20',
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
                                      line[1].split("=")[1],
                                      convert[line[2].split("=")[1]]))
                except KeyError:
                    print(line)
                    print('Improperly formatted meta information')
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


    def distance_based(self):
        pass

    def suspicious_loci(self):
        pass


    def list_samples_with_alternate_allele(self, rsID):
        """ Returns a list of samples that have the alternate non-ref allelex 
        at the variant identified by rsID
        """

        non_ref = np.logical_not(self.vcf.ix[rsID, 9:] == 0)
        nans = np.isnan(np.asarray(self.vcf.ix[rsID, 9:], dtype=np.float))
        return(self.vcf.columns[9:][np.logical_xor(non_ref, nans)])



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
        """ Not all vcfs will have this
        """
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
        allele_counts = [no_counts(i) for i in genotype]
        return(allele_counts)

    # Convienence Functions - Not sure if this is a good thing

    def head(self, **kwargs):
        print(self.vcf.head())

    """
    def __getitem__(self, val):
        self.vcf.ix
    """

def pos_line_convert(line):
    """ Convert chr:pos into a integer
    """
    line = line.rstrip("\n").split("\t")
    try:
        return int(line[1]) + (int(line[0].lstrip("chr"))*200000000)
    except ValueError:
        print("Malformed VCF file")


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



def is_indel(line):
    """ Returns if a line is an indel or not
    """
    line = line.rstrip("\n").split("\t")
    if len(line[3]) > 1 or len(line[4]) > 1:
        return True
    else:
        return False





