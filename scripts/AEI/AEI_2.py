#!/usr/bin/env python

""" Gets the read bases at each position specified by the first argument, which
is a VCF/BED/GFF/GTF file, from an indexed BAMfile/s or pileup/s:
    1. Whether the loci is heterzygous or homozygous if the genotyping is not
    known
    2. The p-value of the loci exhibiting allelic imbalance if the loci is
    heterzygous.

Usage:
    python getReads.py bedfile/vcf bamfile1 [bamfile2 ...] -o OUTPUT

Limtations: Doesn't accurately account paralogous gene expression, which can
manifest itself in allelic imbalance.  One way would be to have an input file
of all paralogous genes and calculated the probability of a read mismatching
and then the liklihood that it is something not valulable.

:TODO Add VCF support. Had it in originally
:TODO Do comparison with reference allele
:TODO Account for mis-alignments of the minor allele
:TODO Need to be able to handle insertion/deletions in the sequence as compared
to reference. This will be hard to do
:TODO add an option that includes a gene annotation file
:TODO add support for pileup files
:TODO also how does PCR duplicates affect?

:TODO Output the results into a VCF file

Written by Jeffrey Hsu
2010
"""


import pysam, optparse, re, sys
from pySeq.formats.VCF import VCFfile
import numpy as np
import pySeq.stats.likelihood_funcs as lf


def callback_fn(pileups, position, samples, phredThreshold):
    position

class lociInformation():
    """ Contains information about a given position in a BAM file
    """

    BASE_INDEX = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
    def __init__(self, region, position, samples, phredThreshold=0):
        """ Location is one-based.  Internally this will get changed to zero
        base by pysam.
        """
        self.region = region
        self.position = int(position)
        self.samples = samples
        # allele_counts is a dictionary of dictionaries
        # the key is the read group ID/sample_ID  Values is a numpy array of
        # the allele counts
        # So allele_counts['sample'] = [0,14,0,15]
        # BASE_INDEX converts index to base, ie C, in 14
        self.allele_counts = {}
        self.strand_bias = {}
        self.phredThreshold = phredThreshold
        #self.genotype = genotype
        for i in self.samples:
            self.allele_counts[i] = np.zeros(4, dtype=np.uint32)

    def __call__(self, pileups, position=None, samples=None,
        phredThreshold=None, genotype=None):
        """ Genotype file contains the infromation the actualy genotyping
        calls.

        """
        if position == None: position = self.position
        else:  pass
        if samples == None: n = self.samples
        else: pass
        if phredThreshold == None: qualT = self.phredThreshold
        if pileups.pos != position-1: pass
        else:
            for read in pileups.pileups:
                if read.alignment.is_duplicate or read.alignment.mapq <= 50:
                    pass
                else:
                    qpos = read.qpos
                    base_quality = read.alignment.qual[qpos]
                    if ord(base_quality)-33 > qualT:
                        try:
                            read_group = read.alignment.opt('RG')
                        except:
                            # Case where bamFile only has 1 sample
                            read_group = self.samples[0]
                        base = read.alignment.seq[qpos]
                        base = self.BASE_INDEX[base]
                        #strand = read.alignment.is_reverse
                        self.allele_counts[read_group][base] += 1
                    else: pass



def threshold_counts(counts, threshold=30, number=1):
    """ Returns True

    Makes sure that at least one of the samples meet a read count threshold.

    """
    counter = 0
    for i in counts:
        if sum(i) > threshold: counter += 1
    if counter >= number:return True
    else: return False


def main():
    """ Main loop that iterates over the desired loci to calculate Aellic
    Imbalance

    """


    #################################################################
    # Argument and Options Parsing
    #################################################################

    p = optparse.OptionParser(__doc__)
    p.add_option("-o", "--output", dest="filename", help="write \
            report to FILE")
    p.add_option("-G", "--genotype", dest="G", help=\
            "Use imputed/genotypes if available, should be in VCF file format")
    p.add_option("-v",  "--vcf_file", action="store_true", dest="inputisvcfile",
                 help="the input is a VCF file")
    p.add_option("-q", "--quality_threshold", type="int", dest="qual",
                 help="base quality threshold to take allele counts from")
    p.add_option("-p", "--pileup", action="store_true", dest="p",
                 help= "Input files are pileup files")
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-c", "--count-threshold", action="store", type="int",
                 dest="c", help="Set the count threshold for making AEI calls")
    p.add_option("-V", "--output_vcf", action="store_true", dest="outputVCF",
                 help="Output the results to a VCF file")

    options, args = p.parse_args()
    if options.qual: pass
    else: options.qual = 20


    # For testing purposes
    debug = 1

    # Open the bedfile/vcf file
    # file_a = csv.reader(open(args[0], "rU"), delimiter="\t")
    # Right now defaulting to VCF file, since speed is an issue

    """
    if options.inputisvcfile: file_a = VCFfile(args[0])
    else: file_a = BEDfile(args[0])
    """
    file_a = open(args[0],"rb")

    # Handling of multiple BAM/SAM inputs
    bam_Names = args[1:]
    bam_files = []
    for filename in bam_Names:
        bam_files.append(pysam.Samfile(filename,"rb"))
    # Creates a dictionary with the bam_file name as the key, and the samples
    # by readgroup as the value i.e. {"bamfile":[RG1, RG2, RG3]
    # "bamfile2":[RG4,RG5]"

    # Also creates a read group sample dictionary
    readGroup_sample = {}
    bam_ReadGroups = {}
    for bam_file, bamName in map(None, bam_files, bam_Names):
        samples = []
        # Grab only the header information
        m = re.compile('@RG.*')
        readGroups = m.findall(bam_file.text)
        for r in readGroups:
            r = r.split('\t')
            for i in r:
                if i[0:3] == "ID:": ID = i[3:]
                elif i[0:3] == "SM:": SM = i[3:]
                else : pass
            readGroup_sample[ID] = SM
            samples.append(ID)

        if len(bam_files) == 1 and len(readGroups) == 0:
            readGroup_sample['No_SM'] = 'No_ID'
            samples.append('No_SM')
            bam_ReadGroups[bamName] = samples
            break
        elif len(bam_files) == 1 and len(readGroups) == 0:
            print('If you have more than 1 bam file, all the bam files \
                    need sample information from the read group')
            sys.exit()
        else:
            bam_ReadGroups[bamName] == samples

        bam_ReadGroups[bamName] = samples


    # Print the header
    header = ["chr", "pos", "rsID"]
    for i in bam_Names:
        ReadGroupsinBam = bam_ReadGroups[i]
        for t in ReadGroupsinBam:
            header.append(readGroup_sample[t])
            header.append("Genotype(Maj/Min)")
            header.append("Ratio")
        print("\t".join(header))

    INDEX_BASE = ['A', 'C', 'G', 'T']
    if options.c:
        count_threshold = options.c
    else:
        count_threshold = 30
    if options.G != None:
        geno = open(options.G, "rU")
        # :TODO ensure allow different well known file formats and also check
        geno_samples = geno.next().strip('\n').split('\t')[4:]
        # Initialize first row
        geno_line = geno.next().strip('\n').split('\t')
        # match geno_samples to bam_samples
        geno_to_bam_ind = []
    else: pass

    for line in file_a:
        line = line.strip('\n').split('\t')
        counts = []
        # Counts is a list of numpy arrays

        ##################################################################
        #
        # Grab the information for bam files or :TODO pileups.  Seems like
        # something that Hadoop will be very good at.
        #
        ##################################################################
        if options.inputisvcfile:
            region = str(line[0])
            position = int(line[1])
        else:
            region = str(line[0])
            position = int(line[2])

        for bamfile, bamNames in map(None, bam_files, bam_Names):
            # :TODO in the VCF and bed files make sure to type the attributes
            variant = lociInformation(region, position,
                                      samples=bam_ReadGroups[bamNames],
                                      phredThreshold=options.qual)

            bamfile.pileup(variant.region, variant.position,
                           variant.position+1, callback=variant)
            for i in variant.samples:
                counts.append(variant.allele_counts[i])
        # First determines if any of the samples meet the read threshold
        # Secondly determines if there are any heterozygotes in the sample (if
        # there aren't any it skips printing that line.
        # There are several ways it calculates this, if imputed genotypes are
        # given it will use use that,
        # otherwise the posterior probability of being a heterozygote
        # given the data is calculated.
        # Need to map genotype names in the files with the bamfiles or sample
        # Need to move to the same position as where the SNP is.
        # Ideally all this information would be loaded into a database.
        # Maybe load into memory?
        # convert into VCF file
        if options.G != None:
            while geno_line[2] < position:
            # Need to handle edge cases at the end of chromosomes
                geno_line = geno.next().strip('\n').split('\t')
                if region == "Stuff":
                    pass
            if geno_line[2] == position: pass
            # Reorder line to match samples

        if threshold_counts(counts, threshold=count_threshold):
            p_values = []
            for sample_c in counts:
                ind = sample_c.argsort()[::-1]
            # Should it also return the value of the estimate?
                any_hets = []
                if sample_c.sum() >= count_threshold:
                    if lf.isHet(sample_c):
                        p_values.append(lf.ratioLik(sample_c))
                        p_values.append("%s:%s" % (INDEX_BASE[ind[0]],
                                                   INDEX_BASE[ind[1]]))
                        any_hets.append(True)
                    else:
                        p_values.append("HOMO")
                        p_values.append(str(INDEX_BASE[ind[0]]))
                        any_hets.append(False)
                    p_values.append("%s:%s" % (sample_c[ind[0]],
                                               sample_c[ind[1]]))
                else:
                    p_values.append("NA")
                    p_values.append(str(INDEX_BASE[ind[0]]))
                    p_values.append("%s:%s" % (sample_c[ind[0]],
                                               sample_c[ind[1]]))
                    any_hets.append(False)
                if any(any_hets):
                # Only print if there are heterozygotes
                    print("\t".join([variant.region,
                                     str(variant.position),"\t"]) +
                        "\t".join([str(i) for i in list(p_values)]))
                else: pass

        # For testing purposes
        if options.D:
            if debug > 2000: break
            else: debug += 1
        else:pass

if __name__ == '__main__':
    main()
