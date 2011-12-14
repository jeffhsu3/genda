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
