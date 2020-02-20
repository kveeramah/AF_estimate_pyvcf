# AF_estimate_pyvcf

This python script takes a vcf file and calculates allele frequencies using genotype likelihoods for a subset of individuals at specific loci using the method of Kim et al. BMC Bioinformatics 2011, 12:231 Estimation of allele frequency and associationmapping using next-generation sequencing data

Only biallelic loci with a reference allele and one alternate allele are considered. 

numpy, scipy and pyvcf must be installed in the version of python used.

vcf files must be bgzipped and tabix indexed.

The sample_pop_file file must have the tab separated fields in the following order: sample_name\tpopulation.
Each sample should be on a different line, samples can appear in more than one population.

The snp_position_file must have the tab separated fields in the following order: chromosome\tsnp_position. Each snp should be on a different line.

95% confidence intervals are calculated the using the likelihood profile method.

Usage is ./AF_estimate_pyvcf.py <in.vcf.gz> <sample_pop_file> <snp_position_file> <min_nb_samps>

Output will be to screen and also written to the sample_pop_file name appended with '.AF;


