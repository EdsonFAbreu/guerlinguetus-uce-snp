## Converting a .vcf (Variant Call Format) file into a .str (STRUCTURE format) file using plink software:

plink --vcf input.vcf --recode structure --out output

# PGDSpider can also be used for this.

## STRUCTURE ANALYSIS
# STRUCTURE requirements:

1. A genotype file (.str format).
	1st column = Sample ID
	Subsequent columns = Genotypes (-9 for missing data)
	* A column with assigment of samples to populations can be added after the 1st column and before the genotype columns.

2. Two parameter files (mainparams and extraparams) to set the run.

# To run STRUCTURE:

./structure -i input.str -m mainparams -e extraparams -o output.txt -K 3

# Set the value of K (number of clusters) according to your hypothesized number of putative populations. The recommended value is always the number of putative populations plus one.

# STRUCTURE results:

STRUCTURE generates:

output_f: Summary statistics
output_results: Cluster assignments

For visualizing results use:

STRUCTURE Harvester (http://taylor0.biology.ucla.edu/structureHarvester/) to determine the best K.

CLUMPAK (http://clumpak.tau.ac.il/) for clustering consistency.