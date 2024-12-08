##################################
#### Preparing the input file ####
##################################

#### Installing stacks using conda
#conda init
#conda creat -n stacks
#conda activate stacks
#conda install stacks -c bioconda -c conda-forge

#### Preparing the vcf file to stacks, this code will remove SNPs that have more than two alleles
awk '/^#/ || $5 !~ /,/' final_25%_missing_snps.vcf > final_25%_missing_snps_clean.vcf

#### Running stacks
populations --in_vcf final_25%_missing_snps_clean.vcf -O . --treemix -M 6species_map.txt

#### Removing stacks header
echo "$(tail -n +2 final_25%_missing_snps_clean.p.treemix)" > final_25%_missing_snps_clean_2.p.treemix

# Gzip input file for treemix
gzip final_25%_missing_snps_clean_2.p.treemix


##################################
######### Running treemix ########
##################################

#### Install treemix using conda
#conda init
#conda create -n treemix
#conda activate treemix
#conda install treemix -c bioconda -c conda-forge

#### Using treemix on the conda environment
conda init
conda activate treemix

#### Running treemix
cd /desktop/treemix_analysis/

#run treemix with m0
treemix -i input.p.treemix.gz -root aestuansA -o treem0

#add 1 migration edge
treemix -i input.p.treemix.gz -m 1 -g treem0.vertices.gz treem0.edges.gz -o treem1

#add 2 migration edges
treemix -i input.p.treemix.gz -m 1 -g treem1.vertices.gz treem1.edges.gz -o treem2

#add 3 migration edges
treemix -i input.p.treemix.gz -m 1 -g treem2.vertices.gz treem2.edges.gz -o treem3


#### Bootstraps over 200s snps with migration edges
for i in {1..100}; do
    treemix -i input.p.treemix.gz -m 1 -g treem1.vertices.gz treem1.edges.gz -bootstrap -k 100 -o $i.treemix
done;

# unzip the tree files
for i in *treeout.gz; do
    gzip -d $i
done;


#### In R
# summarize bootstraps
x <- list.files(pattern="*treeout")
for(a in 1:length(x)) {
	if (a==1) {
		output <- scan(x[a], what="character")[1]
	} else {
		output <- c(output, scan(x[a], what="character")[1])
	}
}
write(output, file="final.100.bootstraps.trees", ncolumns=1)

#### In bash 
# summarize bootstraps to tree
sumtrees.py --output=final.100.boots.summed.tre --min-clade-freq=0.05 final.100.bootstraps.trees




#### Below is an example run of treemix from a UNIX terminal for M={1-5} and 0 iterations per M:
for m in {1..5}
   do
   treemix \
      -i input.p.treemix.gz \
      -o test.${m} \
      -root aestuansA \
      -m ${m} \
      -k 100
done


#### Below is an example run of treemix from a UNIX terminal for M={1-5} and 5 iterations per M:
for m in {1..5}
   do
   for i in {1..5}
      do
      treemix \
         -i input.p.treemix.gz \
         -o test.${i}.${m} \
         -root aestuansA \
         -m ${m} \
         -k 100
      done 
done
