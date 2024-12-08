#################################################################################################################

#install packages
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("vegan" %in% rownames(installed.packages()) == FALSE){install.packages("vegan")
} else {print (paste0("'vegan' has already been installed in library"))}
if("ade4" %in% rownames(installed.packages()) == FALSE){install.packages("ade4")
} else {print (paste0("'ade4' has already been installed in library"))}
if("AssocTests" %in% rownames(installed.packages()) == FALSE){install.packages("AssocTests")
} else {print (paste0("'AssocTests' has already been installed in library"))}
if("fossil" %in% rownames(installed.packages()) == FALSE){install.packages("fossil")
} else {print (paste0("'fossil' has already been installed in library"))}
if("ecodist" %in% rownames(installed.packages()) == FALSE){install.packages("ecodist")
} else {print (paste0("'ecodist' has already been installed in library"))}
if("reshape2" %in% rownames(installed.packages()) == FALSE){install.packages("reshape2")
} else {print (paste0("'reshape2' has already been installed in library"))}
if("dartR" %in% rownames(installed.packages()) == FALSE){install.packages("dartR")
} else {print (paste0("'dartR' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}
if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("LEA")

#load packages
pacman::p_load("adegenet", "dartR", "vegan","ade4", "AssocTests", 'r2vcftools', 'vcfR', 'dplyr', "tibble", "LEA", "reshape2", "ggplot2", "fossil")

#################################################################################################################
############################################# INPUT FILES #######################################################
# 1- structure file with the genomic data 
# 2- cvs file with the geographic coordinates of each individual

setwd("~/") #set working directory

#################################################################################################################
################################################# PCA ###########################################################

setwd("./PCA/")

#Load .str file:
str <- read.structure('final_25%_missing_snps_edit.str', n.ind = 64, n.loc = 1169, col.lab = 1, 
				col.pop = 2, onerowperind = TRUE, row.marknames = 0, NA.char= 0, 
				ask = FALSE)
str2 <- scaleGen (str, cent=TRUE,scale=TRUE,NA.method = "mean")

pca <- prcomp(str2,center = TRUE, scale. =TRUE)
summary(pca)

pca$x

write.csv(pca$x, file = "scores_pca.csv") # save pca scores
plot(pca$x[1,], pca$x[2,])

###Plot PCA 1n2 #### 
tiff("pca1n2_moleculas_guerlinguetus.tiff", width=20, height=15, unit="cm", res=300)
quartz.options(height=10, width=12, dpi=72);
plot.new();
par(oma=c(1,1,1,1));
par(mar=c(5,5,1,1));
plot.window(xlim=c(-40,70), ylim=c(-40, 40));

points(pca$x[64:64,1],pca$x[64:64,2], col = 'peachpuff', bg = 'peachpuff', cex = 1.5, pch=21) #Rio Grande do Sul 
points(pca$x[62:63,1],pca$x[62:63,2], col = 'pink4', bg = 'pink4', cex = 1.5, pch=21) #Santa Catarina 
points(pca$x[50:61,1],pca$x[50:61,2], col = 'salmon2', bg = 'salmon2', cex = 1.5, pch=21) #Sao Paulo 
points(pca$x[46:49,1],pca$x[46:49,2], col = 'rosybrown1', bg = 'rosybrown1', cex = 1.5, pch=21) #Minas Gerais 
points(pca$x[41:45,1],pca$x[41:45,2], col = 'plum2', bg = 'plum2', cex = 1.5, pch=21) #Espírito Santo 
points(pca$x[38:40,1],pca$x[38:40,2], col = 'purple', bg = 'purple', cex = 1.5, pch=21) #Bahia 
points(pca$x[37:37,1],pca$x[37:37,2], col = 'goldenrod', bg = 'goldenrod', cex = 1.5, pch=21) #Mato Grosso 
points(pca$x[30:36,1],pca$x[30:36,2], col = 'darkolivegreen3', bg ='darkolivegreen3', cex = 1.5, pch=21) #Para brasiliensis
points(pca$x[10:29,1],pca$x[10:29,2], col = 'forestgreen', bg = 'forestgreen', cex = 1.5, pch=21) #Para aestuans 
points(pca$x[8:9,1],pca$x[8:9,2], col = 'darkseagreen4', bg = 'darkseagreen4', cex = 1.5, pch=21) #Amazonas
points(pca$x[7:7,1],pca$x[7:7,2], col = 'cyan', bg = 'cyan', cex = 1.5, pch=21) #Amapá
points(pca$x[3:6,1],pca$x[3:6,2], col = 'darkblue', bg = 'darkblue', cex = 1.5, pch=21) #Guyanas
points(pca$x[1:2,1],pca$x[1:2,2], col = 'blue', bg = 'blue', cex = 1.5, pch=21) #Venezuela

axis(1, at=seq(-40, 40, by=20), cex.axis=1.15);
axis(2, at=seq(-40, 40, by=20), cex.axis=1.15, las=1);

mtext(side=1, text='PC1 (25.89%)',line=2.5, cex=1)
mtext(side=2, text='PC2 (10.43%)', line=2.8, cex=1)
legend<- c("Venezuela", "Guyanas", "Amapá", "Amazonas", "Para aestuans", "Para brasiliensis", "Mato Grosso",
           "Bahia", "Espírito Santo", "Minas Gerais", "Sao Paulo", "Santa Catarina", "Rio Grande do Sul" )
col=c('blue', 'darkblue','cyan','darkseagreen4','forestgreen','darkolivegreen3','goldenrod','purple',
      'plum2','rosybrown1','salmon2','pink4','peachpuff')
legend("topright", legend = legend, cex=0.7, bty="n",  col = col, pch= 16)
dev.off()

#################################################################################################################
################################### TRACY-WIDOM TEST FOR EIGENVALUES ############################################

#Tracy CA and Widom H. (1994). Level spacing distributions and the bessel kernel. Commun Math Phys. 161 :289--309.
#Patterson N, Price AL and Reich D. (2006). Population structure and eigenanalysis. PLoS Genet. 2 :20.

eigenvalues<-pca$sdev^2

eigenL<- length(eigenvalues)

#criticalpoint: a numeric value corresponding to the significance level. If the significance level is 0.05, 0.01, 
#0.005, or 0.001,the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, accordingly. The default
# is 2.0234

tw<- tw(eigenvalues, eigenL, criticalpoint = 0.9793)

tw$SigntEigenL #the number of significant eigenvalues

setwd('..') #move directory backward

#################################################################################################################
################################################# DAPC ##########################################################


setwd("./DAPC/")

#Load .vcf file after removing outlier SNPs:
vcf = read.vcfR("final_25%_missing_snps.vcf", verbose = FALSE)

#Convert "VCF" to "GENIND"
input = vcfR2genind(vcf)
input

#Perform a PCA to choose the number of PC in the DAPC:
input_scaled = scaleGen (input, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE)

# % of PC variation
pc_pca = as.data.frame(pca_input$eig)
pc_pca[,2] = (pc_pca/sum(pc_pca))*100

# Rule of 100% of variance:
index_100 = length(rownames(input@tab))-1
index_100 # number of PC to reach to 100% of explained variance

# Rule of at least 95% of variance:
#index_95 = length(which(cumsum(pc_pca[,2]) <= 95))
#index_95

# Rule of at least 70% of variance:
#index_70 = length(which(cumsum(pc_pca[,2]) <= 70))
#index_70 

#Rule of minimum variance per PCs:
#variance_pc = 100/(nrow(input@tab)-1)
#variance_pc #PCs that increase the explained variance bellow this threshold will be removed
#calculate number of PCs
#index_min = length(pc_pca[,2][pc_pca[,2] >= variance_pc])
#index_min

#Identification of the clusters (We specify that we want to evaluate up to k = 40 groups (max.n.clust = 40)
#If you see a plateau on the graph, you can choose a number of PCs in the beginning of this plateau. If there is no plateau test different number of PC's
set.seed(1) #set a seed
grp = find.clusters(input, max.n.clust=15, scale = TRUE) # center=T by default.
#Digit the number of PCs to retain. Complete by the results (index_100 = 63)

#grp_index_95 = find.clusters(index_95, max.n.clust=40, scale = TRUE)

#grp_index_70 = find.clusters(index_70, max.n.clust=40, scale = TRUE)

#grp_index_min = find.clusters(index_min, max.n.clust=40, scale = TRUE)


head(grp$grp, 10)
grp$size
#[1]  7  8  3  5  8  3  7 19  4

#Select the ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
#save best k graph
dapc_df = as.data.frame(grp$Kstat) %>%
  mutate(K = c(1:15)) %>%
  setNames(c("BIC", "K"))

p = ggplot(data=dapc_df, mapping = aes(x=K,y= BIC)) +
  geom_line(size = 1) +
  geom_point(size = 5) +
  scale_x_continuous(breaks = c(1:40)) +
  theme_bw() +
  xlab("Best K") + ylab("Bayesian Information Criterion (BIC)") #set labels

p

pdf("./Bestk_DAPC.pdf", onefile = T)
p
dev.off()

#Choose the best PCs number to recover correctly the clusters. 

pdf("./Best_PCs_Number_DAPC.pdf", onefile = T)
number_PCs = xvalDapc(tab(input, NA.method = "mean"), grp$grp, scale = T, n.pca.max = (nrow(input@tab)-1))
dev.off()


#Verify the number of PCs and DA used and summary of DAPC, and percentage explained by the three first dapc
number_PCs$DAPC$n.pca
number_PCs$DAPC$n.da
summary(number_PCs$DAPC)
(number_PCs$DAPC$eig[1]/sum(number_PCs$DAPC$eig))*100
(number_PCs$DAPC$eig[2]/sum(number_PCs$DAPC$eig))*100
(number_PCs$DAPC$eig[3]/sum(number_PCs$DAPC$eig))*100

#Verify plots and define the group colors and names:
#plot graph
table(grp$grp)
#[1]  7  8  3  5  8  3  7 19  4


snps = vcfLink("final_25%_missing_snps.vcf", overwriteID = T)
snps@meta[which(grp$grp == 1),] 
snps@meta[which(grp$grp == 2),] 
snps@meta[which(grp$grp == 3),] 
snps@meta[which(grp$grp == 4),] 
snps@meta[which(grp$grp == 5),] 
snps@meta[which(grp$grp == 6),] 
snps@meta[which(grp$grp == 7),] 
snps@meta[which(grp$grp == 8),] 
snps@meta[which(grp$grp == 9),] 

#define colors and name for populations
col=c('#CC6633', '#FF99FF','#FF3399','#FF9966','#9900CC','#33CC99','#003399','#006633','#00CCFF')
       #color for the genetic cluster

labels_dapc = c("POP A", "POP B", "POP C", "POP D","POP E","POP F","POP G","POP H","POP I") #name for the genetic cluster

#ggplot graph:
df_dapc = number_PCs$DAPC$ind.coord %>%
  as.data.frame() %>%
  rownames_to_column("IND") %>%
  mutate(POP = as.data.frame(grp$grp)[,1])


d_gg = ggplot(df_dapc, aes(x=LD1, y=LD2, fill=POP))+
  geom_point(size=4, pch=21)+
  scale_fill_manual(values= col,
                    labels= labels_dapc) +
  theme_bw()+
  xlab("DAPC1 (41.71%)") + ylab("DAPC2 (29.40%)") #set labels

d_gg

pdf("./scatter_DAPC_ggplot.pdf", onefile = T)
d_gg
dev.off()

#classic plot graph
pdf("./scatter_DAPC_classical.pdf", onefile = T)
tiff("./scatter_DAPC_classical.tiff")
scatter(number_PCs$DAPC, cex =3, legend = TRUE, col = col, solid=0.7, txt.leg = labels_dapc,
        clabel = FALSE, posi.leg = "bottomright", pch=20, scree.da=TRUE,scree.pca=TRUE, 
        posi.pca = "bottomleft", posi.da = "topright", cleg = 1, xax = 1, yax = 2, inset.solid = 1)
dev.off()

setwd('..') #move directory backward

#################################################################################################################
################################################# sNMF ##########################################################

setwd("./sNMF")

#Load .vcf file after removing outlier SNPs:
snps = vcfLink("final_25%_missing_snps.vcf", overwriteID = T) 
head(snps@meta)

ord_samples_phylog <- read.table("phylo_samp_ord.txt") #phylo order to organizar the plot
colnames(ord_samples_phylog)<- "phylo"
snps@meta = cbind(snps@meta, ord_samples_phylog)
head(snps@meta)

vcf2geno('final_25%_missing_snps.vcf')


#B. Create folders for alpha values and copy .geno object in each folder:
#the regularization parameter of snmf was increased to α = 1000 (default value α = 10) to account for the relatively small number of loci (Frichot, et al., 2014)

alpha_values = c(10, 50, 100, 500) #as we have less than 10000 snps we choose alpha values greater than 1000
for (i in alpha_values){
  path = paste0("./Results/Alpha", i, "n")
  if (dir.exists(file.path(getwd(), path)) == FALSE)
  {dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))
  file.copy(paste0("final_25%_missing_snps.geno"), path )
}


#Set parameters to run SNMF (LEA) using different alpha.
K = c(1:10) # K to be tested
replications = 10 # number of replication by K
ploidy = 2 # species ploidy
CPU = 2 #Number of cores

#Run sNMF (LEA) using different alpha.
set.seed(136)
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("Results/Alpha", i,"n/","final_25%_missing_snps.geno")
  pro_snmf = snmf(path, K = K, rep = replications, alpha = i, entropy = T, ploidy = ploidy , project = "new", CPU= CPU)
  assign(paste0("project_snmf", loop, "n"), pro_snmf)
}


#To load the SNMF projects in a new R session (after quitting R).
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
loop = loop +1
path = paste0("Results/Alpha", i,"n/","final_25%_missing_snps.snmfProject")
pro_snmf = load.snmfProject(path)
assign(paste0("project", loop, "n"), pro_snmf)
}

#Summary of the project
summary(project1n)
summary(project2n)
summary(project3n)
summary(project4n)
#summary(project5n)
#summary(project6n)
#summary(project7n)
#summary(project8n)
#summary(project9n)

#View Cross-Entropy plots

PlotK(project1n) #5
PlotK(project2n) #6 
PlotK(project3n) #6 
PlotK(project4n) #5
#PlotK(project5n) 
#PlotK(project6n) 
#PlotK(project7n) 
#PlotK(project8n) 
#PlotK(project9n) 
#PlotK(project10n) 
#PlotK(project11n) 
#PlotK(project12n) 
#PlotK(project13n) 


#Save graphs of cross-Entropy plot with standard deviation error bars
for (i in alpha_values){
  pdf(paste0("./Results/Cross_Entropy_sNMF_Alpha2_",  i, "n.pdf"), onefile = F)
  path = paste0("Results/Alpha", i,"n/", "final_25%_missing_snps.snmfProject")
  print(PlotK(load.snmfProject(path)))
  dev.off()
}

#Select optimal K value
optimal_K = 5

#Select best run (lowest cross-entropy)
best_run = Best.run(nrep=10, optimalK=optimal_K, p1=project1n, p2=project2n, p3=project3n,
                    p4=project4n)
#load the best project
best_run_split = scan(text = gsub('[[:punct:] ]+',' ', best_run), what = "")
path_best_run = paste0("Results/Alpha", alpha_values[as.numeric(best_run_split[5])],"n/", "final_25%_missing_snps.snmfProject")
#set the values
project = load.snmfProject(path_best_run)
run=as.numeric(best_run_split[7])


pdf("./Results/sNMF_Guerlgraph2.pdf", onefile = F)
my.colors = rainbow(optimal_K)
LEA::barchart(project, K = optimal_K, run = run, border = NA, space = 0, col = my.colors, xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1, cex.axis = .3)
dev.off()

#Add admixture coefficient and replace the population ID to vcf file
Qmat_k5 =
  Q(project, run=run, K=optimal_K) %>%
  as.data.frame() %>%
  setNames(c("Adx_Coeff_1", "Adx_Coeff_2","Adx_Coeff_3","Adx_Coeff_4","Adx_Coeff_5")) %>% #add here the number of k
  mutate(PopID_snmf = apply(., 1, which.max))
head(Qmat_k5)

snps@meta = cbind(snps@meta, Qmat_k5)
head(snps@meta)

write.csv(snps@meta, "res_sNMF_K5.csv")

#Barplot in ggplot2:
#define colors and name for populations
#colors_guerl = c('#00CCFF', '#9900CC','#33CC99','#003399','#FF99FF')

colors_guerl = c('#CC6633','#9900CC',  '#FF3399', '#FF99FF', '#003399')


#marrom, roxo, pink, rosa, azul escuro

#color for thegenetic cluster
labels_guerl = c("POP A", "POP D", "POP C", "POP E", "POP B") #name for the genetic cluster


#create a dataframe (df) for ggplot2, all localities
df = snps@meta %>%
  dplyr::select(sample_num, phylo, Adx_Coeff_1, Adx_Coeff_2, Adx_Coeff_3, Adx_Coeff_4, Adx_Coeff_5) %>%
  dplyr::arrange(phylo) %>%
  melt(., id.vars=c("sample_num", "phylo"))
head(df)

p = ggplot(data=df, mapping = aes(x=reorder(factor(sample_num), phylo),
                                  y= value*100,
                                  fill = factor(variable))) +
  geom_bar(stat="identity", width = 1, size=0.3, colour="black") +
  scale_fill_manual(values= colors_guerl,
                    labels= labels_guerl) +
  theme_minimal() +
  xlab("") + ylab("")

p

pdf(paste0("./Results/Guerl_ggplot2_snmf_K5", ".pdf"), onefile =F)
plot(p)
dev.off()


#Select optimal K value
optimal_K = 6

#Select best run (lowest cross-entropy)
best_run = Best.run(nrep=10, optimalK=optimal_K, p1=project1n, p2=project2n, p3=project3n,
                    p4=project4n)
#load the best project
best_run_split = scan(text = gsub('[[:punct:] ]+',' ', best_run), what = "")
path_best_run = paste0("Results/Alpha", alpha_values[as.numeric(best_run_split[5])],"n/", "final_25%_missing_snps.snmfProject")
#set the values
project = load.snmfProject(path_best_run)
run=as.numeric(best_run_split[7])


pdf("./Results/sNMF_Guerlgraph2.pdf", onefile = F)
my.colors = rainbow(optimal_K)
LEA::barchart(project, K = optimal_K, run = run, border = NA, space = 0, col = my.colors, xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1, cex.axis = .3)
dev.off()


#Add admixture coefficient and replace the population ID to vcf file
Qmat_k6 =
  Q(project, run=run, K=optimal_K) %>%
  as.data.frame() %>%
  setNames(c("Adx_Coeff_1", "Adx_Coeff_2","Adx_Coeff_3","Adx_Coeff_4","Adx_Coeff_5","Adx_Coeff_6")) %>% #add here the number of k
  mutate(PopID_snmf = apply(., 1, which.max))
head(Qmat_k6)


snps = vcfLink("final_25%_missing_snps.vcf", overwriteID = T) 
snps@meta = cbind(snps@meta, ord_samples_phylog)
snps@meta = cbind(snps@meta, Qmat_k6)
head(snps@meta)

write.csv(snps@meta, "res_sNMF_K6.csv")

#Barplot in ggplot2:
#define colors and name for populations
colors_guerl = c('#CC6633','#9900CC',  '#FF3399', '#FF99FF', '#003399', "yellow")

colors_guerl = c('#CC6633','#33CC99',  '#FF3399',  '#003399', '#FF3399',  "yellow")

#marrom, roxo, pink, rosa, azul escuro

#color for thegenetic cluster
labels_guerl = c("POP A", "POP D", "POP C", "POP E", "POP B") #name for the genetic cluster


#create a dataframe (df) for ggplot2, all localities
df = snps@meta %>%
  dplyr::select(sample_num, phylo, Adx_Coeff_1, Adx_Coeff_2, Adx_Coeff_3, Adx_Coeff_4, Adx_Coeff_5, Adx_Coeff_6) %>%
  dplyr::arrange(phylo) %>%
  melt(., id.vars=c("sample_num", "phylo"))
head(df)

p = ggplot(data=df, mapping = aes(x=reorder(factor(sample_num), phylo),
                                  y= value*100,
                                  fill = factor(variable))) +
  geom_bar(stat="identity", width = 1, size=0.3, colour="black") +
  scale_fill_manual(values= colors_guerl,
                    labels= labels_guerl) +
  theme_minimal() +
  xlab("") + ylab("")

p

pdf(paste0("./Results/Guerl_ggplot2_snmf_K6b", ".pdf"), onefile =F)
plot(p)
dev.off()

setwd('..') #move directory backward

#################################################################################################################
################################################ Mantel #########################################################

setwd("./Mantel/")

### Genetic Data  #####

str <- read.structure('project_data.str', n.ind = 64, n.loc = 1169, col.lab = 1, col.pop = 2, onerowperind = TRUE, 
						row.marknames = 0, NA.char= 0, ask = FALSE)
str2 <- scaleGen (str, cent=TRUE,scale=TRUE,NA.method = "mean")
pca <- prcomp(str2)
summary(pca)

### Geographic data #####
coord <- read.csv("coord_ind.csv", header = T)

### Amazon Samples

dir.create("./Amazon/")
setwd("./Amazon/")


### Genetic Distance #####

PCs<- pca$x[1:37,1:14] ### cobre 70% da variação

PCA_dist<-ecodist::distance(PCs, method = "mahalanobis", sprange=NULL, spweight=NULL) 

PCA_dist2<-as.matrix(dist(PCA_dist))

write.csv(PCA_dist2, "PCA_gen_dist_amazon.csv")

### Geographic Distance #####

coord_amazon <- coord[1:37,2:3]

mDIST <- earth.dist(coord_amazon, dist = TRUE)
mDIST <- as.dist(mDIST)
mDIST_vec <- as.vector(mDIST)
hist(mDIST_vec)

#### Perform Mantel Test #####

mantel <- mantel.rtest(PCA_dist, mDIST, 10000)

#Monte-Carlo test
#Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)

#Observation: 0.1717168 

#Based on 10000 replicates
#Simulated p-value: 0.06579342 
#Alternative hypothesis: greater 

#    Std.Obs Expectation    Variance 
#1.530211125 0.001566613 0.012364093 

plot.new()
tiff("mantel_amazon.tiff", width=50, height=40, unit="cm", res=300)
par(mar=c(5,5,1,1))
plot(mDIST, PCA_dist, pch = 16,  ylab = "Genetic Distance", xlab = "Geographic Distance", cex = 1.5, cex.lab=1.5, axes= T) 
abline(lm(PCA_dist  ~ mDIST), col = "gray30", lwd = 3)
legend('topleft', legend = "r = 0.17, p=value = 0.06", bty = 'n', cex = 1)
dev.off()


setwd('..') #move directory backward

### Atlantic Forest Samples

dir.create("./Atlantic Forest/")
setwd("./Atlantic Forest/")


### Genetic Distance #####

PCs<- pca$x[38:64,1:14] ### cobre 70% da variação

PCA_dist<-ecodist::distance(PCs, method = "mahalanobis", sprange=NULL, spweight=NULL) 

PCA_dist2<-as.matrix(dist(PCA_dist))

write.csv(PCA_dist2, "PCA_gen_dist_AF.csv")

### Geographic Distance #####

coord_AF <- coord[38:64,2:3]

mDIST <- earth.dist(coord_AF, dist = TRUE)
mDIST <- as.dist(mDIST)
mDIST_vec <- as.vector(mDIST)
hist(mDIST_vec)

#### Perform Mantel Test #####

mantel <- mantel.rtest(PCA_dist, mDIST, 10000)

#Monte-Carlo test
#Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)

#Observation: 0.2346292 

#Based on 10000 replicates
#Simulated p-value: 0.00969903 
#Alternative hypothesis: greater 

#     Std.Obs  Expectation     Variance 
# 2.301541769 -0.001393977  0.010516504 

plot.new()
tiff("mantel_AF.tiff", width=50, height=40, unit="cm", res=300)
par(mar=c(5,5,1,1))
plot(mDIST, PCA_dist, pch = 16,  ylab = "Genetic Distance", xlab = "Geographic Distance", cex = 1.5, cex.lab=1.5, axes= T) 
abline(lm(PCA_dist  ~ mDIST), col = "gray30", lwd = 3)
legend('topleft', legend = "r = 0.23, p=value = 0.009", bty = 'n', cex = 1)
dev.off()

setwd('..') #move directory backward

### Mantel per Population

# Pop 1 
# only two individuals
# more variables than observations

# Pop 2
# only two individuals
# more variables than observations

# Pop 3
dir.create("./Pop3/")
setwd("./Pop3/")

### Genetic Distance #####

PCs<- pca$x[6:20,1:14] ### cobre 70% da variação

PCA_dist<-ecodist::distance(PCs, method = "mahalanobis", sprange=NULL, spweight=NULL) 

PCA_dist2<-as.matrix(dist(PCA_dist))

write.csv(PCA_dist2, "PCA_gen_dist_Pop3.csv")

### Geographic Distance #####

coord_pop3 <- coord[6:20,2:3]

mDIST <- earth.dist(coord_pop3, dist = TRUE)
mDIST <- as.dist(mDIST)
mDIST_vec <- as.vector(mDIST)
hist(mDIST_vec)

#### Perform Mantel Test #####

mantel <- mantel.rtest(PCA_dist, mDIST, 10000)

#Monte-Carlo test
#Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)

#Observation: 0.3833144 

#Based on 10000 replicates
#Simulated p-value: 0.00569943 
#Alternative hypothesis: greater 

#Std.Obs   Expectation      Variance 
#2.9390152329 -0.0005360443  0.0170577006  

plot.new()
tiff("mantel_Pop3.tiff", width=50, height=40, unit="cm", res=300)
par(mar=c(5,5,1,1))
plot(mDIST, PCA_dist, pch = 16,  ylab = "Genetic Distance", xlab = "Geographic Distance", cex = 1.5, cex.lab=1.5, axes= T) 
abline(lm(PCA_dist  ~ mDIST), col = "gray30", lwd = 3)
legend('topleft', legend = "r = 0.38, p=value = 0.005", bty = 'n', cex = 1)
dev.off()

setwd('..') #move directory backward

# Pop 4
# 10 individuals
# more variables than observations

# Pop 5
# 7 individuals
# more variables than observations

# Pop 6
dir.create("./Pop6/")
setwd("./Pop6/")


### Genetic Distance #####

PCs<- pca$x[38:64,1:14] ### cobre 70% da variação

PCA_dist<-ecodist::distance(PCs, method = "mahalanobis", sprange=NULL, spweight=NULL) 

PCA_dist2<-as.matrix(dist(PCA_dist))

write.csv(PCA_dist2, "PCA_gen_dist_Pop6.csv")

### Geographic Distance #####

coord_pop6 <- coord[38:64,2:3]

mDIST <- earth.dist(coord_pop6, dist = TRUE)
mDIST <- as.dist(mDIST)
mDIST_vec <- as.vector(mDIST)
hist(mDIST_vec)

#### Perform Mantel Test #####

mantel <- mantel.rtest(PCA_dist, mDIST, 10000)

#Monte-Carlo test
#Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)

#Observation: 0.2346292 

#Based on 10000 replicates
#Simulated p-value: 0.00879912 
#Alternative hypothesis: greater 

#Std.Obs   Expectation      Variance 
#2.330237e+00 -4.591004e-05  1.014224e-02  

plot.new()
tiff("mantel_Pop6.tiff", width=50, height=40, unit="cm", res=300)
par(mar=c(5,5,1,1))
plot(mDIST, PCA_dist, pch = 16,  ylab = "Genetic Distance", xlab = "Geographic Distance", cex = 1.5, cex.lab=1.5, axes= T) 
abline(lm(PCA_dist  ~ mDIST), col = "gray30", lwd = 3)
legend('topleft', legend = "r = 0.23, p=value = 0.008", bty = 'n', cex = 1)
dev.off()

setwd('..') #move directory backward


#################################################################################################################
############################################ Fst par a par ######################################################

setwd("./Fst/")

str <- read.structure('str_w_pops.str', n.ind = 64, n.loc = 1169, col.lab = 1, col.pop = 2, onerowperind = TRUE, 
						row.marknames = 0, NA.char= 0, ask = FALSE)
input_fst <- gi2gl (str) 

### fst based on Weir and Cocker- ham (1984) ###

FST_stat = gl.fst.pop(input_fst, nboots = 1000, percent = 95, nclusters = 1)

FST_stat$Fsts

#       1         2         3         4         5  6
#1        NA        NA        NA        NA        NA NA
#2 0.6995642        NA        NA        NA        NA NA
#3 0.3975508 0.4200179        NA        NA        NA NA
#4 0.4068067 0.4541672 0.3848286        NA        NA NA
#5 0.6616520 0.6839179 0.5324958 0.4852293        NA NA
#6 0.7162045 0.7255604 0.6708395 0.6503799 0.4810281 NA



Pvalues<-FST_stat$Pvalues

#bonferroni correction

#bonferroni<- p.adjust(Pvalues, method = "bonferroni", n = length(Pvalues))

setwd('..') #move directory backward