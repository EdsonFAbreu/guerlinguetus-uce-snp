library(vcfR)

source("HDplot.R")

vcfInput<-read.vcfR("thinned_snps.vcf")

HDplotResults<-HDplot(vcfInput)
head(HDplotResults)

#plot H and D
HDplotResults %>% ggplot()+geom_point(aes(x=H,y=D))

#plot H and ratio
HDplotResults %>% ggplot()+geom_point(aes(x=H,y=ratio))


snp_to_exclude<- subset(HDplotResults, H > 0.75 & D > 10)

#plot H and D
snp_to_exclude %>% ggplot()+geom_point(aes(x=H,y=D))

write.table(snp_to_exclude[,1:2], file="snp_to_exclude.txt", sep = '\t', row.names = F, col.names = F)
