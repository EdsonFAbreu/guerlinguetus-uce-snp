#### Script to plot treemix output ####

#step 1: copy the entire treemix outdirectory from the KU cluster to my local machine

#step 2: source plotting functions that are distributed with treemix
source("plotting_funcs.R")

#step 3: move into the treemix output directory and plot trees
#setwd("~/Desktop/aph.data/treemix/")
#setwd("~/Downloads/7spec/")

#0 edge
plot_tree("treem0")
#1 edge
plot_tree("test.1.1", plus = 0.02, arrow=.1, ybar = 0, scale=F, lwd=1.5)
plot_tree("test.1.6", plus = 0.02, arrow=.1, ybar = 0, scale=F, lwd=1.5)
plot_tree("test.1.11", plus = 0.02, arrow=.1, ybar = 0, scale=F, lwd=1.5)
#2 edges
plot_tree("test.2.1")
plot_tree("test.2.6")
plot_tree("test.2.11")
#3 edges
plot_tree("test.3.1")
plot_tree("test.3.6")
plot_tree("test.3.11")

#plot to see how much variance is explained by each edge
#m=NULL
#for(i in 1:5){
#  m[i+1] <- get_f(paste0("test.1.",i))
#}
#
#m
#plot(seq(0,5),m,pch="*",cex=2,col="blue", type="b",xlab="migration edge number", ylab="% explained variance")

pdf(file="treemix1.pdf", width=4.5, height = 4)
plot_tree("test.1.1", plus = 0.02, arrow=.1, ybar = 0.3, scale=F, lwd=1.2)
dev.off()

pdf(file="treemix2.pdf", width=4, height = 4)
plot_tree("test.2.1", plus = 0.02, arrow=.1, ybar = 0.3, scale=F, lwd=1.2)
dev.off()

pdf(file="treemix3.pdf", width=4, height = 4)
plot_tree("test.3.1", plus = 0.02, arrow=.1, ybar = 0.3, scale=F, lwd=1.2)
dev.off()

#pdf(file="variance_explained.pdf", width=4, height = 4)
#plot(seq(0,5),m,pch="*",cex=2,col="blue", type="b",xlab="migration edge number", ylab="% explained variance")
#dev.off()


####
pop.uniq <- c("aestuansA","aestuansB","aestuansC1","aestuansC2","brasiliensisA","brasiliensisB")
write.table(pop.uniq, "poporder", sep = " ", quote = F, row.names = F, col.names = F)

plot_resid("test.1.1", "poporder")
plot_resid("treem2", "poporder")
plot_resid("treem3", "poporder")


####
install.packages("SiZer")
install.packages("OptM")
library(OptM)
library(SiZer)

#folder <- system.file("extdata", package = "OptM")
#test.optM = optM(folder)
#plot_optM(test.optM, method = "Evanno")

test.optM = optM("~/Desktop/treemix_test_k100+k200+k300/")
plot_optM(test.optM, method = "Evanno")

#pdf(file="plot_optM.pdf", width=4, height = 8)
#plot_optM(test.optM, method = "Evanno")
#dev.off()

test.linear = optM("~/Desktop/treemix_test_k100+k200+k300/", method = "linear")
plot_optM(test.linear, method = "linear")

test.sizer = optM("~/Desktop/treemix_test_k100+k200+k300/", method = "SiZer")
plot_optM(test.sizer, method = "SiZer")
