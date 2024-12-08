###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

###############################################################################
### AUTHORED BY: CAROLINA S. CARVALHO, BRENNA R. FORESTER, SIMONE K. MITRE, ###
########## RONNIE ALVES, VERA L. IMPERATRIZ-FONSECA, SILVIO J. RAMOS, #########
##### LUCIANA C. RESENDE-MOREIRA, JOSÉ 0. SIQUEIRA, LEONARDO C. TREVELIN, #####
############# CECILIO F. CALDEIRA, MARKUS GASTAUER, RODOLFO JAFFÉ #############
###############################################################################


VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

PlotK <- function(snmfproject){
  
  Mk1 <- mean(cross.entropy(snmfproject, K=1))
  SDk1 <- sd(cross.entropy(snmfproject, K=1))
  Mk2 <- mean(cross.entropy(snmfproject, K=2))
  SDk2 <- sd(cross.entropy(snmfproject, K=2))
  Mk3 <- mean(cross.entropy(snmfproject, K=3))
  SDk3 <- sd(cross.entropy(snmfproject, K=3))
  Mk4 <- mean(cross.entropy(snmfproject, K=4))
  SDk4 <- sd(cross.entropy(snmfproject, K=4))
  Mk5 <- mean(cross.entropy(snmfproject, K=5))
  SDk5 <- sd(cross.entropy(snmfproject, K=5))
  Mk6 <- mean(cross.entropy(snmfproject, K=6))
  SDk6 <- sd(cross.entropy(snmfproject, K=6))
  Mk7 <- mean(cross.entropy(snmfproject, K=7))
  SDk7 <- sd(cross.entropy(snmfproject, K=7))
  Mk8 <- mean(cross.entropy(snmfproject, K=8))
  SDk8 <- sd(cross.entropy(snmfproject, K=8))
  Mk9 <- mean(cross.entropy(snmfproject, K=9))
  SDk9 <- sd(cross.entropy(snmfproject, K=9))
  Mk10 <- mean(cross.entropy(snmfproject, K=10))
  SDk10 <- sd(cross.entropy(snmfproject, K=10))
  
  CE <- data.frame(K=c(1:10), Mean = c(Mk1,Mk2,Mk3,Mk4,Mk5,Mk6,Mk7,Mk8,Mk9,Mk10),
                   SD = c(SDk1,SDk2,SDk3,SDk4,SDk5,SDk6,SDk7,SDk8,SDk9,SDk10))
  
  library(ggplot2)
  
  ggplot(CE, aes(x=K, y=Mean)) + 
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2)+
    geom_line() + 
    geom_point(size=4, shape=21, fill="red", color="darkred") + xlab("Number of ancestral populations") + ylab("Cross-entropy")+
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, vjust=0.5), axis.text = element_text(size=15, face ="bold" , color = "black"), axis.title.x = element_text(size=15, face="bold", color="black"),axis.title.y = element_text(size=15, face="bold", color="black")) +
    scale_x_continuous(breaks = seq(0,10, by=2))
}

Best.run <- function(nrep, optimalK, p1, p2, p3, p4){
  ce1 = LEA::cross.entropy(p1, K = optimalK) # get the cross-entropy of each run for optimal K
  ce2 = LEA::cross.entropy(p2, K = optimalK) # get the cross-entropy of each run for optimal K
  ce3 = LEA::cross.entropy(p3, K = optimalK) # get the cross-entropy of each run for optimal K
  ce4 = LEA::cross.entropy(p4, K = optimalK) # get the cross-entropy of each run for optimal K
  AllProjects <- rbind(ce1, ce2, ce3, ce4)
  rownames(AllProjects) <- NULL
  AllProjects <- as.data.frame(AllProjects)
  AllProjects$Project <- c(rep(1, nrep), rep(2,nrep), rep(3, nrep), rep(4, nrep))
  Best_project <- AllProjects[AllProjects[, 1]==min(AllProjects[, 1]), 2]
  Best_runs <- AllProjects[AllProjects$Project==Best_project, ]
  Best_runs$Nrun <- 1:nrow(Best_runs)
  Best_run <- Best_runs[Best_runs[, 1]==min(Best_runs[, 1]), 3]
  print(paste0("Best run is: ", "project = ", Best_project, ", run = ", Best_run)) 
}

barplotK <- function(Qfile, Pop, Run_B){
  Q.matrix_1 <-LEA::Q(Qfile, K = Pop, run = Run_B)
  Q.matrix <- as.qmatrix(Q.matrix_1)
  barplot(Q.matrix, xlab = "Sampled individuals",
          ylab = "Ancestry coefficients",
          main = "", cex.axis = 1.5, cex.lab = 1)
}