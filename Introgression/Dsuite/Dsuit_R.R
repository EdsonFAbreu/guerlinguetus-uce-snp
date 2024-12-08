D_BBAA <- read.table("sets_pops_BBAA.txt",as.is=T,header=T)

plot(D_BBAA$Dstatistic, ylab="D",xlab="trio number") # The D statistic
plot(D_BBAA$p.value, ylab="p-value",xlab="trio number") # The p-value
#plot(D_BBAA$f4.ratio, ylab="f4-ratio",xlab="trio number", ylim=c(0,0.2))
plot(D_BBAA$f4.ratio, ylab="f4-ratio",xlab="trio number")

plot(p.adjust(D_BBAA$p.value,method="BH"), ylab="corected p value",xlab="trio number",ylim=c(0,0.05))


