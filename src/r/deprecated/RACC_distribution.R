rm(list=ls())
dev.off()

  ###################################
### ---- FULL-SYNTH DATASET --- ###
###################################

setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Use Cases/Synthetic Dataset/Results")

original.distro.per.specie <- read.csv("../gastro_high.pa", header=FALSE, sep = '\t')
reads.distro.per.specie <- read.csv("../MEGAN/BLAST_format_gastro_reads-taxon_distro.txt", header=FALSE, sep = ',')
contigs.distro.per.specie <- read.csv("../MEGAN/BLAST_contigs-taxon_distro.txt", header=FALSE, sep = ',')

# Make all species be in all sets
which.reads <- which(!(reads.distro.per.specie$V1 %in% contigs.distro.per.specie$V1)==TRUE)
which.contigs <- which(!(contigs.distro.per.specie$V1 %in% reads.distro.per.specie$V1)==TRUE)
which.reads.info <- reads.distro.per.specie[which.reads,]
which.contigs.info <- contigs.distro.per.specie[which.contigs,]
which.reads.null <- cbind(as.vector(which.contigs.info$V1), rep(0, length(which.contigs.info$V1)))
which.contigs.null <- cbind(as.vector(which.reads.info$V1), rep(0, length(which.reads.info$V1)))

reads.distro.per.specie <- rbind(reads.distro.per.specie[-which.reads,],which.reads.info,which.reads.null)
contigs.distro.per.specie <- rbind(contigs.distro.per.specie[-which.contigs,], which.contigs.null, which.contigs.info)

# CHECK ITS THE SAME
which(reads.distro.per.specie$V1 != contigs.distro.per.specie$V1)

orig.desc <- c()
for(specie in reads.distro.per.specie[,1]){
  orig.desc <- c(orig.desc, sum(original.distro.per.specie[grep(specie, original.distro.per.specie[,1]),2])) 
}

# Prepare Values
orig.desc <- as.numeric(orig.desc)*100
reads.desc <- as.numeric(reads.distro.per.specie[,2])
contigs.desc <- as.numeric(contigs.distro.per.specie[,2])
names(orig.desc) <- names(reads.desc) <- names(contigs.desc) <- reads.distro.per.specie[,1]

#--
distro.info <- rbind(reads.desc, orig.desc, contigs.desc)
distro.info.2 <- rbind((orig.desc-reads.desc), (orig.desc-contigs.desc))

fsd.distro.info <- distro.info
fsd.distro.info.2 <- distro.info.2
distro.info <- fsd.distro.info
distro.info.2 <- fsd.distro.info.2

reads.desc <- distro.info[1,]
orig.desc <- distro.info[2,]
contigs.desc <- distro.info[3,]

# Distribution lines
sort.data <- data.frame(orig=orig.desc,reads=reads.desc,contig=contigs.desc)
sort.data <- sort.data[order(sort.data$orig, decreasing = TRUE),]
dev.off()
#ylim=c(0.001,10)
t <- sapply(1:length(rownames(sort.data)), function(x) "")
par(mar=c(6,6,5,1.5))
plot(sort.data$orig,type="l",col="green",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Species",xaxt="n", main = "Percentage of Sequences per Specie (FSD)", log="y", cex.main=2, cex.lab = 2.5, cex.axis=2)

axis(1, at=c(1:nrow(sort.data)), labels = t)
lines(sort.data$reads,type="l",col="red",lty=1,lwd=2)
lines(sort.data$contig,type="l",col="blue",lty=1,lwd=2)
lines(sort.data$orig, type="l", col="green", lty=1, lwd=2)
axis(1, at=c(1:nrow(sort.data)), labels = t)
#axis(1,at=c(1:nrow(sort.data)),labels=rownames(sort.data), cex.axis = 0.4)
grid()
legend("topright",legend=c("Original","Reads","Contigs"),col=c("green","red","blue"),bg="white",lwd=2, cex = 1)

legend("topright",legend=c("Original"),col=c("green"),bg="white",lwd=2, cex = 1.6)
# Distribucion Normal
library(lattice)

max.values <- c(reads.sd,contigs.sd)
equal.count(reads.discrete)


discrete.values <- seq(min(max.values)-0.1,max(max.values)+0.1, length.out = 100)
discrete.diff <- discrete.values[2]-discrete.values[1]

reads.discrete <- shingle(reads.sd,cbind(discrete.values ,discrete.values+discrete.diff))#table(cut(reads.sd, discrete.values,include.lowest = TRUE))
contigs.discrete <- shingle(contigs.sd,cbind(discrete.values ,discrete.values+discrete.diff))#table(cut(contigs.sd, discrete.values,include.lowest = TRUE))

reads.discrete

plot(x = y = reads.discrete[20:35], type="l", col="red",lwd=2)
lines(y = contigs.discrete[20:35], col="blue", type="l",lwd=2)
lines(reads.density, col="red", type="l",lwd=2)
lines(c(0,0),c(0,250),col="black", lty=2)
#
dev.off()
reads.sd <- reads.desc-orig.desc
contigs.sd <- contigs.desc-orig.desc

half.length <- 5
number.points <- 1000

reads.norm <- dnorm(seq(-half.length,half.length,length.out = number.points),mean=mean(reads.sd), sd=sd(reads.sd))
contigs.norm <- dnorm(seq(-half.length,half.length,length.out = number.points),mean=mean(contigs.sd), sd=sd(contigs.sd))
sd.table <- data.frame(reads=reads.norm, contigs=contigs.norm)

plot(x=seq(-half.length,half.length,length.out = number.points),y=sd.table$reads*100,type="l",col="red",lty=1,ylab="Probability of Error (%)",lwd=2,
     xlab="Error Assigning a Specie", main = "Distribution of errors with Reads & Contigs", xaxt = "n")
lines(x=seq(-half.length,half.length,length.out = number.points),y = sd.table$contigs*100,type="l",col="blue",lty=1,lwd=2)
axis(side=1, at=seq(-half.length,half.length,by=0.5), labels=seq(-half.length,half.length,by=0.5), cex=0.3)
grid()
legend("topleft",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

# Diferencia porcentual

read.p.diff <- (reads.desc-orig.desc)*100/orig.desc
contig.p.diff <- (contigs.desc-orig.desc)*100/orig.desc
p.diff <- data.frame(reads=read.p.diff, contigs=contig.p.diff,abs.diff=(orig.desc-orig.desc))
p.diff <- p.diff[order(p.diff$abs.diff, decreasing = TRUE),]
p.diff <- p.diff[-which(p.diff[,2]=="Inf"),]

#ylim=c(-100,100)
plot(p.diff$reads,type="l",col="red",lty=1,ylab="Percentual Difference of Sequences Assigned per Specie",lwd=2,
     xlab="Specie",xaxt="n", main = "Percentual Difference of Sequences Assigned per Specie")
lines(p.diff$contigs,type="l",col="blue",lty=1,lwd=2)
lines(orig.desc-orig.desc,type="l",col="green",lty=1,lwd=2)
legend("topleft",legend=c("Reads","Contigs","Expected"),col=c("red","blue","green"),bg="white",lwd=2)
par(las=2)
axis(1,at=c(1:nrow(p.diff)),labels=rownames(p.diff), cex.axis = 0.4)
grid()

dev.off()

sum(abs(p.diff$reads))
sum(abs(p.diff$contigs))

# Reads Vs Contigs Difference of Percentual Difference

orig.sort <- order(orig.desc, decreasing = TRUE)
read.sort.p.diff <- read.p.diff[orig.sort]
contig.sort.p.diff <- contig.p.diff[orig.sort]
sort.p.diff <- data.frame(reads=read.sort.p.diff, contigs=contig.sort.p.diff)
read.p.diff <- read.sort.p.diff
contig.p.diff <- contig.sort.p.diff

p.diff.diff <- sapply(1:length(sort.p.diff$reads), function(x) 
  abs(sort.p.diff$reads[x]) - abs(sort.p.diff$contigs[x]))
p.diff.diff <- sort(p.diff.diff, decreasing = TRUE)
p.diff.diff <- p.diff.diff[orig.sort]
plot(p.diff.diff,type="l",col="red", ylim=c(-2*max(p.diff.diff),max(p.diff.diff)), lty=1,ylab="Difference of the Percentual Difference between Reads and Contigs",lwd=2,
     xlab="Specie",xaxt="n", main = "Percentual Difference of Sequences Assigned per Specie")
lines(orig.desc-orig.desc,type="l",col="green",lty=1,lwd=2)

dev.off()

# Sum of absolute differences
distro.sum <- c(sum(abs(reads.desc-orig.desc)), sum(abs(contigs.desc-orig.desc)))
barplot(distro.sum, col = c('red','blue'), ylim = c(min(distro.sum)-1, max(distro.sum)+1), xpd = FALSE, main = "Sum of the Error Percentage of Sequences assigned")
axis(1,at=c(1,2), labels = c("Reads","Contigs"))
#legend("bottomright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

dev.off()

# RMSE
distro.rmse <- c(sum((reads.desc-orig.desc)^2)/length(orig.desc), sum((contigs.desc-orig.desc)^2)/length(orig.desc))
barplot(distro.rmse, col = c('red','blue'), ylim = c(0.1, 5), xpd = FALSE, main = "RMSE of the Sequences Assigned", log="y")
axis(1,at=c(1,2), labels = c("Reads","Contigs"))
#legend("bottomright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

dev.off()
# Relacion diferencia porcentual / distribucion original

#----> MORE READS ARE BETTER THAN CONTIGS

###################################
### ---- SEMI-SYNTH DATASET --- ###
###################################
# Name: Reads, Contigs
# Total: 521334, 52953
# Assigned: 517033, 52680

setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Use Cases/SemiSynth Dataset/results/")

ssd.original.distro.per.specie <- read.csv("../soil_high_mod.pa", header=FALSE, sep = '\t')
ssd.reads.distro.per.specie <- read.csv("../MEGAN/BLAST_ssd_reads-taxon_distro.txt", header=FALSE, sep = ',')
ssd.contigs.distro.per.specie <- read.csv("../MEGAN/BLAST_ssd_contigs-taxon_distro.txt", header=FALSE, sep = ',')

which(ssd.reads.distro.per.specie$V1 != ssd.contigs.distro.per.specie$V1)

orig.desc <- c()
for(specie in ssd.reads.distro.per.specie[,1]){
  orig.desc <- c(orig.desc, ssd.original.distro.per.specie[grep(specie, ssd.original.distro.per.specie[,2]),5]) 
}

reads.distro.per.specie <-  ssd.reads.distro.per.specie
contigs.distro.per.specie <- ssd.contigs.distro.per.specie

# Prepare Values
orig.desc <- as.numeric(orig.desc)
reads.desc <- as.numeric(reads.distro.per.specie[,2])
contigs.desc <- as.numeric(contigs.distro.per.specie[,2])
names(orig.desc) <- names(reads.desc) <- names(contigs.desc) <- reads.distro.per.specie[,1]

distro.info <- rbind(reads.desc, orig.desc, contigs.desc)
distro.info.2 <- rbind(abs(reads.desc-orig.desc), abs(contigs.desc-orig.desc))

# -- Sorted
sort.data <- data.frame(orig=orig.desc,reads=reads.desc,contig=contigs.desc)
sort.data <- sort.data[order(sort.data$orig, decreasing = TRUE),]

#ylim=c(0.001,10)
t <- sapply(1:length(rownames(sort.data)), function(x) "")
par(mar=c(6,6,5,1.5))
plot(sort.data$orig,type="l",col="green",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Species",xaxt="n", main = "Percentage of Sequences per Specie (SSD)", log="y", cex.main=2, cex.lab = 2.5, cex.axis=2)

axis(1, at=c(1:nrow(sort.data)), labels = t)

lines(sort.data$reads,type="l",col="red",lty=1,lwd=2)
lines(sort.data$contig,type="l",col="blue",lty=1,lwd=2)
lines(sort.data$orig, type="l", col="green", lty=1, lwd=2)
par(las=2)
#axis(1,at=c(1:nrow(sort.data)),labels=rownames(sort.data), cex.axis = 0.4)
grid()
legend("topright",legend=c("Original","Reads","Contigs"),col=c("green","red","blue"),bg="white",lwd=2, cex = 1)

legend("topright",legend=c("Original"),col=c("green"),bg="white",lwd=2, cex = 1.6)






# Diferencia porcentual

read.p.diff <- (reads.desc-orig.desc)*100/orig.desc
contig.p.diff <- (contigs.desc-orig.desc)*100/orig.desc
p.diff <- data.frame(reads=read.p.diff, contigs=contig.p.diff,abs.diff=read.p.diff-contig.p.diff)
p.diff <- p.diff[order(p.diff$abs.diff, decreasing = TRUE),]

plot(p.diff$reads,type="l",ylim=c(min(p.diff),max(p.diff)),col="red",lty=1,ylab="Percentual Differences",lwd=2,
     xlab="Species",xaxt="n", main = "Percentual Difference of Sequences Assigned per Specie")
lines(p.diff$contigs,type="l",col="blue",lty=1,lwd=2)
lines(orig.desc-orig.desc,type="l",col="green",lty=1,lwd=2)
legend("topleft",legend=c("Reads","Contigs","Expected"),col=c("red","blue","green"),bg="white",lwd=2)
par(las=2)
axis(1,at=c(1:ncol(distro.info)),labels=colnames(distro.info), cex.axis = 0.6)
grid()

# Standard Deviation
reads.sd <- reads.desc-orig.desc
contigs.sd <- contigs.desc-orig.desc
half.length <- 5
number.points <- 1000

reads.norm <- dnorm(seq(-half.length,half.length,length.out = number.points),mean=mean(reads.sd), sd=sd(reads.sd))
contigs.norm <- dnorm(seq(-half.length,half.length,length.out = number.points),mean=mean(contigs.sd), sd=sd(contigs.sd))
sd.table <- data.frame(reads=reads.norm, contigs=contigs.norm)

plot(x=seq(-half.length,half.length,length.out = number.points),y=sd.table$reads*100,type="l",col="red",lty=1,ylab="Probability of Error (%)",lwd=2,
     xlab="Error Assigning a Specie", main = "Distribution of errors with Reads & Contigs", xaxt = "n")
lines(x=seq(-half.length,half.length,length.out = number.points),y = sd.table$contigs*100,type="l",col="blue",lty=1,lwd=2)
axis(side=1, at=seq(-half.length,half.length,by=0.5), labels=seq(-half.length,half.length,by=0.5), cex=0.3)
grid()
legend("topleft",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

dev.off()

# Sum of absolute differences

distro.sum <- c(sum(abs(reads.desc-orig.desc)), sum(abs(contigs.desc-orig.desc)))
barplot(distro.sum, col = c('red','blue'), ylim = c(min(distro.sum)-1, max(distro.sum)+1), xpd = FALSE, main = "Sum of the Error Percentage of Sequences assigned")
axis(1,at=c(1,2), labels = c("Reads","Contigs"))
#legend("bottomright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

# RMSE
distro.rmse <- c(sum((reads.desc-orig.desc)^2)/length(orig.desc), sum((contigs.desc-orig.desc)^2)/length(orig.desc))
barplot(distro.rmse, col = c('red','blue'), ylim = c(0.1, 5), xpd = FALSE, main = "RMSE of the Sequences Assigned", log="y")
axis(1,at=c(1,2), labels = c("Reads","Contigs"))
#legend("bottomright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)



#############################################################
# ---------------------
# Original Distribution
#,ylim=c(min(distro.info[2,]),max(distro.info[2,]))
plot(distro.info[2,],type="l",col="green",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Species",xaxt="n", main = "Species Distribution Synthetic Dataset (%)", log="y", cex=1.5)
lines(distro.info[1,],type="l",col="red",lty=1,lwd=2)
lines(distro.info[3,],type="l",col="blue",lty=1,lwd=2)
lines(distro.info[2,], type="l", col="green", lty=1, lwd=2)
par(las=2)
axis(1,at=c(1:ncol(distro.info)),labels=colnames(distro.info), cex.axis = 0.4)
grid()
legend("topright",legend=c("Original","Reads","Contigs"),col=c("green","red","blue"),bg="white",lwd=2, cex = 0.7)


legend("topright",legend=c("Original"),col=c("dark green"),bg="white",lwd=2, cex = 0.7)

#############################################################
# ---------------------
# Original Distribution
plot(distro.info[2,],type="o",ylim=c(min(distro.info),max(distro.info)),col="dark green",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Species",xaxt="n", main = "Species Distribution Semi-Synthetic Dataset (%)", log="y")
par(las=2)
axis(1,at=c(1:ncol(distro.info)),labels=colnames(distro.info), cex.axis = 0.4)
grid()

legend("topright",legend=c("Original"),col=c("dark green"),bg="white",lwd=2, cex = 0.7)

# Distribution lines
plot(distro.info[2,],type="o",ylim=c(min(distro.info),max(distro.info)),col="dark green",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Specie",xaxt="n", main = "Percentage of Sequences per Specie", log = "y")
legend("topright",legend=c("Original","Reads","Contigs"),col=c("dark green","red","blue"),bg="white",lwd=2)
lines(distro.info[1,],type="o",col="red",lty=1,lwd=2)
lines(distro.info[3,],type="o",col="blue",lty=1,lwd=2)
par(las=2)
axis(1,at=c(1:ncol(distro.info)),labels=colnames(distro.info), cex.axis = 0.4)
grid()

#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################

#---- DEPRECATED

# Distribution lines
sort.data <- data.frame(orig=orig.desc,reads=reads.desc,contig=contigs.desc)
sort.data <- sort.data[order(sort.data$orig, decreasing = TRUE),]

plot(sort.data$orig,type="l",col="green",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Specie",xaxt="n", main = "Percentage of Sequences per Specie", log="y")

lines(sort.data$reads,type="l",col="red",lty=1,lwd=2)
lines(sort.data$contig,type="l",col="blue",lty=1,lwd=2)
par(las=2)
axis(1,at=c(1:ncol(distro.info)),labels=colnames(distro.info), cex.axis = 0.4)
grid()
legend("topleft",legend=c("Original","Reads","Contigs"),col=c("dark green","red","blue"),bg="white",lwd=2, cex = 0.7)

# idk
plot(orig.desc,type="l",ylim=c(0.001,10),col="green",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Specie",xaxt="n", main = "Percentage of Sequences per Specie", log="y")

lines(distro.info[1,],type="l",col="red",lty=1,lwd=2)
options("scipen"=999, "digits"=4)
plot(sort(distro.info.2[1,]-distro.info.2[2,]),type="l",ylim=c(0.00001,1.5),col="red",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Specie",xaxt="n", main = "Percentage of Sequences per Specie")
mierda <- c()
for(specie in names(sort(distro.info.2[1,]))){
  mierda <- c(mierda, distro.info.2[2,specie]) 
}
lines(mierda,type="l",col="blue",lty=1,lwd=2)

lines(sort(distro.info.2[2,]),type="l",col="blue",lty=1,lwd=2)
par(las=2)
axis(1,at=c(1:ncol(distro.info)),labels=colnames(distro.info), cex.axis = 0.4)
grid()
legend("topleft",legend=c("Original","Reads","Contigs"),col=c("dark green","red","blue"),bg="white",lwd=2, cex = 0.7)


# Absolute Difference Best Performance

better.reads <- distro.info.2[1,which(distro.info.2[1,] > distro.info.2[2,])]
names(better.reads) <- colnames(distro.info.2)[which(distro.info.2[1,] > distro.info.2[2,])]
better.contigs <- distro.info.2[2,which(distro.info.2[2,] > distro.info.2[1,])]
names(better.contigs) <- colnames(distro.info.2)[which(distro.info.2[2,] > distro.info.2[1,])]

equal.reads.contigs <- distro.info.2[which(distro.info.2[1,] == distro.info.2[2,])]
equal.reads.contigs

#PLOTS - Better Reads
plot(orig.desc[order(p.diff$abs.diff, decreasing = TRUE)],type="l",ylim=c(0.001,10),col="green",lty=1,ylab="Percentage Difference of Sequences",lwd=2,
     xlab="Specie",xaxt="n", main = "Difference Percentage of Sequences where Reads are better than Contigs", log="y")

lines(distro.info.2[1,names(better.reads)],type="l",col="red",lty=1,lwd=2)
lines(distro.info.2[2,names(better.reads)],type="l",col="blue",lty=1,lwd=2)
par(las=2)
axis(1,at=c(1:length(names(better.reads))),labels=names(better.reads), cex.axis = 0.4)
grid()
legend("topleft",legend=c("Original","Reads","Contigs"),col=c("green","red","blue"),bg="white",lwd=2, cex = 0.7)
dev.off()

#PLOTS - Better Contigs
plot(orig.desc[names(better.contigs)],type="l",ylim=c(0.001,10),col="green",lty=1,ylab="Percentage Difference of Sequences",lwd=2,
     xlab="Specie",xaxt="n", main = "Difference Percentage of Sequences where Reads are better than Contigs", log="y")

lines(distro.info.2[2,names(better.contigs)],type="l",col="blue",lty=1,lwd=2)
lines(distro.info.2[1,names(better.contigs)],type="l",col="red",lty=1,lwd=2)
par(las=2)
axis(1,at=c(1:length(names(better.contigs))),labels=names(better.contigs), cex.axis = 0.4)
grid()
legend("topleft",legend=c("Original","Reads","Contigs"),col=c("green","red","blue"),bg="white",lwd=2, cex = 0.7)