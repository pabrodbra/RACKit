rm(list=ls())

###################################
### ---- FULL-SYNTH DATASET --- ###
###################################

setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/Results")

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

# All bars - NO!!! DONT!!!
distro.info <- rbind(reads.desc, orig.desc, contigs.desc)
barplot(abs(distro.info), beside = TRUE, col = c('red','black','blue'), ylim = c(0,5))
legend("topleft",legend=c("Original","Reads","Contigs"),col=c("black","red","blue"),bg="white",lwd=2)

# Absolute Difference Best Performance
distro.info.2 <- rbind(abs(reads.desc-orig.desc), abs(contigs.desc-orig.desc))

better.reads <- distro.info.2[1,which(distro.info.2[1,] > distro.info.2[2,])]
names(better.reads) <- colnames(distro.info.2)[which(distro.info.2[1,] > distro.info.2[2,])]
better.contigs <- distro.info.2[2,which(distro.info.2[2,] > distro.info.2[1,])]
names(better.contigs) <- colnames(distro.info.2)[which(distro.info.2[2,] > distro.info.2[1,])]

equal.reads.contigs <- distro.info.2[which(distro.info.2[1,] == distro.info.2[2,])]
equal.reads.contigs

#PLOTS - Better Reads
plot(orig.desc[names(better.reads)],type="l",ylim=c(0.001,10),col="green",lty=1,ylab="Percentage Difference of Sequences",lwd=2,
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

dev.off()

#----> MORE READS ARE BETTER THAN CONTIGS

# Absolute Difference

barplot(distro.info.2, beside = TRUE, col = c('red','blue'), main = "Absolute difference against original distribution", log="y")
legend("topright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

# Chi Squared Difference
distro.info.2 <- rbind(abs(reads.desc-orig.desc), abs(contigs.desc-orig.desc))
barplot(distro.info.2, beside = TRUE, col = c('red','blue'), main = "Chi Squared difference against original distribution")
legend("topright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

# Sum of absolute differences
distro.sum <- c(sum(abs(reads.desc-orig.desc)), sum(abs(contigs.desc-orig.desc)))
barplot(distro.sum, col = c('red','blue'))
legend("topright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

# Distribution lines # BEST? # 
plot(orig.desc,type="l",ylim=c(0.001,10),col="green",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Specie",xaxt="n", main = "Percentage of Sequences per Specie", log="y")

lines(distro.info[1,],type="l",col="red",lty=1,lwd=2)
plot(distro.info[1,],type="l",ylim=c(0.001,10),col="red",lty=1,ylab="Percentage of Sequences",lwd=2,
                                                          xlab="Specie",xaxt="n", main = "Percentage of Sequences per Specie", log="y")
lines(distro.info[3,],type="l",col="blue",lty=1,lwd=2)
par(las=2)
axis(1,at=c(1:ncol(distro.info)),labels=colnames(distro.info), cex.axis = 0.4)
grid()
legend("topleft",legend=c("Original","Reads","Contigs"),col=c("dark green","red","blue"),bg="white",lwd=2, cex = 0.7)

dev.off()


###################################
### ---- SEMI-SYNTH DATASET --- ###
###################################

setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/SemiSynth/results/")

ssd.original.distro.per.specie <- read.csv("../soil_high_mod.pa", header=FALSE, sep = '\t')
ssd.reads.distro.per.specie <- read.csv("../MEGAN/BLAST_ssd_reads-taxon_distro.txt", header=FALSE, sep = ',')
ssd.contigs.distro.per.specie <- read.csv("../MEGAN/BLAST_ssd_contigs-taxon_distro.txt", header=FALSE, sep = ',')

which(ssd.reads.distro.per.specie$V1 != ssd.contigs.distro.per.specie$V1)

orig.desc <- c()
for(specie in ssd.reads.distro.per.specie[,1]){
  orig.desc <- c(orig.desc, ssd.original.distro.per.specie[grep(specie, ssd.original.distro.per.specie[,2]),5]) 
}

# Prepare Values
orig.desc <- as.numeric(orig.desc)
reads.desc <- as.numeric(reads.distro.per.specie[,2])
contigs.desc <- as.numeric(contigs.distro.per.specie[,2])
names(orig.desc) <- names(reads.desc) <- names(contigs.desc) <- reads.distro.per.specie[,1]

# All bars
distro.info <- rbind(reads.desc, orig.desc, contigs.desc)
barplot(distro.info, beside = TRUE, col = c('red','black','blue'))
legend("topleft",legend=c("Original","Reads","Contigs"),col=c("black","red","blue"),bg="white",lwd=2)

# Absolute Difference
distro.info.2 <- rbind(abs(reads.desc-orig.desc), abs(contigs.desc-orig.desc))
barplot(distro.info.2, beside = TRUE, col = c('red','blue'), main = "Absolute difference against original distribution")
legend("topright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

# Chi Squared Difference
distro.info.2 <- rbind(abs(reads.desc-orig.desc), abs(contigs.desc-orig.desc))
barplot(distro.info.2, beside = TRUE, col = c('red','blue'), main = "Chi Squared difference against original distribution")
legend("topright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

# Sum of absolute differences
distro.sum <- c(sum(abs(reads.desc-orig.desc)), sum(abs(contigs.desc-orig.desc)))
barplot(distro.sum, col = c('red','blue'))
legend("topright",legend=c("Reads","Contigs"),col=c("red","blue"),bg="white",lwd=2)

# Distribution lines # BEST? # 
plot(distro.info[2,],type="o",ylim=c(min(distro.info),max(distro.info)),col="dark green",lty=1,ylab="Percentage of Sequences",lwd=2,
     xlab="Specie",xaxt="n", main = "Percentage of Sequences per Specie")
legend("topleft",legend=c("Original","Reads","Contigs"),col=c("dark green","red","blue"),bg="white",lwd=2)
lines(distro.info[1,],type="o",col="red",lty=1,lwd=2)
lines(distro.info[3,],type="o",col="blue",lty=1,lwd=2)
par(las=2)
axis(1,at=c(1:ncol(distro.info)),labels=colnames(distro.info), cex.axis = 0.4)
grid()

dev.off()
