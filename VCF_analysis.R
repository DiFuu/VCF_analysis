
##############
#EXTRACT FILE#
##############

#To unzip the vcf file
library(R.utils)
gunzip("Variants_file.vcf.gz", remove=FALSE)

#To use the vcf file
library(vcfR)
vcf <- read.vcfR("Variants_file.vcf", verbose = FALSE )


#############
#GENOME SIZE#
#############

#Read all the information of the contigs
contigs <- queryMETA(vcf, element = 'contig')

#First we are going to calculate the size of the reference genome
contigs_length <- sapply(contigs, function(x){x[2]}) #Extract the length part
c_length <- as.numeric(gsub("length=", "", contigs_length)) #Extract the number

sum(c_length)/1000000 #calculate the Mb -> 687.4581 Mb

#We are going to look for the reference genome
reference <- queryMETA(vcf, element = 'reference')
reference[[1]] #Xgladius_genome.fasta -> Xiphias gladius


################
#MISSING VALUES#
################

#number of missing values in first sample
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
sum(is.na(dp[,1]))


#missing values in all samples
myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) }) #1 files, 2 columns, sum of all missing values
myMiss <- myMiss/nrow(vcf) #percent of missing values per sample

library(RColorBrewer)
palette(brewer.pal(n=12, name = 'Set3'))

par(mar = c(12,4,4,2))
barplot(myMiss, las = 2, col = 1:12)
title(ylab = "Missing Values (%)")

#To export a table with the information of missing values per sample
mv_table <- as.data.frame(myMiss)
mv_table <- tibble::rownames_to_column(mv_table, "Sample")
colnames(mv_table) <- c("Sample", "Missing Values (%)")
write.table(mv_table,file="Percentage of Missing Values.txt",row.names=F,quote=F,sep="\t")
