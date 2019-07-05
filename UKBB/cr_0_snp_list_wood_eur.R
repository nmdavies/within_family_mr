#Neil Davies 24/04/19
#This select SNPs to use as instruments for height using MR-Base

#Load MR base packages
library(devtools)
install_github("MRCIEU/TwoSampleMR")
library(foreign)
#Load TwoSampleMR package
library(TwoSampleMR)

#Install MR base instruments
devtools::install_github("MRCIEU/MRInstruments")

library(MRInstruments)

#Extract IDs, coefficients and effect alleles for the covariate allele scores
#Extract 
phenotypes <-c(89)
instruments<-extract_instruments(outcomes=phenotypes, p1 = 5e-08,r2 = 0.001)

colnames(instruments)
write.table(instruments,paste(path1,"/snplists/wood_eur_effect.txt",sep=""),quote = FALSE)
write.table(instruments$SNP,paste(path1,"/snplists/wood_eur_snps.txt",sep=""),quote = FALSE, row.names = FALSE,col.names = FALSE)
