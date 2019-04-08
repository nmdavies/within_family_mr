#Neil Davies 06/12/18
#This select SNPs to use as instruments using MR-Base

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

phenotypes <-c(1001)
instruments<-extract_instruments(outcomes=phenotypes, p1 = 5e-08)

colnames(instruments)
write.table(instruments,paste($path1,"/snplists/okbay_snps_effect.txt",sep(""),quote = FALSE)
write.table(instruments$SNP,paste($path1,"/snplists/okbay_snps.txt",sep(""),quote = FALSE, row.names = FALSE,col.names = FALSE)
