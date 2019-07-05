#Neil Davies 06/12/18
#This select SNPs to use as instruments using MR-Base

#Load MR base packages
library(devtools)
#install_github("MRCIEU/TwoSampleMR")
library(foreign)
#Load TwoSampleMR package
library(TwoSampleMR)

#Install MR base instruments
devtools::install_github("MRCIEU/MRInstruments")

library(MRInstruments)

#Load Stata read package
library(readstata13)

ao <- available_outcomes()
setwd(path1)

#Import available UK Biobank SNPs
#Create a file list using cr_0_ukbb_snps.do

snp_list <-read.dta13(paste(path1,"/workingdata/snps.dta",sep=""))

#Extract IDs, coefficients and effect alleles for the covariate allele scores
#Education variants from the 2016 Okbay paper
phenotypes <-c(1001)
instruments<-extract_instruments(outcomes=phenotypes, p1 = 5e-08,r2 = 0.001)

write.table(instruments,paste(path1,"/snplists/okbay_snps_effect.txt",sep=""),quote = FALSE)
save.dta13(instruments,paste(path1,"/snplists/okbay_snps_effect.dta",sep=""))
write.table(instruments$SNP,paste(path1,"/snplists/okbay_snps.txt",sep=""),quote = FALSE, row.names = FALSE,col.names = FALSE)

#BMI variants from the Locke 2015 height paper (Europeans)
phenotypes <-c(835)
instruments<-extract_instruments(outcomes=phenotypes, p1 = 5e-08,r2 = 0.001)

write.table(instruments,paste(path1,"/snplists/locke_eur_effect.txt",sep=""),quote = FALSE)
save.dta13(instruments,paste(path1,"/snplists/locke_eur_effect.dta",sep=""))
write.table(instruments$SNP,paste(path1,"/snplists/locke_eur_snps.txt",sep=""),quote = FALSE, row.names = FALSE,col.names = FALSE)

#386 Height variants from the Wood 2015 height paper (Europeans) 
#1 variant not in UK Biobank rs9825951, so will restrict to UK Biobank SNPs and then clump
phenotypes <-c(89)
instruments<-extract_instruments(outcomes=phenotypes, p1 = 5e-08,clump = FALSE)

#26593 SNPs available, subset to those available in UKBB
instruments<-instruments[instruments$SNP %in% snp_list$rsid,]

#This drops the sample to 26319 SNPs, next clump on r2=0.001 and default parameters
instruments<-clump_data(instruments)

#This results in 386 SNPs, all of which are available in UKBB.
write.table(instruments,paste(path1,"/snplists/wood_eur_effect.txt",sep=""),quote = FALSE)
save.dta13(instruments,paste(path1,"/snplists/wood_eur_effect.dta",sep=""))
write.table(instruments$SNP,paste(path1,"/snplists/wood_eur_snps.txt",sep=""),quote = FALSE, row.names = FALSE,col.names = FALSE)

#Create list of SNPs to extract from UKBB
instruments1<-as.data.frame(instruments$SNP)
colnames(instruments1)<-"V1"
instruments2<-read.table(paste(path1,"/snplists/locke_eur_snps.txt",sep=""))
instruments3<-read.table(paste(path1,"/snplists/okbay_snps.txt",sep=""))
instruments_all<-rbind(instruments1,instruments2,instruments3)

#Output list for extraction
write.table(instruments_all,paste(path1,"/snplists/all_snps.txt",sep=""),quote = FALSE, row.names = FALSE,col.names = FALSE)

