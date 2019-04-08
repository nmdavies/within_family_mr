# Neil Davies 06/12/18
# This runs option C
# within family - PLM method with robust SE for each SNP individually for each phenotype

# For each SNP-phenotype relationship create a file with:

#- trait
#- SNP
#- beta
#- se
#- pvalue 
#- effect allele
#- non-effect allele
#- effect allele frequencies
#- sample size

# 1. Load libraries

#install.packages('plm')
#install.packages("ivpack")
#install.packages("data.table") 

library(data.table)
library("plm")
library(ivpack)
library(data.table)  
library('dplyr')

# 2. Functions
substrright <- function(x, n, m){
  substr(x, nchar(x)-n+1, nchar(x)-n+m)
}

# 3. Delete all objects
rm(list = ls())

# 4. Load data

#Load genotypic and phenotypic data
data<-fread(paste($path1,"/workingdata/analysisdata.csv",sep(""))

#Create the data set for analysis
#Extract siblings
data_sibs<-subset(data,n_sibs!=1 & n_sibs!=336538)

#Display objects
names(data)

#Select siblings
data_sibs<-subset(data,famid!="NA")
dim(data_sibs)

#5. Regressions 
# Run the within family fixed effects with 'plm' and robust standard errors for each SNP seperately

snp_plm <- function(outcome, exposure){
  SNP<-substr(exposure,1,nchar(exposure)-4)
  effect_allele<-substrright(exposure,3,1)
  other_allele<-substrright(exposure,1,1)
  
  #Run plm model with family fixed effects + robust standard errors
  estimates<-plm(get(outcome)~get(exposure)+cov_age+cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20 , data = data_sibs, index = "famid", model = "within", inst.method = "bvk")
  estimates_robust<-coeftest(estimates,vcov=vcovHC(estimates,type="HC0",cluster="group"))
  
  beta<-estimates_robust[1,1]
  se<-estimates_robust[1,2]
  pvalue<-estimates_robust[1,3]
  sample_size<-nrow(model.frame(estimates))
  EAF<-mean(model.frame(estimates)[,2])/2
  
  return(c(outcome,SNP,beta,se,pvalue,effect_allele,other_allele,EAF,sample_size))
}

phenotype_snp_est<-as.data.frame(t(c(0,0,0,0,0,0,0,0,0)))
colnames(phenotype_snp_est)<-c("trait","SNP","beta","se","pvalue","effect_allele","other_allele","EAF","sample_size")

for (outcome in c("eduyears2","out_highbloodpressure","out_diabetes","out_height","out_bmi")){
  snps<-grep(names(data_sibs),pattern = "^rs",value = TRUE)
  for  (exposure in snps){
    message(exposure)
    message(outcome)
    phenotype_snp_est<-rbind(phenotype_snp_est,snp_plm(outcome,exposure))
  }
}

write.table(phenotype_snp_est,paste($path1,"/results/individual_snps.txt",sep(""),quote=FALSE,row.names = FALSE)
