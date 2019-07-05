# Neil Davies 06/12/18
# This runs option C
# within family - PLM method with robust SE for each SNP individually for each phenotype
# source('reg_2_option_c_individual_snps-v5.R')

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
#install.packages("lfe") 

library('data.table')
library('plm')
library('ivpack')
library('data.table')  
library('dplyr')
library('lfe')

# 2. Functions
substrright <- function(x, n, m){
  substr(x, nchar(x)-n+1, nchar(x)-n+m)
}

# 3. Delete all objects
rm(list = ls())

# 4. Load data

#Load genotypic and phenotypic data
load("data/data_sibs_v4.RData")

#Create the data set for analysis
#Extract siblings
dim(data_sibs)

#5. Regressions 
# Run the within family fixed effects with 'plm' and robust standard errors for each SNP seperately

snp_plm <- function(outcome, exposure){

  #Run plm model with family fixed effects + robust standard errors
  estimates<-plm(get(outcome)~get(exposure)+cov_age+cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+batch , data = data_sibs, index = "famid", model = "within", inst.method = "bvk")
  estimates_robust<-coeftest(estimates,vcov=vcovHC(estimates,type="HC0",cluster="group"))
  
  beta<-estimates_robust[1,1]
  se<-estimates_robust[1,2]
  pvalue<-estimates_robust[1,4]
  sample_size<-nrow(model.frame(estimates))
  EAF<-mean(model.frame(estimates)[,2])/2
  
  return(c(outcome,exposure,beta,se,pvalue,EAF,sample_size))
}

phenotype_snp_est<-as.data.frame(t(c(0,0,0,0,0,0,0,0,0)))
colnames(phenotype_snp_est)<-c("trait","SNP","beta","se","pvalue","EAF","sample_size", "outcome", "exposure")

for (outcome in c("eduyears2","out_highbloodpressure","out_diabetes","out_height","out_bmi")){
  snps<-grep(names(data_sibs),pattern = "^rs",value = TRUE)
  for  (exposure in snps){
    message(exposure)
    message(outcome)
    phenotype_snp_est<-rbind(phenotype_snp_est,snp_plm(outcome,exposure))
  }
}

write.table(phenotype_snp_est,"results/individual_snps_v8.txt",quote=FALSE,row.names = FALSE)
