# Neil Davies 06/12/18
# This runs option B 
# within family - PLM method with robust SE
# population (one sample per family) - Standard IV

# 1. Load libraries

#install.packages('plm')
#install.packages("ivpack")
#install.packages("data.table") 

library(data.table)
library("plm")
library(ivpack)
library(data.table)  
library('dplyr')

#Delete all objects
rm(list = ls())

# 2. Load data

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

#Select one sibling for singlton analysis
data_singltons<-subset(data,famid!="NA" & within_fam_id==1)
dim(data_singltons)

#Run the within family fixed effects with 'plm' and robust standard errors
results_nwf<-c(0,0,0,0,0,0)
results_wf<-c(0,0,0,0,0,0)
for (outcome in c("eduyears2","out_highbloodpressure","out_diabetes")){
  for (exposure in c("out_height","out_bmi")){
    if (exposure=="out_height"){
      iv<-"wood_387_wprs"
    } else{
      iv<-"locke_79_wprs"
    }
    estimates<-robust.se(ivreg(get(outcome)~get(exposure)+cov_age+cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20  | get(iv)+cov_age + cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20, x = T, data = data_singltons))
    beta<-estimates[2,1]
    standard_error<-estimates[2,2]
    p_value<-estimates[2,4]
    results_nwf<-rbind(results_nwf,c(outcome,exposure,beta,standard_error,p_value))
    
    estimates<-plm(get(outcome)~get(exposure)+cov_age+cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20 | get(iv) +cov_age + cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20, data = data_sibs, index = "famid", model = "within", inst.method = "bvk")
    estimates_robust<-coeftest(estimates,vcov=vcovHC(estimates,type="HC0",cluster="group"))
    
    beta<-estimates_robust[1,1]
    standard_error<-estimates_robust[1,2]
    p_value<-estimates_robust[1,4]
    results_wf<-rbind(results_wf,c(outcome,exposure,beta,standard_error,p_value))    
  }
}

#Add column names
colnames(results_nwf)<-c("exposure","outcome","N","beta","se","p-value")
colnames(results_wf)<-c("exposure","outcome","N","beta","se","p-value")

write.table(results_nwf[2:7,],paste($path1,"/results/r_non_within.txt",sep(""),quote=FALSE,row.names = FALSE)
write.table(results_wf[2:7,],paste($path1,"results/r_within.txt",sep(""),quote=FALSE,row.names = FALSE)

