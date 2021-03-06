# Neil Davies 06/12/18
# Modified BMB 08/12/18
# This runs option B 
# within family - PLM method with robust SE
# population (one sample per family) - Standard IV
# source('reg_1_option_b-v8.R')

# 1. Load libraries

#install.packages('plm')
#install.packages("ivpack")
#install.packages("data.table") 

library('data.table')
library('plm')
library('ivpack')
library('dplyr')
library('multiwayvcov')
library('lmtest')

#Delete all objects
rm(list = ls())

# 2. Load data

#Load genotypic and phenotypic data
input_file_name <- "data/allhunt_sibs_v4.Rda"
load(input_file_name)

# Number participated HUNT2 and genotyped
table(data$Part.NT2BLQ1)

# Number genotyped with complete data
cd <- data[c("IID", "cov_age", "cov_male", "locke_69_wprs", "wood_386_wprs",
                                     "out_height","out_bmi", "eduyears2","out_highbloodpressure","out_diabetes")]
cd_nomiss  <- na.omit(cd)
dim(cd_nomiss)

# 4 Select one sibling for singlton analysis
snps<-grep(names(df),pattern = "^rs",value = TRUE)
data_sing_nomiss <- df[c("ID1", "fid", "cov_age", "cov_male", "locke_69_wprs", "wood_386_wprs",
                                     "out_height","out_bmi", "eduyears2","out_highbloodpressure","out_diabetes", paste0("pc", 1:20), "batch", snps)]
data_sing_nomiss  <- na.omit(data_sing_nomiss)

# Order siblings
data_sing_nomiss <- data_sing_nomiss %>%
  group_by(fid) %>%
  mutate(ranks = order(order(ID1, decreasing=TRUE)))

data_singltons<-subset(data_sing_nomiss, ranks==1) # Select the first id out of each group
dim(data_singltons)

# Descriptive stats
data_sing_nomiss <- data_singltons
summary(data_sing_nomiss)
summary(data_sing_nomiss$cov_age)
sd(data_sing_nomiss$cov_age)
summary(data_sing_nomiss$eduyears2)
sd(data_sing_nomiss$eduyears2)
summary(data_sing_nomiss$out_height)
sd(data_sing_nomiss$out_height)
summary(data_sing_nomiss$out_height[data_sing_nomiss$cov_male==0])
sd(data_sing_nomiss$out_height[data_sing_nomiss$cov_male==0])
summary(data_sing_nomiss$out_height[data_sing_nomiss$cov_male==1])
sd(data_sing_nomiss$out_height[data_sing_nomiss$cov_male==1])
summary(data_sing_nomiss$out_bmi)
sd(data_sing_nomiss$out_bmi)
sapply(data_sing_nomiss, sd, na.rm=T)
prop.table(table(data_sing_nomiss$out_diabetes))
prop.table(table(data_sing_nomiss$out_highbloodpressure))

# rename data
data_sibs <- df

data_sibs$famid <- data_sibs$fid

# How many individual sibs
length(unique(df$GID))

# 5. Check Z > X
summary(lm(wood_386_wprs ~ out_height, data = data_singltons))
#Residual standard error: 0.3613 on 20369 degrees of freedom
#(5816 observations deleted due to missingness)
#Multiple R-squared:  0.05302,	Adjusted R-squared:  0.05297
#F-statistic:  1140 on 1 and 20369 DF,  p-value: < 2.2e-16

summary(lm(locke_69_wprs ~ out_bmi, data = data_singltons))
#Residual standard error: 0.1523 on 20323 degrees of freedom
#(5862 observations deleted due to missingness)
#Multiple R-squared:  0.02557,	Adjusted R-squared:  0.02552
#F-statistic: 533.2 on 1 and 20323 DF,  p-value: < 2.2e-16

# Remove missing sibs
data_sibs_nomiss <- data_sibs[c("ID1", "cov_age", "cov_male", "locke_69_wprs", "wood_386_wprs",
                                "out_height","out_bmi", "eduyears2","out_highbloodpressure","out_diabetes", "famid", paste0("pc", 1:20), "batch", snps)]

data_sibs_nomiss <- na.omit(data_sibs_nomiss)

a <- as.data.frame(table(data_sibs_nomiss$famid))
table(a$Freq)
b <- a[a$Freq>=2,]
keep <- as.numeric(as.character(b$Var1))

data_sibs_nomiss  <- data_sibs_nomiss[(data_sibs_nomiss$famid %in% keep), ]
length(unique(data_sibs_nomiss$ID1))
check <- as.data.frame(table(data_sibs_nomiss$famid))
table(check$Freq)
length(unique(data_sibs_nomiss$famid))

data_sibs  <- data_sibs_nomiss

# Check numbers
print(dim(data_singltons))
#[1] 13347    12
print(length(unique(data_sibs$ID1)))
#[1] 28823
#[1] 28777

# Descriptives sibs
table(data_sibs$cov_male)
prop.table(table(data_sibs$cov_male))

t <- data.frame(table(data_sibs$famid))
t$n <- as.numeric(t$Freq)
sum(t$n)
t$n2 <- t$n -2
sum(t$n2)
t$n3 <- 1
sum(t$n3)
#> 11448+5881
#[1] 17329

# 6. Run the within family fixed effects with 'plm' and robust standard errors
results_ols<-c(0,0,0,0,0,0)
results_ols_wf<-c(0,0,0,0,0,0)
results_nwf<-c(0,0,0,0,0,0)
results_nwf_sibs<-c(0,0,0,0,0,0)
results_wf<-c(0,0,0,0,0,0)

for (outcome in c("eduyears2","out_highbloodpressure","out_diabetes")){
  for (exposure in c("out_height","out_bmi")){
    if (exposure=="out_height"){
      iv<-"wood_386_wprs"
    } else{
      iv<-"locke_69_wprs"
    }
    # OLS
    lm_mod <- lm(get(outcome) ~ get(exposure), data = data_sibs)
    vcov_firm <- cluster.vcov(lm_mod, data_sibs$famid)
    estimates_robust<-coeftest(lm_mod, vcov_firm)
    
    n_obs <- (nobs(lm_mod)) 
    beta<-estimates_robust[2,1]
    standard_error<-estimates_robust[2,2]
    p_value<-estimates_robust[2,4]
    results_ols<-rbind(results_ols,c(exposure,outcome,n_obs,beta,standard_error,p_value))
    
    # OLS family FE
    estimates<-plm(get(outcome)~get(exposure)+cov_age+cov_male, data = data_sibs, index = "famid", model = "within")
    estimates_robust<-coeftest(estimates,vcov=vcovHC(estimates,type="HC0",cluster="group"))
    
    n_obs <- (nobs(estimates)) # Is plm using sibs when one is missing information?
    beta<-estimates_robust[1,1]
    standard_error<-estimates_robust[1,2]
    p_value<-estimates_robust[1,4]
    results_ols_wf<-rbind(results_ols_wf,c(exposure,outcome,n_obs,beta,standard_error,p_value))
    
    # Std MR
    iv_reg_mod <- ivreg(get(outcome)~get(exposure)+cov_age+cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+batch  | get(iv)+cov_age + cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+batch, x = T, data = data_singltons)
    estimates<-robust.se(iv_reg_mod)
    
    n_obs<-nobs(iv_reg_mod) 
    beta<-estimates[2,1]
    standard_error<-estimates[2,2]
    p_value<-estimates[2,4]
    results_nwf<-rbind(results_nwf,c(exposure,outcome,n_obs,beta,standard_error,p_value))
    
    #Std MR but with sibs
    iv_reg_mod <- ivreg(get(outcome)~get(exposure)+cov_age+cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+batch  | get(iv)+cov_age + cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+batch, x = T, data = data_sibs)
    estimates<-robust.se(iv_reg_mod)
    
    n_obs<-nobs(iv_reg_mod) 
    beta<-estimates[2,1]
    standard_error<-estimates[2,2]
    p_value<-estimates[2,4]
    results_nwf_sibs<-rbind(results_nwf_sibs,c(exposure,outcome,n_obs,beta,standard_error,p_value))
    
    # MR within families
    estimates<-plm(get(outcome)~get(exposure)+cov_age+cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+batch | get(iv) +cov_age + cov_male+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+batch, data = data_sibs, index = "famid", model = "within", inst.method = "bvk")
    estimates_robust<-coeftest(estimates,vcov=vcovHC(estimates,type="HC0",cluster="group"))
    
    n_obs <- (nobs(estimates)) # Is plm using sibs when one is missing information?
    beta<-estimates_robust[1,1]
    standard_error<-estimates_robust[1,2]
    p_value<-estimates_robust[1,4]
    results_wf<-rbind(results_wf,c(exposure,outcome,n_obs,beta,standard_error,p_value))
    print(exposure)
    print(outcome)
  }
}

#Add column names
colnames(results_nwf)<-c("exposure","outcome","N","beta","se","p-value")
colnames(results_wf)<-c("exposure","outcome","N","beta","se","p-value")

write.table(results_nwf,"results/r_non_within_v7.txt",quote=FALSE,row.names = FALSE)
write.table(results_wf,"results/r_within_v7.txt",quote=FALSE,row.names = FALSE)

write.table(results_ols,"results/r_non_within_ols_v7.txt",quote=FALSE,row.names = FALSE)
write.table(results_ols_wf,"results/r_within_ols_v7.txt",quote=FALSE,row.names = FALSE)
write.table(results_nwf_sibs,"results/r_non_within_sibs_v7.txt",quote=FALSE,row.names = FALSE)

#save
save(data_sibs, file = "data/data_sibs_v4.RData")

print("done")
