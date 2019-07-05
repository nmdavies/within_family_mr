# This script creates a sibling dataset with geno and phenotype information
# for analysis
rm(list = ls())

# dependencies
library('dplyr')

# input
# add height snps
load(file="/home/benb/projects/ped/data/harm_geno_grs.Rda")

a <- geno
rm(geno)

load(file="#####/harm_geno_grs.Rda")
load(file="#####/kin0.Rda")
pheno <- read.delim2("####/pheno_genotypeID.txt", header=T)

# merge height and bmi geno
geno <- merge(geno, a, by = "IID")
colnames(geno)[colnames(geno)=="cnt"] <- "cnt.hei.389"
colnames(geno)[colnames(geno)=="grs"] <- "grs.hei.389"

# output
output_file_name <- "#####/allhunt_sibs_v2.Rda"

# Display objects
head(kin[,-1:-4]) # ID columns not displayed
head(pheno[,-1]) # ID column not displayed
names(geno)

# Merge pheno geno
data <- merge(pheno, geno, by = "IID")

# 3. Define variables

# Education
#HUNT
#1 Primary education -> 10 (left at age 14 started at 4, 10 years of education)
#2 Lower secondary -> 12 (left at 16 started at 4, 12 years education)
#3 Upper secondary -> 14 (left at 18 started at 4, 14 years education 
#4 Degree -> 17 (left at 21 started at 4, 17 years education)

#1->10 (second stage of basic education) - this is basic education
#2->12 (Lower) secondary education)) - to be completed for qualifying for a trade
#3->13 (Upper secondary education) - to be completed for qualifying for university
#4->16 (First stage tertiary education) - Bachelor
#5->17 + (Second stage tertiary education) - Master, PhD

levels(data$Educ.NT2BLQ1)
#[1] "High school, intermediate school, vocational school, 1-2 years high school"
#[2] "Primary school 7-10 years, continuation school, folk high school"
#[3] "University/college, 4 years or more"
#[4] "University or other post-secondary education, less than 4 years"
#[5] "University qualifying examination, junior college, A levels"
levels(data$Educ.NT2BLQ1) <- c(2,1,5,4,3)

#1->10 (second stage of basic education) - this is basic education
#2->12 (Lower) secondary education)) - to be completed for qualifying for a trade
#3->13 (Upper secondary education) - to be completed for qualifying for university
#4->16 (First stage tertiary education) - Bachelor, Master, PhD

levels(data$Educ.NT2BLQ1) <- c(12,10,16,16,13)
data$eduyears2 <- as.numeric(as.character(data$Educ.NT2BLQ1))

table(data$eduyears2)
#10    12    13    16
#20722 19036  3790  8074

mean(data$eduyears2)
#[1] 11.89621

sd(data$eduyears2)
#[1] 2.042277

# Blood pressure
# Hypertension was defined as systolic blood pressure (SBP)???140 mmHg or diastolic blood pressure (DBP)???90mmHg or use of anti-hypertensive medication.
#BPMedCu.NT2BLQ1
#BPSystMn23.NT2BLM
#BPDiasMn23.NT2BLM

data$out_highbloodpressure <- NA
data$out_highbloodpressure[!is.na(data$BPMedCu.NT2BLQ1) | !is.na(data$BPSystMn23.NT2BLM) | !is.na(data$BPDiasMn23.NT2BLM)] = 0
data$out_highbloodpressure[data$BPMedCu.NT2BLQ1=="Currently taking medication" | data$BPSystMn23.NT2BLM>=140 | data$BPDiasMn23.NT2BLM>=90] = 1
#table(data$out_highbloodpressure)
#0     1
#29389 22233

# Diabetes
table(data$DiaEv.NT2BLQ1)
#No   Yes
#50163  1397
data$out_diabetes <- NA
data$out_diabetes[data$DiaEv.NT2BLQ1=="No"] = 0
data$out_diabetes[data$DiaEv.NT2BLQ1=="Yes"] = 1
table(data$out_diabetes)

# Height 
data$out_height <- data$Hei.NT2BLM
summary(data$out_height)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#120.0   164.0   170.0   170.6   177.0   206.0

# BMI
data$out_bmi <- as.numeric(as.character(data$BMI.NT2BLM)) #Find where this is converted to factor
summary(data$out_bmi)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#  15.00   23.70   25.90   26.37   28.50   52.50      59

# GRS height
data$wood_387_wprs <- data$grs.hei.389
summary(data$wood_387_wprs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#10.38   11.61   11.86   11.86   12.11   13.40

# GRS BMI
data$locke_79_wprs <- data$grs.bmi.79
summary(data$locke_79_wprs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.352   1.869   1.973   1.972   2.075   2.628

# Sex
summary(data$Hei.NT2BLM[data$Sex==1]) # Males
data$cov_male <- NA
data$cov_male[data$Sex==1] = 1
data$cov_male[data$Sex==2] = 0
table(data$cov_male)
#0     1
#26322 25300

# Age
data$cov_age <- as.numeric(as.character(data$PartAg.NT2BLQ1))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#19.20   37.80   48.00   49.28   60.20   98.6

# PC 1-20
data[,c(416:435)]= apply(data[,c(416:435)], 2, function(x) as.numeric(as.character(x))) # Need to add stringAsfactors=F on read
names(data)[416:435] <- tolower(names(data)[416:435])

# Batch
table(data$batch)
#0     1
#10858 40764

# 4 .Select one unrelated individual for singlton analysis
data_singltons<-subset(data,UNRELATED==T) # 2nd degree, all HUNT individuals, might be better soulution than singltons
dim(data_singltons)

# Create the data set for sibling analysis
# for plm this needs to be one data.frame with a index for sib pairs

# Add famid to to kin
kin$ID1 <- as.character(kin$ID1)
kin$ID2 <- as.character(kin$ID2)
kin$fid <- NA

i=0
for (i in 1:nrow(kin)){
  print(i)
  id1 <- kin$ID1[i]
  id2 <- kin$ID2[i]
  
  kin$fid[(kin$ID1 == id1 | kin$ID1 == id2 | kin$ID2 == id1 | kin$ID2 == id2) & is.na(kin$fid)==T] <- i
}
head(kin)

# Add sib id
kin$sid <- rownames(kin)

# Rbind siblings
sib1 <- kin[c("ID1", "fid")]
sib2 <- kin[c("ID2", "fid")]
colnames(sib2) <- c("ID1", "fid")

sibs <- unique(rbind(sib1, sib2))

# Merge kin with pheno_geno
df <- merge(sibs[c("ID1", "fid")], data, by.x = "ID1", by.y = "IID" )

# data frame to save
dim(df)

# save
save(df, data, file = output_file_name)