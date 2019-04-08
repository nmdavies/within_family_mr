# Neil Davies 15/03/19
# Forest plot figue 2 sib MR
# This uses data from the Stata script reg_1_height_education.do

install.packages("readstata13")
install.packages("forestplot")

library(readstata13)
library(forestplot)

dat <- read.dta13(paste($path1,"results/final_table_r.dta",sep(""))

# EDUCATION AND BMI GRAPH

df <- dat[dat$outcome=="eduyears3" & dat$exposure=="out_bmi",]
coef_ci=c("Beta (95% CI)",paste(format(round(df$coef,2), nsmall = 2)," (",format(round(df$lci,2), nsmall = 2)," to ",format(round(df$uci,2), nsmall = 2),")",sep = ""))
N=c("N",df$N)
N[6]<-61361
N[7]<-61361
pval<-c("P-value",formatC(df$pval,format = "e", digits = 2))


# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef), 
    lower = c(NA,df$lci),
    upper = c(NA,df$uci)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -6L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  c("Method", "Ordinary least squares", "MR-PRS unrelateds", 
    "MR-PRS siblings", "MR-PRS siblings family effects", "2MR-IVW siblings", "2SMR-IVW siblings - split sample"),
  c("Cohort", rep("UKBB",4), "UKBB+HUNT", "UKBB+HUNT"),
  N,
  coef_ci,
  pval)

# Plot

# 2. Create the plot
postscript(paste($path1,"results/graph_edu_bmi.eps",sep(""), onefile=FALSE,width = 20, height = 2)
forestplot(tabletext, 
           all_results,
           graph.pos = 4, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,6)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           #clip=c(0.1,2.5), # crop plot
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Beta (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7)) 
dev.off()

# EDUCATION AND HEIGHT GRAPH

df <- dat[dat$outcome=="eduyears3" & dat$exposure=="out_height",]
coef_ci=c("Beta (95% CI)",paste(format(round(df$coef,2), nsmall = 2)," (",format(round(df$lci,2), nsmall = 2)," to ",format(round(df$uci,2), nsmall = 2),")",sep = ""))
N=c("N",df$N)
N[6]<-61414
N[7]<-61414
pval<-c("P-value",formatC(df$pval2,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef), 
    lower = c(NA,df$lci),
    upper = c(NA,df$uci)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -6L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  c("Method", "Ordinary least squares", "MR-PRS unrelateds", 
    "MR-PRS siblings", "MR-PRS siblings family effects", "2MR-IVW siblings", "2SMR-IVW siblings - split sample"),
  c("Cohort", rep("UKBB",4), "UKBB+HUNT", "UKBB+HUNT"),
  N,
  coef_ci,
  pval)

# Plot

# 2. Create the plot
postscript(paste($path1,"results/graph_edu_height.eps",sep(""), onefile=FALSE,width = 20, height = 2)
forestplot(tabletext, 
           all_results,
           graph.pos = 4, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,6)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           #clip=c(0.1,2.5), # crop plot
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Beta (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7)) 
dev.off()


# BMI AND BLOOD PRESSURE GRAPH

df <- dat[dat$outcome=="out_highbloodpressure2" & dat$exposure=="out_bmi",]
coef_ci=c("Risk difference*100 (95% CI)",paste(format(round(df$coef*100,2), nsmall = 2)," (",format(round(df$lci*100,2), nsmall = 2)," to ",format(round(df$uci*100,2), nsmall = 2),")",sep = ""))
N=c("N",df$N)
N[6]<-61361
N[7]<-61361
pval<-c("P-value",formatC(df$pval2,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef*100), 
    lower = c(NA,df$lci*100),
    upper = c(NA,df$uci*100)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -6L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  c("Method", "Ordinary least squares", "MR-PRS unrelateds", 
    "MR-PRS siblings", "MR-PRS siblings family effects", "2MR-IVW siblings", "2SMR-IVW siblings - split sample"),
  c("Cohort", rep("UKBB",4), "UKBB+HUNT", "UKBB+HUNT"),
  N,
  coef_ci,
  pval)

# Plot

# 2. Create the plot
postscript(paste($path1,"results/graph_bloodpressure_bmi.eps",sep(""), onefile=FALSE,width = 20, height = 2)
forestplot(tabletext, 
           all_results,
           graph.pos = 4, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,6)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           #clip=c(0.1,2.5), # crop plot
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Risk difference*100 (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7)) 
dev.off()


# BMI AND DIABETES GRAPH

df <- dat[dat$outcome=="out_diabetes" & dat$exposure=="out_bmi",]
coef_ci=c("Risk difference*100 (95% CI)",paste(format(round(df$coef*100,2), nsmall = 2)," (",format(round(df$lci*100,2), nsmall = 2)," to ",format(round(df$uci*100,2), nsmall = 2),")",sep = ""))
N=c("N",df$N)
N[6]<-61110
N[7]<-61110
pval<-c("P-value",formatC(df$pval2,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef*100), 
    lower = c(NA,df$lci*100),
    upper = c(NA,df$uci*100)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -6L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  c("Method", "Ordinary least squares", "MR-PRS unrelateds", 
    "MR-PRS siblings", "MR-PRS siblings family effects", "2MR-IVW siblings", "2SMR-IVW siblings - split sample"),
  c("Cohort", rep("UKBB",4), "UKBB+HUNT", "UKBB+HUNT"),
  N,
  coef_ci,
  pval)

# Plot

# 2. Create the plot
postscript(paste($path1,"results/graph_diabetes_bmi.eps",sep(""), onefile=FALSE,width = 20, height = 2)
forestplot(tabletext, 
           all_results,
           graph.pos = 4, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,6)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           #clip=c(0.1,2.5), # crop plot
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Risk difference*100 (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7)) 
dev.off()