# Neil Davies 15/03/19
# Forest plot figue 2 sib MR
# This uses data from the Stata script reg_1_height_education.do

install.packages("readstata13")
install.packages("forestplot")
install.packages("adimpro")

library(readstata13)
library(forestplot)
library(lattice) 
library(grid) 
library(ggplot2)
library(adimpro)

dat <- read.dta13(paste(path1,"results/final_table_r.dta",sep=""))

#Text for figures 

# EDUCATION AND BMI GRAPH
method<-c("Effect of BMI on education", "IPD Ordinary least squares","IPD OLS Family FE", "IPD MR-PRS unrelateds", 
          "IPD MR-PRS siblings", "IPD MR-PRS siblings family FE", "2SMR IVW siblings", "2SMR IVW siblings - split sample")

df <- dat[dat$outcome=="eduyears" & dat$exposure=="out_bmi",]
coef_ci=c("       Mean difference (95% CI)        ",paste(format(round(df$coef,2), nsmall = 2)," (",format(round(df$lci,2), nsmall = 2)," to ",format(round(df$uci,2), nsmall = 2),")",sep = ""))
N=c("N",df$N)
N[7]<-61008 
N[8]<-61008
pval<-c("P-value",formatC(df$pval,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef), 
    lower = c(NA,df$lci),
    upper = c(NA,df$uci)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -8L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  method, 
  N,
  coef_ci,
  pval)

# Plot

# 2. Create the plot
postscript(paste(path1,"results/graph_edu_bmi.eps",sep=""), onefile=FALSE,width = 20, height = 2)
graph1<-forestplot(tabletext, 
           all_results,
           graph.pos = 4, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,7)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           #clip=c(0.1,2.5), # crop plot
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Mean difference (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7))
dev.off()


# EDUCATION AND HEIGHT GRAPH
method<-c("Effect of height on education", "IPD Ordinary least squares","IPD OLS Family FE", "IPD MR-PRS unrelateds", 
          "IPD MR-PRS siblings", "IPD MR-PRS siblings family FE", "2SMR IVW siblings", "2SMR IVW siblings - split sample")

df <- dat[dat$outcome=="eduyears" & dat$exposure=="out_height",]
coef_ci=c("Mean difference (95% CI)",paste(format(round(df$coef,2), nsmall = 2)," (",format(round(df$lci,2), nsmall = 2)," to ",format(round(df$uci,2), nsmall = 2),")",sep = ""))
N=c("N",df$N)
N[7]<-61008 
N[8]<-61008
pval<-c("P-value",formatC(df$pval2,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef), 
    lower = c(NA,df$lci),
    upper = c(NA,df$uci)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -7L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  method,
  N,
  coef_ci,
  pval)

# Plot

# 2. Create the plot
postscript(paste(path1,"results/graph_edu_height.eps",sep=""), onefile=FALSE,width = 20, height = 2)
graph2<-forestplot(tabletext, 
           all_results,
           graph.pos = 4, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,7)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           #clip=c(0.1,2.5), # crop plot
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Mean difference (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7)) 
dev.off()


# BMI AND BLOOD PRESSURE GRAPH
method<-c("Effect of BMI on HBP", "IPD Ordinary least squares","IPD OLS Family FE", "IPD MR-PRS unrelateds", 
          "IPD MR-PRS siblings", "IPD MR-PRS siblings family FE", "2SMR IVW siblings", "2SMR IVW siblings - split sample")

df <- dat[dat$outcome=="out_highbloodpressure2" & dat$exposure=="out_bmi",]
coef_ci=c("Risk difference*100 (95% CI)",paste(format(round(df$coef*100,2), nsmall = 2)," (",format(round(df$lci*100,2), nsmall = 2)," to ",format(round(df$uci*100,2), nsmall = 2),")",sep = ""))
N=c("N",df$N)
N[7]<-61008 
N[8]<-61008
pval<-c("P-value",formatC(df$pval2,format = "e", digits = 2))
pval[2]<-"<1.00e-300"
pval[3]<-"<1.00e-300"

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef*100), 
    lower = c(NA,df$lci*100),
    upper = c(NA,df$uci*100)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -7L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  method,
  N,
  coef_ci,
  pval)

# Plot

# 2. Create the plot
postscript(paste(path1,"results/graph_bloodpressure_bmi.eps",sep=""), onefile=FALSE,width = 20, height = 2)
graph3<-forestplot(tabletext, 
           all_results,
           graph.pos = 4, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,7)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           #clip=c(0.1,2.5), # crop plot
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Risk difference*100 (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7)) 
dev.off()


# BMI AND DIABETES GRAPH
method<-c("Effect of BMI on diabetes", "IPD Ordinary least squares","IPD OLS Family FE", "IPD MR-PRS unrelateds", 
          "IPD MR-PRS siblings", "IPD MR-PRS siblings family FE", "2SMR IVW siblings", "2SMR IVW siblings - split sample")

df <- dat[dat$outcome=="out_diabetes" & dat$exposure=="out_bmi",]
coef_ci=c("Risk difference*100 (95% CI)",paste(format(round(df$coef*100,2), nsmall = 2)," (",format(round(df$lci*100,2), nsmall = 2)," to ",format(round(df$uci*100,2), nsmall = 2),")",sep = ""))
N=c("N",df$N)
N[7]<-61008 
N[8]<-61008
pval<-c("P-value",formatC(df$pval2,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef*100), 
    lower = c(NA,df$lci*100),
    upper = c(NA,df$uci*100)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -7L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  method,
  N,
  coef_ci,
  pval)

# Plot

# 2. Create the plot
postscript(paste(path1,"results/graph_diabetes_bmi.eps",sep=""), onefile=FALSE,20, height = 2)
graph4<-forestplot(tabletext, 
           all_results,
           graph.pos = 4, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,7)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           #clip=c(0.1,2.5), # crop plot
           #graphwidth = unit(2.95,"cm"),
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Risk difference*100 (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7)) 
dev.off()

#Paste figures togeather
plot(graph4)