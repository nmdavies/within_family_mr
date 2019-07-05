# Neil Davies 15/03/19
# Forest plot figue 2 sib MR for the MR-Egger and pleiotropy results
# This uses data from the Stata script reg4_mr_egger_weighted_median.do


install.packages("readstata13")
install.packages("forestplot")

library(readstata13)
library(forestplot)

dat <- read.dta13(paste(path1,"workingdata/mregger_results_for_r_no_cons.dta",sep=""))

# EDUCATION AND BMI GRAPH

df <- dat[dat$outcome=="eduyears3" & dat$exposure=="out_bmi",]
coef_ci=c("Mean difference (95% CI)",paste(format(round(df$coef,2), nsmall = 2)," (",format(round(df$lci,2), nsmall = 2)," to ",format(round(df$uci,2), nsmall = 2),")",sep = ""))
N=c("N",rep(61008,4))

pval<-c("P-value",formatC(df$pval,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef), 
    lower = c(NA,df$lci),
    upper = c(NA,df$uci)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -5L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  c("Effect of BMI on education", "MR-IVW", "MR weighted median", 
    "MR weighted modal", "MR-Egger slope"),
  coef_ci,
  pval)

# Plot

# 2. Create the plot

postscript(paste(path1,"results/graph_edu_bmi_egger.eps",sep=""), onefile=FALSE, width = 20, height = 2)
forestplot(tabletext, 
           all_results,
           graph.pos = 2, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,6)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           #clip=c(0.1,2.5), # crop plot
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Mean difference (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7)) 
dev.off()

# EDUCATION AND HEIGHT GRAPH

df <- dat[dat$outcome=="eduyears3" & dat$exposure=="out_height",]
coef_ci=c("Mean difference (95% CI)",paste(format(round(df$coef,2), nsmall = 2)," (",format(round(df$lci,2), nsmall = 2)," to ",format(round(df$uci,2), nsmall = 2),")",sep = ""))

pval<-c("P-value",formatC(df$pval,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef), 
    lower = c(NA,df$lci),
    upper = c(NA,df$uci)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -5L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  c("Effect of height on education", "MR-IVW", "MR weighted median", 
    "MR weighted modal", "MR-Egger slope"),
  coef_ci,
  pval)


# Plot

# 2. Create the plot
postscript(paste(path1,"results/graph_edu_height_egger.eps",sep=""), onefile=FALSE, width = 20, height = 2)
forestplot(tabletext, 
           all_results,
           graph.pos = 2, # graph position
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,6)), # define heading lines
           hrzl_lines = gpar(col="#444444"),
           clip=c(-0.03,0.05), # crop plot
           xlog=FALSE,
           col=fpColors(box="black",line="red"),
           boxsize = .1, # set the box size
           xlab="Mean difference (95% CI) ", 
           txt_gp = fpTxtGp(cex = 0.7)) 
dev.off()

# DIABETES AND BMI GRAPH

df <- dat[dat$outcome=="diabetes" & dat$exposure=="out_bmi",]
coef_ci=c("Risk difference*100 (95% CI)",paste(format(round(df$coef,2), nsmall = 2)," (",format(round(df$lci,2), nsmall = 2)," to ",format(round(df$uci,2), nsmall = 2),")",sep = ""))

pval<-c("P-value",formatC(df$pval,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef), 
    lower = c(NA,df$lci),
    upper = c(NA,df$uci)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -5L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  c("Effect of BMI on diabetes", "MR-IVW", "MR weighted median", 
    "MR weighted modal", "MR-Egger slope"),
  coef_ci,
  pval)


# Plot

# 2. Create the plot
postscript(paste(path1,"results/graph_diabetes_bmi_egger.eps",sep=""), onefile=FALSE, width = 20, height = 2)
forestplot(tabletext, 
           all_results,
           graph.pos = 2, # graph position
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


# BMI AND HIGHT BLOOD PRESSURE GRAPH

df <- dat[dat$outcome=="hbp2" & dat$exposure=="out_bmi",]
coef_ci=c("Risk difference*100 (95% CI)",paste(format(round(df$coef,2), nsmall = 2)," (",format(round(df$lci,2), nsmall = 2)," to ",format(round(df$uci,2), nsmall = 2),")",sep = ""))

pval<-c("P-value",formatC(df$pval,format = "e", digits = 2))

# Create results table
all_results <- 
  structure(list(
    mean  = c(NA,df$coef), 
    lower = c(NA,df$lci),
    upper = c(NA,df$uci)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -5L),
    class = "data.frame")

# Create text
tabletext<-cbind(
  c("Effect of BMI on high blood pressure", "MR-IVW", "MR weighted median", 
    "MR weighted modal", "MR-Egger slope"),
  coef_ci,
  pval)


# Plot

# 2. Create the plot
postscript(paste(path1,"results/graph_hbp2_bmi_egger.eps",sep=""), onefile=FALSE, width = 20, height = 2)
forestplot(tabletext, 
           all_results,
           graph.pos = 2, # graph position
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