########################################################################################################
# Response to Technical Comment by Abruzzi et al. (2021). Code for:
# Figure 1C
########################################################################################################

# ######################################################################################################
# set the working directory
setwd("~/Downloads/")
rm(list = ls(all = TRUE))
# ######################################################################################################


# ######################################################################################################
# load all the needed packages
# NB: you may need to install some packages before running this script
# ######################################################################################################

library(R.utils)
library(tools)
library(tidyr)
library(dplyr)
library(reshape2)
library(rain)
library(scales)
library(lmtest)
library(gplots)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(heatmap3)
library(circular)
library(eulerr)
library(RColorBrewer)
library(plotrix)
library(UpSetR)

source("~/Downloads/JTK_CYCLEv3.1.R") # downloaded from https://openwetware.org/wiki/HughesLab:JTK_Cycle

options(stringsAsFactors=FALSE)

# ######################################################################################################
# Download the data from GEO and unpack it to data directory
# ######################################################################################################

dir.create("data", showWarnings=FALSE)
setwd("data")

download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111696/suppl/GSE111696%5FBmal1%5FexLiver%5Fgenes%2Efpkm%5Ftracking%2Etxt%2Egz", "GSE111696_Bmal1_exLiver_genes.fpkm_tracking.txt.gz")
gunzip("GSE111696_Bmal1_exLiver_genes.fpkm_tracking.txt.gz", overwrite = T)

setwd("..")

# ######################################################################################################
# reading the unpacked data 
# ######################################################################################################

exLiver_fpkm <- read.delim("data/GSE111696_Bmal1_exLiver_genes.fpkm_tracking.txt", quote = "")

# ######################################################################################################
# subset the data
# ######################################################################################################

exLiver_WT_data <- exLiver_fpkm[, grep("WT_.*_FPKM", colnames(exLiver_fpkm))]
exLiver_KO_data <- exLiver_fpkm[, grep("KO_.*_FPKM", colnames(exLiver_fpkm))]
rownames(exLiver_WT_data) <- exLiver_fpkm$gene_id
rownames(exLiver_KO_data) <- exLiver_fpkm$gene_id

# ######################################################################################################
# check/view the data structure
# ######################################################################################################

colnames(exLiver_WT_data); head(exLiver_WT_data)
colnames(exLiver_KO_data); head(exLiver_KO_data)

# ######################################################################################################
# data cleanup, remove rows with any zero values, Log Transformation
# ######################################################################################################

# What we did 
# exLiver_WT_data <- exLiver_WT_data[apply(exLiver_WT_data,1,function(y) !any(y==0)),]
# exLiver_KO_data <- exLiver_KO_data[apply(exLiver_KO_data,1,function(y) !any(y==0)),]

# exLiver_WT_data_log <- log2(exLiver_WT_data)
# exLiver_KO_data_log <- log2(exLiver_KO_data)

# as done by Abruzzi et al.  
exLiver_WT_data <- subset(exLiver_WT_data,rowMeans(exLiver_WT_data) >3 & apply(exLiver_WT_data, 1, function(x) length(x[x==0])) < 6)
exLiver_WT_data_log <- log2(1+exLiver_WT_data)

exLiver_KO_data <- subset(exLiver_KO_data,rowMeans(exLiver_KO_data) >3 & apply(exLiver_KO_data, 1, function(x) length(x[x==0])) < 6)
exLiver_KO_data_log <- log2(1+exLiver_KO_data)

# ###################################################################################################
# Convert 69 hours data to a 24 hours data with replicates as done by Abruzzi et al. 
# ###################################################################################################
posi=seq(1,24,8)
posi.all=NULL
for(i in 0:7){posi.all=c(posi.all,posi+i)}

exLiver_WT_data_rep <- exLiver_WT_data_log[,posi.all]
exLiver_KO_data_rep <- exLiver_KO_data_log[,posi.all]

# ###################################################################################################
# Running JTK
# ###################################################################################################
# wrapper function for running JTK 
runJTK <- function(dat, sampling.time.sequence, sampling.resolution, replicates, period.low, period.high){

  jtkdist(timepoints=length(sampling.time.sequence),reps=replicates)
  periods <- (period.low/sampling.resolution):(period.high/sampling.resolution)
  jtk.init(periods, sampling.resolution)
  
  res <- apply(dat,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  
  res = as.data.frame(t(res))
  bhq = p.adjust(unlist(res[,1]),"BH")
  res = cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- res
} 

WT_Bmal_Tissue_log_jtk_results <- runJTK(dat=exLiver_WT_data_log, sampling.time.sequence=seq(0,69,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)
KO_Bmal_Tissue_log_jtk_results <- runJTK(dat=exLiver_KO_data_log, sampling.time.sequence=seq(0,69,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)

WT_Bmal_Tissue_rep_jtk_results <- runJTK(dat=exLiver_WT_data_rep, sampling.time.sequence=seq(0,21,by=3), sampling.resolution=3, replicates=3, period.low=24, period.high=24)
KO_Bmal_Tissue_rep_jtk_results <- runJTK(dat=exLiver_KO_data_rep, sampling.time.sequence=seq(0,21,by=3), sampling.resolution=3, replicates=3, period.low=24, period.high=24)

# ###################################################################################################
# Running harmonic regression
# ###################################################################################################
# wrapper function for running harmonic regression
harm_reg <- function(x, t, period){
  n=length(x)
  fit0=lm(x~1)
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)
  fit1=lm(x~c+s)
  a=coef(fit1)[2]
  b=coef(fit1)[3]
  p.val=lrtest(fit1, fit0)$Pr[2]
  amp=2*sqrt(a^2+b^2)
  phase=atan2(b,a)%%(2*pi)
  phase=period*phase/(2*pi)
  
  c(p.val_harm=p.val,phase_harm=phase,amp_harm=amp,period_arm=period)
}

sampling.time.log=rep(seq(0,69,3),each=1)  ### time sequence for what we did

WT_Bmal_Tissue_log_harmonics_results <- as.data.frame(t(apply(exLiver_WT_data_log,1,harm_reg,t=sampling.time.log,period=24)))
WT_Bmal_Tissue_log_harmonics_results$q.val_harm <- p.adjust(WT_Bmal_Tissue_log_harmonics_results$p.val_harm,method="BH")

KO_Bmal_Tissue_log_harmonics_results <- as.data.frame(t(apply(exLiver_KO_data_log,1,harm_reg,t=sampling.time.log,period=24)))
KO_Bmal_Tissue_log_harmonics_results$q.val_harm <- p.adjust(KO_Bmal_Tissue_log_harmonics_results$p.val_harm,method="BH")

sampling.time.rep=rep(seq(0,21,3),each=3)  ### time sequence for replicate format

WT_Bmal_Tissue_rep_harmonics_results <- as.data.frame(t(apply(exLiver_WT_data_rep,1,harm_reg,t=sampling.time.rep,period=24)))
WT_Bmal_Tissue_rep_harmonics_results$q.val_harm <- p.adjust(WT_Bmal_Tissue_rep_harmonics_results$p.val_harm,method="BH")

KO_Bmal_Tissue_rep_harmonics_results <- as.data.frame(t(apply(exLiver_KO_data_rep,1,harm_reg,t=sampling.time.rep,period=24)))
KO_Bmal_Tissue_rep_harmonics_results$q.val_harm <- p.adjust(KO_Bmal_Tissue_rep_harmonics_results$p.val_harm,method="BH")

# ###################################################################################################
# Running RAIN
# ###################################################################################################

# wrapper function for running RAIN with defaults settings
runRAIN <- function(dat, deltat, period.delta = 0, replicates) {
  rain(x = t(dat), deltat = deltat, period = 24, period.delta = period.delta, nr.series = replicates, adjp.method = "ABH", peak.border = c(0.1,0.9)) 
}

# ex-vivo liver
WT_Bmal_Tissue_log_rain_results <- runRAIN(dat = exLiver_WT_data_log, deltat = 3, replicates=1)
WT_Bmal_Tissue_log_rain_results$p.adj <- round(p.adjust(WT_Bmal_Tissue_log_rain_results$pVal,method='BH'), digit = 6)
KO_Bmal_Tissue_log_rain_results <- runRAIN(dat = exLiver_KO_data_log, deltat = 3, replicates=1)
KO_Bmal_Tissue_log_rain_results$p.adj <- round(p.adjust(KO_Bmal_Tissue_log_rain_results$pVal,method='BH'), digit = 6)

WT_Bmal_Tissue_rep_rain_results <- runRAIN(dat = exLiver_WT_data_rep, deltat = 3, replicates=3)
WT_Bmal_Tissue_rep_rain_results$p.adj <- round(p.adjust(WT_Bmal_Tissue_rep_rain_results$pVal,method='BH'), digit = 6)
KO_Bmal_Tissue_rep_rain_results <- runRAIN(dat = exLiver_KO_data_rep, deltat = 3, replicates=3)
KO_Bmal_Tissue_rep_rain_results$p.adj <- round(p.adjust(KO_Bmal_Tissue_rep_rain_results$pVal,method='BH'), digit = 6)

# ###################################################################################################
# save all the output dataframes
# ###################################################################################################
dir.create("output", showWarnings=FALSE)
current_ession <- sessionInfo()
save(WT_Bmal_Tissue_log_jtk_results, KO_Bmal_Tissue_log_jtk_results,
  WT_Bmal_Tissue_log_rain_results, KO_Bmal_Tissue_log_rain_results,
  WT_Bmal_Tissue_log_harmonics_results, KO_Bmal_Tissue_log_harmonics_results,
  WT_Bmal_Tissue_rep_jtk_results, KO_Bmal_Tissue_rep_jtk_results,
  WT_Bmal_Tissue_rep_rain_results, KO_Bmal_Tissue_rep_rain_results,
  WT_Bmal_Tissue_rep_harmonics_results, KO_Bmal_Tissue_rep_harmonics_results,
  file = "output/Results_Abruzzi_1C.RData" )

# ###################################################################################################
# output number of cycling gene
# ###################################################################################################
# load("~/Downloads/output/Results_Abruzzi_1C.RData")

thresh <- c(1e-5,1e-4,1e-3,1e-2,0.05,0.1,0.2) #,0.4,0.6,0.8,1)

output_nb.r <- function(RES.JTK,RES.RAIN,RES.HARM,thresh){ 
  cycler.JTK.pval=NULL
  cycler.JTK.qval=NULL
  cycler.RAIN.pval=NULL
  cycler.RAIN.qval=NULL
  cycler.HARM.pval=NULL
  cycler.HARM.qval=NULL
  for(k in thresh){
    cycler.JTK.pval=c(cycler.JTK.pval,nrow(subset(RES.JTK, ADJ.P < k )))
    cycler.JTK.qval=c(cycler.JTK.qval,nrow(subset(RES.JTK, BH.Q < k)))
    cycler.RAIN.pval=c(cycler.RAIN.pval,nrow(subset(RES.RAIN, pVal < k )))
    cycler.RAIN.qval=c(cycler.RAIN.qval,nrow(subset(RES.RAIN, p.adj < k )))
    cycler.HARM.pval=c(cycler.HARM.pval,nrow(subset(RES.HARM, p.val_harm < k )))
    cycler.HARM.qval=c(cycler.HARM.qval,nrow(subset(RES.HARM, q.val_harm < k )))
  }
  list(JTK.pval=cycler.JTK.pval, JTK.qval=cycler.JTK.qval, RAIN.pval=cycler.RAIN.pval, RAIN.qval=cycler.RAIN.qval, HARM.pval=cycler.HARM.pval, HARM.qval=cycler.HARM.qval )
}


nb.log.WT <- output_nb.r(RES.JTK=WT_Bmal_Tissue_log_jtk_results,RES.RAIN=WT_Bmal_Tissue_log_rain_results,RES.HARM=WT_Bmal_Tissue_log_harmonics_results,thresh)
nb.log.KO <- output_nb.r(RES.JTK=KO_Bmal_Tissue_log_jtk_results,RES.RAIN=KO_Bmal_Tissue_log_rain_results,RES.HARM=KO_Bmal_Tissue_log_harmonics_results,thresh)

DF_log <- data.frame(pval=c(nb.log.WT$JTK.pval,nb.log.KO$JTK.pval,nb.log.WT$RAIN.pval,nb.log.KO$RAIN.pval,nb.log.WT$HARM.pval,nb.log.KO$HARM.pval), 
              qval=c(nb.log.WT$JTK.qval,nb.log.KO$JTK.qval,nb.log.WT$RAIN.qval,nb.log.KO$RAIN.qval,nb.log.WT$HARM.qval,nb.log.KO$HARM.qval), 
              genotype=c(rep('WT JTK',length(nb.log.WT$JTK.pval)),rep('KO JTK',length(nb.log.KO$JTK.pval)),rep('WT Rain',length(nb.log.WT$RAIN.pval)),rep('KO Rain',length(nb.log.KO$RAIN.pval)), rep('WT Harmonic',length(nb.log.WT$HARM.pval)),rep('KO Harmonic',length(nb.log.KO$HARM.pval)) ), 
              threshold=c(thresh,thresh,thresh,thresh,thresh,thresh))


nb.rep.WT <- output_nb.r(RES.JTK=WT_Bmal_Tissue_rep_jtk_results,RES.RAIN=WT_Bmal_Tissue_rep_rain_results,RES.HARM=WT_Bmal_Tissue_rep_harmonics_results,thresh)
nb.rep.KO <- output_nb.r(RES.JTK=KO_Bmal_Tissue_rep_jtk_results,RES.RAIN=KO_Bmal_Tissue_rep_rain_results,RES.HARM=KO_Bmal_Tissue_rep_harmonics_results,thresh)

DF_rep <- data.frame(pval=c(nb.rep.WT$JTK.pval,nb.rep.KO$JTK.pval,nb.rep.WT$RAIN.pval,nb.rep.KO$RAIN.pval,nb.rep.WT$HARM.pval,nb.rep.KO$HARM.pval), 
              qval=c(nb.rep.WT$JTK.qval,nb.rep.KO$JTK.qval,nb.rep.WT$RAIN.qval,nb.rep.KO$RAIN.qval,nb.rep.WT$HARM.qval,nb.rep.KO$HARM.qval), 
              genotype=c(rep('WT JTK',length(nb.rep.WT$JTK.pval)),rep('KO JTK',length(nb.rep.KO$JTK.pval)),rep('WT Rain',length(nb.rep.WT$RAIN.pval)),rep('KO Rain',length(nb.rep.KO$RAIN.pval)), rep('WT Harmonic',length(nb.rep.WT$HARM.pval)),rep('KO Harmonic',length(nb.rep.KO$HARM.pval)) ), 
              threshold=c(thresh,thresh,thresh,thresh,thresh,thresh))

# ###################################################################################################
# Figure 1C 
# Plots on psuedo replicates as in the comments by Abruzzi et al. 
# ###################################################################################################

g3ko <- ggplot(DF_rep,aes(x=threshold,y=pval)) + 
      geom_point(aes(col=genotype)) + 
      geom_line(aes(col=genotype)) + 
      theme_minimal() + 
      xlab("p-value") + 
      ylab("Number of Genes") +
      theme(aspect.ratio=1, axis.title=element_text(size=15), axis.text=element_text(size=12)) + 
      scale_color_manual(values=c('darkred','orangered', 'maroon2', NA,NA,NA)) + 
      geom_text_repel(aes(label=pval,col=genotype),size=3,force=2,nudge_x=-0.3) + scale_x_log10(labels = function(x) format(x, scientific = TRUE))

g4ko <- ggplot(DF_rep,aes(x=threshold,y=qval)) + 
      geom_point(aes(col=genotype)) + 
      geom_line(aes(col=genotype)) + 
      theme_minimal() + 
      xlab("p-value BH corrected") + 
      ylab("Number of Genes") +
      theme(aspect.ratio=1, axis.title=element_text(size=15), axis.text=element_text(size=12)) + 
      scale_color_manual(values=c('darkred','orangered', 'maroon2', NA,NA,NA)) + 
      geom_text_repel(aes(label=qval,col=genotype),size=3,force=2,nudge_x=-0.3) + scale_x_log10(labels = function(x) format(x, scientific = TRUE))

pdf("output/Figure_1C_pvalue_qvalue.pdf", width=10, height=5 )
grid.arrange(g3ko,g4ko,ncol=2)
dev.off()


# ################################################################################################################


