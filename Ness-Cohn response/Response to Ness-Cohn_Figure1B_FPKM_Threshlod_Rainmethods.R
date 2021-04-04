# ##########################################################################################################################
# Response to Technical Comment by Ness-Cohn et al. (2021). Code for:
# Figure 1B
# ##########################################################################################################################

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
library(gplots)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(heatmap3)
library(eulerr)
library(RColorBrewer)
library(plotrix)
library(UpSetR)


# ######################################################################################################
# Download the data from GEO and unpack it to data directory
# ######################################################################################################

dir.create("data", showWarnings=FALSE)
setwd("data")

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE134333&format=file&file=GSE134333%5FDataset%5FAM%5FPM%5FTempComp%5FFPKM%5Fvalues%2Etxt%2Egz", "GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt.gz")
gunzip("GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt.gz", overwrite = T)

setwd("..")

# ######################################################################################################
# reading the unpacked data 
# ######################################################################################################

AM_PM_TempComp_FPKM <- read.delim("data/GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt", quote = "")

# ######################################################################################################
# subset the data
# ######################################################################################################

Bmal_WT_37C <- AM_PM_TempComp_FPKM[, grep("WT.*.37C.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_KO_37C <- AM_PM_TempComp_FPKM[, grep("Bmal1.*.37C.FPKM", colnames(AM_PM_TempComp_FPKM))]

rownames(Bmal_WT_37C) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_KO_37C) <- AM_PM_TempComp_FPKM$gene_id

# ######################################################################################################
# check/view the data structure
# ######################################################################################################

colnames(Bmal_WT_37C); head(Bmal_WT_37C)
colnames(Bmal_KO_37C); head(Bmal_KO_37C)

dir.create("output", showWarnings=FALSE)
# wrapper function for running RAIN 
runRAIN <- function(dat, deltat, period.delta = 0) {
  rain(x = t(dat), deltat = deltat, period = 24, period.delta = period.delta, adjp.method = "ABH", method = "longitudinal", peak.border = c(0.3,0.7)) 
}

# ######################################################################################################
# Loop for mean fpkm and rain calculation
# ######################################################################################################

fpkm_thresh_list <- c("none", 0, 0.05, 0.1, 0.5, 1, 1.5, 2, 3)

for(fpkmtresh in fpkm_thresh_list){

	thresh_fpkm_value=fpkmtresh
	# data cleanup, remove rows with lower fmkm than the thresh 
	if (thresh_fpkm_value == "none" ){
		WT_37C <- Bmal_WT_37C
		KO_37C <- Bmal_KO_37C
	} else {
		WT_37C <- subset(Bmal_WT_37C,rowMeans(Bmal_WT_37C)>thresh_fpkm_value & apply(Bmal_WT_37C, 1, function(x) length(x[x==0])) < 6)
		KO_37C <- subset(Bmal_KO_37C,rowMeans(Bmal_KO_37C)>thresh_fpkm_value & apply(Bmal_KO_37C, 1, function(x) length(x[x==0])) < 6)
	}
	
	# Log2 Transformation
	WT_37C_log <- log2(WT_37C+1)
	KO_37C_log <- log2(KO_37C+1)

	# ######################################################################################################
	# Running RAIN
	# ######################################################################################################

	WT37 <- runRAIN(dat = WT_37C_log, deltat = 2)
	KO37 <- runRAIN(dat = KO_37C_log, deltat = 2)

	WT37$p.adj <- round(p.adjust(WT37$pVal,method='BH'), digit = 6)
	KO37$p.adj <- round(p.adjust(KO37$pVal,method='BH'), digit = 6)

	assign(paste("WT37_fpkm_",thresh_fpkm_value,"_results",sep=""),WT37)
	assign(paste("KO37_fpkm_",thresh_fpkm_value,"_results",sep=""),KO37)

	rm(thresh_fpkm_value,WT_37C,KO_37C,WT_37C_log,KO_37C_log,WT37,KO37)
}

	WT_37C <- Bmal_WT_37C[apply(Bmal_WT_37C,1,function(y) !any(y==0)),]
	KO_37C <- Bmal_KO_37C[apply(Bmal_KO_37C,1,function(y) !any(y==0)),]

	WT_37C_log <- log2(WT_37C)
	KO_37C_log <- log2(KO_37C)

	WT37 <- runRAIN(dat = WT_37C_log, deltat = 2)
	KO37 <- runRAIN(dat = KO_37C_log, deltat = 2)

	WT37$p.adj <- round(p.adjust(WT37$pVal,method='BH'), digit = 6)
	KO37$p.adj <- round(p.adjust(KO37$pVal,method='BH'), digit = 6)

	assign(paste("WT37_Individual0","_results",sep=""),WT37)
	assign(paste("KO37_Individual0","_results",sep=""),KO37)

	rm(WT_37C,KO_37C,WT_37C_log,KO_37C_log,WT37,KO37,runRAIN)


# ######################################################################################################
# wrapper function for running RAIN with defaults settings as in Ray et al. 
# ######################################################################################################

runRAIN <- function(dat, deltat, period.delta = 0) {
  rain(x = t(dat), deltat = deltat, period = 24, period.delta = period.delta, adjp.method = "ABH", method = "independent", peak.border = c(0.1,0.9)) 
}

	WT_37C <- Bmal_WT_37C[apply(Bmal_WT_37C,1,function(y) !any(y==0)),]
	KO_37C <- Bmal_KO_37C[apply(Bmal_KO_37C,1,function(y) !any(y==0)),]

	WT_37C_log <- log2(WT_37C)
	KO_37C_log <- log2(KO_37C)

	WT37 <- runRAIN(dat = WT_37C_log, deltat = 2)
	KO37 <- runRAIN(dat = KO_37C_log, deltat = 2)

	WT37$p.adj <- round(p.adjust(WT37$pVal,method='BH'), digit = 6)
	KO37$p.adj <- round(p.adjust(KO37$pVal,method='BH'), digit = 6)

	assign(paste("WT37_RayEtAl","_results",sep=""),WT37)
	assign(paste("KO37_RayEtAl","_results",sep=""),KO37)

	rm(WT_37C,KO_37C,WT_37C_log,KO_37C_log,WT37,KO37,runRAIN)

# ######################################################################################################
# save all the output dataframes
# ######################################################################################################
rm(AM_PM_TempComp_FPKM)
current_session <- sessionInfo()
save.image(file = "output/FPKM_Threshlod_Rainmethods.RData")  

# ######################################################################################################
# Figure 1B
# ######################################################################################################

KO37_cyclic_genes_list <- list( "KO37 No threshold" = rownames(KO37_fpkm_none_results[which(KO37_fpkm_none_results$p.adj<0.1),]),
								"KO37 Mean FMKM > 0" = rownames(KO37_fpkm_0_results[which(KO37_fpkm_0_results$p.adj<0.1),]), 
								"KO37 Mean FPKM > 0.05" = rownames(KO37_fpkm_0.05_results[which(KO37_fpkm_0.05_results$p.adj<0.1),]), 
								"KO37 Mean FPKM > 0.1" = rownames(KO37_fpkm_0.1_results[which(KO37_fpkm_0.1_results$p.adj<0.1),]), 
								"KO37 Mean FPKM > 0.5" = rownames(KO37_fpkm_0.5_results[which(KO37_fpkm_0.5_results$p.adj<0.1),]), 
								"KO37 Mean FPKM > 1" = rownames(KO37_fpkm_1_results[which(KO37_fpkm_1_results$p.adj<0.1),]), 
								"KO37 Mean FPKM > 1.5" = rownames(KO37_fpkm_1.5_results[which(KO37_fpkm_1.5_results$p.adj<0.1),]), 
								"KO37 Mean FPKM > 2" = rownames(KO37_fpkm_2_results[which(KO37_fpkm_2_results$p.adj<0.1),]), 
								"KO37 Mean FPKM > 3" = rownames(KO37_fpkm_3_results[which(KO37_fpkm_3_results$p.adj<0.1),]), 
								"KO37 Individual FPKM > 0" = rownames(KO37_Individual0_results[which(KO37_Individual0_results$p.adj<0.1),]), 
								"KO37 RayEtAl" = rownames(KO37_RayEtAl_results[which(KO37_RayEtAl_results$p.adj<0.1),]) )

summary(KO37_cyclic_genes_list)

KO37_cyclic_genes_numbers <- as.data.frame(unclass(summary(KO37_cyclic_genes_list)))
KO37_cyclic_genes_numbers[,1] <- as.numeric(as.character(KO37_cyclic_genes_numbers[,1]))
KO37_cyclic_genes_numbers$Names <- rownames(KO37_cyclic_genes_numbers)

plot1 <- ggplot(KO37_cyclic_genes_numbers, aes(x=Names, y=Length)) + geom_bar(stat="identity") + 
			scale_y_continuous(name ="Number of Genes") + scale_x_discrete(limits = rownames(KO37_cyclic_genes_numbers)) + 
			theme(axis.text.x = element_text(angle=90) )


pdf("output/Figure1B_FPKM_Threshlod_Rainmethods.pdf")
print(plot1)
dev.off()


# ######################################################################################################






