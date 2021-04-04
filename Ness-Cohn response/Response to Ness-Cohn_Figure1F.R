# ##########################################################################################################################
# Response to Technical Comment by Ness-Cohn et al. (2021). Code for:
# Figure 1F
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

library(annotate)
library(data.table)
library(foreach)
library(GEOquery)
library(ggplot2)
library(knitr)
library(limma)
library(limorhyde)
library(org.Mm.eg.db)
library(eulerr)

source(system.file('extdata', 'vignette_functions.R', package = 'limorhyde'))

# ######################################################################################################
# Download the data from GEO and unpack it to data directory
# ######################################################################################################

dir.create("output", showWarnings=FALSE)
dir.create("data", showWarnings=FALSE)
setwd("data")

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE134333&format=file&file=GSE134333%5FDataset%5FAM%5FPM%5FTempComp%5FFPKM%5Fvalues%2Etxt%2Egz", "GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt.gz")
gunzip("GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt.gz", overwrite = T)

download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134333/matrix/GSE134333_series_matrix.txt.gz", "GSE134333_series_matrix.txt.gz")
gunzip("GSE134333_series_matrix.txt.gz", overwrite = T)

# ######################################################################################################

period = 24
qvalRhyCutoff = 0.1
qvalDrCutoff = 0.1

eset = getGEO(filename = "GSE134333_series_matrix.txt")

sm_all = data.table(pData(phenoData(eset)))
sm_all = sm_all[, .(title, sample = geo_accession, genotype = `genotype/variation:ch1`)]
sm_all[, time := as.numeric(sub('[[:space:]].*', '', sub('.*[[:space:]]CT', '', title)))]
sm_all[, genocond := factor(sub('[[:space:]]CT..[[:space:]]', '_', title))]
sm_all <- sm_all[1:68, ]
sm_all[, title := sub('[[:space:]]', '.', sub('[[:punct:]][[:punct:]][[:punct:]]', '...', sub('[[:space:]]', '.', title)))]
sm_all[, cond := factor(genocond, c('WT_AM', 'Bmal1-/-_AM', 'WT_PM', 'Bmal1-/-_PM'), c('WildTypeAM', 'KnockOutAM','WildTypePM', 'KnockOutPM'))]

setorderv(sm_all, c('cond', 'time'))

sm_all = cbind(sm_all, limorhyde(sm_all$time, 'time_'))

AM_PM_TempComp_FPKM <- read.delim("GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt", quote = "")
AM_PM_data <- AM_PM_TempComp_FPKM[, grep(".*.M.FPKM", colnames(AM_PM_TempComp_FPKM))]
rownames(AM_PM_data) <- AM_PM_TempComp_FPKM$gene_id
colnames(AM_PM_data) <- sub('.FPKM', '', colnames(AM_PM_data))
AM_PM_data <- AM_PM_data[c(sm_all$title)]
colnames(AM_PM_data) <- sm_all$sample

## removes rows  if all colums have zeros
AM_PM_data <- AM_PM_data[apply(AM_PM_data,1,function(z) any(z!=0)),]

emat_all <- as.matrix(log2(AM_PM_data+1))

# # # WT_AM 1:17
# # # KO_AM 18:34
# # # WT_PM 35:51
# # # KO_PM 52:68

sm_WT <- sm_all[c(1:17,35:51),]
sm_WT[, cond := factor(genocond, c('WT_AM', 'WT_PM'), c('WildTypeAM', 'WildTypePM'))]

sm_KO <- sm_all[c(18:34,52:68),]
sm_KO[, cond := factor(genocond, c('Bmal1-/-_AM', 'Bmal1-/-_PM'), c('KnockOutAM', 'KnockOutPM'))]

emat_WT <- emat_all[ ,c(1:17,35:51)]
emat_KO <- emat_all[ ,c(18:34,52:68)]

emat_list <- c("emat_all", "emat_WT", "emat_KO")

for(myemat in emat_list){

	emat <- get(myemat)
	if(myemat=="emat_all") {sm <- sm_all}
	if(myemat=="emat_WT") {sm <- sm_WT}
	if(myemat=="emat_KO") {sm <- sm_KO}

	####### Identify rhythmic genes

	rhyLimma = foreach(condNow = unique(sm$cond), .combine = rbind) %do% {
	  design = model.matrix(~ time_cos + time_sin, data = sm[cond == condNow])
	  fit = lmFit(emat[, sm$cond == condNow], design)
	  fit = eBayes(fit, trend = TRUE)
	  rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
	  setnames(rhyNow, 'rn', 'geneId')
	  rhyNow[, cond := condNow]
	}

	rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = geneId]
	rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
	setorderv(rhyLimmaSummary, 'adj.P.Val')

	kable(rhyLimmaSummary[1:5, ])

	####### Identify differentially rhythmic genes

	design = model.matrix(~ cond * (time_cos + time_sin), data = sm)

	fit = lmFit(emat, design)
	fit = eBayes(fit, trend = TRUE)
	drLimma = data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
	setnames(drLimma, 'rn', 'geneId')

	drLimma = drLimma[geneId %in% rhyLimmaSummary[adj.P.Val <= qvalRhyCutoff]$geneId]
	drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
	setorderv(drLimma, 'adj.P.Val')

	kable(drLimma[1:5, ])

	rm(design, fit)

	####### Identify differentially expressed genes

	design = model.matrix(~ cond + time_cos + time_sin, data = sm)

	fit = lmFit(emat, design)
	fit = eBayes(fit, trend = TRUE)
	deLimma = data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = TRUE)
	setnames(deLimma, 'rn', 'geneId')

	deLimma = deLimma[!(geneId %in% drLimma[adj.P.Val <= qvalDrCutoff]$geneId)]
	deLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
	setorderv(deLimma, 'adj.P.Val')

	kable(deLimma[1:5, ])


	assign(paste(myemat,"_","rhyLimma",sep=""),rhyLimma)
	assign(paste(myemat,"_","rhyLimmaSummary",sep=""),rhyLimmaSummary)
	assign(paste(myemat,"_","drLimma",sep=""),drLimma)
	assign(paste(myemat,"_","deLimma",sep=""),deLimma)
			
	rm(myemat, emat, sm, design, fit, rhyLimma, rhyLimmaSummary, drLimma, deLimma)
}

setwd("..")

current_session <- sessionInfo()
save(current_session, AM_PM_data, emat_all, sm_all, sm_KO, sm_WT,
	emat_all_deLimma, emat_all_drLimma, emat_all_rhyLimma, emat_all_rhyLimmaSummary,
	emat_KO, emat_KO_deLimma, emat_KO_drLimma, emat_KO_rhyLimma, emat_KO_rhyLimmaSummary,
	emat_WT, emat_WT_deLimma, emat_WT_drLimma, emat_WT_rhyLimma, emat_WT_rhyLimmaSummary,
    file = "output/limorhyde_results.RData" )

# load('~/Downloads/output/limorhyde_results.RData')

###### differentially rhythmic analysis

drLimma = emat_WT_drLimma
emat = emat_WT
sm = sm_WT

geneIdsNow = c( drLimma$geneId[1:5] )
df = data.table(t(emat[geneIdsNow, ]))
df[, sample := colnames(emat[geneIdsNow, ])]

df = merge(df, sm[, .(sample, cond, time)], by = 'sample')
df = melt(df, measure.vars = geneIdsNow, variable.name = 'geneID',
          value.name = 'expr')
levels(df$geneID) = geneIdsNow
WT_drLimma_significant <- drLimma[which(drLimma$adj.P.Val<0.1),]
nrow(WT_drLimma_significant)
rm(drLimma,emat,sm)


drLimma = emat_KO_drLimma
emat = emat_KO
sm = sm_KO

geneIdsNow = c( drLimma$geneId[1:5] )
df = data.table(t(emat[geneIdsNow, ]))
df[, sample := colnames(emat[geneIdsNow, ])]

df = merge(df, sm[, .(sample, cond, time)], by = 'sample')
df = melt(df, measure.vars = geneIdsNow, variable.name = 'geneID',
          value.name = 'expr')
levels(df$geneID) = geneIdsNow
KO_drLimma_significant <- drLimma[which(drLimma$adj.P.Val<0.1),]
nrow(KO_drLimma_significant)
rm(drLimma,emat,sm)




WT_KO_AM_PM_drGenes <- list(  WT_drGenes = WT_drLimma_significant$geneId,
  KO_drGenes = KO_drLimma_significant$geneId )
summary(WT_KO_AM_PM_drGenes)


# ######################################################################################################

