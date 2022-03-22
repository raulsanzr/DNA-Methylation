#
# Josep Carreras Leukaemia Research Institute (IJC)                            
# Authors: Raúl Sanz, Angelika Merkel | Bioinformatics Unit                    
# Project: PIK3Ca mutation study                                                     
# Collaborators: Mariona Graupera, Damiana Álvarez, Sandra del Castillo                           
#
# Description: Workflow to analyze mouse methylation array data from Illumina's Infinium Mouse Methylation BeadChip.
#
# Note: minfi has not already implemented functions to analyze mouse methylation data. For this reason, in this script, the procedure aims to
#       achieve the same results in a different way.
#
# * files produced with the script mouse_annotation.R
#

# REQUIRED PACKAGES

library(readxl)
library(ENmix)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(tidyverse)
library(limma)
library(bumphunter)

# LOADING THE DATA

# reading the metadata file
metadata <- as.data.frame(read_excel("data/mouse/PIK3Ca_mouse_samples.xlsx"))
names(metadata) <- c("Sample", "Organism", "Tissue", "Type", "Condition", "Preservation", "DNA_quantity", 
                     "EPIC_ID", "EPIC_position", "EPIC_barcode")
## renaming some columns to adapt the metadata
metadata$Cond <- gsub(" ", "_", as.character(metadata$Condition))
metadata$Observation <- substr(metadata$Sample, 0, 4)
metadata <- metadata %>% separate(Condition, c("Group", "Time"), " ")

# reading all the methylation array files (.idat) from a folder
rgSet <- readidat(path="data/mouse/", manifestfile="data/Infinium_Mouse_Methylation_v1.0_A1_GS_Manifest.csv")
colnames(rgSet) <- metadata$Sample[which(colnames(rgSet) %in% metadata$EPIC_barcode)]

# reading the manifest file *
manifest <- read.csv("data/manifest_mouse.csv")

# PREPROCESSING

# applying noob normalization and obtaining the beta values
beta_values <- mpreprocess(rgSet, qnorm=T, bgParaEst="oob", nCores=6)

# removing sex chromosomes
cg_autosomes <- manifest$Name_long[which(manifest$CHR != "X" & manifest$CHR != "Y")]
beta_values <- beta_values[which(row.names(beta_values) %in% cg_autosomes),]

# removing probes with known SNPs
cg_probes <-  manifest$Name_long[which(manifest$Probe_Type == "cg" | manifest$Probe_Type == "ch")]
beta_values <- beta_values[which(row.names(beta_values) %in% cg_probes),]

# removing outliers (in case there are)
metadata <- metadata[metadata$Observation != "OB42",]
beta_values <- beta_values[, grep("OB42", colnames(beta_values), invert=T)]

# PRINCIPAL COMPONENT ANALYSIS

# selecting the top 100 most variable CpGs to reduce the number of inputs
sdv  <- apply(beta_values, 1, sd)
keep <- names(head(sort(sdv, decreasing=T), 100))
beta_top100 <- beta_values[keep, ]

# PCA on the top 100 sites
pca_res <- prcomp(t(beta_top100), scale=T, center=T)

## plotting PC1 and PC2 by group
autoplot(pca_res, data=metadata, colour="Group", shape="Time") +
  geom_text_repel(aes(label=Sample, color=Group), hjust=-0.3, vjust=0, show.legend=F, size=3.5) +
  labs(colour="Condition")+
  theme_bw() + 
  xlim(c(-0.5, 0.4)) + ylim(c(-0.5, 0.6))

# DIFFERENTIAL METHYLATION ANALYSIS

# differentially methylated positions (DMPs)
## building the design matrix
designMat <- model.matrix(~0 + Observation, data=metadata)
designMat <- cbind(designMat, X=metadata$Group == "Vehicle" )
designMat <- cbind(designMat, Y=metadata$Time == "96h")

## building the contrast matrix
contMat <- makeContrasts(X-Y, levels=colnames(designMat))

## fitting linear model with limma
fit <- lmFit(beta_values, designMat)
fit2 <- contrasts.fit(fit, contMat)
fit2 <- eBayes(fit2)

DMP_summary <- summary(decideTests(fit2, p.value=0.01))
rownames(DMP_summary) <- c("Hyper", "NotSig", "Hypo")
colnames(DMP_summary) <- gsub(" - ","_", colnames(DMP_summary))

## extracting the DMPs
### getting the significant positions
CpG_sig <- topTable(fit2, p.value=0.01, number=Inf)
### adding the beta values
DMPs <- data.frame(beta_values[which(rownames(beta_values) %in% rownames(CpG_sig)),])
### adding the p-values and the adjusted p-values
DMPs$adj.pval <- as.numeric(CpG_sig$adj.P.Val[which(row.names(CpG_sig) %in% row.names(DMPs))])
DMPs$pval <- as.numeric(CpG_sig$P.Value[which(row.names(CpG_sig) %in% row.names(DMPs))])

### difference in methylation between controls and mutation samples at the two time points
DMPs$mean_diff_16h <- rowMeans(DMPs[c("OB39_et_16h", "OB41_et_16h")]) - rowMeans(DMPs[c("OB39_4OH_16h", "OB41_4OH_16h")])
DMPs$mean_diff_96h <- rowMeans(DMPs[c("OB39_et_96h", "OB41_et_96h")]) - rowMeans(DMPs[c("OB39_4OH_96h", "OB41_4OH_96h")])

### getting the positions with higher difference
DMPs[order(DMPs$mean_diff_16h),]

## annotating the DMPs
annNCBIshort <- read.csv("data/annNCBIshort_mouse.csv") # NCBI annotation file *

### joining the DMPs with the manifest
DMPs$Name_long  <- row.names(DMPs)
DMPs_ann <- left_join(DMPs, manifest, by="Name_long")

### Adding the NCBI annotation
DMPs_ann$name <- unlist(lapply(strsplit(DMPs_ann$Name_long,"_"), "[[", 1))
DMPs_ann_NCBI <- left_join(DMPs_ann, annNCBIshort, by="name")

## selecting the DMPs above a threshold (mean diff >= 0.2) in the difference of methylation
DMPs_ann_short <- DMPs_ann_NCBI[(abs(DMPs_ann_NCBI$mean_diff_16h)>= 0.2) | (abs(DMPs_ann_NCBI$mean_diff_96h)>= 0.2), c(1:13, 16, 17, 20:22, 26)]

## saving the results
write.csv(DMPs_ann_NCBI, file="results/DMPs/DMP_ann.csv")

# differentially methylated regions (DMRs)
## finding and annotating the DMPs
DMR_Control_Mutation <- bumphunter(beta_values, design=designMat, chr=manifest$CHR, pos=manifest$MAPINFO, cutoff=0.5, coef=3)

## saving the results
write.csv(DMR_Control_Mutation$table, file="results/DMR/DMR_Control_Mutation.csv")

# Note: Those DMRs found can be plotted using the plot_DMR.R script.