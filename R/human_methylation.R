#
# Josep Carreras Leukaemia Research Institute (IJC)                            
# Authors: Raúl Sanz, Angelika Merkel | Bioinformatics Unit                    
# Project: PIK3Ca mutation study                                                     
# Collaborators: Mariona Graupera, Sandra del Castillo                           
#
# Description: Workflow to analyze human methylation array data from Illumina's Infinium MethylationEPIC BeadChip.
#

# REQUIRED PACKAGES

library(minfi)
library(readxl)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggrepel)
library(RColorBrewer)
library(ggfortify)
library(maxprobes)
library(limma)
library(DMRcate)

# LOADING THE DATA

# reading all the methylation array files (.idat) from a folder
rgSet <- read.metharray.exp("data/human/")

# reading the metadata file
metadata <- as.data.frame(read_excel("data/human/PIK3CA_samples_SC.xlsx"))
names(metadata)[1:10] <- c("Sample", "Organism", "Tissue", "Type", "Condition", "Preservation", "DNA_quantity", 
                           "EPIC_ID", "EPIC_position", "EPIC_barcode")

## renaming some columns to adapt the metadata
metadata$CellType <- paste(metadata$CellType1, metadata$CellType2, sep="_" )
metadata$CellType <- gsub("thelial","", gsub("_NA","", as.character(metadata$CellType)))
metadata$Condition <- gsub("pik3ca ","", as.character(metadata$Condition))
metadata$Type <- gsub("Vascular malformation", "Mutation", metadata$Type)
metadata <- metadata[,-c(2, 3, 6, 8, 9, 12, 13, 14)]

# QUALITY CONTROL

# Generating the quality control report for the raw data 
qcReport(rgSet, pdf="results/human/qcReport.pdf", sampGroups=metadata$Condition, sampNames=metadata$Sample)

# detection p-values
p_values <- detectionP(rgSet, type = "m+u")

# mean detection p-values
mean_p <- data.frame(p_values=colMeans(p_values), Sample=metadata$Sample)

## plotting the resulting mean detection p-values
ggplot(mean_p, aes(x=Sample, y=p_values, fill=p_values))+
  geom_col()+
  theme_bw()+
  ggtitle("Mean detection p-value per sample")+
  theme(legend.position='none')+
  geom_hline(yintercept=0.01, linetype="dashed", color="red")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))

# PREPROCESSING

mSet <- preprocessNoob(rgSet)

# removing low-quality probes by detection p-value
keep <- rowSums(p_values<0.01)==ncol(mSet)
mSet <- mSet[keep,]

# removing probes with known SNPs
gmSet <- dropLociWithSnps(mapToGenome(mSet))

# removing cross reactive probes
xreactive_probes <- xreactive_probes(array_type="EPIC")
keep <- !(featureNames(gmSet) %in% xreactive_probes)
gmSet <- gmSet[keep,]

# removing sex chromosomes
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(gmSet) %in% ann$Name[ann$chr %in% c("chrX","chrY")])
gmSet <- gmSet[keep,]

# obtaining the beta values
beta_values <- getBeta(gmSet)
colnames(beta_values) <- metadata$Sample
head(beta_values)

# PRINCIPAL COMPONENT ANALYSIS

# selecting the top 100 most variable CpG's to reduce the number of inputs
sdv <- apply(beta_values, 1, sd)
keep <- names(head(sort(sdv,decreasing=T), 100))
beta_top100 <- beta_values[keep,]

# PCA on the top 100 sites
pca_res <- prcomp(t(beta_top100), scale=T, center=T)

## plotting PC1 and PC2 by condition
autoplot(pca_res, x=1, y=2, data=metadata, colour="Condition")+
  geom_text_repel(aes(label=Sample, color=Condition),hjust=-0.2, vjust=0, show.legend=F, size=3.5)+
  labs(colour="Mutation")+
  xlim(c(-0.5,0.3))+
  theme_bw()

## plotting PC1 and PC2 by cell type
autoplot(pca_res, x=1, y=2,data=metadata, colour="CellType", shape = "Type")+
  geom_text_repel(aes(label=Sample, color=CellType),hjust=-0.2, vjust=0, show.legend=F, size=3.5)+
  scale_color_brewer(palette = "Set1")+ 
  labs(colour="Cell Type") + xlim(c(-0.5,0.3))+
  theme_bw()

# DIFFERENTIAL METHYLATION ANALYSIS

# differentially methylated positions (DMPs)
## building the design matrix
mutation <- factor(metadata$Condition)
designMat <- model.matrix(~0+mutation, data=metadata)
colnames(designMat) <- levels(mutation)

## building the contrast matrix
contMat <- makeContrasts(Control-E542K, 
                         Control-E545K, 
                         Control-H1047R, 
                         E542K-E545K, 
                         E542K-H1047R, 
                         H1047R-E545K, levels=designMat)

## fitting linear model with limma
fit <- lmFit(beta_values, designMat)
fit2 <- contrasts.fit(fit, contMat)
fit2 <- eBayes(fit2)

DMP_summary <- summary(decideTests(fit2, p.value=0.01))
rownames(DMP_summary) <- c("Hyper", "NotSig", "Hypo")
colnames(DMP_summary) <- gsub(" - ","_", colnames(DMP_summary))

## annotating the DMPs
contrasts <- colnames(contMat)
### matching beta values with annotation by the probe names
annEPICSub <- ann[match(rownames(beta_values),ann$Name), c(1:4,12:19,22:ncol(ann))]

DMP.list <- list() # list to store the DMPs
cg <- list() # list to store the names of significant probes

for (i in 1:length(contrasts)){
  # extract the significant DMPs and annotate
  DMP <- topTable(fit2, num=Inf, coef=i, genelist=annEPICSub, p.value=0.01)
  DMP$Contrast <- contrasts[i]
  # saving the results
  DMP.list[[i]] <- DMP
  cg[[i]] <- row.names(DMP)
  # write.csv(DMP_1, file=paste0("results/DMP/DMP_", contrasts[i], ".csv"))
}

### merging all the DMP's
DMP_ann <- do.call(rbind, DMP.list)

### classifying DMP's according to its change in methylation
DMP_ann$Type <- "Hypermethylated"
DMP_ann$Type[which(DMP_ann$logFC > 0)] <- "Hypomethylated"

### gene feature annotation
DMP_ann$UCSC_RefGene_Group[which(DMP_ann$UCSC_RefGene_Group == "")] <- "."
DMP_ann$UCSC_RefGene_Group_short <- unlist(lapply(strsplit(DMP_ann$UCSC_RefGene_Group, ";"),'[[', 1))

## barplot relative to CpG Islands
DMP_annCGI <- DMP_ann[, c("Contrast", "Relation_to_Island", "Type")]

ggplot(DMP_annCGI, aes(Contrast, fill=Relation_to_Island))+
  facet_wrap(.~Type, scales="free_x")+
  geom_histogram(stat="count", width=0.75)+
  theme_bw()+ 
  theme(axis.text.x=element_text(angle=30, hjust=1, size=8))+
  ylab("Number of DMPss")+ 
  xlab("")

## barplot relative to genetic features
DMP_Gene_Group <- DMP_ann[,c("Contrast","UCSC_RefGene_Group_short", "Type")]

ggplot(DMP_Gene_Group, aes(Contrast, fill=UCSC_RefGene_Group_short))+
  facet_wrap(.~Type, scales="free_x")+
  geom_histogram(stat="count", width=0.75)+ 
  theme_bw()+ 
  scale_fill_brewer(palette="Set2")+
  theme(axis.text.x=element_text(angle=30, hjust=1, size=8))+
  ylab("Number of DMPs")+
  xlab("")+ 
  labs(fill="UCSC_RefGene")

## heatmap
DMP_beta <- data.frame(beta_values[unlist(cg), ])
colnames(DMP_beta) <- paste0(metadata$Sample,"_",metadata$Condition)

colors <- colorRampPalette(c("royalblue", "white", "red"))(n=100) # range of colors
heatmap.2(as.matrix(t(unique(DMP_beta))), trace="none", density.inf="none", margins=c(7,10), col=colors, cexRow = 1, lwid = c(5,15), lhei = c(5,15))

# differentially methylated regions (DMRs)
## finding and annotating the DMPs
DMR.list <- list()
for(i in 1:length(contrasts)){
  # annotating the CpG's for each contrast + finding significant CpG's
  myAnnotation <- cpg.annotate(object=beta_values, datatype="array", what="Beta", analysis.type="differential", design=designMat, contrasts=T, 
                               cont.matrix=contMat, coef=gsub("_", " - ", as.character(contrasts[i])), arraytype="EPIC")
  
  if (sum(myAnnotation@ranges@elementMetadata@listData$is.sig)!=0){ # if there are significant CpG's
    # test for DMRs
    DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
    # extract genomic ranges
    results.ranges <- extractRanges(DMRs)
    # saving the results
    DMR.list[[i]] <- data.frame(Contrast=contrasts[i], results.ranges)
    # write.csv(results.ranges, paste0("results/DMR/DMR_", contrasts[i], ".csv"))
  }
}

## plotting the DMRs
pal <- brewer.pal(8, "Dark2") # palette

### colour by condition
groups <- pal[1:length(unique(metadata$Condition))]
names(groups) <- levels(factor(metadata$Condition))
cols <- groups[as.character(factor(metadata$Condition))]

### is needed to convert the data frame into a genomic ranges object first
DMR.GR <- makeGRangesFromDataFrame(DMR.list[3])

DMR.plot(ranges=DMR.GR, dmr=10, CpGs=beta_values, phen.col=cols, what="Beta", arraytype="EPIC", genome="hg19")