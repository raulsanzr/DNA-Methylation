#
# Author: Raul Sanz               
# Description: Pipeline to analyze human methylation array data from Illumina's Infinium MethylationEPIC array.
#

# DEPENDENCIES

library(minfi) # BiocManager::install("minfi")
library(limma) # BiocManager::install("limma")
library(DMRcate) # BiocManager::install("DMRcate")
library(maxprobes) # remotes::install_github("markgene/maxprobes")
library(readxl) # install.packages("readxl")
library(dplyr) # install.packages("dplyr")
library(ggplot2) # install.packages("ggplot2")
library(ggrepel) # install.packages("ggrepel")
library(ggfortify) # install.packages("ggfortify")
library(gplots) # install.packages("gplots")
library(RColorBrewer) # install.packages("RColorBrewer")

# LOADING THE DATA

# reading all the DNA methylation array files (idat) from a folder
rgSet <- read.metharray.exp("data/")

# reading the metadata
metadata <- as.data.frame(read_excel("data/metadata.xlsx"))
names(metadata)[1:10] <- c("Sample", "Organism", "Tissue", "Type", "Condition", "Preservation", "DNA_quantity", 
                           "EPIC_ID", "EPIC_position", "EPIC_barcode")

## formatting the metadata
metadata$CellType <- paste(metadata$CellType1, metadata$CellType2, sep="_" )
metadata$CellType <- gsub("thelial","", gsub("_NA","", as.character(metadata$CellType)))
metadata$Condition <- gsub("pik3ca ","", as.character(metadata$Condition))
metadata$Type <- gsub("Vascular malformation", "Mutation", metadata$Type)
metadata$Group <- paste0(metadata$Type,"_", metadata$CellType)

# QUALITY CONTROL

# generating the QC report
qcReport(rgSet, pdf="results/qcReport.pdf", sampGroups=metadata$Condition, sampNames=metadata$Sample)

# detection p-values
p_values <- detectionP(rgSet, type = "m+u")

## mean detection p-values
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

# ssNoob normalization
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

# removing probes mapping sex chromosomes
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(gmSet) %in% ann$Name[ann$chr %in% c("chrX","chrY")])
gmSet <- gmSet[keep,]

# obtaining the beta values
beta_values <- getBeta(gmSet)
colnames(beta_values) <- metadata$Sample

# PRINCIPAL COMPONENT ANALYSIS

# selecting the top 100 most variable CpG sites
sdv <- apply(beta_values, 1, sd)
top100 <- names(head(sort(sdv,decreasing=T), 100))
beta_top100 <- beta_values[top100,]

# PCA on the top 100 sites
pca_res <- prcomp(t(beta_top100), scale=T, center=T)

## plotting PC1 and PC2 by condition
autoplot(pca_res, x=1, y=2, data=metadata, colour="Condition")+
  geom_text_repel(aes(label=Sample, color=Condition),hjust=-0.2, vjust=0, show.legend=F, size=3.5)+
  labs(colour="Mutation")+
  xlim(c(-0.5,0.3))+
  theme_bw()

# DIFFERENTIAL METHYLATION ANALYSIS

## DMPs

# design matrix
designMat <- model.matrix(~0+Group, data=metadata)
colnames(designMat) <- gsub("group","",colnames(designMat))

# contrast matrix
contMat <- makeContrasts(Normal_blood-Normal_lymphatic,
                         Mutation_blood-Mutation_lymphatic,
                         levels=designMat)

# fitting linear model
fit <- lmFit(beta_values, designMat)
fit2 <- contrasts.fit(fit, contMat)
fit2 <- eBayes(fit2)

DMP_summary <- summary(decideTests(fit2, p.value=0.01))
rownames(DMP_summary) <- c("Hyper", "NotSig", "Hypo")
colnames(DMP_summary) <- gsub(" - ","_", colnames(DMP_summary))

# annotating the DMPs
## matching beta values with annotation by the probe names
annEPICSub <- ann[match(rownames(beta_values),ann$Name), c(1:4,12:19,22:ncol(ann))]

DMP.list <- data.frame() # df to store the DMPs
cg <- list() # list to store the names of significant probes

contrasts <- colnames(contMat)
for (i in 1:length(contrasts)){ # for every comparison...
  # extract the significant DMPs and annotate
  DMP <- topTable(fit2, num=Inf, coef=i, genelist=annEPICSub, p.value=0.01)
  DMP$Contrast <- contrasts[i]
  DMP.list <- rbind(DMP.list, DMP)
  cg[[i]] <- row.names(DMP)
}

## saving the results
write.csv(DMP.list, "results/DMP_list.csv")

## merging all the DMPs
DMP_ann <- DMP.list

## classifying DMPs according to its change in methylation
DMP_ann$Type <- "Hypermethylated"
DMP_ann$Type[which(DMP_ann$logFC > 0)] <- "Hypomethylated"

## gene feature annotation
DMP_ann$UCSC_RefGene_Group[which(DMP_ann$UCSC_RefGene_Group == "")] <- "."
DMP_ann$UCSC_RefGene_Group_short <- unlist(lapply(strsplit(DMP_ann$UCSC_RefGene_Group, ";"),'[[', 1))

## DMPs relative to the nearest CGI
DMP_annCGI <- DMP_ann[, c("Contrast", "Relation_to_Island", "Type")]

ggplot(DMP_annCGI, aes(Contrast, fill=Relation_to_Island))+
  facet_wrap(.~Type, scales="free_x")+
  geom_histogram(stat="count", width=0.75)+
  theme_bw()+ 
  theme(axis.text.x=element_text(angle=30, hjust=1, size=8))+
  ylab("Number of DMPss")+ 
  xlab("")

## DMPs relative to gene elements
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

## DMPs heatmap
DMP_beta <- data.frame(beta_values[unlist(cg), ])
colnames(DMP_beta) <- paste0(metadata$Sample,"_",metadata$Condition)
colors <- colorRampPalette(c("royalblue", "white", "red"))(n=100) # range of colors
heatmap.2(as.matrix(t(unique(DMP_beta))), trace="none", density.inf="none", margins=c(7,10), col=colors, cexRow = 1, lwid = c(5,15), lhei = c(5,15))

## DMRs

DMR.list <- data.frame() # df to store the DMRs

for(i in 1:length(contrasts)){ # for every comparison...
  # annotating and finding significant CpG sites
  myAnnotation <- cpg.annotate(object=beta_values, datatype="array", what="Beta", analysis.type="differential", design=designMat, contrasts=T, 
                               cont.matrix=contMat, coef=gsub("_", " - ", as.character(contrasts[i])), arraytype="EPIC")
  
  if(sum(myAnnotation@ranges@elementMetadata@listData$is.sig)!=0){ # if there are significant sites...
    # search for DMRs
    DMR_raw <- dmrcate(myAnnotation, lambda=1000, C=2)
    # extract the genomic ranges
    results.ranges <- extractRanges(DMR_raw)
    DMR <- data.frame(Contrast=contrasts[i], results.ranges)
    DMR.list <- rbind(DMR.list, DMR)
  }
}

## saving the results
write.csv(DMR.list, "results/DMR_list.csv")
