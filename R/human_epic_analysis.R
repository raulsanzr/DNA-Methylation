# Workflow to analyze human methylation array data from Illumina's Infinium MethylationEPIC array. #

# 0.LIBRARIES
library(minfi) # BiocManager::install("minfi")
library(limma) # BiocManager::install("limma")
library(DMRcate) # BiocManager::install("DMRcate")
library(maxprobes) # remotes::install_github("markgene/maxprobes")
library(readxl) # install.packages("readxl")
library(dplyr) # install.packages("dplyr")

# 1.LOADING THE DATA
# 1.1.Read all the idat from a folder
rgSet <- read.metharray.exp("data/")

# 1.2.Read the metadata
metadata <- as.data.frame(read_excel("metadata.xlsx"))

# 2.QUALITY CONTROL
# 2.1.QC report
qcReport(rgSet, pdf="results/qcReport.pdf", sampGroups=metadata$Condition, sampNames=metadata$Sample)

# 2.2.Detection p-values
p_values <- detectionP(rgSet, type = "m+u")

# 3.PREPROCESSING
# 3.1.Normalization (Check the minfi documentation)
mSet <- preprocessNoob(rgSet) 

# 3.2.Remove low-quality probes by detection p-value
keep <- rowSums(p_values<0.01)==ncol(mSet)
mSet <- mSet[keep,]

# 3.3.Remove probes with known SNPs
gmSet <- dropLociWithSnps(mapToGenome(mSet))

# 3.4.Remove cross reactive probes
xreactive_probes <- xreactive_probes(array_type="EPIC")
keep <- !(featureNames(gmSet) %in% xreactive_probes)
gmSet <- gmSet[keep,]

# 3.5.Remove probes mapping sex chromosomes
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(gmSet) %in% ann$Name[ann$chr %in% c("chrX","chrY")])
gmSet <- gmSet[keep,]

# 3.6.Beta values
beta_values <- getBeta(gmSet)
colnames(beta_values) <- metadata$Sample

# 4.DIFFERENTIAL METHYLATION ANALYSIS
# 4.1.Differentially methylated positions
# 4.1.1.Design matrix
designMat <- model.matrix(~0+Group, data=metadata)
colnames(designMat) <- gsub("group","",colnames(designMat))

# 4.1.2.Contrast matrix
contMat <- makeContrasts(Normal_blood-Normal_lymphatic,
                         Mutation_blood-Mutation_lymphatic,
                         levels=designMat)

# 4.1.3.Fit the linear model
fit <- lmFit(beta_values, designMat)
fit2 <- contrasts.fit(fit, contMat)
fit2 <- eBayes(fit2)

dmps <- summary(decideTests(fit2, p.value=0.01))
rownames(dmps) <- c("Hyper", "NotSig", "Hypo")
colnames(dmps) <- gsub(" - ","_", colnames(dmps))

# 4.1.4. Annotate the DMPs
## Match beta values with their annotation by the probe names
annEPICSub <- ann[match(rownames(beta_values),ann$Name), c(1:4,12:19,22:ncol(ann))]

dmplist <- data.frame() # df to store the DMPs
cg <- list() # list to store the names of significant probes

contrasts <- colnames(contMat)
for (i in 1:length(contrasts)){ # for every comparison...
  # extract the significant DMPs and annotate
  dmp <- topTable(fit2, num=Inf, coef=i, genelist=annEPICSub, p.value=0.01)
  dmp$Contrast <- contrasts[i]
  dmplist <- rbind(dmplist, dmp)
  cg[[i]] <- row.names(dmp)
}

## Save the results
write.csv(dmplist, "results/dmplist.csv")

## Classify DMPs according to its change in methylation
dmplist$Type <- "Hypermethylated"
dmplist$Type[which(dmplist$logFC > 0)] <- "Hypomethylated"

## Heatmap
dmpbeta <- data.frame(beta_values[unlist(cg), ])
colnames(dmpbeta) <- paste0(metadata$Sample,"_",metadata$Condition)
colors <- colorRampPalette(c("royalblue", "white", "red"))(n=100) # range of colors
heatmap.2(as.matrix(t(unique(dmpbeta))), trace="none", density.inf="none", margins=c(7,10), col=colors, cexRow = 1, lwid = c(5,15), lhei = c(5,15))

# 4.2.Differentially methylated regions
dmrlist <- data.frame() # df to store the DMRs

for(i in 1:length(contrasts)){ # for every comparison...
  # Annotate and and find significant CpG sites
  myAnnotation <- cpg.annotate(object=beta_values, datatype="array", what="Beta", analysis.type="differential", design=designMat, contrasts=T, 
                               cont.matrix=contMat, coef=gsub("_", " - ", as.character(contrasts[i])), arraytype="EPIC")
  
  if(sum(myAnnotation@ranges@elementMetadata@listData$is.sig)!=0){ # if there are significant sites...
    # Search for DMRs
    dmr <- dmrcate(myAnnotation, lambda=1000, C=2)
    # Extract the genomic ranges
    results.ranges <- extractRanges(dmr)
    dmr <- data.frame(Contrast=contrasts[i], results.ranges)
    dmrlist <- rbind(dmrlist, dmr)
  }
}

## Save the results
write.csv(dmrlist, "results/dmrlist.csv")