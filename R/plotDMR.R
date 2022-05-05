#
# Author: Raúl Sanz
#
# Description: Script to plot the methylation and the detected DMRs in a defined range of the genome.
#

# PACKAGES INSTALLATION

# BiocManager::install("Gviz")
# install.packages("dplyr")
# BiocManager::install("coMET")
# BiocManager::install("GenomicFeatures")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm39.refGene")

# REQUIRED PACKAGES

library(Gviz)
library(dplyr)
library(coMET)
library(GenomicFeatures)
library(org.Hs.eg.db)
# library(org.Mm.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(TxDb.Mmusculus.UCSC.mm39.refGene)

# FUNCTION

plot.DMR <- function(genome, chr, start, end, CpGs, DMRs, enhancers, feature){
  # genomic coordinates
  gtrack <- GenomeAxisTrack()
  
  # representation of the chromosome
  itrack <- IdeogramTrack(genome=genome, chromosome=chr)
  
  # gene model
  txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene # human model
  ### mouse model: txdb_mm39 <- TxDb.Mmusculus.UCSC.mm39.refGene
  gene_model <- GeneRegionTrack(txdb_hg19, genome=genome, chromosome=chr, showId=TRUE, geneSymbol=TRUE, name="UCSC")
  ## show symbol ID instead of transcript ID
  symbols <- unlist(mapIds(org.Hs.eg.db, gene(gene_model), "SYMBOL", "ENTREZID", multiVals="first"))
  ### mouse model: symbols <- unlist(mapIds(org.Mm.eg.db, gene(gene_model), "SYMBOL", "ENTREZID", multiVals="first"))
  symbol(gene_model) <- symbols[gene(gene_model)]
  
  # enhancers
  enh <- AnnotationTrack(start=c(enhancers$start), end=c(enhancers$end), chromosome=chr, name="enh", fill="darkgreen", col="darkgreen")
  
  # detected DMRs
  DMR_track <- AnnotationTrack(start=c(DMRs$start), end=c(DMRs$end), chromosome=chr, name="DMRs", col="purple4", fill="purple4")
  
  # heatmap of the methylation values at every CpG site
  heatmap <- DataTrack(CpGs, name=" ",chromosome = chr, type="heatmap", showSampleNames=T, cex.sampleNames=0.7, 
                      gradient=c(colorRampPalette(c("blue", "white", "red"))(n = 299)), separator=2)
  
  # average methylation level per group
  if(missing(feature)){ # if feature to group is not specified
    methylation <- DataTrack(CpGs, name="Methylation", chromosome=chr, type="a", groups=colnames(CpGs@elementMetadata))
  } else{
    methylation <- DataTrack(CpGs, name="Methylation", chromosome=chr, type="a", groups=feature)
  }
  
  plotTracks(list(itrack, gtrack, gene_model, enh, DMR_track, heatmap, methylation), from=start, to=end, 
             extend.left=0.1, extend.right=0.1, sizes=c(2,2,5,2,2,10,5))
}

# USAGE

# reference genome of the specie to plot
genome <- "hg19"
# chromosome to represent
chr <- "chr3"
# range to plot
start <- 37493945
end <- 37498950

# annotated enhancers from FANTOM5 project
enh_FANTOM <- DNaseI_FANTOM(gen=genome, chr=chr, start=start, end=end, bedFilePath="/home/rsanz/Documents/FDP/data/human/human_permissive_enhancers_phase_1_and_2.bed",
                           featureDisplay='enhancer', stacking_type="full")
enh_FANTOM.df <- as.data.frame(enh_FANTOM@range@ranges)
enhancers <- enh_FANTOM.df

# GenomicRanges object including the position of the CpG site and its beta values associated
CpGs <- gmSet@rowRanges
values(CpGs) <- beta_values

# data frame containing the detected DMRs with the human_methylation.R script 
DMRs <- read.csv("/home/rsanz/Documents/FDP/supplementary material/human/DMR_list.csv") 
### looking for DMRs at a specific gene: DMRs <- DMRs %>% filter(overlapping.genes=="ITGA9")

# feature to group the samples (optional)
feature <- metadata$Condition

plot.DMR(genome, chr, start, end, CpGs, DMRs, enhancers, feature)


# DMR.list <- data.frame()
# for(i in 1:length(contrasts)){
#   # annotating the CpG's for each contrast + finding significant CpG's
#   myAnnotation <- cpg.annotate(object=beta_values, datatype="array", what="Beta", analysis.type="differential", design=designMat, contrasts=T, 
#                                cont.matrix=contMat, coef=gsub("_", " - ", as.character(contrasts[i])), arraytype="EPIC", fdr = 0.5)
#   
#   if(sum(myAnnotation@ranges@elementMetadata@listData$is.sig)!=0){ # if there are significant CpG's
#     # test for DMRs
#     DMR_raw <- dmrcate(myAnnotation, lambda=1000, C=2)
#     # extract genomic ranges
#     results.ranges <- extractRanges(DMR_raw)
#     # saving the results
#     DMR <- data.frame(Contrast=contrasts[i], results.ranges)
#     DMR.list <- rbind(DMR.list, DMR)
#   }
# }
# 
# DMR.list$overlapping.genes
# DMR.list %>% filter(overlapping.genes=="ITGA9")
# 
# DMR.list$overlapping.genes <- unlist(DMR.list$overlapping.genes)
# # write.csv(DMR.list, "results/human/DMR_list.csv")
# 
# ## plotting the DMRs
# pal <- brewer.pal(8, "Dark2") # palette
# 
# ### colour by condition
# groups <- pal[1:length(unique(metadata$Condition))]
# names(groups) <- levels(factor(metadata$Condition))
# cols <- groups[as.character(factor(metadata$Condition))]
# 
# ### is needed to convert the data frame into a genomic ranges object first
# DMR.GR <- makeGRangesFromDataFrame(DMR.list)
# 
# DMR.plot(ranges=DMR.GR, dmr=5, CpGs=beta_values, phen.col=cols, what="Beta", arraytype="EPIC", genome="hg19")
# 
DMRs <- DMR.list
DMRs <- DMR.list %>% filter(overlapping.genes=="ITGA9")
DMRs <- DMRs[c(1,2,3),]

colSums(table(DMR.list$Contrast, DMR.list$start))==1
colSums(table(DMR.list$Contrast, DMR.list$end))==1

DMR.list$uniqueid <- paste(DMR.list$seqnames,DMR.list$start, DMR.list$end,sep="_")

table(DMR.list$Contrast, DMR.list$uniqueid)
