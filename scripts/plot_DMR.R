#
# Author: Raúl Sanz
#
# Description: Script to create plots for the differential methylation and the detected DMRs in a specified range of the genome with Gviz. Aimed
#              to be used with the results produced with the script human_methylation.R or mouse_methylation.R.
#

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

# DEFINING ARGUMENTS

# reference genome of the specie to plot
genome <- "hg19"
# chromosome to represent
chr <- "chr11"
# range to plot
start <- 70670000
end <- 70675000

# annotated enhancers from FANTOM5 project
enh_FANTOM <- DNaseI_FANTOM(gen=genome, chr=chr, start=start, end=end, bedFilePath="data/human_permissive_enhancers_phase_1_and_2.bed",
                           featureDisplay='enhancer', stacking_type="full")
enh_FANTOM.df <- as.data.frame(enh_FANTOM@range@ranges)
enhancers <- enh_FANTOM.df

# GenomicRanges object including the position of the CpG site and its beta values associated
CpGs <- gmSet@rowRanges
values(CpGs) <- beta_values

# data frame containing the detected DMRs with the human_methylation.R script 
DMRs <- read.csv("results/DMRs/DMR_Control_E542K.csv") 
### looking for DMRs at a specific gene: DMRs <- DMRs %>% filter(overlapping.genes=="ITGA9")

# feature to group the samples (optional)
feature <- metadata$Condition

# USAGE

plot.DMR(genome, chr, start, end, CpGs, DMRs, enhancers, feature)
