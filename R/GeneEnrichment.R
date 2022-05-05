#
# Author: Raúl Sanz
# 
# Description: Script to extract the set of genes overlapping the detected DMRs with R/MouseMethylation.R or R/MethylationEPIC.R
#

# PACKAGES INSTALLATION

# BiocManager::install("biomaRt")

# REQUIRED PACKAGES

library(biomaRt)

# OVERLAPPING GENES

# getting the human overlapping genes from the DMR list
DMR_human <- read.csv("results/human/DMR_list.csv")
genes_human <- DMR_human$overlapping.genes

## clean the gene list (delete NAs and repeated)
genes_human_flt <- unique(unlist(sapply(genes_human, function(i) strsplit(i, ", "))))
genes_human_flt <- na.omit(genes_human_flt)

## saving the results
# write(paste(genes_human_flt, collapse="\n"), file="results/human/genelist.txt")

# getting the mouse overlapping genes from the DMR set
DMR_mouse <- read.csv("results/mouse/DMR_list.csv")

## annotating the mouse DMRs
ann <- read.csv("data/annNCBIshort_mouse.csv")
rownames(ann) <- ann$name
positions <- list()
for(i in c(1:6,8:nrow(DMR_mouse))){
  positions <- append(positions, c(DMR_mouse[i,]$indexStart:DMR_mouse[i,]$indexEnd))
}
cg_names <- substr(rownames(beta_values[unlist(positions),]),1,10)
genes_mouse <- unique(ann[cg_names,]$NCBI_Gene)

## clean the gene list (delete NAs and repeated)
genes_mouse_flt <- unique(unlist(sapply(genes_mouse, function(i) strsplit(i, ";"))))
genes_mouse_flt <- na.omit(genes_mouse_flt)

## saving the results
# write(paste(genes_mouse_flt, collapse="\n"), file="results/mouse/genelist.txt")

# COMPARING GENE LISTS

# finding the homologous genes in human from the mouse genes
human_ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
mouse_ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
genes_homologous <- getLDS(attributes=c("hgnc_symbol","chromosome_name", "start_position"), 
                           filters="hgnc_symbol", values=genes_mouse_flt, mart=human_ensembl,
                           attributesL=c("chromosome_name","start_position"), martL=mouse_ensembl)
## saving the results
# write(paste(genes_homologous$HGNC.symbol, collapse="\n"), file="results/mouse/genelist_homologous.txt")

# comparing the human set with the (homologous) mouse set of genes 
intersect(genes_homologous$HGNC.symbol, genes_human_flt)