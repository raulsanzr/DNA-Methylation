#
# Author: Raúl Sanz (adapted from an Angelika Merkel script)
#
# Description: Script to generate the manifest and the annotation objects for the Infinium Mouse Methylation BeadChip from Illumina.
#

library(tidyverse)

annot <- read.csv("data/MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv")

# NCBI annotation
annotNCBI<- annot[annot$Source=="NCBI",]
annNCBIshort <- annotNCBI %>% group_by(name) %>% summarise(NCBI_Gene=paste(unique(Gene),collapse=";"),
                                                           NCBI_Transcript=paste(unique(Transcript),collapse=";"),
                                                           NCBI_Feature=paste(unique(Feature),collapse=";"))
write.csv(annNCBIshort, file="data/annNCBIshort_mouse.csv")

# UCSC annotation
annotUCSC<- annot[annot$Source=="UCSC",]
annUCSCshort <- annotUCSC %>% group_by(name) %>% summarise(UCSC_Gene= paste(unique(Gene),collapse=";"),
                                                           UCSC_Transcript=paste(unique(Transcript),collapse=";"),
                                                           UCSC_Feature=paste(unique(Feature),collapse=";"))
write.csv(annUCSCshort, file="data/annUCSCshort_mouse.csv")

# Manifest file
manifest <- read_csv("data/Infinium_Mouse_Methylation_v1.0_A1_GS_Manifest.csv", skip=7)
manifest$Name_long <- paste(manifest$Name, manifest$AddressA_ID, sep="_")
manifest <- manifest[, c("Name_long", "Probe_Type", "CHR", "MAPINFO", "Species",
                         "Genome_Build", "N_Shelf", "N_Shore", "CpG_Island",
                         "S_Shore", "S_Shelf")]
write.csv(manifest, file = "data/manifest_mouse.csv")
