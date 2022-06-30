#
# Author: Raul Sanz
# Description: Script to generate and format the manifest and the annotation files for the Infinium Mouse Methylation array.
#

# DEPENDENCIES

library(tidyverse) # install.packages("tidyverse")

# ANNOTATION

# reading the annotation
ann <- read.csv("data/MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv")

# NCBI annotation
annNCBI<- ann[ann$Source=="NCBI",]
annNCBIshort <- annNCBI %>% group_by(name) %>% summarise(NCBI_Gene=paste(unique(Gene),collapse=";"),
                                                         NCBI_Transcript=paste(unique(Transcript),collapse=";"),
                                                         NCBI_Feature=paste(unique(Feature),collapse=";"))
## saving the results
write.csv(annNCBIshort, file="data/ann_NCBI_mouse.csv")

# UCSC annotation
annUCSC<- ann[ann$Source=="UCSC",]
annUCSCshort <- annUCSC %>% group_by(name) %>% summarise(UCSC_Gene= paste(unique(Gene),collapse=";"),
                                                         UCSC_Transcript=paste(unique(Transcript),collapse=";"),
                                                         UCSC_Feature=paste(unique(Feature),collapse=";"))
## saving the results
write.csv(annUCSCshort, file="data/ann_UCSC_mouse.csv")

# MANIFEST

# generating the manifest file
manifest <- read_csv("data/Infinium_Mouse_Methylation_v1.0_A1_GS_Manifest.csv", skip=7)
manifest$Name_long <- paste(manifest$Name, manifest$AddressA_ID, sep="_")
manifest <- manifest[, c("Name_long", "Probe_Type", "CHR", "MAPINFO", "Species",
                         "Genome_Build", "N_Shelf", "N_Shore", "CpG_Island", "S_Shore", "S_Shelf")]
## saving the results
write.csv(manifest, file = "data/manifest_mouse.csv")
