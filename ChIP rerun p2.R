
##### Data Visualisation
## Install Bioconductor core packages and Load CHIPseeker packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.13")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("clusterProfiler")
BiocManager::install("ChIPseeker", force = TRUE)
BiocManager::install("ReactomePA", force = TRUE)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

## Load libraries

library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## txdb contains database generated genomic annotation data from the UCSC
## split into 2 objects specific to chromosomes 2 and 8

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

library(clusterProfiler)

## Load in data

setwd("C:/Users/Shiven/Documents/R/Honours Project/Data/BED Files")




# CHECKPOINT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PXDNTFs = read.table("PXDNFiltered2.bed", header = TRUE, quote = '"', sep = " ", comment.char = "")
PXDNLTFs = read.table("TFsPXDNL2.bed", header = TRUE, quote = '"', sep = " ", comment.char = "")

PXDNTFs
PXDNLTFs

## Isolate columns for chrom, start, end to make a peak file

peakPXDNTFs = PXDNTFs %>% select(1:3)
peakPXDNLTFs = PXDNLTFs %>% select(1:3)

## Write peak into new bed file
## make sure separators are "\t", and file ends with an emptyline

## !!! access file and make sure there is no header line

write.table(peakPXDNTFs, "peakPXDNTab.bed", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(peakPXDNLTFs, "peakPXDNLTab.bed", sep = "\t", row.names = FALSE, col.names = FALSE)

## Annotate Peak files in 5kb radius of TSS

peakAnnoPXDN = annotatePeak(
  "peakPXDNTab.bed", 
  tssRegion=c(-5000, 5000), 
  TxDb=txdb, 
  level = "gene", 
  verbose = FALSE)

peakAnnoPXDNL = annotatePeak(
  "peakPXDNLTab.bed", 
  tssRegion=c(-5000, 5000), 
  TxDb=txdb, 
  level = "gene", 
  verbose = FALSE)

## Visualise Peak Data as pie

plotAnnoPie(peakAnnoPXDN, col = blues9)
plotAnnoPie(peakAnnoPXDNL, col = blues9)

## Visualise Peak Data as bar

plotAnnoBar(peakAnnoPXDN)

plotAnnoBar(peakAnnoPXDNL)

## Visualise data relative to TSS

plotDistToTSS(peakAnnoPXDN, title="")

plotDistToTSS(peakAnnoPXDNL, title = "")

## Store annotation data

PXDNanno = as.data.frame(peakAnnoPXDN@anno)

PXDNLanno = as.data.frame(peakAnnoPXDNL@anno)

write.table(PXDNanno, file="PXDNanno.txt", sep = "\t", row.names = FALSE)

write.table(PXDNLanno, file="PXDNLanno.txt", sep = "\t", row.names = FALSE)

## Load in annotated data

PXDNanno <- read.table("PXDNanno.txt", sep = "\t")

PXDNLanno <- read.table("PXDNLanno.txt", sep = "\t")

## Attach annotation data to tf data

library(tidyverse)

## Load missing columns

missingDataPXDN = read.table("PXDNFiltered2.bed", header = FALSE, quote = '"', sep = " ", comment.char = "")
missingDataPXDNL = read.table("TFsPXDNL2.bed", header = FALSE, quote = '"', sep = " ", comment.char = "")

missingDataPXDN
missingDataPXDNL

Macs_SRX_TF_CellType = select(missingDataPXDN, -c(V1, V2, V3))
Macs_SRX_TF_CellTypePXDNL = select(missingDataPXDNL, -c(V1, V2, V3))

Macs_SRX_TF_CellType

fullAnnoPXDN = cbind(PXDNanno, Macs_SRX_TF_CellType)
fullAnnoPXDNL = cbind(PXDNLanno, Macs_SRX_TF_CellTypePXDNL)

fullAnnoPXDN

colnames(fullAnnoPXDN) = c("chr", "start", "end", "width", "strand", "annotation", "geneCh", "geneStart", "geneEnd", "geneLen", "geneStr", "geneID", "distanceToTSS", "MACSqValue", "SRX_ID", "TF", "cell_type")
colnames(fullAnnoPXDNL) = c("chr", "start", "end", "width", "strand", "annotation", "geneCh", "geneStart", "geneEnd", "geneLen", "geneStr", "geneID", "distanceToTSS", "MACSqValue", "SRX_ID", "TF", "cell_type")

write.table(fullAnnoPXDN, file="PXDNFullAnno.bed", sep = "\t", row.names = FALSE)
write.table(fullAnnoPXDNL, file="PXDNLFullAnno.bed", sep = "\t", row.names = FALSE)

## Test for .bed files

peak = readPeakFile("Test2.bed")

peak = readPeakFile("peakPXDNTab.bed")

##### Data Filtering

library(ChIPseeker)
library(dplyr)

## Load in annotated files if needed

#FullAnnoPXDN = read.table("PXDNFullAnno.bed", header = TRUE, quote = "", sep = "\t", comment.char = "")
#FullAnnoPXDNL = read.table("PXDNLFullAnno.bed", header = TRUE, quote = "", sep = "\t", comment.char = "")

## Filter annotated data within 5 kb of TSS (focus on promoter region, first intron and distal intergenic)

## Filter using gene IDs for pxdn and PXDNL
## PXDN: "7837"
## PXDNL: "137902"

filtAnnoPXDN_1 = fullAnnoPXDN %>% 
  filter(geneID == c("7837"))

filtAnnoPXDNL_1 = fullAnnoPXDNL %>% 
  filter(geneID == c("137902"))


## Filter for regions of interest

filtAnnoPXDN = filtAnnoPXDN_1 %>% 
  filter(annotation %in% c("Promoter (1-2kb)", "Promoter (<=1kb)"))

filtAnnoPXDNL = filtAnnoPXDNL_1 %>% 
  filter(annotation %in% c("Promoter (4-5kb)", "Promoter (<=1kb)"))

## Filter for intronic regions

filtAnnoPXDN_intron = filtAnnoPXDN_1 %>% 
  filter(annotation %in% c("Intron (ENST00000252804.9/7837, intron 1 of 22)"))

filtAnnoPXDNL_intron = filtAnnoPXDNL_1 %>% 
  filter(annotation %in% c("Intron (ENST00000356297.5/137902, intron 1 of 22)"))

## Filter for high Macsq values

filtMacsAnnoPXDN = filtAnnoPXDN %>% 
  filter(MACSqValue > 200)

filtMacsAnnoPXDNL = filtAnnoPXDNL %>% 
  filter(MACSqValue > 200)

filtMacsAnnoPXDN_intron = filtAnnoPXDN_intron %>% 
  filter(MACSqValue > 200)

filtMacsAnnoPXDNL_intron = filtAnnoPXDNL_intron %>% 
  filter(MACSqValue > 200)

## remove any unwanted data, like epitope tags

filtMacsAnnoPXDN = filtAnnoPXDN %>% 
  filter(!(TF == "Epitope"))

filtMacsAnnoPXDN_intron = filtMacsAnnoPXDN_intron %>% 
  filter(!(TF == "Epitope"))

## Write into files

write.table(filtMacsAnnoPXDN, file="filtMacsAnnoPXDN.txt", sep = "\t", row.names = FALSE)
write.table(filtMacsAnnoPXDNL, file="filtMacsAnnoPXDNL.txt", sep = "\t", row.names = FALSE)



write.table(filtMacsAnnoPXDN_intron, file="filtMacsAnnoPXDN_intron.txt", sep = "\t", row.names = FALSE)
write.table(filtMacsAnnoPXDNL_intron, file="filtMacsAnnoPXDNL_intron.txt", sep = "\t", row.names = FALSE)




















##### Graph Plotting

library(ggplot2)

# TFs binding to PXDNs promoter region

par(mfrow=c(4,1))

filtMacsAnnoPXDN$n = 1
filtMacsAnnoPXDN$colorCategory <- ifelse(filtMacsAnnoPXDN$MACSqValue > 500, "MACSqscore > 500", "MACSqscore < 500")

TFPXDN = filtMacsAnnoPXDN %>%
  arrange(MACSqValue) %>%
  ggplot(aes(x = TF, n, fill = colorCategory)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("MACSqscore > 500" = "blue", "MACSqscore < 500" = "deepskyblue"))

TFPXDN + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of studies identifying significant transcription factors \nbinding to PXDN in Cardiovascular Cells") +
  labs(y = "Number of studies", x = "Transcription Factors", fill = "Significance Level") + theme(plot.title = element_text(hjust = 0.5))


# TFs binding to PXDNLs promoter region

filtMacsAnnoPXDNL$n = 1
filtMacsAnnoPXDNL$colorCategory <- ifelse(filtMacsAnnoPXDNL$MACSqValue > 500, "MACSqscore > 500", "MACSqscore < 500")

TFPXDNL = filtMacsAnnoPXDNL %>%
  arrange(MACSqValue) %>%
  ggplot(aes(x = TF, n, fill = colorCategory)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("MACSqscore > 500" = "blue", "MACSqscore < 500" = "deepskyblue"))

TFPXDNL + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of studies identifying significant transcription factors \nbinding to PXDNL in Cardiovascular Cells") +
  labs(y = "Number of studies", x = "Transcription Factors", fill = "Significance Level") + theme(plot.title = element_text(hjust = 0.5))


# TFs binding to PXDN's first intron

filtMacsAnnoPXDN_intron$n = 1
filtMacsAnnoPXDN_intron$colorCategory <- ifelse(filtMacsAnnoPXDN_intron$MACSqValue > 500, "MACSqscore > 500", "MACSqscore < 500")

TFPXDN = filtMacsAnnoPXDN_intron %>%
  arrange(MACSqValue) %>%
  ggplot(aes(x = TF, n, fill = colorCategory)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("MACSqscore > 500" = "blue", "MACSqscore < 500" = "deepskyblue"))

TFPXDN + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of studies identifying significant transcription factors \nbinding to PXDN's first Intron in Cardiovascular Cells") +
  labs(y = "Number of studies", x = "Transcription Factors", fill = "Significance Level") + theme(plot.title = element_text(hjust = 0.5))

# TFs in PXDNLs 1st intron


filtMacsAnnoPXDNL_intron$n = 1
filtMacsAnnoPXDNL_intron$colorCategory <- ifelse(filtMacsAnnoPXDNL_intron$MACSqValue > 500, "MACSqscore > 500", "MACSqscore < 500")

TFPXDNL = filtMacsAnnoPXDNL_intron %>%
  arrange(MACSqValue) %>%
  ggplot(aes(x = TF, n, fill = colorCategory)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("MACSqscore > 500" = "blue", "MACSqscore < 500" = "deepskyblue"))

TFPXDNL + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of studies identifying significant transcription factors \nbinding to PXDNL in Cardiovascular Cells") +
  labs(y = "Number of studies", x = "Transcription Factors", fill = "Significance Level") + theme(plot.title = element_text(hjust = 0.5))


colnames(filtMacsAnnoPXDN)


















## SCRAP CODE FOR FUTURE PROJECTS

# checking annotations in txdb object for 5' UTR
library(AnnotationDbi)
library(GenomicFeatures)




# check valid keytypes of current reference genome

keytypes(txdb)

# check columns of txdb

columns(txdb)

# Get 5' UTRs for all transcripts
five_utrs <- fiveUTRsByTranscript(txdb, use.names = TRUE)

# Check for 3' as well

three_utrs <- threeUTRsByTranscript(txdb, use.names = TRUE)

# Check for the Gene ID type used for UTRs

names(five_utrs)

# create gene_identifiers object using Ensembl gene IDs for PXDN/L

gene_identifiers = c("ENST00000252804.9", "ENST00000356297.5")

# Check if 5' UTRs exist for PXDN and PXDNL

filtered_five_utrs <- five_utrs[names(five_utrs) %in% gene_identifiers]

# repeat for three

filtered_three_utrs <- three_utrs[names(three_utrs) %in% gene_identifiers]


lengths(filtered_five_utrs)

lengths(filtered_three_utrs)








# Find TSS sequence for PXDN and PXDNL

# create object to store gene of interest
gene_numbers = c("7837", "137902")   # gene id for pxdn and pxdnl from NCBI


# get identifiers for PXDN/l

PXDN_name <- select(txdb, keys="7837", keytype="GENEID", columns="GENEID")


all_transcripts <- transcriptsBy(txdb, by="gene")

# Find the TSSs PXDN and PXDNL gene_numbers
# Assuming the gene_id for PXDN is in the returned table and it's the first one (it may return multiple IDs)
all_transcripts
all_transcripts[names(all_transcripts) %in% gene_identifiers]

pxdn_transcripts <- all_transcripts["7837"]

# Get the TSS sites
pxdn_transcripts

# To view the range
pxdn_tss
