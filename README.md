# Bioinformatics-Honours-Project

# Project premise

PXDN (encoding the peroxidasin protein) has been shown to be dysregulated/mutated in a variety of diseases. Research into the regulation of the genes PXDN, and its homolog PXDNL, is lacking. This study aimed to use ChIPseq data mining to identify the regulatory TFs of both genes in cardiovascular cells. From these findings, the study will postulated the roles these genes may play in cardiovascular cells and how their dysregulation may be linked to cardiovascular disease pathogenesis.

# Workflow

All data processing and analysis was done using a combination of ChIP Atlas' built in tools, RStudios and Microsoft Excel. ChIP seq binding data was sourced from ChIP Atlas, a publicly accessible repository for ChIPseq data. The ChIPseq data in ChIP Atlas is derived from sequencing ChIPseq experiments using the MACS2 and Bowtie2 peak calling software. The peak call data is stored in a modified BED9 + GFF3 format. The data relevant to this study was extracted and cleaned using R Studios and the dplyr package. The data was annotated, filtered and analysed in RStudios according to the ChIPseeker workflow as described here (G Yu, LG Wang, QY He. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics 2015, 31(14):2382-2383. doi:[10.1093/bioinformatics/btv145](http://dx.doi.org/10.1093/bioinformatics/btv145)). The code used by this project is stored in this repository as Code_CardioBed and Code_Annotation_V3.

In addition, the project used the Peak Browser and Peak Enrichment Analysis functions in ChIP Atlas.

This was my first foray into bioinformatics, and I am looking to improve my coding and analysis skills in RStudio. As such, any and all comments are welcome.
