# EBV_IDO1
Volcano and Heatmap generation from Transcriptome data 2 days post infection.

Part of publication titles:
A metabolic vulnerability of EBV can be targeted to prevent viral pathology.

In brief, RNA sequencing was done in paired-end mode on an Illumina NextSeq 500 with 38 cycles. The raw reads were aligned against the NCBI human GRCh38 assembly using STAR version 2.7.9a in paired-end mode.
The read counts of each sample were summarized using featureCounts from the Subread package using strand information and the NCBI GRCh38 annotations. 
Differential gene expression was analyzed using DESeq2 version 1.38.2. Gene biotype and gene symbol annotation was fetched from Ensemble using biomaRt version 2.54.1. 
Protein coding genes and genes with more than 10 reads in at least 3 samples were kept. Principal Components Analysis (PCA) were plotted using the 500 most variable genes. 
Transcription changes due to batch effects were corrected for using the removeBatchEffect function of limma version 3.54.2.
Volcano plot was drawn with a significanse cut-off at 0.05 adjusted p-value.

