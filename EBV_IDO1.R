library(BiocManager)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(limma)
library(pathfindR)
library(msigdbr)
library(stringr)
library(pheatmap)

# Load sample mapping
sample_mapping = read_csv("sample_mapping.csv")
#view(sample_mapping)


# Read raw file to get column names
raw_counts = read_delim("counts_ebv_combined.out",
                        delim="\t",
                        skip=1)

# Number of genes
nrow(raw_counts)

# Process column names to shorter labels
col_names = gsub(".*/", "", names(raw_counts[-c(1:6)])) 
col_names_short = gsub("Aligned.sortedByCoord.out.bam", "", col_names) 

col_names_short

# Create Column name string 
col_name_string = list('GeneId', 'Chr', 'Start', 'End', 'Strand', 'Length')

for(i in col_names_short){
  col_name_string[[length(col_name_string)+1]] <- i
}
col_name_string = unlist(col_name_string)
print(col_name_string)

# Read raw file
raw_counts = read_delim("counts_ebv_combined.out",
                        delim="\t",
                        skip=2,
                        col_names = col_name_string)
head(raw_counts)



# Number of genes with no counts
sum(rowSums(raw_counts[-c(1:6)]) == 0)
# Drop these rows
raw_counts = raw_counts[!rowSums(raw_counts[-c(1:6)]) == 0,]
sum(rowSums(raw_counts[-c(1:6)]) == 0)

# Read Human cDNA file for meta data
# ensembl = useEnsembl(biomart="ensembl")
# listDatasets(ensembl)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# listAttributes(ensembl)
gene_info <- getBM(attributes=c('ensembl_gene_id',
                                'external_gene_name', 'gene_biotype'),
                   values = raw_counts$GeneId,
                   mart = ensembl)
# Join counts table with gene information
df_raw_counts = data.frame(raw_counts)
# remove chromosome information columns
df_raw_counts = df_raw_counts[,-c(2:6)]
# view(df_raw_counts)
annot_counts = right_join(x = gene_info, y = df_raw_counts, by=c('ensembl_gene_id' = 'GeneId'))
ncol(annot_counts)
#view(annot_counts[,-c(1)])
df_count = data.frame(annot_counts[,-c(1)], row.names = annot_counts$ensembl_gene_id)
head(df_count)
ncol(df_count)

# Create factors
head(sample_mapping)
sample_id = sample_mapping$SampleName[order(sample_mapping$SampleName)]
sample_group = sample_mapping$SampleGroup[order(sample_mapping$SampleName)]
external_name = sample_mapping$ExternalSampleName[order(sample_mapping$SampleName)]
replicate = sample_mapping$SampleReplicate[order(sample_mapping$SampleName)]
dataset = sample_mapping$Dataset[order(sample_mapping$SampleName)]
replicate_dataset = sample_mapping$Replicate_dataset[order(sample_mapping$SampleName)]
lin = sample_mapping$Lin[order(sample_mapping$SampleName)]

# Sample information (columns)
myfactors = data.frame(sample_group = sample_group,
                       replicate = replicate,
                       id = sample_id,
                       external_name = external_name,
                       dataset = dataset,
                       lin = lin)

# Seperate datasets
df_count_dpi2 = df_count[,sample_mapping$SampleName[sample_mapping$Dataset == 'dpi2']]
df_count_dpi4 = df_count[,sample_mapping$SampleName[sample_mapping$Dataset == 'dpi4']]
ncol(df_count_dpi2)

myfactors_dpi2 = myfactors[myfactors$dataset == 'dpi2',]
myfactors_dpi4 = myfactors[myfactors$dataset == 'dpi4',]
nrow(myfactors_dpi2)

####### DPI-2

# Drop hiEBV-Lin-HD2 - technical problem with this sample
df_count_dpi2_subset = subset(df_count_dpi2, select= -c(BSSE_QGF_212846))
ncol(df_count_dpi2_subset)
myfactors_dpi2_subset = myfactors_dpi2[!(myfactors_dpi2$id %in% c("BSSE_QGF_212846")),]
nrow(myfactors_dpi2_subset)
#view(myfactors_dpi2_subset)


### Read data into DESeq2
dds_2dpi <- DESeqDataSetFromMatrix(countData = df_count_dpi2_subset,
                                   colData = myfactors_dpi2_subset,
                                   design = ~ replicate + sample_group)


# Filter data to keep only genes with 10 or more reads in at least 3 samples
keep <- rowSums(counts(dds_2dpi, normalized=FALSE) >= 10) >=3
dds_2dpi <- dds_2dpi[keep,]

# Access normalized values, normalized using VST
norm_vst_dpi2 <- vst(dds_2dpi, blind=FALSE)
df_norm_dpi2 = data.frame(assay(norm_vst_dpi2))
ncol(df_norm_dpi2)


# PCA before correction, sample_groups
pcaData <- plotPCA(norm_vst_dpi2, intgroup=c("replicate", 'sample_group'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=sample_group, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle('EBV 2 dpi, no correction') +
  coord_fixed()

# Batch correction
batch_corr_dpi2 = norm_vst_dpi2
assay(batch_corr_dpi2) = limma::removeBatchEffect(assay(batch_corr_dpi2), batch_corr_dpi2$replicate)
df_dpi2_corr = data.frame(assay(batch_corr_dpi2))
colnames(df_dpi2_corr) = myfactors_dpi2_subset$external_name

# PCA after correction, sample_groups
pcaData <- plotPCA(batch_corr_dpi2, intgroup=c("replicate", 'sample_group'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=sample_group, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle('EBV 2 dpi, batch correction') +
  coord_fixed()

# PCA after correction, sample_groups simplified
pcaData <- plotPCA(batch_corr_dpi2, intgroup=c('sample_group'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=sample_group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle('EBV 2 dpi, batch correction') +
  coord_fixed()


############# EBV vs EBV - Lin

# Get DEseq results
dds_2dpi_ebv = DESeq(dds_2dpi)
res_EBV = results(dds_2dpi_ebv, contrast=c("sample_group", "EBV-Lin__2dpi", "EBV__2dpi"))
res_EBV = na.omit(res_EBV)

df_EBV = data.frame(res_EBV)

# add annotation data
df_annot_EBV = merge(df_EBV, gene_info, by.x=0, by.y='ensembl_gene_id', all.x=TRUE)
row.names(df_annot_EBV) = df_annot_EBV$Row.names
# remove duplicate column
df_annot_EBV = df_annot_EBV[-1]
# Take only protein coding genes
df_annot_EBV = df_annot_EBV[df_annot_EBV$gene_biotype == 'protein_coding', ]
# Add significance annotation column
df_annot_EBV$sig = ifelse(df_annot_EBV$padj <= 0.05, "yes", "no")


### Volcano plot to get an overview
# Filter protein coding genes
ggplot(df_annot_EBV, aes(x = log2FoldChange, y = -log10(padj), color = sig,)) +
  geom_point() +
  geom_text_repel(aes(label=ifelse((padj <= 0.05 & (log2FoldChange > 2 | log2FoldChange < -2)) | (padj <= 0.000000001) ,
                                   as.character(external_gene_name),'')),color='black',hjust=1,vjust=-3, max.overlaps = 120)


# Filter for only significant hits
df_annot_EBV_sig = df_annot_EBV[df_annot_EBV$sig == 'yes', ]
nrow(df_annot_EBV_sig)

# Take batch corrected read counts of hits significant in EBV-Lin vs EBV
df_dpi2_corr_EBV_sig = df_dpi2_corr[row.names(df_annot_EBV_sig),]
nrow(df_dpi2_corr_EBV_sig)
df_dpi2_corr_EBV_sig = na.omit(df_dpi2_corr_EBV_sig)
nrow(df_dpi2_corr_EBV_sig)

## Annotation column for heatmap
myfactors_dpi2_sorted = myfactors_dpi2_subset[order(myfactors_dpi2_subset$external_name),]
annotation_col = data.frame(as.character(myfactors_dpi2_sorted$sample_group),
                            row.names = myfactors_dpi2_sorted$external_name)
colnames(annotation_col) = "Condition"
head(annotation_col)

# Legend 
legend_breaks =  seq(-4, 4, 0.16)
legend_breaks

# Color map
color_scale = colorRampPalette(c("navy", "white", "firebrick3"))(50)

# Draw heatmap
pheatmap(df_dpi2_corr_EBV_sig, annotation_col = annotation_col,
         show_rownames = FALSE, main = 'Genes Sig. EBV vs EBV-Lin', scale = 'row',
         color = color_scale, breaks = legend_breaks, border_color = NA)

# Write table
write.csv(df_dpi2_corr_EBV_sig, "sig_EBV_vs_EBV-lin.csv")

