#### Libraries ####
library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(pheatmap)
library(RColorBrewer)
library(sceasy)
library(reticulate)
use_condaenv("/home/jacquelyne/miniconda3/envs/sc_env", required=TRUE)
library(presto)
library(tibble)
library(fgsea)
library(DESeq2)
library(edgeR)
library(scRNAtoolVis)

#### Conversion of h5ad to seurat ####
setwd("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/FinalisedData/FinalMerge/H5AD")
sceasy::convertFormat("ADMicroglia_scANVI.h5ad", from = "anndata", to="seurat", outFile = "ADMicrogliascANVI_Seurat.RDS", main_layer="counts")
ADMic = readRDS("ADMicrogliascANVI_Seurat.RDS")
seaborn_palette = c('#1f77b5', '#fa7f11', '#279f68', '#d62829', '#ab40fc', '#8c564a', '#e377c2', '#b4bd60','#31becf', '#afc7e8','#fcba79')

#### Identification of microglial identity of each cluster ####
setwd("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/FinalisedData/FinalMerge/8_DEG")
# Find markers for every cluster compared to all remaining cells, reporting only the positive ones
Idents(ADMic) = 'leiden0_7'

pdf("1_ADMic_scANVI.pdf", height=6, width=7) #double check the latent representation preservation
DimPlot(ADMic, reduction='scANVI_umap', pt.size=2)+scale_color_manual(values=seaborn_palette)
dev.off()

Idents(ADMic) = 'Dataset'
pdf("2_ADMic_Datasets.pdf", height=6, width=8)
DimPlot(ADMic, reduction='scANVI_umap', pt.size=2)+scale_color_manual(values=seaborn_palette)
dev.off()

library(dittoSeq)
pdf('3_ADMic_DatasetDistribution.pdf', height = 2, width=4)
ggplot(ADMic@meta.data, aes(x=leiden0_7, fill=Dataset)) + geom_bar(position='fill')+scale_fill_manual(values=seaborn_palette)
dev.off()

pdf('3a_ADMic_DatasetDistribution.pdf', height = 2, width=4)
ggplot(ADMic@meta.data, aes(x=leiden0_7, fill=Dataset)) + geom_bar(position='fill')+scale_fill_manual(values=seaborn_palette) + coord_flip()
dev.off()

Idents(ADMic) = 'leiden0_7'
pdf("4_ADMic_Diagnosis.pdf",width=12, height=6)
DimPlot(ADMic, reduction = 'scANVI_umap', split.by='Disease', pt.size=2)+scale_color_manual(values=seaborn_palette)
dev.off()

pdf('5_ADMic_DistributionAPOE.pdf')
ggplot(ADMic@meta.data, aes(x=leiden0_7, fill=APOEStatus)) + geom_bar(position='fill')
dev.off()

pdf('6_ADMic_DiseaseDistribution.pdf')
ggplot(ADMic@meta.data, aes(x=leiden0_7, fill=Disease)) + geom_bar(position='fill')
dev.off()

pdf('6a_ADMic_DiseaseDistribution.pdf', height=2, width=4)
ggplot(ADMic@meta.data, aes(x=leiden0_7, fill=Disease)) + geom_bar(position='fill') + coord_flip()
dev.off()

pdf('7_ADMic_SexDistribution.pdf', height=2, width=4)
ggplot(ADMic@meta.data, aes(x=leiden0_7, fill=Sex)) + 
  geom_bar(position='fill') + 
  coord_flip() +
  scale_fill_manual(values = c("Male" = "lightblue3", "Female" = "lightpink3"))
dev.off()

#### muscat DEG analysis ####
library(SingleCellExperiment)
ADMic@meta.data
ADMic@meta.data$Age = gsub("90\\+", "90", ADMic@meta.data$Age)
ADMic@meta.data$Age = gsub("95 or above", "95", ADMic@meta.data$Age)
DefaultAssay(ADMic) = "RNA"
ADMic_sce = as.SingleCellExperiment(ADMic)
ADMic_sce

# removing undetected genes
ADMic_sce = ADMic_sce[rowSums(counts(ADMic_sce) > 0) > 0,]
dim(ADMic_sce) #[1]  43844 167441

# removing lowly expressed genes
library(scater)
qc = perCellQCMetrics(ADMic_sce)
ADMic_sce = ADMic_sce[rowSums(counts(ADMic_sce) > 1) >=10,]
dim(ADMic_sce) #1]  19847 167441
ADMic_sce$id = paste0(ADMic_sce$Disease, ADMic_sce$Sample)
(sce <- prepSCE(ADMic_sce,
  kid='leiden0_7',
  gid= "Disease",
  sid='id',
  drop=FALSE))

nk = length(kids <- levels(sce$cluster_id))
ns = length(sids <- levels(sce$sample_id))
names(kids) = kids; names(sids) = sids
t(table(sce$cluster_id, sce$sample_id))
dim(sce)

pdf('11_muscat_umap.pdf')
plotReducedDim(sce, dimred='SCANVI_UMAP', colour_by = 'cluster_id') + guides(fill=guide_legend(override.aes = list(alpha=1, size=3))) + theme_minimal() + theme(aspect.ratio=1)+ scale_color_manual(values=seaborn_palette)
dev.off()

# Differential state analysis (DS)
# aggregation of single-cell to pseudobulk data
pb = aggregateData(
  sce,
  assay = 'counts',
  fun='sum',
  by=c('cluster_id', 'sample_id')
)
assayNames(pb)
sample_metadata <- unique(as.data.frame(colData(sce)[, c("sample_id", "Age")]))
colData(pb)$sample_id = rownames(colData(pb))
pb_colData = as.data.frame(colData(pb))
age_vector <- sample_metadata$Age[match(rownames(colData(pb)), sample_metadata$sample_id)]
colData(pb)$Age <- age_vector
t(head(assay(pb)))

# pseudobulk level MDS plot
pdf("12_muscatpb_mds.pdf")
pbMDS(pb) + scale_shape_manual(values=c(17,4)) + scale_color_manual(values=seaborn_palette)
dev.off()

## limma analysis
res_limma = pbDS(pb, method='limma-trend', min_cells=2)
tbl_limma = res_limma$table[[1]]
tbl_limma <- lapply(tbl_limma, function(u){
  dplyr::arrange(u, p_adj.glb)
})
tbl_limma_fil <- lapply(tbl_limma, function(u) {
  u = dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 0.25)
  dplyr::arrange(u, p_adj.loc)
})
write_degs <- function(df, name) {
  filename <- paste0(name, "_deg.txt")
  write.table(df, file = filename, sep="\t", row.names=F, quote=F)
  return(filename)
}

setwd("~/13_muscat_deg/limma")
mapply(write_degs, tbl_limma, names(tbl_limma))

#### DEG analysis in seurat object ####
setwd("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/FinalisedData/FinalMerge/8_DEG_AllGenes/14_ADvNC_Clusters")
cluster_numbers = as.matrix(table(ADMic$Disease, ADMic$leiden0_7))
write.table(cluster_numbers, file = '14_leiden07_diagnosisnumber.txt', quote=F, sep="\t")

ADMic$cluster_disease = paste(ADMic$leiden0_7, ADMic$Disease, sep="_")
Idents(ADMic) = "cluster_disease"
Clus0_DEG = FindMarkers(ADMic, ident.1 = "0_AD", ident.2 = "0_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus0_DEG, file="1_Clus0_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus1_DEG = FindMarkers(ADMic, ident.1 = "1_AD", ident.2 = "1_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus1_DEG, file="2_Clus1_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus2_DEG = FindMarkers(ADMic, ident.1 = "2_AD", ident.2 = "2_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus2_DEG, file="3_Clus2_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus3_DEG = FindMarkers(ADMic, ident.1 = "3_AD", ident.2 = "3_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus3_DEG, file="4_Clus3_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus4_DEG = FindMarkers(ADMic, ident.1 = "4_AD", ident.2 = "4_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus4DEG, file="5_Clus4_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus5_DEG = FindMarkers(ADMic, ident.1 = "5_AD", ident.2 = "5_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus5_DEG, file="6_Clus5_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus6_DEG = FindMarkers(ADMic, ident.1 = "6_AD", ident.2 = "6_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus6_DEG, file="7_Clus6_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus7_DEG = FindMarkers(ADMic, ident.1 = "7_AD", ident.2 = "7_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus7_DEG, file="8_Clus7_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus8_DEG = FindMarkers(ADMic, ident.1 = "8_AD", ident.2 = "8_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus8_DEG, file="9_Clus8_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus9_DEG = FindMarkers(ADMic, ident.1 = "9_AD", ident.2 = "9_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus9_DEG, file="10_Clus9_ADvNC.txt", quote=F, sep="\t", row.names=T)

Clus10_DEG = FindMarkers(ADMic, ident.1 = "10_AD", ident.2 = "10_NC", verbose=FALSE, only.pos=FALSE, min.pct=0.25, test.use="wilcox")
write.table(Clus10_DEG, file="11_Clus10_ADvNC.txt", quote=F, sep="\t", row.names=T)

#### Visualisations ####
library(EnhancedVolcano)
setwd("/Users/jacquelyne/Documents/Data/IronMicroglia/FinalMerge/9_DEG_ProteinCodingGenes/13_muscat_deg/limma_sig/clus3")
data = read.table("Clus3.txt", sep="\t", header=T)

pdf(file="1_EnhancedVolcano_Clus3.pdf")
EnhancedVolcano(
  data,
  lab=data$gene,
  ylab=bquote(-Log[10]~'P-adj'),
  x='logFC',
  y="p_adj.glb",
  selectLab = c('*FASN*', '*FTH1*','*FTL*','*ACSL1*','*MDM2*','*P4HA1*','*SLC11A1*','*ATG7*','*PICALM*','*BTBD9*','*FOXP1*','*TTC7A*','*BCL2*','*SLC25A37*'),
  pCutoff=0.69,
  xlim=c(-100,100),
  ylim=c(0,5),
  drawConnectors=TRUE,
  parseLabels=TRUE,
  boxedLabels=TRUE,
  legendLabels=c('Not sig.', 'logFC','p-adj','padj & logFC')
)
dev.off()

pdf(file="1a_EnhancedVolcano_Clus3.pdf")
EnhancedVolcano(
  data,
  lab=NA,
  ylab=bquote(-Log[10]~'P-adj'),
  x='logFC',
  y="p_adj.glb",
  pCutoff=0.69,
  xlim=c(-100,100),
  ylim=c(0,5),
  drawConnectors=TRUE,
  parseLabels=TRUE,
  boxedLabels=TRUE,
  legendLabels=c('Not sig.', 'logFC','p-adj','padj & logFC')
)
dev.off()

#### Gene Expression of genes of interest against Cogdx Score ####
setwd("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/FinalisedData/FinalMerge/9_DEG_ProteinCodingGenes/13_muscat_deg/geneexpression")
ADMic
Idents(ADMic) = "leiden0_7"
MG3 = subset(ADMic, idents="3")
pseudo_MG3 = AggregateExpression(MG3, assays = "RNA", return.seurat = T, group.by = c("Disease", "Sample"))
MG3$sample_disease = paste(MG3$Disease, MG3$Sample, sep="_")
tail(Cells(pseudo_MG3))
MG3_pseudocounts = pseudo_MG3@assays$RNA@layers$counts
MG3_pseudocounts = as.data.frame(MG3_pseudocounts)
dim(MG3)
colnames(MG3_pseudocounts) = Cells(pseudo_MG3)
rownames(MG3_pseudocounts) <- rownames(MG3@assays$RNA)
write.table(MG3_pseudocounts, file="1_MG3_pseudocounts.txt", quote=F, sep="\t", row.names=T)
MG3meta = MG3@meta.data
unique_MG3meta = MG3meta %>% distinct(Sample, .keep_all = TRUE)
write.table(unique_MG3meta, file = "2_MG3_meta.txt", quote=F, sep="\t", row.names=F)

GeneExpr = read.table("1_MG3_pseudocounts.txt", header=TRUE, row.names=1, sep="\t", check.names=F)
SampleInfo = read.table("2_MG3_meta.txt", header=TRUE, row.names=1, sep="\t")
rownames(SampleInfo) = SampleInfo$sample_disease
Samples=intersect(colnames(GeneExpr), rownames(SampleInfo))
GeneExpr = GeneExpr[rowSums(GeneExpr)>0,Samples]
SampleInfo = SampleInfo[Samples,]
all(colnames(GeneExpr)==rownames(SampleInfo))
BigData=cbind(SampleInfo, t(GeneExpr))

gene=c("ACSL1","FASN","MDM2","FTH1","FTL","P4HA1", "ATG7","SLC11A1","ZFP36L1","ZFP36L2")

unique(SampleInfo$cogdx)
BigDataCogdx=BigData[BigData$cogdx %in% c(1,2,3,4,5,6),]
tmp = BigDataCogdx[,c("Disease","cogdx",gene)]
library(reshape2)
tmp1 = melt(tmp, id=c(1,2))
colnames(tmp1) = c("Disease","cogdx","Gene","Expression")
tmp1$Disease = factor(tmp1$Disease, levels=c("NC", "AD"))
tmp1$cogdx_group <- factor(paste0(tmp1$cogdx, "_", tmp1$Disease),
                           levels = c("1_NC", "1_AD", "2_NC", "2_AD", "3_NC", "3_AD", 
                                      "4_NC", "4_AD", "5_NC", "5_AD", "6_NC", "6_AD"))

tmp1$Disease <- factor(tmp1$Disease, levels = c("NC", "AD"))
tmp1$cogdx_merged <- case_when(
  tmp1$cogdx %in% c("1", "2") ~ "1-2",
  tmp1$cogdx %in% c("3", "4") ~ "3-4",
  tmp1$cogdx %in% c("5", "6") ~ "5-6",  # Assuming you want this too
  TRUE ~ as.character(tmp1$cogdx)
)
tmp1$cogdx_merged <- factor(tmp1$cogdx_merged, levels = c("1-2", "3-4", "5-6"))

# Verify
table(tmp1$cogdx, tmp1$cogdx_merged, tmp1$Disease)

# Ensure factors for ordering
tmp1$Disease <- factor(tmp1$Disease, levels = c("NC", "AD"))
tmp1$cogdx <- factor(tmp1$cogdx, levels = c("1", "2", "3", "4", "5", "6"), ordered = TRUE)

g5 <- ggplot(tmp1, aes(x = cogdx_merged, y = log2(Expression + 2), color = Disease)) +  
  geom_point(position = position_jitter(width = 0.2, height = 0), 
             size = 1, alpha = 0.5) +  # Smaller, 50% transparent dots
  stat_summary(fun = "median", geom = "line", color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.3, color = "black") +
  facet_wrap(~ Gene, scales = "free_y", ncol = 1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +  # Remove all grid lines
  scale_color_manual(values = c("NC" = "blue", "AD" = "red")) +
  ggpubr::stat_compare_means(aes(group = cogdx_merged), method = "kruskal.test", label = "p.signif") +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = rel(1.0), face = "italic"),
        axis.text.x = element_text(size = rel(1.0)),
        legend.position = "top") + 
  labs(x = "Cogdx Score", y = "log2(Expression + 2)", color = "Disease") +
  scale_y_continuous(expand = c(0.1, 0.5))

pdf("Cogdx_Expr_scatter_merged.pdf", height = 16, width = 2)
print(g5)
dev.off()

g6<- ggplot(tmp1, aes(x = cogdx_merged, y = log2(Expression + 2), color = Disease)) +  
  geom_point(position = position_jitter(width = 0.2, height = 0), 
             size = 1, alpha = 0.5) +  # Smaller, 50% transparent dots
  stat_summary(fun = "median", geom = "line", color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.3, color = "black") +
  facet_wrap(~ Gene, scales = "free_y", ncol = 1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +  # Remove all grid lines
  scale_color_manual(values = c("NC" = "blue", "AD" = "red")) +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = rel(1.0), face = "italic"),
        axis.text.x = element_text(size = rel(1.0)),
        legend.position = "top") + 
  labs(x = "Cogdx Score", y = "log2(Expression + 2)", color = "Disease") +
  scale_y_continuous(expand = c(0.1, 0.5))

pdf("Cogdx_Expr_scatter_no_significance_merged.pdf", height = 16, width = 2)
print(g6)
dev.off()

kw_merged3 <- tmp1 %>%
  group_by(Gene) %>%
  kruskal_test(log2(Expression + 2) ~ cogdx_merged) %>%
  adjust_pvalue(p.col="p", output.col = "p_adj", method="fdr")
print(kw_merged3)
write.csv(kw_merged3, "cogdx_disease_pvalues_merged.csv", row.names = FALSE)
