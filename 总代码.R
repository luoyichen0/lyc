
##############################
#COPD
library(Seurat)
samples=list.files("./")
dir =samples
names(dir) = samples
scRNAlist <- list()
for(i in 1:length(dir)){
counts <- Read10X(data.dir = dir[i])
scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells=1)
}
scRNA<- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], 
                scRNAlist[[4]]))

for (i in 1:length(scRNAlist)) {
    scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
    scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
}
##以VariableFeatures为基础寻找锚点，运行时间较长
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
##利用锚点整合数据，运行时间较长
scRNA <- IntegrateData(anchorset = scRNA.anchors)


immune.combined <- scRNA
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
pbmc <- immune.combined
DefaultAssay(pbmc) <- "RNA"
pbmc[["percent.mt"]] <-PercentageFeatureSet(pbmc,pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 20)
DimPlot(object = pbmc, reduction = "umap", group.by ='orig.ident')

saveRDS(pbmc, file = "COPD.rds")

levels(pbmc@meta.data$orig.ident) <- c("COPD", "health")
Idents(pbmc) <- as.character(pbmc@meta.data$orig.ident)
DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")

pbmc1 <- readRDS("/SGRNJ03/pipline_test/luoyichen/COPD/COPD/COPD.rds")
pbmc1@meta.data$orig.ident <- gsub("COPD[1-3]", "COPD",pbmc1@meta.data$orig.ident)
pbmc1@meta.data$orig.ident <- gsub("healthy[1-3]", "Health",pbmc1@meta.data$orig.ident)
library(ggplot2)

pdf("COPD_3X3_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc1, reduction = "umap")+ labs(title = "UMAP of 6 samples")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()


pdf("COPD_1X1_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc1, reduction = "umap", group.by = "orig.ident")+ labs(title = "UMAP of Disease")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

labels <- read.csv("/SGRNJ03/pipline_test/luoyichen/COPD/outs3/n4-test_cell.annot.csv", header = TRUE, row.names = 1)
pbmc1@meta.data$labels <- labels$majority_voting

pdf("prediction_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc1, reduction = "umap", group.by = "labels")+ labs(title = "UMAP after prediction")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

pbmc1@meta.data$labels2 <- pbmc1@meta.data$labels
pbmc1@meta.data$labels2[pbmc1@meta.data$labels2 != "Basal"] <- "Other"

pdf("PDF/COPD/Basal_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc1, reduction = "umap", group.by = "labels2",cols =c("red", "grey"))+ labs(title = "Basal Cell in UMAP(COPD)")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000)
pbmc1 <- ScaleData(pbmc1)
pbmc1 <- RunPCA(pbmc1,features = VariableFeatures(object = pbmc1))
pbmc1 <- FindNeighbors(pbmc1, dims = 1:10)
pbmc1 <- FindClusters(pbmc1, resolution = 0.5)

pdf("start_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc2, reduction = "umap")+ labs(title = "UMAP before prediction")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

diff3 <- FindMarkers(pbmc1, min.pct = 0.25, 
                    logfc.threshold = 0.25,
                    group = 'labels3',
                    ident.1 ="COPDBasal",
                    ident.2="HealthBasal")

degdf <- diff3
degdf$symbol <- rownames(diff3)
degdf <- subset(degdf, p_val_adj != 0)
logFC_t=0.5
degdf$Sig = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < -0.5,"down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 0.5,"up","stable"))
library(ggrepel)
library(dplyr)

top10 <- degdf %>%
  filter(avg_log2FC > 0.5 | avg_log2FC < -0.5) %>%
  arrange(p_val_adj) %>%
  head(10)

ggplot(degdf, aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point(alpha = 0.4, size = 2.5, aes(color = Sig)) +
  ylab("-log10(Pvalue)") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(logFC_t, -logFC_t), lty = 4, col = "black", lwd = 0.8) +
  theme_bw() +
  geom_text_repel(data = top10, aes(label = symbol))+
  ggtitle("Volcanic map of basal differential genes(COPD)")+
  theme(plot.title = element_text(hjust = 0.5))


VlnPlot(pbmc1, group.by = 'labels', features = 'KRT14') + ggtitle("KRT14(COPD)")


VlnPlot(pbmc1, group.by = 'labels', features = 'COL17A1') + ggtitle("COL17A1(COPD)")


gene=rownames(diff[diff$avg_log2FC > 0.5&diff$p_val_adj <0.05,])
gene=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                   keys = gene_up1,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))


## 把SYMBOL改为ENTREZID
gene_up1=rownames(diff3[diff3$avg_log2FC > 0,])
gene_down1=rownames(diff3[diff3$avg_log2FC < 0,])
library(org.Hs.eg.db)
gene_up1=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                   keys = gene_up1,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))
gene_down1=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                     keys = gene_down1,
                                                     columns = 'ENTREZID',
                                                     keytype = 'SYMBOL')[,2]))
library(clusterProfiler)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
## 以上调基因为例，下调基因同理
## KEGG
gene_up1 <- unique(gene_up1)
kk.up <- enrichKEGG(gene = gene_up,
                    organism = "hsa",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
dotplot(kk.up)

gene_down <- unique(gene_down)
kk.down <- enrichKEGG(gene = gene_down,
                    organism = "hsa",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
dotplot(kk.dow)



go.up <- enrichGO(gene = gene_up,
                OrgDb = org.Hs.eg.db,
                ont = "ALL" ,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.99,
                qvalueCutoff = 0.99,
                readabl = TRUE)
barplot(go.up,showCategory=10,drop=T)


go.down <- enrichGO(gene = gene_down,
                OrgDb = org.Hs.eg.db,
                ont = "ALL" ,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.99,
                qvalueCutoff = 0.99,
                readabl = TRUE)
barplot(go.down,showCategory=13,drop=T)
################################################
#IPF
setwd("/SGRNJ03/pipline_test/luoyichen/IPF")
library(Seurat)
samples=list.files("./")
dir =samples
names(dir) = samples
scRNAlist <- list()
for(i in 1:length(dir)){
counts <- Read10X(data.dir = dir[i])
scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells=1)
}
scRNA<- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], 
                scRNAlist[[4]]，scRNAlist[[5]]，scRNAlist[[6]]))

for (i in 1:length(scRNAlist)) {
    scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
    scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
}
##以VariableFeatures为基础寻找锚点，运行时间较长
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
##利用锚点整合数据，运行时间较长
scRNA <- IntegrateData(anchorset = scRNA.anchors)


immune.combined <- scRNA
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
pbmc <- immune.combined
DefaultAssay(pbmc) <- "RNA"
pbmc[["percent.mt"]] <-PercentageFeatureSet(pbmc,pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 20)
DimPlot(object = pbmc, reduction = "umap", group.by ='orig.ident')
saveRDS(pbmc, file = "IPF.rds")

pbmc <- readRDS("IPF2.rds")
pbmc@meta.data$orig.ident <- gsub("IPF[1-3]", "IPF",pbmc@meta.data$orig.ident)
pbmc@meta.data$orig.ident <- gsub("healthy[1-3]", "Health",pbmc@meta.data$orig.ident)
library(ggplot2)

pdf("IPF_3X3_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc, reduction = "umap")+ labs(title = "UMAP of 6 samples")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("IPF_1X1_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc, reduction = "umap", group.by = "orig.ident")+ labs(title = "UMAP of Disease")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


labels <- read.csv("outs_ipf2/n4-test_cell.annot.csv", header = TRUE, row.names = 1)
pbmc@meta.data$labels <- labels$majority_voting

pdf("PDF/IPF/1X1_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc, reduction = "umap", group.by = "orig.ident")+ labs(title = "UMAP for samples(IPF)")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("PDF/IPF/prediction_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc, reduction = "umap", group.by = "labels")+ labs(title = "UMAP after prediction(IPF)")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

pbmc@meta.data$labels2 <- pbmc@meta.data$labels
pbmc@meta.data$labels2[pbmc@meta.data$labels2 != "Basal"] <- "Other"

pdf("PDF/IPF/Basal Cell in UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc, reduction = "umap", group.by = "labels2",cols =c("red", "grey"))+ labs(title = "Basal Cell in UMAP(IPF)")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("PDF/IPF/start_UMAP.pdf",width = 7, height = 7)
DimPlot(object = pbmc2, reduction = "umap")+ labs(title = "UMAP before prediction(IPF)")+
theme(plot.title = element_text(hjust = 0.5))
dev.off()

pbmc@meta.data$labels3 <- paste0(pbmc@meta.data$orig.ident,pbmc@meta.data$labels2)
diff <- FindMarkers(pbmc, min.pct = 0.25, 
                    logfc.threshold = 0.25,
                    group = 'labels2',
                    ident.1 ="Basal",
                    ident.2="Other")

degdf <- diff
degdf$symbol <- rownames(diff)
logFC_t=1
P.Value_t = 0.05
degdf$Sig = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < -1,"down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 1,"up","stable"))
library(ggrepel)
library(dplyr)

top20 <- degdf %>%
  filter(avg_log2FC > 1 | avg_log2FC < -1) %>%
  arrange(p_val_adj) %>%
  head(20)

ggplot(degdf, aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point(alpha = 0.4, size = 2.5, aes(color = Sig)) +
  ylab("-log10(Pvalue)") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.5) +
  geom_vline(xintercept = c(logFC_t, -logFC_t), lty = 4, col = "black", lwd = 0.5) +
  theme_bw() +
  geom_text_repel(data = top20, aes(label = symbol))+
  ggtitle("Volcanic map of basal differential genes(IPF)")+
  theme(plot.title = element_text(hjust = 0.5))

features <- c("FHL2","IL33","PLAU","KRT5","HCAR3","CYP24A1","FBLN1")


VlnPlot(pbmc, group.by = 'labels', features = 'KRT5') + ggtitle("KRT5(IPF)")


## 获取上下调基因
gene_up=rownames(diff[diff$avg_log2FC > 0.5&diff$p_val_adj <0.05,])
gene_down=rownames(diff[diff$avg_log2FC < -0.5&diff$p_val_adj <0.05,])

library(org.Hs.eg.db)
gene_up=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                   keys = gene_up,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))
gene_down=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                     keys = gene_down,
                                                     columns = 'ENTREZID',
                                                     keytype = 'SYMBOL')[,2]))

library(clusterProfiler)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
## 以上调基因为例，下调基因同理
## KEGG
gene_up <- unique(gene_up)
kk.up <- enrichKEGG(gene = gene_up,
                    organism = "hsa",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
dotplot(kk.up)
gene_down <- unique(gene_down)
kk.down <- enrichKEGG(gene = gene_down,
                    organism = "hsa",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
dotplot(kk.down,,showCategory=5)



go.up <- enrichGO(gene = gene_up,
                OrgDb = org.Hs.eg.db,
                ont = "ALL" ,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.99,
                qvalueCutoff = 0.99,
                readabl = TRUE)
barplot(go,showCategory=13,drop=T)


go.down <- enrichGO(gene = gene_down,
                OrgDb = org.Hs.eg.db,
                ont = "ALL" ,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.99,
                qvalueCutoff = 0.99,
                readabl = TRUE)
barplot(go.down,showCategory=13,drop=T)


##GSEA
nrDEG = diff[,c('avg_log2FC', 'p_val')]
colnames(nrDEG)=c('log2FoldChange','pvalue') ##更改列名
library(org.Hs.eg.db)
library(clusterProfiler)
## 把SYMBOL转换为ENTREZID，可能有部分丢失
gene <- bitr(rownames(nrDEG),     
             fromType = "SYMBOL",     
             toType =  "ENTREZID",    
             OrgDb = org.Hs.eg.db)
## 基因名、ENTREZID、logFC一一对应起来
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))]
## 构建genelist
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
geneList=sort(geneList,decreasing = T) # 降序，按照logFC的值来排序
## GSEA分析
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 10,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
kk_gse=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
sortkk<-kk_gse[order(kk_gse$enrichmentScore, decreasing = T),]
library(enrichplot)
gseaplot2(kk_gse, 
          "hsa04510", 
          color = "firebrick",
          rel_heights=c(1, .2, .6))