#!/usr/bin/env Rscript

library(dplyr)
library(tibble)
library(Seurat)
library(openxlsx)
library(biomaRt)	
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)


workdir <- "~/CLUSTER/rnaseq_course"

counts_file <- file.path(workdir, "counts", "toxoplasma.finalCounts.txt")
expDesign_file <- file.path(workdir, "counts", "toxoplasma.expDesign.txt")

count_data <- read.table(counts_file,
                         header = TRUE,
                         sep = "\t",
                         row.names = 1,
                         check.names = FALSE)

exp_design <- read.table(expDesign_file,
                         header = TRUE,
                         colClasses = c("character", "character"))
colnames(exp_design) <- c("sample", "group")
exp_design <- exp_design %>% column_to_rownames("sample")
exp_design$group <- as.factor(exp_design$group)

colnames(count_data) <- sub("_$", "", colnames(count_data))
exp_design <- exp_design[colnames(count_data),,drop=FALSE]

splits <- unlist(strsplit(as.character(exp_design$group), "_"))
tissues <- splits[seq(1, length(splits), by=3)]
conditions <- splits[seq(3, length(splits), by=3)]
exp_design2 <- data.frame(row.names=rownames(exp_design), tissue=tissues, condition=conditions)
exp_design2 <- exp_design2[colnames(count_data),,drop=FALSE]
diff2 <- DESeqDataSetFromMatrix(count_data, exp_design2, ~ tissue + condition + tissue:condition)
dds2 <- DESeq(diff2, betaPrior=FALSE)


diff <- DESeqDataSetFromMatrix(count_data, exp_design, ~ group)
diff$group <- relevel(diff$group, ref="Lung_WT_Control")

diff <- estimateSizeFactors(diff)

dds <- DESeq(diff, betaPrior=TRUE)

res <- results(dds, contrast=c("group", "Lung_WT_Case", "Lung_WT_Control"), alpha=0.05)
#res_shrunk <- lfcShrink(dds, coef="group_Lung_WT_Case_vs_Lung_WT_Control", res=res, type='apeglm')

res_table <- results(dds, contrast=c("group", "Lung_WT_Case", "Lung_WT_Control"), alpha=0.05, tidy=TRUE) %>%
                column_to_rownames(var="row")

deg <- res_table %>%
       arrange(padj) %>%
       filter(padj <= 0.05)

symbols <- mapIds(org.Mm.eg.db, keys=rownames(res_table), column=c('SYMBOL'), keytype='ENSEMBL')
res$symbol <- symbols[match(rownames(res), names(symbols))]
res <- res_shrunk[!is.na(res_shrunk$symbol),]

EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'padj')


vsd <- vst(dds)
plotPCA(vsd, intgroup=c("group"))


enrich <- enrichGO(gene          = rownames(deg),
                   universe      = rownames(res_table),
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

barplot(enrich, showCategory=20)
dotplot(enrich, showCategory=30) + ggtitle("dotplot for GO enrichment")

edox <- setReadable(enrich, 'org.Mm.eg.db', 'ENSEMBL')
p1 <- cnetplot(enrich, foldChange=res_table$log2FoldChange)
cowplot::plot_grid(p1, ncol=1)



interactionplot <- function(dds, geneid, xaxis) {
  assay(dds) %>%
    as_tibble(rownames = "gene") %>%
    filter(gene == geneid) %>%
    gather(Label, value, -gene) %>%
    select(-gene) -> 
    genedat
  
  colData(dds) %>%
    as.data.frame %>%
    as_tibble(rownames = "sample") %>%
    full_join(genedat, by=c("sample"="Label")) -> genedat
  
  mygeom  <-  geom_point()
  mypal   <- scale_colour_manual(name="", values=brewer.pal(3, "Set1"))
  mytheme <- theme_bw()
  mytitle <- ggtitle(geneid)
  
  xaxis <- rlang::sym(xaxis)
  gp <- ggplot(genedat, aes(x=!!xaxis, y=value, color=tissue)) + mygeom + mytheme + mypal + mytitle
  
  return(gp)   
}




