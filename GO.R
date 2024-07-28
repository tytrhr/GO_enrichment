library(clusterProfiler)
library(dplyr)
library(ggpubr)
library(tidyr)

dt <- read.table("go_annotation.txt",header = T)
go_bac <- data.frame(term = dt$term, gene = dt$id)
go_ter <- go2term(go_bac$term)
gene <- read.table("gene_id.txt", header = T)


go_enrich <- function(x){
  enrichment = enricher(x, TERM2GENE = go_bac, TERM2NAME = go_ter, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
  write.table(enrichment, quote = F, row.names = F, sep = "\t",file = "go_tmp.txt")
  
  enrichment_fst <- read.table("go_tmp.txt", header = T, sep = "\t", stringsAsFactors = FALSE,quote = "")
  data <- subset(enrichment_fst,enrichment_fst$Description!="NA")
  enrich_fst <- subset(data, pvalue < 1)
  
  go_des <- go2ont(enrich_fst$ID)
  
  dt_go <- na.omit(data.frame(id=enrich_fst$ID, des=enrich_fst$Description, pval=-log10(enrich_fst$pvalue), count=enrich_fst$Count, ont=go_des$Ontology[match(enrich_fst$ID, go_des$go_id)]))
  
  dt_go2 <- na.omit(data.frame(id=enrich_fst$ID, des=enrich_fst$Description, GeneRatio=enrich_fst$GeneRatio, BgRatio=enrich_fst$BgRatio, pvalue=enrich_fst$pvalue, padj=enrich_fst$p.adjust, qvalue=enrich_fst$qvalue, gene=enrich_fst$geneID, count=enrich_fst$Count, ont=go_des$Ontology[match(enrich_fst$ID, go_des$go_id)]))
  write.table(dt_go2, quote = F, row.names = F, sep = "\t",file = "dt_go2.txt")
  
  dt_go$ont <- factor(dt_go$ont, levels = c("MF", "CC", "BP"))
  dt_sort <- arrange(dt_go, ont, pval)
  p = ggbarplot(dt_go, y= "pval", x = "des", 
                fill = "ont",
                color = "white",
                palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                sort.val = "asc",
                sort.by.groups = TRUE,
                rotate = TRUE,
                xlab = "",
                ylab = "-log10(P-value)",
                label = dt_sort$count,
                lab.pos = "out",
                lab.col = "black",
                lab.hjust = -2,
                lab.vjust = .5
  )
  return(p)
}
go_enrich(gene$gene)



