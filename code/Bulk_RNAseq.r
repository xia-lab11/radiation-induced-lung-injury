#### Differential gene analysis
## limma
library(limma)
library(tidyverse)
load("./data/expr_sva.Rdata")

counts <- as.data.frame(expr_sva)
group_list <- metadata$group
identical(colnames(counts),metadata$sample)

design <- model.matrix(~0+factor(group_list)) 
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(counts)

cont.matrix <- makeContrasts(contrasts=paste0("radiation",'-',"control"), 
                             levels = design)

fit <- lmFit(counts,design)
fit2 <- contrasts.fit(fit, cont.matrix) 
fit2 <- eBayes(fit2) 
DiffEG<-topTable(fit2, coef=1, n=Inf) %>% na.omit()


#### Visualization
library(tidyverse)
library(pheatmap)
library(viridisLite)
library(ggplot2)
library(ggrepel) 

exprset2 <- expr_sva
metadata <- metadata
Group <- factor(metadata$group)

diffSig <- DiffEG[with(DiffEG, (abs(logFC)>0.5 & P.Value < 0.05 )), ]
diffUp <- DiffEG[with(DiffEG, (logFC>0.5 & P.Value < 0.05 )), ]
diffDown <- DiffEG[with(DiffEG, (logFC<(-0.5) & P.Value < 0.05 )), ]

## pheatmap
heatdata <- exprset2[rownames(diffSig),]
annotation_col <- data.frame(Group)
rownames(annotation_col) <- colnames(heatdata)

pheatmap(heatdata, 
         cluster_rows = T,
         treeheight_row = 0,
         treeheight_col = 30,
         cluster_cols = T,
         annotation_col =annotation_col, 
         annotation_legend=TRUE, 
         show_rownames = F,
         scale = "row", 
         color = colorRampPalette(c("#007E58", "white", "#BF3512"))(100),  
         cellwidth = 18,
         cellheight = 0.5,
         fontsize = 10)

## volcano plot
data <- DiffEG
log2_FC = 0.5 
PValue = 0.05 

data <- data %>%  mutate(Expression = case_when(logFC >= log2_FC & P.Value <= PValue ~ "Up-regulated",
                                                logFC <= -log2_FC & P.Value <= PValue ~ "Down-regulated",
                                                TRUE ~ "No significant"))
data <- data %>% 
  mutate(
    Significance = case_when(
      abs(logFC) >= log2_FC & P.Value <= 0.05 & P.Value > 0.01 ~ "P.Value 0.05", 
      abs(logFC) >= log2_FC & P.Value <= 0.01 & P.Value > 0.001 ~ "P.Value 0.01",
      abs(logFC) >= log2_FC & P.Value <= 0.001 ~ "P.Value 0.001", 
      TRUE ~ "No significant")
  )

v2 <- ggplot(data, aes(logFC, -log(P.Value,10))) +
  geom_point(aes(color = Expression), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"P.Value")) +
  scale_color_manual(values=c("#09A372", "#d2dae2","#EF601F"))+
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  geom_vline(xintercept = c(0.5,-0.5),lty = "dashed",size = 0.8) +
  geom_hline(yintercept = -log(0.05,10),lty = "dashed",size = 0.8) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=15,face="bold"),
        axis.title.y=element_text(size=15,angle=90,face="bold"),
        axis.text.y=element_text(size=12,face="bold"),
        axis.text.x=element_text(size=12,face="bold"),
        panel.grid=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size =12, face = "bold"),
        legend.key = element_rect(fill = "white"));v2


#### Enrichment
library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(enrichplot)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)
library(org.Mm.eg.db)

GO_database <- 'org.Mm.eg.db' 
KEGG_database <- 'mmu'

DEG <- DiffEG %>% filter(abs(logFC)>0.5) %>% filter(P.Value<0.05)
DEG$change <- ifelse(DEG$logFC > 0,"up","down")
table(DEG$change)

DEG_all_mouse <- bitr(rownames(DEG),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
DEG_up_mouse <-bitr(rownames(DEG[DEG$change == "up",]),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
DEG_down_mouse <- bitr(rownames(DEG[DEG$change == "down",]),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)


## GO
GO_ORA<-enrichGO(DEG_all_mouse$ENTREZID,
              OrgDb = GO_database,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable = T)
GO_ORA_results <- as.data.frame(GO_ORA)

## KEGG
KEGG_ORA<-enrichKEGG(DEG_all_mouse$ENTREZID,
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
KEGG_ORA = setReadable(KEGG_ORA,
                       OrgDb = "org.Mm.eg.db",
                       keyType = "ENTREZID")
KEGG_ORA_results <- as.data.frame(KEGG_ORA)

## Chord diagram
KEGG_ORA_results$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", KEGG_ORA_results$Description)
pathways <- c("Neutrophil extracellular trap formation","Cytokine-cytokine receptor interaction","Chemokine signaling pathway","C-type lectin receptor signaling pathway",
              "IL-17 signaling pathway","Leukocyte transendothelial migration","ECM-receptor interaction")
KEGGplotIn<-KEGG_ORA_results[KEGG_ORA_results$Description %in% pathways,c(3,4,11,13)] 
KEGGplotIn$geneID <-str_replace_all(KEGGplotIn$geneID,'/',',') 
names(KEGGplotIn)<-c('ID','Term','adj_pval','Genes')
genedata <- data.frame(genes = rownames(DEG),logFC = DEG$logFC)

df_go_clean <- KEGGplotIn %>% 
  separate_rows(Genes,sep = ",") %>% 
  left_join(genedata,by = c("Genes" = "genes"))
dfPlot = chord_dat(data.frame(df_go_clean),genes = genedata,  process=unique(df_go_clean$Term))

GOChord(dfPlot,
        title="",      
        space = 0.01,              
        gene.order = "logFC",      
        gene.space = 0.25,         
        gene.size = 6,             
        lfc.col=c('firebrick3', 'white','royalblue3'), 
        border.size = 0.05,               
        ribbon.col=brewer.pal(ncol(dfPlot)-1, "Set2"),   
        process.label = 10)

## GSEA
hallmarks_mouse <- read.gmt("MSigDB/mouse/mh.all.v2024.1.Mm.symbols.gmt") 

genelist <- DiffEG$logFC
names(genelist) <-  rownames(DiffEG)
genelist <- sort(genelist,decreasing = T)
head(genelist)

enrich_hallmarks <- GSEA(genelist,TERM2GENE = hallmarks_mouse)
enrich_hallmarks_results <- as.data.frame(enrich_hallmarks)

p1 <- ridgeplot(enrich_hallmarks,
          showCategory = 16,
          fill = "p.adjust", 
          core_enrichment = TRUE,
          label_format = 60,
          orderBy = "NES",
          decreasing = FALSE)+
  scale_fill_gradient(low = "#EB5D5B", high = "#2CB8AE") + 
  theme_bw() + 
  theme(axis.text.y=element_text(size=12,face="bold"),
        axis.text.x=element_text(size=12,face="bold"),
        panel.grid=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size =12, face = "bold"),
        legend.key = element_rect(fill = "white"))+
  xlim(c(-0.2, 2))


#### Immune infiltration
## load data
load("./data/expr_sva.Rdata")
sig_matrix <- read.csv("./data/ImmuCC.csv") %>% 
  column_to_rownames("X")

## CIBERSORT
library(tidyverse)
library(pheatmap)
library(CIBERSORT)
library(IOBR)

exp <- as.data.frame(expr_sva)
table(is.na(exp))
table(exp<0)
mixture_file <- exp[!apply(exp,1,function(x) any(x < 0)),]

res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE,absolute = F)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
cellnum <- res_cibersort[,1:ncol(sig_matrix)]
cell.prop<- apply(cellnum, 1, function(x){x/sum(x)})
data4plot <- data.frame()
for (i in 1:ncol(cell.prop)) {
  data4plot <- rbind(
    data4plot,
    cbind(cell.prop[,i],rownames(cell.prop),
          rep(colnames(cell.prop)[i],nrow(cell.prop)
          )
    )
  )
}
colnames(data4plot)<-c('proportion','celltype','sample')
data4plot$proportion <- as.numeric(data4plot$proportion)
data4plot$sample <- factor(data4plot$sample,levels = c("GSM622201","GSM622203","GSM1024339","GSM1024340","GSM622214","GSM622215","GSM622216","GSM1024344","GSM1024345","GSM1024346"))
ggplot(data4plot,aes(sample,proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=my36colors)+
  ggtitle("cell portation")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.title.x=element_text(size=1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
  guides(fill=guide_legend(title=NULL))

## mMCP-counter
library(mMCPcounter)
library(data.table)
library(tidyverse)

data <- mMCPcounter.estimate(exp, features = c("Gene.Symbol","ENSEMBL.ID","Probes")[1],
                             genomeVersion = c("GCRm38","GCRm39")[1])

rowscale <- data
rowscale <- rowscale[,apply(rowscale, 2, function(x){sum(x)>0})]

annol_col<-data.frame(Group=factor(metadata$group))
row.names(annol_col)=colnames(rowscale)
annol_color<-list(Group=c(control='#66C2A5',radiation="#FC8D62"))

pheatmap(rowscale,
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         show_colnames = F,
         treeheight_row = 0, 
         annotation_col = annol_col,
         annotation_colors = annol_color,
         gaps_col = 4,
         border_color = "black",
         fontsize =10,
         cellwidth =20,
         cellheight = 15,
         color = colorRampPalette(colors = c("#66C2A5","white","#FC8D62"))(100))


#### GSVA analysis of 16 PCD genesets
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
library(ggsignif)
library(ggpubr)

load("./data/expr_sva.Rdata")
load("./data/rcd_genesets.Rdata")

genes_mouse <- rownames(expr_sva)
library(homologene)
genes_human <- homologene(genes_mouse, inTax = 10090, outTax = 9606)

genes_human <- genes_human %>% 
  distinct(`9606`,.keep_all = T)

expr_sva_human <- expr_sva %>% 
  as.data.frame() %>% 
  rownames_to_column("10090") %>% 
  inner_join(genes_human,by = "10090") %>% 
  dplyr::select(-c("10090","9606_ID","10090_ID")) %>% 
  column_to_rownames("9606")

dat <- as.matrix(expr_sva_human)
geneset <- rcd_genesets

gsvaPar <- gsvaParam(dat,geneset)
gsva.es <- gsva(gsvaPar, verbose=FALSE)

diff_path <- data.frame()
for (i in 1:nrow(gsva.es)) {
  diff_path[i,1] <- rownames(gsva.es)[i]
  diff_path[i,2] <- wilcox.test(as.numeric(gsva.es[i,1:4]),as.numeric(gsva.es[i,5:10]))$p.value
  
}

diff_plot <- gsva.es %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(group = c(rep("control",4),rep("radiation",6))) %>% 
  dplyr::select(c("group",everything())) %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(cols = 3:ncol(.),names_to = "death_type",values_to = "value")
diff_plot$death_type <- str_to_title(diff_plot$death_type)

ggplot(diff_plot,aes(death_type, value, fill = group)) + 
  geom_boxplot() + 
  xlab("") + 
  ylab("Score") +
  theme_bw() + 
  theme(axis.title.y=element_text(size=15,angle=90,face="bold"),
        axis.text.y=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1,size=12,face="bold"),
        panel.grid=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.3),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size =12, face = "bold"),
        legend.key = element_rect(fill = "white")) + 
  stat_compare_means(aes(group = group),method = "wilcox.test",label = "p.signif")
















