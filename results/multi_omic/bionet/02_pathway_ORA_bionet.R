load("~/Dropbox/kidney_cancer_multiomic_pipeline/results/multi_omic/bionet/module_list.Rdata")

library(readr)
library(piano)
library(biomaRt)
library(omicToolsTest)
library(pheatmap)

all_interactions <- as.data.frame(do.call(rbind,module_list))
all_interactions$edgeID <- paste(all_interactions$V1, all_interactions$V2, sep = "___")

edge_counts <- table(all_interactions$edgeID)

all_interactions$count <- 0

for(i in 1:length(all_interactions[,1]))
{
  all_interactions[i,"count"] <- as.numeric(edge_counts[all_interactions[i,"edgeID"]])/length(module_list)
}

unique_interactions <- unique(all_interactions)
plot(hist(unique_interactions$count))

top_confident <- unique_interactions[unique_interactions$count >= 0.95,]

confident_genes <- unique(c(as.character(top_confident$V1), as.character(top_confident$V2)))
confident_genes <- confident_genes[!grepl("Metab",confident_genes)]
confident_genes <- gsub("Gene.+__","",confident_genes)
confident_genes <- unique(confident_genes)

######

meta_network <- as.data.frame(read_csv("Dropbox/Meta_PKN/result/meta_network.csv"))
meta_network <- meta_network[complete.cases(meta_network),]

background <- unique(c(as.character(meta_network$source),as.character(meta_network$target)))
background <- background[!grepl("Metab",background)]
background <- gsub("Gene.+__","",background)
background <- unique(background)
background <- background[grepl("^[0-9]+$",background)]

######

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

G_list <- getBM(filters = "entrezgene", 
                attributes = c("ensembl_peptide_id",'hgnc_symbol','entrezgene', "description"),
                values = background, mart = ensembl)

G_list <- unique(G_list[,c(2,3)])

mapping_vec <- G_list[,1]
names(mapping_vec) <- G_list[,2]

for(i in 1:length(confident_genes))
{
  confident_genes[i] <- mapping_vec[confident_genes[i]]
}
for(i in 1:length(background))
{
  background[i] <- mapping_vec[background[i]]
}

confident_genes <- unique(confident_genes)
background <- unique(background)

gene_set <- gmt_to_csv("~/Dropbox/kidney_cancer_multiomic_pipeline/support/c2.cp.v6.2.symbols.gmt")
gene_set <- gene_set[grepl("KEGG",gene_set[,2]),]
gene_set <- loadGSC(gene_set)

pathway_scores <- runGSAhyper(confident_genes, universe = background, gsc = gene_set, adjMethod = "fdr")
pathway_scores_df <- as.data.frame(pathway_scores$resTab) 
pathway_scores_df_top <- pathway_scores_df[pathway_scores_df$`Adjusted p-value` < 0.05,]

to_heatmap <- pathway_scores_df_top
to_heatmap<- as.data.frame(to_heatmap$`Adjusted p-value`)
row.names(to_heatmap) <- gsub("KEGG_","",row.names(pathway_scores_df_top))
names(to_heatmap) <- "tumor"

pheatmap(to_heatmap, cluster_cols = F, display_numbers = T)
