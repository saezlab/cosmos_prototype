library(readr)
library(vsn)
library(limma)
library(viper)
library(pheatmap)

psites_1 <- as.data.frame(read_csv("data/phospho/final_phoshpo_raw.csv")) #from 01_cleaning_limma_cell_comparison.R
# write_csv(psites_1,"~/Dropbox/COSMOS/data/phosphoproteomic_raw.csv")
row.names(psites_1) <- psites_1$ID
psites_1 <- psites_1[,-1]

targets_1 <- as.data.frame(matrix(NA,length(psites_1[1,]),2))
names(targets_1) <- c("sample","condition")
targets_1$sample <- names(psites_1)
targets_1$condition <- gsub("[0-9]*","",targets_1$sample)


# write_csv(targets_1, "~/Dropbox/COSMOS/support/phosphoproteomic_targets.csv")

sub_df1 <- psites_1[,c("38KI","38TU","15KI","15TU","29KI","29TU","16KI","16TU","32KI","32TU","35KI","35TU")]
sub_df2 <- psites_1[,c("40KI","40TU","24KI","24TU","11KI","11TU")]
sub_df1 <- sub_df1[rowSums(is.na(sub_df1)) < (length(sub_df1[1,])-1),]
sub_df2 <- sub_df2[rowSums(is.na(sub_df2)) < (length(sub_df2[1,])-1),]

sub_df_list <- list(sub_df1,sub_df2)

for(i in 1:length(sub_df_list))
{
  df <- sub_df_list[[i]]
  fit <- vsnMatrix(as.matrix(df))
  meanSdPlot(fit)
  df <- as.data.frame(vsn::predict(fit,as.matrix(df)))
  sub_df_list[[i]] <- df
}

psites_1_vsn <- as.data.frame(do.call(merge,c(sub_df_list, all = T, by = "row.names")))
row.names(psites_1_vsn) <- psites_1_vsn$Row.names
psites_1_vsn <- psites_1_vsn[,-1]
psites_1_vsn <- psites_1_vsn[rowSums(is.na(psites_1_vsn)) < (length(psites_1_vsn[1,])-3),]

psites_1_vsn_bcorrect <- as.data.frame(removeBatchEffect(psites_1_vsn, batch = c(rep("A",2),rep("B",6), rep("C",4), rep("D",6)))) #Based on correlation heatmap and PCA

write.csv(psites_1_vsn_bcorrect, "data/phospho/final_subset_phoshpo_processed.csv")

targets_1 <- targets_1[targets_1$sample %in% names(psites_1),]
write_csv(targets_1, "support/phospho_targets_final.csv")

limmaRes <- runLimma(psites_1_vsn_bcorrect, targets_1, comparisons = list(c(2,-1)))

ttop <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 14243, adjust.method = "fdr"))
write_csv(ttop, "results/phospho/final_samples/ttop_tumorVsHealthy.csv")

dev.off()
volcano_nice(ttop, FCIndex = 2,pValIndex = 5, IDIndex = 1, nlabels = 30, label = T)



############

url <- paste0(
  'http://omnipathdb.org/ptms?',
  'fields=sources,references&genesymbols=1'
)

download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

omnipath_ptm <- download_omnipath()
write_csv(omnipath_ptm, "support/omnipath_ptm_20200205.csv")
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

KSN_viper <- df_to_viper_regulon(KSN)

#Prepare the measurments for viper
viper_expression <- ttop$t
names(viper_expression) <- ttop$ID
###This is just to handle a bug in viper
viper_expression <- as.data.frame(viper_expression)
viper_expression$dummy <- viper_expression$viper_expression
###COmpute TF activity scores
Kin_activity <- as.data.frame(viper(eset = viper_expression, regulon = KSN_viper, minsize = 5, adaptive.size = F, eset.filter = F))
kinases <- row.names(Kin_activity) 
Kin_activity <- as.data.frame(Kin_activity[,-1])
row.names(Kin_activity) <- kinases
names(Kin_activity) <- "NES"

Kin_pvals <- pnorm(Kin_activity[,1])
Kin_pvals <- ifelse(Kin_pvals > 0.5, 1-Kin_pvals, Kin_pvals)
Kin_pvals <- Kin_pvals * 2
Kin_pvals <- as.data.frame(Kin_pvals)
row.names(Kin_pvals) <- kinases
names(Kin_pvals) <- "pvalue"

top_Kin_Activity <- Kin_activity
top_Kin_Activity$kinase <- row.names(top_Kin_Activity)
top_Kin_Activity <- top_Kin_Activity[order(abs(top_Kin_Activity[,1]), decreasing = T),]
top_Kin_Activity <- as.data.frame(top_Kin_Activity[1:31,])
top_kinase <- top_Kin_Activity[,2]
top_Kin_Activity <- as.data.frame(top_Kin_Activity[,-2])
row.names(top_Kin_Activity) <- top_kinase
names(top_Kin_Activity) <- "NES"

t <- as.vector(t(top_Kin_Activity))
palette1 <- createLinearColors(-1*(t[t > 0]),withZero = F , maximum = 37)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = 63)
palette <- c(palette1,palette2)

pheatmap(top_Kin_Activity, display_numbers = T, cluster_cols = F, color = palette, cluster_rows = F, border_color = F)


top_Kin_pvals <- Kin_pvals
top_Kin_pvals$kinase <- row.names(top_Kin_pvals)
top_Kin_pvals <- top_Kin_pvals[order(abs(top_Kin_pvals[,1]), decreasing = F),]
top_Kin_pvals <- as.data.frame(top_Kin_pvals[1:20,])
top_kinase <- top_Kin_pvals[,2]
top_Kin_pvals <- as.data.frame(top_Kin_pvals[,-2])
row.names(top_Kin_pvals) <- top_kinase
names(top_Kin_pvals) <- "pval"

pheatmap(top_Kin_pvals, display_numbers = T, cluster_cols = F, cluster_rows = F, border_color = F)

write.csv(Kin_activity, "results/phospho/final_samples/kinase_activities.csv")
