library(vsn)
library(limma)
library(readr)

source("scripts/support_functions.R")

load("data/transcripto/counts_final_subset.RData")

count_df <- as.data.frame(fc$counts)

diffs <- c()
for(i in seq(1,length(names(count_df)),2))
{
  diffs <- c(diffs,abs((sum(count_df[,i]) - sum(count_df[,i+1]))))
}
plot(hist(diffs)) #cutoff 1.500.000

count_df[count_df == 0] <- NA
column_list <- list()
j <- 1
for(i in seq(1,length(names(count_df)),2))
{
  column_list[[j]] <- rowMeans(count_df[,c(i,i+1)], na.rm = T)
  j <- j+1
}


count_df_integrated <- as.data.frame(do.call(cbind,column_list))
column_names <- gsub("X.Volumes.Kramman_for_mac.RNASeq.bam.RNAseq.","",fc$targets)
column_names <- gsub("_.*","",column_names)
column_names <- unique(column_names)

names(count_df_integrated) <- column_names

to_write <- count_df_integrated
to_write$gene <- row.names(to_write)

# write_csv(to_write[,c(23,1:22)],"~/Dropbox/COSMOS/data/transcriptomic_raw.csv")

count_df_integrated[is.na(count_df_integrated)] <- 0
count_df_integrated <- count_df_integrated[rowMeans(count_df_integrated) > 50,]

count_df_integrated[count_df_integrated == 0] <- 0.5

targets<- as.data.frame(matrix(NA,length(count_df_integrated[1,]),2))
names(targets) <- c("sample","condition")
targets$sample <- names(count_df_integrated)
targets$condition <- gsub("[0-9]*","",targets$sample)

# write_csv(targets,"~/Dropbox/COSMOS/support/transcriptomic_targets.csv")


fit <- vsnMatrix(as.matrix(count_df_integrated))
meanSdPlot(fit)
count_df_integrated_vsn <- as.data.frame(vsn::predict(fit,as.matrix(count_df_integrated)))

limmaRes <- runLimma(count_df_integrated_vsn, targets, comparisons = list(c(2,-1)))

ttop <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 15919, adjust.method = "fdr"))

write_csv(ttop, "results/transcriptomic/final_subset/ttop_tumorvshealthy.csv")

counts_to_write <- count_df_integrated_vsn

RNAseq_entrez_to_symbol <- as.data.frame(read_delim("support/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127

mapping_vec <- gsub(" .*","",RNAseq_entrez_to_symbol$`Gene names`)
names(mapping_vec) <- RNAseq_entrez_to_symbol$`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL`
counts_to_write <- counts_to_write[row.names(counts_to_write) %in% names(mapping_vec),]
for(i in 1:length(counts_to_write[,1]))
{
  row.names(counts_to_write)[i] <- mapping_vec[row.names(counts_to_write)[i]] 
}

counts_to_write$ID <- row.names(counts_to_write)
counts_to_write <- counts_to_write[,c(23,1:22)]

write_csv(counts_to_write, "data/transcripto/counts_vsn.csv")
