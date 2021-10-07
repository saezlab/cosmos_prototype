library(readr)
library(limma)
library(vsn)
library(readxl)
library(ggplot2)
library(ggrepel)

source("scripts/support_functions.R")
#Import data

raw_metabolomic <- as.data.frame(read_excel("data/metabolomic/raw_metabolomic.xlsx", 
                              na = "NF"))

metab_id_correspondance <- as.data.frame(read_delim("support/metab_id_correspondance", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE))

names(raw_metabolomic) <- c("ID", metab_id_correspondance$ID)
names(raw_metabolomic) <- gsub(" T","TU", names(raw_metabolomic))
names(raw_metabolomic) <- gsub(" H","KI", names(raw_metabolomic))

row.names(raw_metabolomic) <- raw_metabolomic[,1]
raw_metabolomic <- raw_metabolomic[,-1]

samples_to_keep <- c(paste(c("11", "15", "16", "24", "29", "32", "38", "40", "35"), "TU", sep = ""), paste(c("11", "15", "16", "24", "29", "32", "38", "40", "35"), "KI", sep = ""))
samples_to_keep <- samples_to_keep[samples_to_keep %in% names(raw_metabolomic)]

raw_metabolomic <- raw_metabolomic[,samples_to_keep]
min(raw_metabolomic, na.rm = T)
max(raw_metabolomic, na.rm = T)

targets_1 <- as.data.frame(matrix(NA,length(raw_metabolomic[1,]),2))
names(targets_1) <- c("sample","condition")
targets_1$sample <- names(raw_metabolomic)
targets_1$condition <- gsub("[0-9]*","",targets_1$sample)

to_write <- raw_metabolomic
to_write$metabolite <- row.names(to_write)
# write_csv(to_write[,c(17,1:16)],"data/raw_metabolomic.csv")
# write_csv(targets_1,"data/metabolomic_targets.csv")

fit <- vsnMatrix(as.matrix(raw_metabolomic))
meanSdPlot(fit)
raw_metabolomic_vsn <- as.data.frame(vsn::predict(fit,as.matrix(raw_metabolomic)))

limmaRes <- runLimma(raw_metabolomic_vsn, targets_1, comparisons = list(c(1,-2)))
ttop_tumour_vs_healthy <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(raw_metabolomic_vsn[,1]), adjust.method = "fdr"))
volcano_nice(ttop_tumour_vs_healthy, FCIndex = 2,pValIndex = 5, IDIndex = 1, nlabels = 15, label = T)

write_csv(ttop_tumour_vs_healthy, "results/metabolomic/ttop_tumour_vs_healthy.csv")

ttop_metabolomic_OvHK2 <- as.data.frame(read_csv("data/metabolomic/ttop_metabolomic_OvHK2.csv"))

ttop_tumour_vs_healthy$ID <- tolower(ttop_tumour_vs_healthy$ID)
ttop_metabolomic_OvHK2$ID <- tolower(ttop_metabolomic_OvHK2$ID)

t_table <- merge(ttop_tumour_vs_healthy[,c(1,4)], ttop_metabolomic_OvHK2[,c(1,4)], by = "ID")

cor.test(t_table$t.x, t_table$t.y)
