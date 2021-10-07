plot_enzyme_target_network <- function(ttop, sif, enzyme_scores, enzyme, ntarget = 10)
{
  names(sif) <- c("target","source","sign")
  enzyme_df <- sif[sif[,2] == enzyme,]
  names(enzyme_df)[1] <- "ID" 
  enzyme_df <- merge(enzyme_df, ttop[,c(1,4)])
  enzyme_df$t <- enzyme_df$t * enzyme_df$sign
  enzyme_df <- enzyme_df[order(abs(enzyme_df$t), decreasing = T),]
  enzyme_df <- enzyme_df[c(1:ntarget),]
  
  
  
  nodes <- enzyme_df[,c(1,4)]
  names(nodes) <- c("id","t")
  nodes$t <- as.numeric(as.character(nodes$t))
  nodes$id <- as.character(nodes$id)
  nodes <- as.data.frame(rbind(nodes,c("id" = enzyme,"t" = enzyme_scores[enzyme_scores[,1] == enzyme,2])))
  nodes$label <- nodes$id
  nodes$color.background <- ifelse(nodes$t > 0, "red","blue")
  nodes$color.border <- "black"
  nodes$font.size = "20"
  
  edges <- enzyme_df[,c(2,1,3)]
  names(edges) <- c("from","to","sign")
  edges$arrows <- "to"
  
  network <- visNetwork(nodes = nodes, edges = edges, width="100%", height="1600px")
  
  enzyme_df <- sif[sif[,2] == enzyme,]
  names(enzyme_df)[1] <- "ID" 
  enzyme_df <- merge(enzyme_df, ttop)
  enzyme_df[,4] <- enzyme_df[,4] * enzyme_df$sign
  volcano <- volcano_nice(enzyme_df, FCIndex = 4, pValIndex = 7, IDIndex = 1)
  
  return(list(network, volcano))
}

volcano_nice <- function(df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL)
{
  df <- df[complete.cases(df), ]
  names(df)[1] <- "X1"
  hAssOri <- hAss
  hAss <- -log10(hAss)
  names(df) <- gsub("adj.P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- ceiling(max(abs(df[, FCIndex])))
    ylimAbs <- ceiling(max(abs(-log10(df[, pValIndex]))))
  }
  else {
    xlimAbs <- max(abs(df[, FCIndex]))
    ylimAbs <- max(abs(-log10(df[, pValIndex])))
  }
  if (is.null(vAss)) {
    vAss <- xlimAbs/10
  }
  xneg <- function(x) abs(hAss - 1 + x/(x + vAss))
  xpos <- function(x) abs(hAss - 1 + x/(x - vAss))
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log10(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss) {
        if (xpos(x) < -log10(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  
  df$couleur <- "0"
  df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),
                                                    as.numeric(x[pValIndex]), vAss))
  df <- df[order(df$adj.P.Val, decreasing = F), ]
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2",
                       df$couleur)
  
  a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val),
                      color = couleur)) + geom_point(alpha = 0.5, size = 10) +
    stat_function(fun = xneg, xlim = c(-xlimAbs,
                                       -vAss), color = "black", alpha = 0.7) + ylim(c(0,
                                                                                      ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                                                                                                             xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
    scale_colour_manual(values = c("grey30", "red",
                                   "royalblue3")) + 
    theme_minimal() + 
    theme(legend.position = "none",axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + 
    ylab("-log10(pvalue)") + 
    xlab("logFC * sign")
  
  return(a)
}

library(readr)
library(viper)

url <- paste0(
  'http://omnipathdb.org/interactions?',
  'datasets=tfregulons&tfregulons_levels=A,B,C&genesymbols=1&fields=sources,tfregulons_level'
)

download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

##Dorothea/viper
dorothea <- download_omnipath()
dorothea <- dorothea[,c(4,3,6,7)]
dorothea$sign <- dorothea$is_stimulation - dorothea$is_inhibition
dorothea$sign <- ifelse(dorothea$sign == 0, 1, dorothea$sign)
dorothea <- dorothea[,c(1,2,5)]
dorothea <- unique(dorothea)


url <- paste0(
  'http://omnipathdb.org/ptms?',
  'fields=sources,references&genesymbols=1'
)

omnipath_ptm <- download_omnipath()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)
KSN <- unique(KSN)

ttop_tumorvshealthy <- as.data.frame(
  read_csv("results/transcriptomic/final_subset/ttop_tumorvshealthy.csv"))
RNAseq_entrez_to_symbol <- as.data.frame(read_delim("support/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127

names(RNAseq_entrez_to_symbol)[1] <- "ID"

ttop_tumorvshealthy <- merge(ttop_tumorvshealthy, RNAseq_entrez_to_symbol[,c(1,6)])
ttop_tumorvshealthy <- ttop_tumorvshealthy[,c(8,2:7)]
names(ttop_tumorvshealthy)[1] <- "ID"
ttop_tumorvshealthy$ID <- gsub(" .*","",ttop_tumorvshealthy$ID)
ttop_tumorvshealthy <- unique(ttop_tumorvshealthy)
ttop_tumorvshealthy_RNA <- ttop_tumorvshealthy
rm(ttop_tumorvshealthy)

ttop_tumorVsHealthy_phospho <- as.data.frame(
  read_csv("results/phospho/final_samples/ttop_tumorVsHealthy.csv"))

kinase_activities <- as.data.frame(
  read_csv("results/phospho/final_samples/kinase_activities.csv"))
names(kinase_activities)[1] <- "ID"

TF_scores <- as.data.frame(
  read_csv("results/transcriptomic/final_subset/TF_scores.csv"))

enzyme_scores <- as.data.frame(rbind(
  kinase_activities,
  TF_scores
))

library(visNetwork)



test <- plot_enzyme_target_network(ttop = ttop_tumorvshealthy_RNA
                           , sif = dorothea
                           , enzyme_scores = enzyme_scores
                           , enzyme = "STAT2"
                           , ntarget = 10)

test2 <- plot_enzyme_target_network(ttop = ttop_tumorVsHealthy_phospho
                           , sif = KSN
                           , enzyme_scores = enzyme_scores
                           , enzyme = "CDK7"
                           , ntarget = 10)
