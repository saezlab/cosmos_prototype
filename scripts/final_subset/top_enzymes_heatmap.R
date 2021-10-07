#These override function of Pheatmap just allow to have tilted labels.
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

library(readr)
library(pheatmap)
library(grid)
library(ggplot2)

source("scripts/support_functions.R")

##KINASE
kinase_activities <- as.data.frame(
  read_csv("results/phospho/final_samples/kinase_activities.csv"))


##TF
TF_scores <- as.data.frame(
  read_csv("results/transcriptomic/final_subset/TF_scores.csv"))

names(TF_scores) <- c("X1","NES")
kinase_activities <- as.data.frame(rbind(kinase_activities, TF_scores))

top_enzymes <- kinase_activities[order(abs(kinase_activities$NES), decreasing = T),]
top_enzymes <- top_enzymes[1:30,]
row.names(top_enzymes) <- top_enzymes[,1]

top_enzymes[] <- top_enzymes[,2]
top_enzymes <- as.data.frame(t(top_enzymes))
top_enzymes <- top_enzymes[-1,]

t <- as.vector(t(top_enzymes))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = 45)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = 55)
palette <- c(palette1,palette2)

pheatmap(top_enzymes, display_numbers = T, cluster_cols = F, color = palette, cluster_rows = F, border_color = "white")

#######
top_enzymes_long <- as.data.frame(t(top_enzymes))
top_enzymes_long$enzyme <- row.names(top_enzymes_long)
top_enzymes_long$enzyme <- factor(top_enzymes_long$enzyme, levels = unique(top_enzymes_long$enzyme))
top_enzymes_long$sign <- ifelse(top_enzymes_long$NES > 0, "pos", "neg")
top_enzymes_long$NES <- abs(top_enzymes_long$NES)

ggplot(top_enzymes_long, aes(x = enzyme, y = NES, fill = sign)) + 
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_manual(values=c("pos" = "#FF5733", "neg" = "#5C90FF")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45)) + ylab("abs(NES)") + xlab("TF/Kinase/Phosphotase")
