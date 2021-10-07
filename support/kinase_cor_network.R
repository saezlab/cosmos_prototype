kinase_cor_network <- function(kinase_activity_df, KSN_filtered)
{
  cor_kinase_activities <- as.data.frame(cor(t(kinase_activity_df), method = "spearman"))
  
  cor_kinase_network <- list()
  k <- 1
  for(i in 1:length(cor_kinase_activities[,1]))
  {
    for(j in 1:i)
    {
      # if(abs(cor_kinase_activities[i,j]) > 0.8 & abs(cor_kinase_activities[i,j]) != 1)
      if(abs(cor_kinase_activities[i,j]) != 1)
      {
        cor_kinase_network[[k]] <- c(row.names(cor_kinase_activities)[i],names(cor_kinase_activities)[j],cor_kinase_activities[i,j])
        # cor_kinase_network[k,1] <- row.names(cor_kinase_activities)[i] 
        # cor_kinase_network[k,2] <- names(cor_kinase_activities)[j]
        # cor_kinase_network[k,3] <- cor_kinase_activities[i,j]
        k <- k+1
      }
    }
  }
  
  cor_kinase_network <- as.data.frame(do.call(rbind,cor_kinase_network))
  names(cor_kinase_network) <- c("kinase_A","kinase_B","correlation")  
  cor_kinase_network$correlation <- as.numeric(as.character(cor_kinase_network$correlation))
  
  enzymes <- as.character(unique(KSN_filtered$enzyme))
  
  kinase_targets_normalised_intersesction <- list()
  k <- 1
  for(i in 1:(length(enzymes)-1))
  {
    enzyme_A <- enzymes[i]
    for(j in (i+1):(length(enzymes)))
    {
      enzyme_B <- enzymes[j]
      targets_A <- KSN_filtered[KSN_filtered$enzyme == enzyme_A,1]
      targets_B <- KSN_filtered[KSN_filtered$enzyme == enzyme_B,1]
      
      normalised_intersection <- sum(targets_A %in% targets_B)/min(length(targets_A), length(targets_B))
      
      kinase_targets_normalised_intersesction[[k]] <- c(enzyme_A,enzyme_B, normalised_intersection)
      k <- k+1
    }
  }
  
  kinase_targets_normalised_intersesction_df <- as.data.frame(do.call(rbind,kinase_targets_normalised_intersesction))
  names(kinase_targets_normalised_intersesction_df) <- c("kinase_A","kinase_B","normalised_intersection")  
  kinase_targets_normalised_intersesction_df$normalised_intersection <- as.numeric(as.character(kinase_targets_normalised_intersesction_df$normalised_intersection))
  
  plot(hist(kinase_targets_normalised_intersesction_df$normalised_intersection, breaks = 100))
  
  kinase_targets_normalised_intersesction_df_inverse <- kinase_targets_normalised_intersesction_df
  names(kinase_targets_normalised_intersesction_df_inverse) <- c("kinase_B","kinase_A","normalised_intersection")
  
  kinase_targets_normalised_intersesction_df <- rbind(kinase_targets_normalised_intersesction_df, kinase_targets_normalised_intersesction_df_inverse)
  
  cor_kinase_network_normalised <- merge(cor_kinase_network, kinase_targets_normalised_intersesction_df, by = c("kinase_A","kinase_B"))
  cor_kinase_network_normalised$normalised_correlation <- cor_kinase_network_normalised$correlation*(1 - cor_kinase_network_normalised$normalised_intersection)
  
  return(cor_kinase_network_normalised)
}