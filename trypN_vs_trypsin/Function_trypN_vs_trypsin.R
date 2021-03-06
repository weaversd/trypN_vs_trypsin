#load libraries
library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)

cognate_pairs <- function(peptide_file, evidence_file, residue1, residue2, output_file = NULL, verbose=FALSE){
  peptides <- read_excel(peptide_file)
  evidence <- read_excel(evidence_file, col_types = c("text", 
                                                      "numeric", "text", "text", "text", 
                                                      "text", "text", "text", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", "numeric", 
                                                      "numeric", "text", "text", "text", "text", 
                                                      "text", "text", "numeric", "numeric", 
                                                      "numeric", "numeric", "text", "text", 
                                                      "text", "text", "text", "text", "text", 
                                                      "numeric", "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", "numeric", 
                                                      "text", "numeric", "numeric", "numeric", 
                                                      "numeric", "text", "text", "text", "text", 
                                                      "numeric", "numeric", "text", "text", 
                                                      "numeric", "numeric", "numeric", "text", 
                                                      "numeric", "text", "numeric", "numeric", 
                                                      "text", "numeric", "text", "numeric", 
                                                      "numeric"))
  if (verbose !=FALSE){
    print("loaded files")
  }
  
  #divide peptides into two dataframes
  #if it starts with K and the next AA after is K, then it is a trypN peptide
  #if it ends with K and the firsta AA before is K, then it is a trypsin peptide
  trypN_peptides <- peptides[peptides$`First amino acid` == residue1 & peptides$`Amino acid after` == residue2,]
  trypsin_peptides <- peptides[peptides$`Last amino acid` == residue2 & peptides$`Amino acid before` == residue1,]
  
  if (verbose !=FALSE){
    print("made two dataframes")
  }
  
  #find the "base" peptide, which is the sequence without the leading K (for trypN) or the trailing K (for trypsin)
  trypN_peptides$base_sequence <- substr(trypN_peptides$Sequence, 2, nchar(trypN_peptides$Sequence))
  trypsin_peptides$base_sequence <- substr(trypsin_peptides$Sequence, 1, nchar(trypsin_peptides$Sequence)-1)

  #create new dataframe that joins the previous two based on "base" sequence. Retains only base sequences that were in both tables
  valid_peptide_combination <- inner_join(trypsin_peptides, trypN_peptides, by = "base_sequence",
                                    suffix = c(".trypsin", ".trypN"))
  if (verbose !=FALSE){
    print("joined dataframes")
  }
  
  #create columns for calculating average Retention time and std dev of retention time
  valid_peptide_combination$trypN_Retention_time <- 0.0
  valid_peptide_combination$trypsin_Retention_time <- 0.0
  
  valid_peptide_combination$trypN_RT_stdev <- 0.0
  valid_peptide_combination$trypsin_RT_stdev <- 0.0
  
  if (verbose !=FALSE){
    print("created columns")
  }
  
  #loop to calculate average and stddev of RT for each peptide
  for (i in 1:nrow(valid_peptide_combination)){
    
    #for trypsin, set the sequence and subset the DF to only include that sequence
    seq_trypsin <- valid_peptide_combination$Sequence.trypsin[i]
    evidence_trypsin <- evidence[evidence$Sequence == seq_trypsin,]
    #calculate the average and add it to the original DF
    trypsin_RT <- mean(evidence_trypsin$`Retention time`)
    valid_peptide_combination$trypsin_Retention_time[i] <- trypsin_RT
    #calcualte the stdev and add it to the original DF
    trypsin_RT_stdev <- sd(evidence_trypsin$`Retention time`)
    valid_peptide_combination$trypsin_RT_stdev[i] <- trypsin_RT_stdev
  
    
    #same as above but for trypN
    seq_trypN <- valid_peptide_combination$Sequence.trypN[i]
    evidence_trypN <- evidence[evidence$Sequence == seq_trypN,]
    trypN_RT <- mean(evidence_trypN$`Retention time`)
    valid_peptide_combination$trypN_Retention_time[i] <- trypN_RT
    trypN_RT_stdev <- sd(evidence_trypN$`Retention time`)
    valid_peptide_combination$trypN_RT_stdev[i] <- trypN_RT_stdev
  }

  
  if (verbose !=FALSE){
    print("finished loop")
  }
  
  minimum_trypsin_RT <- min(valid_peptide_combination$trypsin_Retention_time)
  minimum_trypN_RT <- min(valid_peptide_combination$trypN_Retention_time)
  minimum_RT <- min(minimum_trypN_RT, minimum_trypsin_RT)
  
  max_trypsin_RT <- max(valid_peptide_combination$trypsin_Retention_time)
  max_trypN_RT <- max(valid_peptide_combination$trypN_Retention_time)
  max_RT <- max(max_trypN_RT, max_trypsin_RT)
  
  #plot. Each point is a base peptide, with the x axis being the RT for the trypsin version, and the y axis being the RT for the TrypN version
  plot <- ggplot(valid_peptide_combination, aes(x = trypsin_Retention_time, y = trypN_Retention_time)) +
    geom_point(size = 3) + 
    geom_abline(linetype = "dashed") +
    xlim(minimum_RT - (0.05*max_RT), max_RT + (0.05*max_RT)) +
    ylim(minimum_RT - (0.05*max_RT), max_RT + (0.05*max_RT)) +
    theme_bw(base_size = 20) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    #geom_smooth(method = "lm", se = FALSE, formula=y~x-1, fullrange = TRUE) +
    labs(x = paste0("Retention time of 'peptide+", residue2), y = paste0("Retention time of ", residue1, "+peptide"))
  show(plot)
  if (!is.null(output_file)){
    write_xlsx(valid_peptide_combination, output_file)
    print(paste0("Saved table to excel file: ", output_file))
  }
  
  if (verbose !=FALSE){
    print("completed")
  }
  
  return(valid_peptide_combination)
}


RK_cleavage_pairs <- function(peptide_file, evidence_file, point_col="residues", output_file=NULL){
  k_pairs <- cognate_pairs("peptides.xlsx", "evidence.xlsx", "K", "K")
  r_pairs <- cognate_pairs("peptides.xlsx", "evidence.xlsx", "R", "R")
  rk_pairs <- cognate_pairs("peptides.xlsx", "evidence.xlsx", "R", "K")
  kr_pairs <- cognate_pairs("peptides.xlsx", "evidence.xlsx", "K", "R")
  
  
  k_pairs$cleave_residues <- "KK"
  r_pairs$cleave_residues <- "RR"
  rk_pairs$cleave_residues <- "RK"
  kr_pairs$cleave_residues <- "KR"
  
  k_pairs$cognate <- "cognates"
  r_pairs$cognate <- "cognates"
  rk_pairs$cognate <- "false cognates"
  kr_pairs$cognate <- "false cognates"
  
  all_peptides <- rbind(k_pairs, r_pairs, rk_pairs, kr_pairs)
  
  min_tryp_rt <- min(all_peptides$trypsin_Retention_time)
  min_trypN_rt <- min(all_peptides$trypN_Retention_time)
  min_rt <- min(min_trypN_rt, min_tryp_rt)
  
  max_tryp_rt <- max(all_peptides$trypsin_Retention_time)
  max_trypN_rt <- max(all_peptides$trypN_Retention_time)
  max_rt <- max(max_trypN_rt, max_tryp_rt)
  
  if (point_col == "residues"){
    p <- ggplot(all_peptides, aes(trypsin_Retention_time, trypN_Retention_time)) + 
      geom_point(aes(color = cleave_residues), size = 3) +
      geom_abline(linetype = "dashed") +
      xlim(min_rt - (0.05*max_rt), max_rt + (0.05*max_rt)) +
      ylim(min_rt - (0.05*max_rt), max_rt + (0.05*max_rt)) +
      theme_bw(base_size = 20) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "RT of Trypsin Peptide", y = "RT of TrypN Peptide", color = "Cleavage\nResidues")
    show(p)
  } else {
    p <- ggplot(all_peptides, aes(trypsin_Retention_time, trypN_Retention_time)) + 
      geom_point(aes(color = cognate), size = 3) +
      geom_abline(linetype = "dashed") +
      xlim(min_rt - (0.05*max_rt), max_rt + (0.05*max_rt)) +
      ylim(min_rt - (0.05*max_rt), max_rt + (0.05*max_rt)) +
      theme_bw(base_size = 20) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "RT of Trypsin Peptide", y = "RT of TrypN Peptide", color = "Peptide\nRelationship")
    show(p)
  }
  
  if (!is.null(output_file)){
    write_xlsx(valid_peptide_combination, output_file)
    print(paste0("Saved table to excel file: ", output_file))
  }
  
  return(all_peptides)
}
