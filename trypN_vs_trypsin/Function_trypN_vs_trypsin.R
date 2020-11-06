#load libraries
library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)

cognate_pairs <- function(peptide_file, evidence_file, residue, output_file = NULL){
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
  
  print("loaded files")
  
  
  #divide peptides into two dataframes
  #if it starts with K and the next AA after is K, then it is a trypN peptide
  #if it ends with K and the firsta AA before is K, then it is a trypsin peptide
  trypN_peptides <- peptides[peptides$`First amino acid` == residue & peptides$`Amino acid after` == residue,]
  trypsin_peptides <- peptides[peptides$`Last amino acid` == residue & peptides$`Amino acid before` == residue,]
  
  print("made two dataframes")
  
  #find the "base" peptide, which is the sequence without the leading K (for trypN) or the trailing K (for trypsin)
  trypN_peptides$base_sequence <- substr(trypN_peptides$Sequence, 2, nchar(trypN_peptides$Sequence))
  trypsin_peptides$base_sequence <- substr(trypsin_peptides$Sequence, 1, nchar(trypsin_peptides$Sequence)-1)

  #create new dataframe that joins the previous two based on "base" sequence. Retains only base sequences that were in both tables
  valid_peptide_combination <- inner_join(trypsin_peptides, trypN_peptides, by = "base_sequence",
                                    suffix = c(".trypsin", ".trypN"))
  
  print("joined dataframes")
  
  #create columns for calculating average Retention time and std dev of retention time
  valid_peptide_combination$trypN_Retention_time <- 0.0
  valid_peptide_combination$trypsin_Retention_time <- 0.0
  
  valid_peptide_combination$trypN_RT_stdev <- 0.0
  valid_peptide_combination$trypsin_RT_stdev <- 0.0
  
  print("created columns")

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

  print("finished loop")
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
    labs(x = paste0("Retention time of 'peptide+", residue), y = paste0("Retention time of ", residue, "+peptide"))
  show(plot)
  if (!is.null(output_file)){
    write_xlsx(valid_peptide_combination, output_file)
    print(paste0("Saved table to excel file: ", output_file))
  }
  print("completed")
  return(valid_peptide_combination)
}

