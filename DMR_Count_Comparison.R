# DMR_Count_Comparison
# Script to read in result files from DiMmer from runs on different datasets 
# using different methods (Linear Regression, Repeated Measures ANOVA,
# Friedman Test), write corresponding number of detected DMRs to a 
# file (DMR_statistics_BA.csv) and plot a comparative bar plot for the
# DMR counts

require(data.table)
require(ggplot2)
require(ggpattern)

# which pvalue type is used, saved and plotted
pval_type <- "original" #original fdr
dmr_statistics_path <- "/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/R_Dimmer/Thesis_Results/DMR_statistics_BA.csv"
plot_path <- "/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/R_Dimmer/Thesis_Results/"
write_file <- FALSE # if new dmr statistics file should be written

### Functions ####
write_to_file <- function(dataset, method, dmrfile){
  dmr = data.table::fread(dmrfile, header = TRUE, sep = ",")
  dmr_statistics <<- rbindlist(list(dmr_statistics, data.table("dataset" = dataset, method = method, "p-value-type"= pval_type, "count" = nrow(dmr))))         
}

# adapt files when pvalue type fdr (or the other way round original)
if(write_file){
  # Read in dmr tables and write/update dmr statistic file
  dmr_statistics <- data.table("dataset" = character(),    # Create an empty data.table
                               "method" = factor(),
                               "p-value-type" = factor(),
                               "count" = double())
  
  d <- "1-EX"
  write_to_file(dataset=d, method="Regression", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE220928_SkeletalMuscleTraining/output_HL_Reg/230731094728728_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="ANOVA", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE220928_SkeletalMuscleTraining/output_HL_ANOVA/230731100721721_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="Friedman", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE220928_SkeletalMuscleTraining/output_HL_Friedman/23073110460979_DMRSearch_org/DMRs.csv")
  
  d <- "2-AGE"
  write_to_file(dataset=d, method="Regression", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE62219_5yearsAfterBirth/output6_regr/230821121652852_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="ANOVA", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE62219_5yearsAfterBirth/output5_anova/230821122336836_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="Friedman", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE62219_5yearsAfterBirth/output4_friedman/230821123217817_DMRSearch_org/DMRs.csv")
  
  d <- "3-CA"
  write_to_file(dataset=d, method="Regression", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE157273_cfDNAMethylomePatterns_ProstateCancer/output_Reg/230821124557857_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="ANOVA", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE157273_cfDNAMethylomePatterns_ProstateCancer/output_Anova/230821010345845_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="Friedman", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE157273_cfDNAMethylomePatterns_ProstateCancer/output_Friedman/23082101160585_DMRSearch_org/DMRs.csv")
  
  d <- "4-GEN"
  write_to_file(dataset=d, method="Regression", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE176394_GenderAffirmingHT/output3/output_Reg/230821014717817_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="ANOVA", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE176394_GenderAffirmingHT/output3/output_Anova/230821021655855_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="Friedman", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE176394_GenderAffirmingHT/output3/output_Friedman/23082102330888_DMRSearch_org/DMRs.csv")
  
  d <- "5-LUP"
  write_to_file(dataset=d, method="Regression", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE161476_lupusPatients/output4/output_Reg/23082104160686_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="ANOVA", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE161476_lupusPatients/output4/output_Anova/23082104380282_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="Friedman", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE161476_lupusPatients/output4/output_Friedman/230821124451851_DMRSearch_org/DMRs.csv")
  
  d <- "6-T1D"
  write_to_file(dataset=d, method="Regression", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE142512_methylationDifferencesPrecedeType1diabetes/output_Reg/230801045342842_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="ANOVA", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE142512_methylationDifferencesPrecedeType1diabetes/output_Anova/230801051455855_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="Friedman", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE142512_methylationDifferencesPrecedeType1diabetes/output_Friedman/230801061054854_DMRSearch_org/DMRs.csv")
  
  d <- "7-TNF"
  write_to_file(dataset=d, method="Regression", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE191297_TherapyResponse_InflammatoryBowelDisease/output_Reg/230905025854954_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="ANOVA", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE191297_TherapyResponse_InflammatoryBowelDisease/output_Anova/230821032223823_DMRSearch_org/DMRs.csv")
  write_to_file(dataset=d, method="Friedman", dmrfile="/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE191297_TherapyResponse_InflammatoryBowelDisease/output_Friedman/230821041132832_DMRSearch_org/DMRs.csv")
  
  write.table(
    dmr_statistics,
    file = dmr_statistics_path,
    row.names = FALSE,
    col.names = TRUE,
    sep = ","
  )
}

### Or read in dmr statistics file ####
dmr_statistics <- data.table::fread(dmr_statistics_path)

### Plot Bar plot counts of dmrs between all methods and dataset

# Select datasets that should be plotted (in case you don't wan't to plot all)
dmr_statistics_selection = dmr_statistics[dmr_statistics$dataset %in% c("1-EX","2-AGE","3-CA","4-GEN","5-LUP","6-T1D","7-TNF"),]
# only selected pvalue type is plotted
dmr_statistics_selection = dmr_statistics_selection[dmr_statistics_selection$`p-value-type`==pval_type,]

# Plot Barplot
ggplot(dmr_statistics_selection,                                     
       aes(x = dataset,
           y = as.numeric(count),
           fill = method)) +
  geom_bar(stat = "identity", colour = "white",
           position= position_dodge(width = 0.5)) +
  ylab("DMR Count") +
  xlab("Dataset ID") +
  scale_x_discrete(limits = c("1-EX","2-AGE","3-CA","4-GEN","5-LUP","6-T1D","7-TNF")) +
  scale_fill_manual(name ="Methods", values= c(Regression="#619CFF", ANOVA="#F8766D",Friedman="#00BA38")) +
  ggtitle(paste("DMR Count Comparison (P-value type for CpGs: ",pval_type,")",sep=""))
ggsave(paste0(plot_path,"DMRs_comparison_allDatasets_",pval_type,"_BA.png"))

# to avoid transformation error when count = 0
dmr_statistics_selection[dmr_statistics_selection$count==0,]$count <- NA 

# Plot Logscale Barplot
ggplot(dmr_statistics_selection,                                     
       aes(x = dataset,
           y = as.numeric(count),
           fill = method)) +
  geom_bar(stat = "identity", colour = "white", 
           position= position_dodge(width = 0.5)) +
  ylab("DMR Count") +
  xlab("Dataset ID") +
  scale_y_continuous(trans='log10') +
  scale_x_discrete(limits = c("1-EX","2-AGE","3-CA","4-GEN","5-LUP","6-T1D","7-TNF")) +
  scale_fill_manual(name ="Methods", values= c(Regression="#619CFF", ANOVA="#F8766D",Friedman="#00BA38")) + 
  ggtitle(paste("DMR Count Comparison (P-value type for CpGs: ",pval_type,")",sep=""))
ggsave(paste0(plot_path,"DMRs_comparison_allDatasets_logscale_",pval_type,"_BA.png"))
