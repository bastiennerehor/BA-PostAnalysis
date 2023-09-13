# LOLA_plus_Overlap_Comparison
# Script to perform Genomic Regions Enrichment with LOLA to compare enrichments
# for detected DMRs between Methods (Linear Regression, Repeated Measures ANOVA,
# Friedman Test). If plot is set to TRUE it also produces Venn, Euler and Upset
# plots to determine which detected DMRs overlap between methods and which are
# uniquely found

require(forcats)
require(data.table)
require(dplyr)
require(ggplot2)
require(LOLA)
require(simpleCache)
require(GenomicRanges)
require(seqsetvis)
require(ComplexHeatmap)
require(ggsci)

path_plots <- "/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/R_Dimmer/Thesis_Results/"
dataset <- '1-EX'
plot <- FALSE # determines if overlap plots should be produced

### Read in dmr tables and Preprocess it to get LOLA input ####
regfile <- "/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE220928_SkeletalMuscleTraining/output_HL_Reg/230731094728728_DMRSearch_org/merged_table.csv"
dmr_reg <- data.table::fread(regfile, header = TRUE, sep = ",")
anovafile <- "/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE220928_SkeletalMuscleTraining/output_HL_ANOVA/230731100721721_DMRSearch_org/merged_table.csv"
dmr_anova <- data.table::fread(anovafile, header = TRUE, sep = ",")
friedTfile <- "/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE220928_SkeletalMuscleTraining/output_HL_Friedman/23073110460979_DMRSearch_org/merged_table.csv"
dmr_friedman <- data.table::fread(friedTfile, header = TRUE, sep = ",")

# Build Universe from all DMRs found in Dimmer Search
dmr_reg$Chr <- paste("chr", dmr_reg$Chr, sep="")
dmr_reg_ranges <- makeGRangesFromDataFrame(dmr_reg[,1:3],start.field = "Begin", end.field = "End") # .bed format
dmr_anova$Chr <- paste("chr", dmr_anova$Chr, sep="")
dmr_anova_ranges <- makeGRangesFromDataFrame(dmr_anova[,1:3],start.field = "Begin", end.field = "End") # .bed format
dmr_friedman$Chr <- paste("chr", dmr_friedman$Chr, sep="")
dmr_friedman_ranges <- makeGRangesFromDataFrame(dmr_friedman[,1:3],start.field = "Begin", end.field = "End") # .bed format

combined_DMR_Ranges <- GRangesList(dmr_reg_ranges, dmr_anova_ranges, dmr_friedman_ranges)
names(combined_DMR_Ranges) <- c('dmr_reg','dmr_anova','dmr_friedman')
userUniverse <- buildRestrictedUniverse(combined_DMR_Ranges)

### Test for overlaps (with the found DMRs of the study, of specific DMRs (TimeCourse etc.) ####
# file <- "/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/R_Dimmer/Thesis_Results/6-T1D/studyDMRs.csv"
# dmr_study <- data.table::fread(file, header = TRUE ,sep = ",")
# rownames(dmr_study) <- dmr_study$DMR
# dmr_study$DMR <- NULL
# dmr_study_ranges <- makeGRangesFromDataFrame(dmr_study[,1:3],start.field = "start", end.field = "end") # .bed format
# 
# hits <- overlapsAny(dmr_reg_ranges, dmr_study_ranges)
# hits
# dmr_reg_ranges[hits]
# 
# input <- GRangesList(Anova = dmr_anova_ranges,
#                      Friedman = dmr_friedman_ranges,
#                      Study = dmr_study_ranges)
# 
# olaps = ssvOverlapIntervalSets(input)
# 
# ssvFeatureVenn(olaps, circle_colors=c("#F8766D","#00BA38","#CC79A7")) + 
#   guides(fill = "none", color = "none")
# ggsave(paste0(path_plots,dataset,"/OverlapComparison/VennDiagram_withStudy.png"))
#####
# Build UserSets from all - significant - DMRs found in Dimmer Search
dmr_reg <- dmr_reg[dmr_reg$`p-value`<=0.01,]
dmr_reg_ranges <- makeGRangesFromDataFrame(dmr_reg[,1:3],start.field = "Begin", end.field = "End")
dmr_anova <- dmr_anova[dmr_anova$`p-value`<=0.01,]
dmr_anova_ranges <- makeGRangesFromDataFrame(dmr_anova[,1:3],start.field = "Begin", end.field = "End")
dmr_friedman <- dmr_friedman[dmr_friedman$`p-value`<=0.01,]
dmr_friedman_ranges <- makeGRangesFromDataFrame(dmr_friedman[,1:3],start.field = "Begin", end.field = "End")

userSets <- GRangesList(dmr_reg_ranges, dmr_anova_ranges, dmr_friedman_ranges)
userSets
names(userSets) <- c('dmr_reg','dmr_anova','dmr_friedman')

## Plot the Overlap of significant DMRs ####
if(plot){
  # create combined GRanges Object with logical membership table 
  # ("ranges withing 2 * ext (=0) of one another will be joined during the merge" 
  # meaning they count as an overlap/ shared DMR)
  olaps = ssvOverlapIntervalSets(userSets)
  # Plot Venn Diagram
  ssvFeatureVenn(olaps, circle_colors=c("#619CFF", "#F8766D","#00BA38")) + 
    guides(fill = "none", color = "none")
  ggsave(paste0(path_plots,dataset,"/OverlapComparison/VennDiagram.png"))
  # Plot Euler Diagram
  ssvFeatureEuler(olaps, circle_colors=c("#619CFF", "#F8766D","#00BA38")) + 
    guides(fill = "none", color = "none")
  ggsave(paste0(path_plots,dataset,"/OverlapComparison/EulerDiagram.png"))
  
  # GRangesList with pretty names
  input <- GRangesList(Regression = dmr_reg_ranges,
                       Anova = dmr_anova_ranges,
                       Friedman = dmr_friedman_ranges)
  # create Combination Matrix that indicates how many BasePairs (in total/ sum) 
  # are included for all intersection & non-intersection possibilities
  m = make_comb_mat(input)
  upset <- ComplexHeatmap::UpSet(m, pt_size = unit(5, "mm"), lwd = 3,
                                 comb_col = c("orange", "DeepSkyBlue", "black")[comb_degree(m)],
                                 comb_order = order(comb_size(m)),right_annotation = upset_right_annotation(m,
                                                                                                            gp = gpar(fill = c("#619CFF", "#F8766D","#00BA38"),col=c("#619CFF", "#F8766D","#00BA38")),
                                                                                                            width = unit(4, "cm")))
  png(paste0(path_plots,dataset,"/OverlapComparison/UpSetDiagram.png"), width=600, height=300)
  upset
  dev.off()
  pdf(paste0(path_plots,dataset,"/OverlapComparison/UpSetDiagram.pdf"), width=8, height=4)
  upset
  dev.off()
}
### Continue LOLA ####

# necessary to perform LOLArun iff user Sets have overlaps between sets to split 
# in overlapping and non overlapping range
userSets <- redefineUserSets(userSets, userUniverse, cores = 1)

# Get regionDB (hg19)
regionDB <- loadRegionDB("/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/R_Dimmer/LOLACore_data/LOLACore/hg19")

### Run LOLA Enrichment ####
locResults <- runLOLA(userSets, userUniverse, regionDB, cores=1)
#unique(locResults[locResults$collection=="ucsc_features"]$antibody)
locResults[locResults$userSet=="dmr_anova"][1:20,15:18] #dmr_friedman dmr_anova dmr_reg

# Check which column causes Fisher’s exact test to assign an Inf Odds Ratio
# by being 0 
# -> its the b column meaning all universe sets overlapping with the TFBS/ 
#    reference database can also be found in the query set so there is no region 
#    of interest not covered by the query set
# => therefor it seems reasonable to set the Odds ration to the max value of all
#    Odds ratios (+1 for distinction)
checkInf_oddsr <- locResults[!is.finite(locResults$oddsRatio)]
colSums(checkInf_oddsr==0)

# Get Significant enrichment results
pval <- 0.01

locResults_Significant <- locResults[locResults$pValueLog>=-log10(pval),]

### Visualize LOLA Results ####

locResults_Significant_Copy <- locResults_Significant[,1:5] # get relevant columns
# Set Inf Odds ratio values to the max value of all Odds ratios (+1 for distinction)
locResults_Significant_Copy$oddsRatio[!is.finite(locResults_Significant_Copy$oddsRatio)] <- (max(locResults_Significant_Copy$oddsRatio[is.finite(locResults_Significant_Copy$oddsRatio)]) + 1)
# make max and min pValue columns to plot ranges for pvalues iff a set overlaps 
# with multiple database reference sets that correspond with the same antibody/TFBS
locResults_Significant_Copy$maxpValueLog <- locResults_Significant$pValueLog
locResults_Significant_Copy$minpValueLog <- locResults_Significant$pValueLog
locResults_Significant_Copy$antibody <- locResults_Significant$antibody # get antibody column
# used to safe how many overlaps with database reference sets that 
# correspond with the same antibody/TFBS are detected and summarized in the next step
locResults_Significant_Copy$count <- rep(1,nrow(locResults_Significant_Copy))
locResults_Significant_Copy

# summarize the data to get only one occurrence of an antibody per user set
# in the process the max, min and mean pvalue as well as the mean Odds Ratio
# for that antibody is saved and its counted how many occurences were summarized
locResults_Significant_ALL <- locResults_Significant_Copy %>%
  group_by(userSet) %>%
  group_by(antibody,collection,userSet) %>%
  summarise(across(maxpValueLog, max), across(minpValueLog, min),across(count, sum),across(oddsRatio, mean), across(pValueLog, mean)) 
locResults_Significant_ALL <- na.omit(locResults_Significant_ALL)
locResults_Significant_ALL

# Plot enrichment comparison between methods; log scale x-axis corresponding 
# to Pvalue, Odds Ratio is captured in point/ dot size, grouped by
# collection
p <- ggplot(locResults_Significant_ALL, aes(y = fct_reorder(antibody, collection), x = pValueLog)) +
  geom_errorbarh(aes(xmin = minpValueLog, xmax = maxpValueLog), height = 0.25) +
  geom_point(aes(fill= collection,size = oddsRatio), pch = 21) +
  labs(size="Odds Ratio", fill= "Collection") +
  facet_grid(collection ~ userSet,scales = "free_y", space= "free", 
             labeller = as_labeller(c(dmr_reg = 'Regression', dmr_anova='ANOVA',dmr_friedman = 'Friedman', cistrome_cistrome='',
                                      cistrome_epigenome='',
                                      codex='',
                                      encode_tfbs=''))) + 
  guides(fill = guide_legend(override.aes = list(size=4))) +
  #scale_fill_viridis(option= "inferno", discrete = TRUE) +
  scale_fill_futurama(labels = c(cistrome_cistrome="Cistrome database", cistrome_epigenome="Epigenome database (from Cistrome)", codex="Codex database", encode_tfbs="TFBS from ENCODE")) +  # fill statt color _d for discrete
  theme_minimal() +
  xlab(expression('-log'[10]*'(p-value)')) + 
  ylab("TFBS/ Antibody") +
  ggtitle(paste("LOLA Enrichment - Dataset ",dataset," - Pvalue Cutoff: ",signif(pval, digits=3), sep="")) +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=16),
        axis.text.y = element_text(size=8)) # , legend.text=element_text(size=12),

# Change vertical spacing between facets (optional)
p <- p + theme(panel.spacing.x = unit(1, "lines"))

p
# adapt size if need be
ggsave(paste0(path_plots,dataset,'/LOLA/LOLAEnrichmentResultsCompared_BA.png'), width = 12, height = 12, units = "in")
p
# adapt size if need be
ggsave(paste0(path_plots,dataset,'/LOLA/LOLAEnrichmentResultsCompared_BA.pdf'), width = 12, height = 12, units = "in")
