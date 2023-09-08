# DMR_Time_Course
# Script to plot the time course of CpGs contained in highly significant DMRs
# A section at the end is used to plot the beta-value density per region type of
# a dataset

require(tidyr)
require(forcats)
require(data.table)
require(GenomicRanges)
require(dplyr)
require(ggplot2)
# color palettes
require(ggsci)
require(viridis)

### Functions ####
# returns CpG name lists for CpGs overlapping with given DMR ranges per DMR 
# using get_CpGs_from_manifest function
getOverlappingCpGs <- function(this_DMR_ranges){
  overlaps <- findOverlaps(this_DMR_ranges,manifest_ranges)
  this_CpGs <- as.data.table(overlaps) %>% 
    group_by(queryHits) %>% 
    group_map(~ get_CpGs_from_manifest(.x$subjectHits))
  return(this_CpGs)
}
# Build merged datatable from betamatrix and annotation file for given CpGs
build_dt <- function(this_CpGs){
  beta_matrix_signCpGs <- beta_matrix[beta_matrix$CpG_ID %in% this_CpGs,]
  rownames(beta_matrix_signCpGs) <- beta_matrix_signCpGs$CpG_ID
  beta_matrix_signCpGs$CpG_ID <- NULL
  
  annotation_data$x <- paste(annotation_data$Sentrix_ID, '_', annotation_data$Sentrix_Position, sep = "")
  rownames(annotation_data) <- annotation_data$x
  genomic_idx <- match(colnames(beta_matrix_signCpGs), rownames(annotation_data))
  annotation_data  <- annotation_data[genomic_idx,] # just in case its a diff ordering
  
  # Transpose matrix to match annotation rotation and keep row/colnames
  transpose_matrix <- transpose(beta_matrix_signCpGs)
  rownames(transpose_matrix) <-  colnames(beta_matrix_signCpGs)
  colnames(transpose_matrix) <- rownames(beta_matrix_signCpGs)
  transpose_matrix$x <- rownames(transpose_matrix)
  
  # Build merged datatable from betamatrix and annotation file and preprocess
  dt <- merge.data.table(annotation_data, transpose_matrix, by = 'x', all = TRUE) # , by='row.names', all=TRUE
  rownames(transpose_matrix) <- dt$x
  dt$x <- NULL
  dt$Sentrix_ID <- NULL
  dt$Sentrix_Position <- NULL
  
  return(dt)
}
# returns CpG names corredponing to idxs in manifest
get_CpGs_from_manifest <- function(idxs){
  return(manifest_ranges$IlmnID[idxs])
}

### Load Data ####
path <- "/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/R_Dimmer/Thesis_Results/"
dataset <- "6-T1D"
method <- "Friedman" # Regression ANOVA Friedman

beta_matrix_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE142512_methylationDifferencesPrecedeType1diabetes/output_Friedman/beta_matrix.csv'
annotation_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE142512_methylationDifferencesPrecedeType1diabetes/sample_annotation4_case.csv'
dmr_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE142512_methylationDifferencesPrecedeType1diabetes/output_Friedman/230801061054854_DMRSearch_org/merged_table.csv'
manifestEPIC_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/DimmerJava/src/main/resources/epic_manifest.csv' # to get the CpGs contained in/ overlapping with a DMR
manifest450k_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/DimmerJava/src/main/resources/manifest_summary.csv' # to get the CpGs contained in/ overlapping with a DMR

beta_matrix <- data.table::fread(beta_matrix_file)
annotation_data <- data.table::fread(annotation_file)
# change this according to the array type of the dataset
manifest <- data.table::fread(manifestEPIC_file) 

# Get the CpGs contained in the highest scoring DMR (or selected DMR)
manifest_pos <- data.table(Chr=manifest$CHR, POS=manifest$MAPINFO, IlmnID=manifest$IlmnID, regionType=manifest$Relation_to_UCSC_CpG_Island)
manifest_pos$Chr <- paste("chr", manifest_pos$Chr, sep="")
manifest_ranges <- makeGRangesFromDataFrame(manifest_pos,start.field = "POS", end.field = "POS", keep.extra.columns=TRUE) # .bed format
manifest_ranges

# produce Granges Object for ANOVA and Friedman DMRs
dmr <- data.table::fread(dmr_file)
dmr$Chr <- paste("chr", dmr$Chr, sep="")

# change this according to the dmrs you want to plot
dmr_sign <- head(dmr[order(`p-value`, decreasing=FALSE),],1)

dmr_range <- makeGRangesFromDataFrame(dmr_sign[,1:3],start.field = "Begin", end.field = "End") # .bed format

dmrs_CpGs <- unlist(getOverlappingCpGs(dmr_range))
dmrs_CpGs_filtered <- dmrs_CpGs[dmrs_CpGs %in% beta_matrix$CpG_ID]

# Get methylation and annotation values for those CpGs
dt <- build_dt(unlist(dmrs_CpGs_filtered))

long_DT <- dt %>% gather("CpG", "value", unlist(dmrs_CpGs))

#sum(is.na(long_DT))

breaks_start = min(dt$timestamp)
breaks_end = max(dt$timestamp)
step = (breaks_end - breaks_start) / (length(unique(dt$timestam)) - 1)

if (method == "Regression"){
  my_colors <- rev(pal_material(palette = "blue")(length(unique(long_DT$PatientNr))))
}else if(method == "ANOVA"){
  my_colors <- rev(pal_material(palette = "red")(length(unique(long_DT$PatientNr))))
}else{
  my_colors <- rev(pal_material(palette = "light-green")(length(unique(long_DT$PatientNr))))
}

g <- ggplot(long_DT,aes(x = timestamp, y = value, colour = factor(PatientNr))) +
  geom_line() + ylim(0.0,1.0) + scale_x_continuous(breaks=seq(breaks_start, breaks_end, step)) + scale_color_manual(values = my_colors) + #seq(breaks_start, breaks_end, step) c(3,12,24,48,60)
  geom_point() + xlab("Timestamps") + ylab("Beta Value") + guides(color=guide_legend("Patient No.")) +
  #scale_colour_viridis(option= "inferno", discrete = TRUE) + #option= "turbo", discrete = TRUE
  #scale_color_futurama() + 
  ggtitle(paste("Time Course of Top",method,"DMR - Dataset",dataset)) +
  facet_wrap(~CpG, nrow = 2, ncol = 7) #change according to length of DMR
g
ggsave(paste0(path,dataset,"/TimeCourse_",method,"_TopDMR.png"),
       width=9,height=5,units = "in") #change according to length of DMR

# 14 CpGs: width=9,height=6
# 8 CpGs: width=10,height=3
# 6 CpGs: width=8,height=3

### Plot Denisity per Region for this Dataset ####
manifest_dt <- data.table("CpGs" = manifest_ranges$IlmnID, "regionType" = manifest_ranges$regionType)

dt_allCpGs[1:10,1:10]
dt_allCpGs <- build_dt(beta_matrix$CpG_ID)
dt_allCpGs <- dt_allCpGs[dt_allCpGs$timestamp==3] # only take the basline beta_values
dt_allCpGs$Group_ID <- NULL
dt_allCpGs$Gender_ID <- NULL
dt_allCpGs$PatientNr <- NULL
dt_allCpGs$timestamp <- NULL
dt_allCpGs_withoutPatient <- dt_allCpGs[, parallel::mclapply(.SD, FUN= function(x){
  mean(x,na.rm = TRUE)}, mc.cores = as.numeric(4))] # average over the patients
dt_allCpGs_withoutPatient[1:3,1]
dt_allCpGs_withoutPatient <- transpose(dt_allCpGs_withoutPatient)
dt_allCpGs_withoutPatient$CpGs <- colnames(dt_allCpGs)
joined <- inner_join(x = dt_allCpGs_withoutPatient, y = manifest_dt, by = "CpGs")
joined$group <- joined$regionType
joined$group[joined$group == 'N_Shore'] <- 'Shore'
joined$group[joined$group == 'S_Shore'] <- 'Shore'
joined$group[joined$group == 'N_Shelf'] <- 'Shelf'
joined$group[joined$group == 'S_Shelf'] <- 'Shelf'
joined$group[joined$group == ''] <- 'Other'

ggplot(joined, aes(x = V1, y = after_stat(density), colour = group)) +
  geom_density(lwd = 1.2, linetype = 1) + scale_color_npg() + xlab("Beta Values") + 
  ylab("Density") + labs(colour = "Region Type") +ggtitle("Beta Value Density per Region Type") + 
  theme(text=element_text(size=16),axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(paste0(path,dataset,"/BetaValueDensity.png"))

