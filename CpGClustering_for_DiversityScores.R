# CpGClustering_for_Diversity_Scores
# Script to perform CpG clustering for both the CpGs detected by Dimmer using 
# ANOVA and Friedman Test (for comparison). Afterwards the clustering results
# are used to score DMRs based on the diversity of Clusters assigned to CpGs
# they contain (alpha diversity). These scores can be used to filter DMRs with 
# inhomogeneous trends which might not be as relevant or might be functionally 
# split.
# Beta diversity between DMRs is also calculated based on the clustering tree.
# Every big Section (readin, clustering, diversity scores) is spilt to do 
# everything ANOVA first and the same procedure for Friedman afterwards

require(data.table)
require(dplyr)
require(stats)
require(dtwclust) 
require(ggpubr)
require(ggplot2)
require(ComplexHeatmap)
require(gridtext)
require(GenomicRanges)
require(vegan)
require(picante)
require(phyloseq)
require(ggrepel)
# color palettes
require(adegenet)
require(colorRamp2)
require(ggsci)
require(viridis)

### Functions ####
# function for euclidean distance
euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))
# returns CpG names corredponing to idxs in manifest
get_CpGs_from_manifest <- function(idxs){
  return(manifest_ranges$IlmnID[idxs])
}
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
# Plot CVI scores; input list of chosen scores from your clustering and give the 
# amount of clusters you want to plot; choose   Scorename, Scorenameshort,  
# Methodname accordingly; decide if geom_point() should be used
plotScore <- function(scores,Scorename, Scorenameshort,clusterList,Methodname, points){
  g <- ggplot(data.frame(), aes(x=clusterList, y=scores[clusterList], group= 1)) +
    geom_line() +
    xlab("Number of Clusters") +
    ylab(Scorename) +
    theme(text=element_text(size=16))
    #ggtitle(paste("Internal CVI:",Scorename))
  if(points){
    g <- g + geom_point()
  }
  g
  ggsave(paste0(path,Methodname,"/",Scorenameshort,"_",as.character(head(clusterList, n=1)),"_",as.character(tail(clusterList, n=1)),"_",type_cluster,"_euclidean.png"),width = 10,
         height = 6)
}
# input: list dmrs containing list of contained CpGs, output: list dmrs containing 
# list of contained CpGs and their corresponding cluster (if they have none NA)
get_Clusters <- function(dmr_cpgs){
  return(sapply(unlist(dmr_cpgs), FUN=function(cpg){
    if(cpg %in% colnames(cluster_memb)){
      return(cluster_memb[[cpg]])
    }
    else{
      return(NA)
    }
  }))
}
# Functions to return the Frequencies, Percentages of CpGs of each Clusteror 
# simply the CpGs with their Cluster  for each DMR (to build data tables)
get_Percentage <- function(CpGs_clusters){
  dt <- as.data.table(prop.table(table(CpGs_clusters, useNA = "always")))
  cnames <- dt$CpGs_clusters
  cnames[is.na(cnames)] <- "NoCluster"
  dt$CpGs_clusters <- NULL
  dt <- as.data.frame(t(dt))
  colnames(dt)<- as.character(cnames)
  return(dt)
}
get_Frequencies <- function(CpGs_clusters){
  dt <- as.data.table(table(CpGs_clusters, useNA = "always"))
  cnames <- dt$CpGs_clusters
  cnames[is.na(cnames)] <- "NoCluster"
  dt$CpGs_clusters <- NULL
  dt <- as.data.frame(t(dt))
  colnames(dt)<- cnames
  return(dt)
}
just_Return <- function(CpGs_clusters){
  dt <- transpose(as.data.table(CpGs_clusters))
  colnames(dt) <- as.character(c(1:ncol(dt)))
  return(as.data.table(dt))
}
onehot_encoding_DMRs_per_CpGs <- function(CpGs){
  dt <- as.data.table(table(CpGs, useNA = "always"))
  cnames <- dt$CpGs
  cnames[is.na(cnames)] <- "NoCpG"
  dt$CpGs <- NULL
  dt <- as.data.frame(t(dt))
  colnames(dt)<- as.character(cnames)
  return(dt)
}
# Corresponding make data table function that expects one of the above function 
# as one input along with the fitting dmr list of lists the right amount of cols
# as an integer, the colnames and the rownames
makeDT <- function(dmr_list, collength, colnm, func, rownm){
  DT_temp <- as.data.table(matrix(nrow = 1, ncol = collength))
  colnames(DT_temp) <- colnm
  all_dmrs <- bind_rows(DT_temp, lapply(dmr_list, FUN=func))
  all_dmrs <- all_dmrs[-1,]
  rownames(all_dmrs) <- rownm
  return(all_dmrs)
}

### Load Data ####
path <- "/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/R_Dimmer/Thesis_Results/2-AGE/CpGClustering/"
beta_matrix_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE62219_5yearsAfterBirth/output4_friedman/beta_matrix.csv'
annotation_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE62219_5yearsAfterBirth/sample_annotation2.csv'
dimmer_anova_project_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE62219_5yearsAfterBirth/output5_anova/dimmer_project.csv'
dimmer_friedman_project_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE62219_5yearsAfterBirth/output4_friedman/dimmer_project.csv'
dmr_anova_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE62219_5yearsAfterBirth/output5_anova/230821122336836_DMRSearch_org/merged_table.csv'
dmr_friedman_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/Datensets/GSE62219_5yearsAfterBirth/output4_friedman/230821123217817_DMRSearch_org/merged_table.csv'
manifest450k_file <- '/Users/basti/Documents/Uni/Bioinformatik/Bachelorarbeit/DimMeR/DimmerJava/src/main/resources/manifest_summary.csv' # to get the CpGs contained in/ overlapping with a DMR

beta_matrix <- data.table::fread(beta_matrix_file)
annotation_data <- data.table::fread(annotation_file)
dimmer_project_anova <- data.table::fread(dimmer_anova_project_file)
dimmer_project_friedman <- data.table::fread(dimmer_friedman_project_file)
manifest <- data.table::fread(manifest450k_file)

### Get DMRs and corresponding CpGs ####
All_CpGs_order <- beta_matrix$CpG_ID

# produce GRanges Object with positions of those CpGs (and their names & regions)
manifest_pos <- data.table(Chr=manifest$CHR, POS=manifest$MAPINFO, IlmnID=manifest$IlmnID, regionType=manifest$regionType)
manifest_pos$Chr <- paste("chr", manifest_pos$Chr, sep="")
manifest_ranges <- makeGRangesFromDataFrame(manifest_pos,start.field = "POS", end.field = "POS", keep.extra.columns=TRUE) # .bed format
manifest_ranges

# produce Granges Object for ANOVA and Friedman DMRs
dmr_anova <- data.table::fread(dmr_anova_file)
dmr_anova$Chr <- paste("chr", dmr_anova$Chr, sep="")
dmr_anova_ranges <- makeGRangesFromDataFrame(dmr_anova[,1:3],start.field = "Begin", end.field = "End") # .bed format

dmr_friedman <- data.table::fread(dmr_friedman_file)
dmr_friedman$Chr <- paste("chr", dmr_friedman$Chr, sep="")
dmr_friedman_ranges <- makeGRangesFromDataFrame(dmr_friedman[,1:3],start.field = "Begin", end.field = "End") # .bed format

anova_dmrs_CpGs <- getOverlappingCpGs(dmr_anova_ranges)
friedman_dmrs_CpGs <- getOverlappingCpGs(dmr_friedman_ranges)

# Unique CpG List from all CpGs contained in DMRs (there should be no double 
# mentioning of CpGs so this is basically just making one linst of all CpGs 
# from all DMRs)
anova_CpGs <- unique(unlist(anova_dmrs_CpGs))
length(anova_CpGs)
friedman_CpGs <- unique(unlist(friedman_dmrs_CpGs))
length(friedman_CpGs)

# Filter out all low quality or non significant CpGs that occur in the DMRs 
# (all those CpGs are either not mentioned in the dimmer project file (that 
# saves calculated CpG p-values) or have a p-value <= 0.05)
anova_CpGs_filtered <- anova_CpGs[anova_CpGs %in% dimmer_project_anova[dimmer_project_anova$ORG <= 0.05]$CPG]
friedman_CpGs_filtered <- friedman_CpGs[friedman_CpGs %in% dimmer_project_friedman[dimmer_project_friedman$ORG <= 0.05]$CPG]

#### Get Clustering Data/ TimeSeries Beta values from beta_matrix and annotationfile for those CpGs ####
anova_dt <- build_dt(anova_CpGs_filtered)
friedman_dt <- build_dt(friedman_CpGs_filtered)


# individually delete columns of annotation that are not needed/ used here
colnames(anova_dt)[1:15]
anova_dt <- anova_dt %>% select(-c("PatientNr","Gender_ID", "Date_of_birth", "HLA-DR-DQ_haplotype", "T1D_RiskClass", "Delivery", "smoking_during_pregnancy", "exclusive_end_breast-feeding(m)", "total_end_breast-feeding(m)"))
friedman_dt <- friedman_dt %>% select(-c("PatientNr","Gender_ID", "Date_of_birth", "HLA-DR-DQ_haplotype", "T1D_RiskClass", "Delivery", "smoking_during_pregnancy", "exclusive_end_breast-feeding(m)", "total_end_breast-feeding(m)"))

#### Reduce Dimensionality by substracting out the Patient dimension: ####
# Calculate mean values over all Patients grouped by timestamp
ncores = 4
anova_dt_red <- anova_dt[, parallel::mclapply(.SD, mean, mc.cores = as.numeric(ncores)), by = timestamp]
friedman_dt_red <- friedman_dt[, parallel::mclapply(.SD, mean, mc.cores = as.numeric(ncores)), by = timestamp]

# Get and transpose Matrix for Clustering
m_anova <- data.matrix(anova_dt_red[,-1])
m_anova <- t(m_anova)
m_friedman <- data.matrix(friedman_dt_red[,-1])
m_friedman <- t(m_friedman)

# delete rows with NA values
sum(is.na(m_anova))
m_anova <- na.omit(m_anova)
sum(is.na(m_friedman))
m_friedman <- na.omit(m_friedman)

# CpG Clustering ----------------------------------------------------------
## Cluster type and centroid type for all following ####
type_cluster <-"hierarchical" # partitional tadpole hierarchical

## ANOVA First get the best value for k (the amount of clusters you want) ####
method <- "ANOVA"

print(Sys.time())
clust.all <- tsclust(m_anova, type=type_cluster, k=2L:200L, distance="Euclidean")
print(Sys.time())
names(clust.all) <- c(2L:200L)
# calculate CVI scored for all clusterings depending on k
scores <- lapply(clust.all,FUN=function(x){
  return(cvi(x,type = c("Sil","DB","DBstar","CH","SF")))
})
# calculate average within-cluster distance to centroid for all clusterings 
# depending on k
sum_average_clustdist  <- lapply(clust.all,FUN=function(x){
  return(mean(x@clusinfo$av_dist))
})

# Plot respective score
setDT(scores)
scores <- transpose(scores)
colnames(scores) <- c("Sil","DB","DBstar","CH","SF")

plotScore(scores$Sil, "Silhouette index", "Sil", c(2L:200L), method,FALSE) # to be maximized
plotScore(scores$Sil, "Silhouette index", "Sil", c(2L:20L), method,TRUE)
plotScore(scores$DB, "Davies-Bouldin index", "DB", c(2L:200L), method,FALSE) # to be minimized
plotScore(scores$DB, "Davies-Bouldin index", "DB", c(2L:20L), method,TRUE)
plotScore(scores$DBstar, "Modified Davies-Bouldin index", "DBstar", c(2L:200L), method,FALSE) # to be minimized
plotScore(scores$DBstar, "Modified Davies-Bouldin index", "DBstar", c(2L:20L), method,TRUE)
plotScore(scores$CH, "Calinski-Harabasz index", "CH", c(2L:200L), method,FALSE) # to be maximized
plotScore(scores$CH, "Calinski-Harabasz index", "CH", c(2L:20L), method,TRUE)
plotScore(scores$SF, "Score Function", "SF", c(2L:200L), method,FALSE) # to be minimized
plotScore(scores$SF, "Score Function", "SF", c(2L:20L), method,TRUE)

# Plot the average within-cluster distance to Centroid
ggplot(data.frame(), aes(x=c(2L:20L), y=unlist(sum_average_clustdist)[1:19], group= 1)) +
  geom_line() +
  geom_point() +
  xlab("Number of Clusters") +
  ylab("value") +
  ggtitle("Average Within-Cluster Distance to Centroid")
ggsave(paste0(path,"ANOVA/WithinClusterDistance_2_20_",type_cluster,"_euclidean.png"),width = 10,
       height = 6)

## ANOVA Clustering with selected k ####
method <- "ANOVA"
k_cluster <- 12L # 12L 15L 

print(Sys.time())
clust.hist <- tsclust(m_anova, type=type_cluster, k=k_cluster , distance="Euclidean")
print(Sys.time())

height <- clust.hist$height[length(clust.hist$height) - k_cluster + 1]

# Plot Clustering results
png(file=paste0(path, method, "/HierarchicalTree_",k_cluster,"_euclidean.png"), # _euclidean_
    width=5000, height=800)
plot(clust.hist)# + ylim(0.0,1.0)
rect.hclust(clust.hist , k = k_cluster, border = 2:6)
abline(h = height, col = "orange")
dev.off()
png(file=paste0(path, method, "/All_",k_cluster,"_Cluster_TimeCourses_",type_cluster,"_euclidean.png"), # _euclidean_
    width=1300, height=800)
plot(clust.hist, type = "sc") + ylim(0.0,1.0) +  scale_colour_viridis(option= "inferno", discrete = TRUE, direction = -1) +  #TODO Test color and x axis labels
   scale_x_continuous(breaks=c(1,2,3,4,5), labels= c(3,12,24,48,60)) + xlab("") + ylab("Beta Value") +   ggtitle("")
dev.off()
png(file=paste0(path, method, "/All_",k_cluster,"_Cluster_Centroids_",type_cluster,"_euclidean.png"), # _euclidean_
    width=1300, height=800)
plot(clust.hist, type = "centroids")  + ylim(0.0,1.0) + scale_x_discrete(labels= c(3,12,24,48,60))
dev.off()

## Friedman First get the best value for k (the amount of clusters you want) ####
method <- "Friedman"

print(Sys.time())
clust.all2 <- tsclust(m_friedman, type=type_cluster, k=2L:200L, distance="Euclidean")
print(Sys.time())
names(clust.all2) <- c(2L:200L)
scores2 <- lapply(clust.all2,FUN=function(x){
  return(cvi(x,type = c("Sil","DB","DBstar","CH","SF")))
})
sum_average_clustdist2  <- lapply(clust.all2,FUN=function(x){
  return(mean(x@clusinfo$av_dist))
})

# Plot respective score
setDT(scores2)
scores2 <- transpose(scores2)
colnames(scores2) <- c("Sil","DB","DBstar","CH","SF")

plotScore(scores2$Sil, "Silhouette index", "Sil", c(2L:200L), method,FALSE) # to be maximized
plotScore(scores2$Sil, "Silhouette index", "Sil", c(2L:20L), method,TRUE)
plotScore(scores2$DB, "Davies-Bouldin index", "DB", c(2L:200L), method,FALSE) # to be minimized
plotScore(scores2$DB, "Davies-Bouldin index", "DB", c(2L:20L), method,TRUE)
plotScore(scores2$DBstar, "Modified Davies-Bouldin index", "DBstar", c(2L:200L), method,FALSE) # to be minimized
plotScore(scores2$DBstar, "Modified Davies-Bouldin index", "DBstar", c(2L:20L), method,TRUE)
plotScore(scores2$CH, "Calinski-Harabasz index", "CH", c(2L:200L), method,FALSE) # to be maximized
plotScore(scores2$CH, "Calinski-Harabasz index", "CH", c(2L:20L), method,TRUE)
plotScore(scores2$SF, "Score Function", "SF", c(2L:200L), method,FALSE) # to be minimized
plotScore(scores2$SF, "Score Function", "SF", c(2L:20L), method,TRUE)

# Plot the average within-cluster distance to Centroid
ggplot(data.frame(), aes(x=c(2L:20L), y=unlist(sum_average_clustdist2)[1:19], group= 1)) +
  geom_line() +
  geom_point() +
  xlab("Number of Clusters") +
  ylab("value") +
  ggtitle("Average Within-Cluster Distance to Centroid")
ggsave(paste0(path,"Friedman/WithinClusterDistance_2_20_",type_cluster,"_euclidean.png"),width = 10,
       height = 6)

## Friedman Clustering with selected k ####
method <- "Friedman"
k_cluster <- 12L # 12L 14L 96L

print(Sys.time())
clust.hist2 <- tsclust(m_friedman, type=type_cluster, k=k_cluster , distance="Euclidean")
print(Sys.time())

height <- clust.hist2$height[length(clust.hist2$height) - k_cluster + 1]

# Plot Clustering results
png(file=paste0(path, method, "/HierarchicalTree_",k_cluster,"_euclidean.png"), # _euclidean_
    width=5000, height=800)
plot(clust.hist2)# + ylim(0.0,1.0)
rect.hclust(clust.hist2 , k = k_cluster, border = 2:6)
abline(h = height, col = "orange")
dev.off()
png(file=paste0(path, method, "/All_",k_cluster,"_Cluster_TimeCourses_",type_cluster,"_euclidean.png"), # _euclidean_
    width=1300, height=800)
plot(clust.hist2, type = "sc") + ylim(0.0,1.0) + scale_colour_viridis(option= "viridis", discrete = TRUE, direction = -1) + #TODO Test color and x axis labels
  scale_x_continuous(breaks=c(1,2,3,4,5), labels= c(3,12,24,48,60)) + xlab("") + ylab("Beta Value") +   ggtitle("")
dev.off()
png(file=paste0(path, method, "/All_",k_cluster,"_Cluster_Centroids_",type_cluster,"_euclidean.png"), # _euclidean_
    width=1300, height=800)
plot(clust.hist2, type = "centroids")+ ylim(0.0,1.0) + scale_x_discrete(labels= c(3,12,24,48,60))
dev.off()

# Diversity Score ---------------------------------------------------------
## ANOVA Alpha Diversity ####
method <- "ANOVA"

# Taxonomic Diversity with vegan package

# Get the Cluster Membership for each CpG
cluster_memb <- transpose(as.data.table(clust.hist@cluster))
rownames(cluster_memb) <- c("Assigned_Clusters")
colnames(cluster_memb) <- colnames(anova_dt)[-1]

#calculate distance between centroids using euclidean distance
centroids_table <- t(as.data.table(clust.hist@centroids))
rownames(centroids_table) <- c(1:k_cluster)
# Get the Distance Matrix for Clusters
dismat_cluster = as.matrix(dist(centroids_table, method = "euclidean", diag = TRUE, upper = TRUE))
png(file=paste0(path, method, "/Heatmap_ClusterDistance.png"),
    width=800, height=800)
Heatmap(dismat_cluster, show_column_dend = TRUE, show_row_dend = FALSE, 
        row_order= c(1:k_cluster),
        col =  viridis(200, option="inferno", direction = -1), column_names_rot = 0, 
        column_dend_height=unit(5,"cm"),
        row_names_gp= gpar(fontsize = 20), column_names_gp= gpar(fontsize = 20),
        heatmap_legend_param = 
          list(title="Euclidean\nDistance\n", legend_height = unit(4, "cm"),
               title_gp = gpar(fontsize = 20),
               labels_gp = gpar(fontsize = 15)))
dev.off()

# Get the Clusters of the CpGs in a DMR
dmrs_CpGs_clusters <- sapply(anova_dmrs_CpGs, FUN=get_Clusters)

# Create a temp data.table with all clusters as column names and a "NoCluster"
# Column for low quality CpG that have been filtered out by Dimmer or 
# non-significant CpGs and are therefore not included in the clustering
all_dmr_clust_freq <- makeDT(dmr_list=dmrs_CpGs_clusters, collength=(k_cluster + 1),
                             colnm=as.character(c(1:k_cluster,"NoCluster")), 
                             func=get_Frequencies, rownm=rownames(dmr_anova))
all_dmr_clust_freq[is.na(all_dmr_clust_freq)] <- 0

# We're not interested in the NoCluster CpGs for the score
all_dmr_clust_freq_noNC <- all_dmr_clust_freq
all_dmr_clust_freq_noNC$NoCluster <- NULL

diversity_dmr_cluster_temp <- all_dmr_clust_freq

# Calculate Taxonomic Diversity
taxonomic_diversity <- taxondive(as.matrix(all_dmr_clust_freq_noNC), dis = dismat_cluster, match.force = TRUE)
diversity_dmr_cluster_temp$taxonomic_diversity <- taxonomic_diversity$D
diversity_dmr_cluster_temp$amountDiffCluster <- taxonomic_diversity$Species

# inspect Scores
summary(diversity_dmr_cluster_temp$taxonomic_diversity)

## Visualize Diversity/ Heterogenity Score
# Make Heatmap of DMRs (Position, Cluster) 

## creat continous color function representing diversity score ##
mi <- min(diversity_dmr_cluster_temp$taxonomic_diversity)
mx <- max(diversity_dmr_cluster_temp$taxonomic_diversity)
mid <- (mi + mx) /2
RowSideCols <- colorRamp2(c(mi,(mid-mi)/2,mid+(mid-mi)/2,mx), c("white","blue","blue", "black"))


# colors for the Clusters
new_color_palatte <- as.character(num2col(c(-1:12), col.pal=redpal))
names(new_color_palatte) <- c("-1","0","1","2","3","4","5","6","7","8","9","10","11","12")
new_color_palatte["-1"] <- "gainsboro" # Insignificant CpGs (or those that are not even in the beta matrix)
new_color_palatte["0"] <- "white" # No CpGs/ empty Position 
new_color_palatte <- as.factor(new_color_palatte)

# get max DMR length
max_dmr_length <- max(unlist(lapply(dmrs_CpGs_clusters,FUN=length)))

# Change this according to how the Insignificant CpGs are supposed to look
dmrs_CpGs_clusters <- lapply(dmrs_CpGs_clusters, function(x) replace(x,is.na(x),-1))
# Create a temp data.table of the max DMR length and with DMR positions as columns
all_dmr_all_clusters_with_position <- makeDT(dmr_list=dmrs_CpGs_clusters, 
                                              collength=max_dmr_length, 
                                              colnm=as.character(c(1:max_dmr_length)), 
                                              func=just_Return, rownm=rownames(dmr_anova))

# Change this to change how the empty postions appear for dmrs that are shorter 
# than the max length dmr
all_dmr_all_clusters_with_position[is.na(all_dmr_all_clusters_with_position)] <- 0
all_dmr_all_clusters_with_position$DMR <- rownames(dmr_anova)
all_dmr_all_clusters_with_position$taxonomic_diversity <- diversity_dmr_cluster_temp$taxonomic_diversity

##### BEST Get all the best scored (low score) DMRs
temp <- all_dmr_all_clusters_with_position[all_dmr_all_clusters_with_position$taxonomic_diversity <= 0.15]
temp <- temp[order(temp$taxonomic_diversity, decreasing = FALSE), ]

# What is the max length among them, convert sized dt to matrix
max_dmr_temp_length <- max(unlist(lapply(dmrs_CpGs_clusters[as.numeric(temp$DMR)],FUN=length)))
input <- as.matrix(temp[,1:max_dmr_temp_length])
rownames(input) <- temp$DMR

# Include Diversity Score as extra Heatmap Annotation
ha = rowAnnotation(df= data.frame(divers =temp$taxonomic_diversity), col = list(divers = RowSideCols),
                   show_annotation_name = FALSE, annotation_legend_param = list(divers =list(
                     title = "Diversity Score", at = c(0.0, 0.07, 0.15), labels = c("0", "0.07", "0.15"))),
                   border = FALSE,gp = gpar(col = "gainsboro"))

# Create Heatmap corresponding to the CpGs Clusters and their Position on the DMR
h1 <- Heatmap(input, col=new_color_palatte, name = "Clusters",
        column_title = "Positions on DMR", column_title_side = "bottom", column_names_rot = 0, 
        row_title = gt_render("DMRs", padding = unit(c(-4, -4, -4, -4), "pt")),
        cluster_columns = FALSE,
        row_title_side = "left", row_names_side = "left", row_split= factor(c(1:nrow(input)), levels = c(1:nrow(input))), gap = unit(0.1,"mm"),
        cluster_rows = FALSE, row_names_gp = gpar(fontsize = 5),
        border_gp = gpar(col = "white"), left_annotation = ha,
        heatmap_legend_param = list(
          title = "Clusters", at= c(-1:12),
          labels =  c("Insignificant CpG", "", 1:12)
        ))

# plot and save figure
pdf(file=paste0(path, method, "/Heatmap_DiversityScore_DMRClusterPositions_LOWDiversity.pdf"),
    width=(max_dmr_temp_length), height=(0.1* nrow(input)))
draw(h1)
dev.off()

#### WORST Get all the worst scored (high score) DMRs
temp <- all_dmr_all_clusters_with_position[all_dmr_all_clusters_with_position$taxonomic_diversity >= 0.8]
temp <- temp[order(temp$taxonomic_diversity,), ]

# What is the max length among them, convert sized dt to matrix
max_dmr_temp_length <- max(unlist(lapply(dmrs_CpGs_clusters[as.numeric(temp$DMR)],FUN=length)))
input <- as.matrix(temp[,1:max_dmr_temp_length])
rownames(input) <- temp$DMR

# Include Diversity Score as extra Heatmap Annotation
#RowSideCols <- setNames(as.list(unique(temp$taxonomic_diversity)), num2col(unique(temp$taxonomic_diversity), col.pal=redpal))
ha = rowAnnotation(df = data.frame(divers =temp$taxonomic_diversity), col = list(divers = RowSideCols),   
                   show_annotation_name = FALSE, annotation_legend_param = list(divers = list(
                     title = "Diversity Score", at = c(0.8, 1.0, 1.18), labels = c("0.8", "1.0", "1.18"))),
                   border = FALSE,gp = gpar(col = "gainsboro"))

# Create Heatmap corresponding to the CpGs Clusters and their Position on the DMR
h1 <- Heatmap(input, col=new_color_palatte, name = "Clusters",
              column_title = "Positions on DMR", column_title_side = "bottom", column_names_rot = 0, 
              row_title = gt_render("DMRs", padding = unit(c(-4, -4, -4, -4), "pt")),
              cluster_columns = FALSE,
              row_title_side = "left", row_names_side = "left", row_split= factor(c(1:nrow(input)), levels = c(1:nrow(input))), gap = unit(0.1,"mm"),
              cluster_rows = FALSE, row_names_gp = gpar(fontsize = 5),
              border_gp = gpar(col = "white"), left_annotation = ha,
              heatmap_legend_param = list(
                title = "Clusters", at= c(-1:12),
                labels =  c("Insignificant CpG", "", 1:12)
              ))

# plot and save figure
pdf(file=paste0(path, method,"/Heatmap_DiversityScore_DMRClusterPositions_HIGHDiversity.pdf"),
      width=(0.8*max_dmr_temp_length), height=(0.1* nrow(input)))
draw(h1)
dev.off()

# Plot Heterogenity against DMRScore or Pvalue
diversity_dmr_cluster_temp$Pvalue <- dmr_anova$`p-value`
diversity_dmr_cluster_temp$NumCpGs <- dmr_anova$`#CpG`
diversity_dmr_cluster_temp$color <- "Insign"
diversity_dmr_cluster_temp$color[diversity_dmr_cluster_temp$Pvalue < 0.01 & diversity_dmr_cluster_temp$taxonomic_diversity < 0.8] <- "Sign"

ggplot(diversity_dmr_cluster_temp, aes(x=taxonomic_diversity, y=-log10(Pvalue), color = color)) +
  geom_point() + geom_vline(xintercept = 0.8, colour="darkgray", linewidth=1) + 
  geom_hline(yintercept = -log10(0.01), colour="darkgray", linewidth=1) +
  theme(legend.position = "none") + 
  ggtitle(paste("Scatterplot",method,"DMRs - P-value vs. Diversity Score")) +
  xlab("Heterogenity Score") +
  scale_color_manual(values = c(Sign="#F8766D", Insign="gray")) 
ggsave(paste0(path,method,"/ScatterPlot_TaxonomicDiversity_DMRPvalue.png", sep = ""))

# Have a look how much Taxonomic Diversity correlates with Number of CpGs that 
# are accounted for in the score
ggplot(diversity_dmr_cluster_temp, aes(x=taxonomic_diversity, y=NumCpGs)) +
  geom_point() + 
  theme(legend.position = "none") + 
  ggtitle(paste("Scatterplot",method,"DMRs - Number of CpGs vs. Diversity Score")) +
  ylab("Number of CpGs (accounted for in the Score)") +
  xlab("Heterogenity Score") +
  scale_color_manual(values = c("gray"))
ggsave(paste0(path,method,"/ScatterPlot_TaxonomicDiversity_DMR_Size.png", sep = ""))

## ANOVA Beta Diversity ####
# Create a temp data.table with all CpGs as column names / an OTU table with CpGs as taxa
OTU_temp <- makeDT(dmr_list=anova_dmrs_CpGs, 
                   collength=length(colnames(anova_dt_red)[-1])+1, 
                   colnm=as.character(c(colnames(anova_dt_red)[-1],"NoCpG")), 
                   func=onehot_encoding_DMRs_per_CpGs, rownm=rownames(dmr_anova))
OTU_temp[is.na(OTU_temp)] <- 0
# check if one CpG was not consisten with the names otherwise delete obsolete col
unique(OTU_temp$NoCpG)
OTU_temp$NoCpG <- NULL
#which(colSums(OTU_temp==0) == nrow(OTU_temp))

abundance <- as.matrix(OTU_temp)
rownames(abundance) <- rownames(OTU_temp)
colnames(abundance) <- colnames(OTU_temp)
OTU = otu_table(abundance, taxa_are_rows = FALSE)

#which(colSums(df==1) == nrow(df))

myTree <- as.phylo(clust.hist)
phy <- phyloseq(OTU, myTree)

unifrac_dist <- UniFrac(phy, weighted=TRUE)

unifrac_dist_mat <- data.matrix(unifrac_dist)   

colnames(unifrac_dist_mat) <- rownames(dmr_anova)
rownames(unifrac_dist_mat) <- rownames(dmr_anova)

#unifrac_dist_mat[c(40,59,265,275,360,409,416),275]

#all_dmr_all_clusters_with_position[310,]
#all_dmr_all_clusters_with_position[275,]
#all_dmr_all_clusters_with_position[c(130, 247, 263, 291, 304, 347, 377),]

dmr_nr = 130

which(unifrac_dist_mat[,dmr_nr] <= 0.17, arr.ind=TRUE)
co_modified_dmrs <- as.vector(which(unifrac_dist_mat[,dmr_nr] <= 0.17, arr.ind=TRUE))

# create dfs for plotting
df <- data.frame (DMRs  = rownames(unifrac_dist_mat),
                  unifr = unifrac_dist_mat[,dmr_nr]
)
df_marked <- data.frame(DMRs  = as.vector(co_modified_dmrs),
                        unifr = as.vector(unifrac_dist_mat[co_modified_dmrs,dmr_nr])
)

ggplot(df, aes(x=as.numeric(DMRs), y=unifr, label = DMRs)) +
  geom_point(color="darkgrey") +
  ggtitle(paste("Unifrac Distance between DMR no. ",dmr_nr," and all other DMRs",sep="")) +
  theme(legend.position = "none",axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_point(data=df_marked,
             aes(x=as.numeric(DMRs),y=unifr, 
                 size=5),color='violet') +
  geom_label_repel(data=df_marked, box.padding = 1.5) + 
  ylab("Unifrac Distance") + xlab("DMRs")
ggsave(paste0(path,method,"/ScatterPlot_UnifracScore_DMR_",dmr_nr,".png", sep = ""))

sum(is.na(as.vector(unifrac_dist_mat)))
# Make Histogram of unifrac distances (not really reflective since every value 
# except the diagonal contributes once to offen)
ggplot() + aes(x=as.vector(unifrac_dist_mat)) + geom_histogram(binwidth=0.01,color="black", fill="#F8766D") +
  xlab("Unifrac Distance") + xlim(0,1)
ggsave(paste0(path,method,"/UnifracDistanceHistogram.png", sep = ""),width = 1200,
       height = 1200,units = "px")

colnames(unifrac_dist_mat) <- NULL
rownames(unifrac_dist_mat) <- NULL
# Make Heatmap of unifrac distances
hb <- Heatmap(unifrac_dist_mat,show_column_dend = FALSE, col = viridis(200, option="inferno"),
              name = "Unifrac Distance")
png(file=paste0(path,method,"/Heatmap_UnifracDistances.png", sep = ""))
draw(hb)
dev.off()

### Friedman Alpha Diversity ####
method <- "Friedman"

# Taxonomic Diversity with vegan package
# Get the Cluster Membership for each CpG
cluster_memb <- transpose(as.data.table(clust.hist2@cluster))
rownames(cluster_memb) <- c("Assigned_Clusters")
colnames(cluster_memb) <- colnames(friedman_dt)[-1]

#calculate distance between centroids using euclidean distance
centroids_table <- t(as.data.table(clust.hist2@centroids))
rownames(centroids_table) <- c(1:k_cluster)
# Get the Distance Matrix for Clusters
dismat_cluster = as.matrix(dist(centroids_table, method = "euclidean", diag = TRUE, upper = TRUE))
png(file=paste0(path, method, "/Heatmap_ClusterDistance.png"),
    width=800, height=800)
Heatmap(dismat_cluster, show_column_dend = TRUE, show_row_dend = FALSE, 
        row_order= c(1:k_cluster),
        col =  viridis(200, direction = -1), column_names_rot = 0, 
        column_dend_height=unit(5,"cm"),
        row_names_gp= gpar(fontsize = 20), column_names_gp= gpar(fontsize = 20),
        heatmap_legend_param = 
          list(title="Euclidean\nDistance\n", legend_height = unit(4, "cm"),
               title_gp = gpar(fontsize = 20),
        labels_gp = gpar(fontsize = 15)))
dev.off()

# Get the Clusters of the CpGs in a DMR
dmrs_CpGs_clusters <- sapply(friedman_dmrs_CpGs, FUN=get_Clusters)

# Create a temp data.table with all clusters as column names and a "NoCluster"
# Column for low quality CpG that have been filtered out by Dimmer or 
# non-significant CpGs and are therefore not included in the clustering
all_dmr_clust_freq <- makeDT(dmr_list=dmrs_CpGs_clusters, collength=(k_cluster + 1),  #TODO change clustername if friedman has different clusters
                             colnm=as.character(c(1:k_cluster,"NoCluster")), 
                             func=get_Frequencies, rownm=rownames(dmr_friedman))
all_dmr_clust_freq[is.na(all_dmr_clust_freq)] <- 0

# We're not interested in the NoCluster CpGs for the score
all_dmr_clust_freq_noNC <- all_dmr_clust_freq
all_dmr_clust_freq_noNC$NoCluster <- NULL

diversity_dmr_cluster_temp <- all_dmr_clust_freq

# Calculate Taxonomic Diversity
taxonomic_diversity <- taxondive(as.matrix(all_dmr_clust_freq_noNC), dis = dismat_cluster, match.force = TRUE)
diversity_dmr_cluster_temp$taxonomic_diversity <- taxonomic_diversity$D
diversity_dmr_cluster_temp$amountDiffCluster <- taxonomic_diversity$Species

# inspect Scores
summary(diversity_dmr_cluster_temp$taxonomic_diversity)

## Visualize Diversity/ Heterogenity Score
# Make Heatmap of DMRs (Position, Cluster) 

## creat continous color function representing diversity score ##
mi <- min(diversity_dmr_cluster_temp$taxonomic_diversity)
mx <- max(diversity_dmr_cluster_temp$taxonomic_diversity)
mid <- (mi + mx) /2
RowSideCols <- colorRamp2(c(mi,(mid-mi)/2,mid+(mid-mi)/2,mx), c("white","blue","blue", "black"))

# colors for the Clusters
new_color_palatte <- as.character(num2col(c(-1:12), col.pal=greenpal))
names(new_color_palatte) <- c("-1","0","1","2","3","4","5","6","7","8","9","10","11","12")
new_color_palatte["-1"] <- "gainsboro" # Insignificant CpGs (or those that are not even in the beta matrix)
new_color_palatte["0"] <- "white" # No CpGs/ empty Position 
new_color_palatte <- as.factor(new_color_palatte)

# get max DMR length
max_dmr_length <- max(unlist(lapply(dmrs_CpGs_clusters,FUN=length)))

# Change this according to how the Insignificant CpGs are supposed to look
dmrs_CpGs_clusters <- lapply(dmrs_CpGs_clusters, function(x) replace(x,is.na(x),-1))
# Create a temp data.table of the max DMR length and with DMR positions as columns
all_dmr_all_clusters_with_position <- makeDT(dmr_list=dmrs_CpGs_clusters, 
                                             collength=max_dmr_length, 
                                             colnm=as.character(c(1:max_dmr_length)), 
                                             func=just_Return, rownm=rownames(dmr_friedman))

# Change this to change how the empty postions appear for dmrs that are shorter 
# than the max length dmr
all_dmr_all_clusters_with_position[is.na(all_dmr_all_clusters_with_position)] <- 0
all_dmr_all_clusters_with_position$DMR <- rownames(dmr_friedman)
all_dmr_all_clusters_with_position$taxonomic_diversity <- diversity_dmr_cluster_temp$taxonomic_diversity

##### BEST Get all the best scored (low score) DMRs
temp <- all_dmr_all_clusters_with_position[all_dmr_all_clusters_with_position$taxonomic_diversity <= 0.15]
temp <- temp[order(temp$taxonomic_diversity, decreasing = FALSE), ]

# What is the max length among them, convert sized dt to matrix
max_dmr_temp_length <- max(unlist(lapply(dmrs_CpGs_clusters[as.numeric(temp$DMR)],FUN=length)))
input <- as.matrix(temp[,1:max_dmr_temp_length])
rownames(input) <- temp$DMR

# Include Diversity Score as extra Heatmap Annotation
ha = rowAnnotation(df= data.frame(divers =temp$taxonomic_diversity), col = list(divers = RowSideCols),
                   show_annotation_name = FALSE, annotation_legend_param = list(divers =list(
                     title = "Diversity\nScore", at = c(0.0, 0.07, 0.15), labels = c("0", "0.07", "0.15"))),
                   border = FALSE,gp = gpar(col = "gainsboro"))

# Create Heatmap corresponding to the CpGs Clusters and their Position on the DMR
h1 <- Heatmap(input, col=new_color_palatte, name = "Clusters",
              column_title = "Positions on DMR", column_title_side = "bottom", column_names_rot = 0, 
              row_title = gt_render("DMRs", padding = unit(c(-4, -4, -4, -4), "pt")),
              cluster_columns = FALSE,
              row_title_side = "left", row_names_side = "left", row_split= factor(c(1:nrow(input)), levels = c(1:nrow(input))), gap = unit(0.1,"mm"),
              cluster_rows = FALSE, row_names_gp = gpar(fontsize = 5),
              border_gp = gpar(col = "white"), left_annotation = ha,
              heatmap_legend_param = list(
                title = "Clusters", at= c(-1:12),
                labels =  c("Insignificant CpG", "", 1:12)
              ))

# plot and save figure
pdf(file=paste0(path, method, "/Heatmap_DiversityScore_DMRClusterPositions_LOWDiversity.pdf"),
    width=(max_dmr_temp_length), height=(0.1* nrow(input)))
draw(h1)
dev.off()

#### WORST Get all the worst scored (high score) DMRs
temp <- all_dmr_all_clusters_with_position[all_dmr_all_clusters_with_position$taxonomic_diversity >= 0.8]
temp <- temp[order(temp$taxonomic_diversity,), ]

# What is the max length among them, convert sized dt to matrix
max_dmr_temp_length <- max(unlist(lapply(dmrs_CpGs_clusters[as.numeric(temp$DMR)],FUN=length)))
input <- as.matrix(temp[,1:max_dmr_temp_length])
rownames(input) <- temp$DMR

# Include Diversity Score as extra Heatmap Annotation
#RowSideCols <- setNames(as.list(unique(temp$taxonomic_diversity)), num2col(unique(temp$taxonomic_diversity), col.pal=redpal))
ha = rowAnnotation(df = data.frame(divers =temp$taxonomic_diversity), col = list(divers = RowSideCols),   
                   show_annotation_name = FALSE, annotation_legend_param = list(divers = list(
                     title = "Diversity Score", at = c(0.8, 1.0, 1.18), labels = c("0.8", "1.0", "1.18"))),
                   border = FALSE,gp = gpar(col = "gainsboro"))

# Create Heatmap corresponding to the CpGs Clusters and their Position on the DMR
h1 <- Heatmap(input, col=new_color_palatte, name = "Clusters",
              column_title = "Positions on DMR", column_title_side = "bottom", column_names_rot = 0, 
              row_title = gt_render("DMRs", padding = unit(c(-4, -4, -4, -4), "pt")),
              cluster_columns = FALSE,
              row_title_side = "left", row_names_side = "left", row_split= factor(c(1:nrow(input)), levels = c(1:nrow(input))), gap = unit(0.1,"mm"),
              cluster_rows = FALSE, row_names_gp = gpar(fontsize = 5),
              border_gp = gpar(col = "white"), left_annotation = ha,
              heatmap_legend_param = list(
                title = "Clusters", at= c(-1:12),
                labels =  c("Insignificant CpG", "", 1:12)
              ))

# plot and save figure
pdf(file=paste0(path, method,"/Heatmap_DiversityScore_DMRClusterPositions_HIGHDiversity.pdf"),
    width=(0.8*max_dmr_temp_length), height=(0.1* nrow(input)))
draw(h1)
dev.off()

# Plot Heterogenity against DMRScore or Pvalue
diversity_dmr_cluster_temp$Pvalue <- dmr_friedman$`p-value`
diversity_dmr_cluster_temp$NumCpGs <- dmr_friedman$`#CpG`
diversity_dmr_cluster_temp$color <- "Insign"
diversity_dmr_cluster_temp$color[diversity_dmr_cluster_temp$Pvalue < 0.01 & diversity_dmr_cluster_temp$taxonomic_diversity < 0.8] <- "Sign"

ggplot(diversity_dmr_cluster_temp, aes(x=taxonomic_diversity, y=-log10(Pvalue), color = color)) +
  geom_point() + geom_vline(xintercept = 0.8, colour="darkgray", linewidth=1) + 
  geom_hline(yintercept = -log10(0.01), colour="darkgray", linewidth=1) +
  theme(legend.position = "none") + 
  ggtitle(paste("Scatterplot",method,"DMRs - P-value vs. Diversity Score")) +
  xlab("Heterogenity Score") +
  scale_color_manual(values = c(Sign="#00BA38", Insign="gray")) 
ggsave(paste0(path,method,"/ScatterPlot_TaxonomicDiversity_DMRPvalue.png", sep = ""))

# Have a look how much Taxonomic Diversity correlates with Number of CpGs that 
# are accounted for in the score
ggplot(diversity_dmr_cluster_temp, aes(x=taxonomic_diversity, y=NumCpGs)) +
  geom_point() + 
  theme(legend.position = "none") + 
  ggtitle(paste("Scatterplot",method,"DMRs - Number of CpGs vs. Diversity Score")) +
  ylab("Number of CpGs (accounted for in the Score)") +
  xlab("Heterogenity Score") +
  scale_color_manual(values = c("gray"))
ggsave(paste0(path,method,"/ScatterPlot_TaxonomicDiversity_DMR_Size.png", sep = ""))

## Friedman Beta Diversity ####
# Create a temp data.table with all CpGs as column names / an OTU table with CpGs as taxa
OTU_temp <- makeDT(dmr_list=friedman_dmrs_CpGs, 
                   collength=length(colnames(friedman_dt_red)[-1])+1, 
                   colnm=as.character(c(colnames(friedman_dt_red)[-1],"NoCpG")), 
                   func=onehot_encoding_DMRs_per_CpGs, rownm=rownames(dmr_friedman))
OTU_temp[is.na(OTU_temp)] <- 0
# check if one CpG was not consisten with the names otherwise delete obsolete col
unique(OTU_temp$NoCpG)
OTU_temp$NoCpG <- NULL
#which(colSums(OTU_temp==0) == nrow(OTU_temp))

abundance <- as.matrix(OTU_temp)
rownames(abundance) <- rownames(OTU_temp)
colnames(abundance) <- colnames(OTU_temp)
OTU = otu_table(abundance, taxa_are_rows = FALSE)

#which(colSums(df==1) == nrow(df))

myTree <- as.phylo(clust.hist2)
phy <- phyloseq(OTU, myTree)

unifrac_dist <- UniFrac(phy, weighted=TRUE)

unifrac_dist_mat <- data.matrix(unifrac_dist)   

colnames(unifrac_dist_mat) <- rownames(dmr_friedman)
rownames(unifrac_dist_mat) <- rownames(dmr_friedman)

#unifrac_dist_mat[c(40,59,265,275,360,409,416),275]

#all_dmr_all_clusters_with_position[310,]
#all_dmr_all_clusters_with_position[275,]

dmr_nr = 130
threshold = 0.2

a <- which(unifrac_dist_mat[,dmr_nr] <= threshold, arr.ind=TRUE)
all_dmr_all_clusters_with_position[a,]
co_modified_dmrs <- as.vector(which(unifrac_dist_mat[,dmr_nr] <= threshold, arr.ind=TRUE))

# create dfs for plotting
df <- data.frame (DMRs  = rownames(unifrac_dist_mat),
                  unifr = unifrac_dist_mat[,dmr_nr]
                  )
df_marked <- data.frame(DMRs  = as.vector(co_modified_dmrs),
                           unifr = as.vector(unifrac_dist_mat[co_modified_dmrs,dmr_nr])
)

ggplot(df, aes(x=as.numeric(DMRs), y=unifr, label = DMRs)) +
  geom_point(color="darkgrey") +
  ggtitle(paste("Unifrac Distance between DMR no. ",dmr_nr," and all other DMRs",sep="")) +
  theme(legend.position = "none",axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept = threshold, colour="darkgray", linewidth=0.8, alpha=0.3) +
  geom_point(data=df_marked,
             aes(x=as.numeric(DMRs),y=unifr, 
                  size=5),color='violet') +
  geom_label_repel(data=df_marked, box.padding = 1.5) + 
  ylab("Unifrac Distance") + xlab("DMRs")
ggsave(paste0(path,method,"/ScatterPlot_UnifracScore_DMR_",dmr_nr,".png", sep = ""))

sum(is.na(as.vector(unifrac_dist_mat)))
# Make Histogram of unifrac distances (not really reflective since every value 
# except the diagonal contributes once to offen)
ggplot() + aes(x=as.vector(unifrac_dist_mat)) + geom_histogram(binwidth=0.01,color="black", fill="#00BA38") +
  xlab("Unifrac Distance") + xlim(0,1)
ggsave(paste0(path,method,"/UnifracDistanceHistogram.png", sep = ""),width = 1200,
       height = 1200,units = "px")

colnames(unifrac_dist_mat) <- NULL
rownames(unifrac_dist_mat) <- NULL
# Make Heatmap of unifrac distances
hb <- Heatmap(unifrac_dist_mat,show_column_dend = FALSE, col = viridis(200),
              name = "Unifrac Distance")
png(file=paste0(path,method,"/Heatmap_UnifracDistances.png", sep = ""))
draw(hb)
dev.off()
