# intro -------------------------------------------------------------------
rm(list=ls())
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
os="mac"

{
  if (os == "mac"){
    source("/Users/Max/OneDrive/Amplicon/r_scripts/2_phyloseq_start_qiime.R", echo = T, spaced = T)  
    source("/Users/Max/OneDrive/Amplicon/r_scripts/ancom_v2.1.R", echo = T, spaced = T)
    p.colors <- read.table("/Users/Max/OneDrive/Amplicon/r_scripts/phyla_colors.tsv", header=T)
    setwd("/Users/Max/OneDrive/Amplicon/") #Set to your 'phyloseq' folder!
  } else if (os == "ubuntu") {
    source("/home/mfbeuq/Amplicon/r_scripts/2_phyloseq_start_qiime.R", echo = T, spaced = T) 
    source("/home/mfbeuq/Amplicon/r_scripts/ancom_v2.1.R", echo = T, spaced = T)
    p.colors <- read.table("/home/mfbeuq/Amplicon/r_scripts/phyla_colors.tsv", header=T)
    setwd("/home/mfbeuq/Amplicon/") #Set to your 'phyloseq' folder!
  } else {
    print("Error")
  }
  library(dplyr)
  library(ANCOMBC)
  library(nlme)
  library(tidyverse)
  library(compositions)
  library(readr)
  library(DT)
  library(matrixStats)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(ComplexHeatmap) #Heatmap function
  library(circlize)
  library(round)
  library(readr)
  }

# Read data -----------------------------------------------------------------

temp <- readRDS("V2b_mh/phyloseq/ps.RDS")
conc <- readRDS("Marion_A1/phyloseq/ps.RDS")
straw <- readRDS("Marion_E1/phyloseq/ps.RDS")

# ANCOMBC -----------------------------------------------------------------

## Temporal: -----------------------------------------------------------------
ps <- temp
pseq <- ps 
genus_data <- tax_glom(pseq, taxrank <- rank_names(pseq)[4], NArm = FALSE)
var1 <- "compartment"


# Bulk Soil
sample_data(genus_data)$compartment <- factor(sample_data(genus_data)$compartment, levels = 
                                         c("B" ,"L", "T"))

### Run ancombc function
out = ancombc(phyloseq = genus_data, formula = var1,
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,
              group = var1, struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = FALSE)
res = out$res

### Coefficients
tab_coef = res$W
colnames(tab_coef) 
col_name = c("B - L", "B - T")
colnames(tab_coef) = col_name
source("r_scripts/ancom_bc_rename_variables.R", echo = T, spaced = T)
tmp <- tab_diff
bulk <- merge(tmp, tab_w, by=0) 
rownames(bulk) <- bulk$Row.names
bulk$Row.names = NULL


# Loosely
sample_data(genus_data)$compartment <- factor(sample_data(genus_data)$compartment, levels = 
                                                c("L" ,"B", "T"))

### Run ancombc function
out = ancombc(phyloseq = genus_data, formula = var1,
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,
              group = var1, struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = FALSE)
res = out$res

### Coefficients
tab_coef = res$W
colnames(tab_coef) 
col_name = c("L - B", "L - T")
colnames(tab_coef) = col_name
source("r_scripts/ancom_bc_rename_variables.R", echo = T, spaced = T)
tmp <- tab_diff
loosely <- merge(tmp, tab_w, by=0) 
rownames(loosely) <- loosely$Row.names
loosely$Row.names = NULL
loosely <- loosely[,-c(1,3)]


tmp <- merge(bulk, loosely, by=0); rownames(tmp) <- tmp$Row.names; tmp$Row.names = NULL
tmp3 <- as.data.frame(tax_table(genus_data))
tmp4 <- merge(tmp3, tmp, by=0); rownames(tmp4) <- tmp4$Row.names; tmp4$Row.names = NULL
tmp4 <- tmp4[,-c(1,5:6)]
head(tmp4)


PGroup <- transform_sample_counts(genus_data, function(x)100* x / sum(x))
OTUg <- otu_table(PGroup)
TAXg <- tax_table(PGroup)
AverageD <- as.data.frame(rowMeans(OTUg))
names(AverageD) <- c("Mean")
SD <- as.data.frame(rowSds(OTUg),na.rm = T)
names(SD) <- c("SD")
tmp <- cbind(AverageD,SD)
GTable <- merge(TAXg, tmp, by=0, all=TRUE); rownames(GTable) <- GTable$Row.names; GTable$Row.names = NULL
GTable <- GTable[,-c(1,5:6)]
head(GTable)


tmp5 <- merge(tmp4, GTable, by=0); row.names(tmp5) <- tmp5$Row.names; tmp5$Row.names = NULL
tmp6 <- tmp5[,-c(10:12)]
temp <- tmp6[order(tmp6$Mean, decreasing = TRUE),]
View(temp)
write_csv(temp, file = "PHPP Application effects/Data/temporal_top_families.csv")


## Concentration: -----------------------------------------------------------------
ps <- conc
pseq <- ps 
genus_data <- tax_glom(pseq, taxrank <- rank_names(pseq)[4], NArm = FALSE)
var1 <- "Compartment"

# Bulk Soil
sample_data(genus_data)$Compartment <- factor(sample_data(genus_data)$Compartment, levels = 
                                                c("B" ,"L", "T"))

### Run ancombc function
out = ancombc(phyloseq = genus_data, formula = var1,
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,
              group = var1, struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = FALSE)
res = out$res

### Coefficients
tab_coef = res$W
colnames(tab_coef) 
col_name = c("B - L", "B - T")
colnames(tab_coef) = col_name
source("r_scripts/ancom_bc_rename_variables.R", echo = T, spaced = T)
tmp <- tab_diff
bulk <- merge(tmp, tab_w, by=0) 
rownames(bulk) <- bulk$Row.names
bulk$Row.names = NULL


# Loosely
sample_data(genus_data)$Compartment <- factor(sample_data(genus_data)$Compartment, levels = 
                                                c("L" ,"B", "T"))

### Run ancombc function
out = ancombc(phyloseq = genus_data, formula = var1,
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,
              group = var1, struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = FALSE)
res = out$res

### Coefficients
tab_coef = res$W
colnames(tab_coef) 
col_name = c("L - B", "L - T")
colnames(tab_coef) = col_name
source("r_scripts/ancom_bc_rename_variables.R", echo = T, spaced = T)
tmp <- tab_diff
loosely <- merge(tmp, tab_w, by=0) 
rownames(loosely) <- loosely$Row.names
loosely$Row.names = NULL
loosely <- loosely[,-c(1,3)]


tmp <- merge(bulk, loosely, by=0); rownames(tmp) <- tmp$Row.names; tmp$Row.names = NULL
tmp3 <- as.data.frame(tax_table(genus_data))
tmp4 <- merge(tmp3, tmp, by=0); rownames(tmp4) <- tmp4$Row.names; tmp4$Row.names = NULL
tmp4 <- tmp4[,-c(1,5:6)]
head(tmp4)


PGroup <- transform_sample_counts(genus_data, function(x)100* x / sum(x))
OTUg <- otu_table(PGroup)
TAXg <- tax_table(PGroup)
AverageD <- as.data.frame(rowMeans(OTUg))
names(AverageD) <- c("Mean")
SD <- as.data.frame(rowSds(OTUg),na.rm = T)
names(SD) <- c("SD")
tmp <- cbind(AverageD,SD)
GTable <- merge(TAXg, tmp, by=0, all=TRUE); rownames(GTable) <- GTable$Row.names; GTable$Row.names = NULL
GTable <- GTable[,-c(1,5:6)]
head(GTable)


tmp5 <- merge(tmp4, GTable, by=0); row.names(tmp5) <- tmp5$Row.names; tmp5$Row.names = NULL
tmp6 <- tmp5[,-c(10:12)]
conc <- tmp6[order(tmp6$Mean, decreasing = TRUE),]
View(conc)
write_csv(conc, file = "PHPP Application effects/Data/concentration_top_families.csv")



## Strawberry: -----------------------------------------------------------------
ps <- straw
pseq <- ps 
genus_data <- tax_glom(pseq, taxrank <- rank_names(pseq)[4], NArm = FALSE)
var1 <- "Compartment"

# Bulk Soil
sample_data(genus_data)$Compartment <- factor(sample_data(genus_data)$Compartment, levels = 
                                                c("B" ,"L", "T"))

### Run ancombc function
out = ancombc(phyloseq = genus_data, formula = var1,
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,
              group = var1, struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = FALSE)
res = out$res

### Coefficients
tab_coef = res$W
colnames(tab_coef) 
col_name = c("B - L", "B - T")
colnames(tab_coef) = col_name
source("r_scripts/ancom_bc_rename_variables.R", echo = T, spaced = T)
tmp <- tab_diff
bulk <- merge(tmp, tab_w, by=0) 
rownames(bulk) <- bulk$Row.names
bulk$Row.names = NULL


# Loosely
sample_data(genus_data)$Compartment <- factor(sample_data(genus_data)$Compartment, levels = 
                                                c("L" ,"B", "T"))

### Run ancombc function
out = ancombc(phyloseq = genus_data, formula = var1,
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 10000,
              group = var1, struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = FALSE)
res = out$res

### Coefficients
tab_coef = res$W
colnames(tab_coef) 
col_name = c("L - B", "L - T")
colnames(tab_coef) = col_name
source("r_scripts/ancom_bc_rename_variables.R", echo = T, spaced = T)
tmp <- tab_diff
loosely <- merge(tmp, tab_w, by=0) 
rownames(loosely) <- loosely$Row.names
loosely$Row.names = NULL
loosely <- loosely[,-c(1,3)]


tmp <- merge(bulk, loosely, by=0); rownames(tmp) <- tmp$Row.names; tmp$Row.names = NULL
tmp3 <- as.data.frame(tax_table(genus_data))
tmp4 <- merge(tmp3, tmp, by=0); rownames(tmp4) <- tmp4$Row.names; tmp4$Row.names = NULL
tmp4 <- tmp4[,-c(1,5:6)]
head(tmp4)


PGroup <- transform_sample_counts(genus_data, function(x)100* x / sum(x))
OTUg <- otu_table(PGroup)
TAXg <- tax_table(PGroup)
AverageD <- as.data.frame(rowMeans(OTUg))
names(AverageD) <- c("Mean")
SD <- as.data.frame(rowSds(OTUg),na.rm = T)
names(SD) <- c("SD")
tmp <- cbind(AverageD,SD)
GTable <- merge(TAXg, tmp, by=0, all=TRUE); rownames(GTable) <- GTable$Row.names; GTable$Row.names = NULL
GTable <- GTable[,-c(1,5:6)]
head(GTable)


tmp5 <- merge(tmp4, GTable, by=0); row.names(tmp5) <- tmp5$Row.names; tmp5$Row.names = NULL
tmp6 <- tmp5[,-c(10:12)]
straw <- tmp6[order(tmp6$Mean, decreasing = TRUE),]
View(straw)
write_csv(straw, file = "PHPP Application effects/Data/strawberry_top_families.csv")


# MERGE: -----------------------------------------------------------------
temp <- read.csv("PHPP Application effects/Data/temporal_top_families.csv")
conc <- read.csv("PHPP Application effects/Data/concentration_top_families.csv")
straw <- read.csv("PHPP Application effects/Data/strawberry_top_families.csv")


for (trial in c("temporal", "concentration","strawberry")){
  
  tax <- "Order"
  ###ANCOM results
  ancom <- read.csv(file = paste0("PHPP Application effects/Data/",trial,"_top_families.csv"))
  ancom$Phylum <- ancom$P.x
  ancom <- merge(ancom, p.colors, by="Phylum", all.y = T, )
  ancom$O.x[ancom$O.x=="uncultured"] <- "Unclassified"
  ancom$O.x[is.na(ancom$O.x)] <- "Unclassified"
  ancom <- add_column(ancom, Name = ancom$O.x, .after = "O.x")
  #select only root section columns
  write_csv(ancom[,-2], file = paste0("PHPP Application effects/Data/",trial,"_top_families_mod.csv"))
}

#merge all and then redo in excel!!
temp <- read.csv("PHPP Application effects/Data/temporal_top_families_mod.csv")[,-c(2,3)]
conc <- read.csv("PHPP Application effects/Data/concentration_top_families_mod.csv")[,-c(2,3)]
straw <- read.csv("PHPP Application effects/Data/strawberry_top_families_mod.csv")[,-c(2,3)]

colnames(temp) <- c("Phylum","Name","B:L_sig","B:T_sig","B:L_W","B:T_W","L:T_sig","L:T_W","Temporal Mean","Temporal SD", "P.color")
colnames(conc) <- c("Phylum","Name","B:L_sig","B:T_sig","B:L_W","B:T_W","L:T_sig","L:T_W","Concentration Mean","Concentration SD", "P.color")
colnames(straw) <- c("Phylum","Name","B:L_sig","B:T_sig","B:L_W","B:T_W","L:T_sig","L:T_W","Strawberry Mean","Strawberry SD", "P.color")

tmp <- merge(temp, conc, by="Name", all.y=T, all.x =T)
tmp <- merge(tmp, straw, by="Name", all.y=T, all.x =T)
tmp$Average <- as.numeric(rowMeans(tmp[c(9,19,29)],na.rm = TRUE))
tmp <- na.omit(tmp[tmp$Average>0.5,])
write_csv(tmp, file = "PHPP Application effects/Data/top_families_mod.csv")
#redo in excel!!!!

# Load final data ---------------------------------------------------------

#change in excel
tmp2 <- read.csv("PHPP Application effects/Data/top_families_mod.csv")
tmp3 <- merge(tmp2, p.colors, by="Phylum"); row.names(tmp3) <- tmp3$Name; tmp3$Name = NULL;
ancom <- tmp3[(rowSums(tmp3[,c(2,3,6,8,9,12,14,15,18)])!=0),]

# Heatmap: -----------------------------------------------------------------
mat.df <- ancom
str(mat.df)
cat(paste(shQuote(colnames(mat.df), type="cmd"), collapse=", "))

colnames(mat.df) <- c("Phylum", 
                      "B.L_sig.x", "B.T_sig.x", "Temporal: B->L", "Temporal: B->T", "L.T_sig.x", "Temporal: L->T", 
                      "B.L_sig.y", "B.T_sig.y", "Concentration: B->L", "Concentration: B->T", "L.T_sig.y", "Concentration: L->T", 
                      "B.L_sig", "B.T_sig", "Strawberry: B->L", "Strawberry: B->T", "L.T_sig", "Strawberry: L->T", 
                      "Average", "P.color")

#IMPORTANT: specify ONLY the columns with the differentials
col.order <- c("Temporal: B->L", "Concentration: B->L", "Strawberry: B->L",
               "Temporal: L->T", "Concentration: L->T", "Strawberry: L->T",
               "Temporal: B->T", "Concentration: B->T", "Strawberry: B->T")

mat.diff <- as.matrix(mat.df[,c(4,5,7,10,11,13,16,17,19)])
sig_mat <- as.matrix(mat.df[,c(2,3,6,8,9,12,14,15,18)])
sig_mat[is.na(sig_mat)] <- FALSE
min_lc <- min(mat.diff, na.rm = T)
max_lc <- max(mat.diff, na.rm = T)
colors <- colorRamp2(c(min_lc, 0, max_lc), c("blue", "white", "red"))

hb1 = rowAnnotation("Mean\n Abundance [%]" = anno_barplot(as.vector(mat.df$Average), gp = gpar(fill = "black"), bar_width = 0.9, extend=0.1, 
                                                       axis_param = list(side = "bottom", labels_rot = 90, gp = gpar(fontsize = 24)), 
                                                       width = unit(10, "cm")), annotation_name_gp = gpar(fontsize = 24, fontface = "bold"),
                   annotation_name_side = "top", annotation_name_offset = unit(0.3, "cm"), annotation_name_rot = 0)

lgd = Legend(at = c(round(min_lc+min_lc*0.05,digits = 1), 0, round(max_lc+max_lc*0.05, digits = 1)),  col_fun = colors, 
             title = "W-value", title_gp = gpar(fontsize = 32, fontface = "bold"), title_gap = unit(0.5, "cm"),
             labels_gp = gpar(col = "black", fontsize = 28, fontface = "bold"), title_position = "topcenter", 
             grid_height = unit(3, "cm"), legend_width = unit(10, "cm"), direction = "horizontal") 

ht = Heatmap(mat.diff , name = "logfold change",
             column_gap = unit(5, "mm"),
             row_gap = unit(5, "mm"), 
             row_split = mat.df$Phylum,
             row_names_gp = gpar(fontsize = 32, col = mat.df$P.color),
             na_col = "grey",
             column_order = col.order,
             col = colors,
             rect_gp = gpar(col = "white", lwd = 2), #white spaces in between
             border = TRUE, 
             cluster_rows = FALSE, #remove cluster         
             show_column_dend = FALSE, #remove cluster 
             show_heatmap_legend = FALSE, #remove legend
             row_names_side = "right",
             row_names_max_width = unit(2, "cm"),
             row_names_rot = 0,
             row_names_centered = FALSE,
             row_title_gp = gpar(fontsize = 36, fontface = "bold"),
             row_title_rot = 0,
             column_title = paste0("ANCOM-BC"), 
             column_title_gp = gpar(fill = "black", col = "white", fontsize = 28, fontface = "bold", border = "black"), 
             column_title_side = "top",
             column_names_max_height = unit(6, "cm"),
             column_names_gp = gpar(fontsize = 32, fontface = "bold", col = "black"),
             column_names_rot = 90,
             column_names_centered = TRUE,
             show_parent_dend_line = FALSE, 
             cell_fun = function(j, i, x, y, width, height, fill) {
                if(sig_mat[i, j] == "TRUE")
                  grid.text("*", x, y, gp = gpar(fontsize = 30, fontface = "bold"))
             }, 
             right_annotation = c(hb1))
png("PHPP Application effects/Plots/BLT-comparison_order.png", width = 2400, height = 600+45*(nrow(mat.diff)))
draw(ht, padding = unit(c(5, 1, 1, 30), "cm")) # add space for titles
draw(lgd, x = unit(55, "cm"), y = unit(4, "cm"))
dev.off()

tiff("PHPP Application effects/Plots/BLT-comparison_order.png", units="in", width=30, height=5+0.5*(nrow(mat.diff)), pointsize = 12, res=300)
draw(ht, padding = unit(c(1.5,1,1,10.5), "in")) # add space for titles
draw(lgd, x = unit(14, "in"), y = unit(1.9, "in"))
dev.off()
#



