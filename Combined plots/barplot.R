rm(list=ls())
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

{
  source("/Users/Max/Science/Amplicon/r_scripts/4_stackbar_start.R", echo = T, spaced = T)  
  library(microbiome)
  library(viridis)
  library(forcats)
}

#made for macos
# Phylum  -----------------------------
temp <- readRDS("/Users/Max/Science/Amplicon/PHPP_temp/phyloseq/ps.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "P", NArm = FALSE)
conc <- readRDS("/Users/Max/Science/Amplicon/PHPP_conc/phyloseq/ps.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("Compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "P", NArm = FALSE)
straw <- readRDS("/Users/Max/Science/Amplicon/PHPP_strawberry/phyloseq/ps.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("Compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "P", NArm = FALSE)


dat_taxa <- data.table(psmelt(temp))
dat_taxa$P <- as.character(dat_taxa$P)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "P"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:17)]
dat_taxa$trial <- as.character("temp")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder1 <- dat_taxa[(mean <= 1), P := "Other"]
remainder1 <- remainder1[,-3]

dat_taxa <- data.table(psmelt(conc))
dat_taxa$P <- as.character(dat_taxa$P)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "P"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:12)]
dat_taxa$trial <- as.character("conc")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder2 <- dat_taxa[(mean <= 1), P := "Other"]

dat_taxa <- data.table(psmelt(straw))
dat_taxa$P <- as.character(dat_taxa$P)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "P"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:16)]
dat_taxa$trial <- as.character("straw")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder3 <- dat_taxa[(mean <= 1), P := "Other"]

#merge
remove(remainder)
remainder <- rbind(remainder1, remainder2, remainder3)

remainder[remainder$P=="Proteobacteria",]

####################### Order -sorted ####################### 
temp <- readRDS("/Users/Max/Science/Amplicon/PHPP_temp//phyloseq/ps.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "O", NArm = FALSE)
conc <- readRDS("/Users/Max/Science/Amplicon/PHPP_conc//phyloseq/ps.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("Compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "O", NArm = FALSE)
straw <- readRDS("/Users/Max/Science/Amplicon/PHPP_strawberry//phyloseq/ps.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("Compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "O", NArm = FALSE)


dat_taxa <- data.table(psmelt(temp))
dat_taxa$O <- as.character(dat_taxa$O)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "O"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:17)]
dat_taxa$trial <- as.character("temp")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder1 <- dat_taxa[(mean <= 1), O := "Other"]
remainder1 <- remainder1[,-3]

dat_taxa <- data.table(psmelt(conc))
dat_taxa$O <- as.character(dat_taxa$O)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "O"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:12)]
dat_taxa$trial <- as.character("conc")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder2 <- dat_taxa[(mean <= 1), O := "Other"]

dat_taxa <- data.table(psmelt(straw))
dat_taxa$O <- as.character(dat_taxa$O)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "O"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:16)]
dat_taxa$trial <- as.character("straw")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder3 <- dat_taxa[(mean <= 1), O := "Other"]

#merge
remove(remainder)
remainder <- rbind(remainder1, remainder2, remainder3)
#remainder <-remainder[!(remainder$axis=="Spat-B"),]
#remainder <-remainder[!(remainder$axis=="DLR-B"),]
#remainder <-remainder[!(remainder$axis=="Spat-B"),]
remainder$axis <- mapvalues(remainder$axis, from = c("temp-L", "temp-T", "temp-B", "conc-T", "conc-L", "conc-B", "straw-L", "straw-T", "straw-B"),
                            to = c("Temporal-L", "Temporal-T", "Temporal-B", "Concentration-T", "Concentration-L", "Concentration-B", "Strawberry-L", "Strawberry-T", "Strawberry-B"))
level_order1 <- factor(remainder$axis, level = c("Temporal-B","Concentration-B","Strawberry-B", "Temporal-L","Concentration-L","Strawberry-L", "Temporal-T","Concentration-T","Strawberry-T"))
level_order2 <- factor(remainder$axis, level = c("Temporal-B","Temporal-L","Temporal-T",
                                                 "Concentration-B","Concentration-L","Concentration-T",
                                                 "Strawberry-B","Strawberry-L","Strawberry-T"))

## 1 percent ######
remove(df)
remove(tmp1)
remove(tmp2)
tmp1 <- remainder
tmp2 <- tmp1#[(mean <= 1), O := "Other"]

df <- tmp2
names(df)[names(df) == "P"] <- "Phylum"
names(df)[names(df) == "O"] <- "Order"
df$group <- paste0(df$Phylum, "-", df$Order, sep = "")
sort(unique(df$group))
cat(paste(shQuote(sort(unique(df$Phylum)), type="cmd"), collapse=", "))

#test
phylums <- c("Acidobacteriota","Actinobacteriota","Bacteroidota","Dependentiae","Firmicutes",
             "Myxococcota","Proteobacteria")


df$Order[df$Phylum=="Acidobacteriota" & df$Order=='Acidobacteriota-Other'] <- "Other Acidobacteria"

df$Order[df$Phylum=="Actinobacteriota" & df$Order=='Actinobacteriota-Other'] <- "Other Actinobacteriota"

df$Order[df$Phylum=="Bacteroidota" & df$Order=='Bacteroidota-Other'] <- "Other Bacteroidota"

df$Order[df$Phylum=="Dependentiae" & df$Order=='Dependentiae-Other'] <- "Other Dependentiae"

df$Order[df$Phylum=="Firmicutes" & df$Order=='Firmicutes-Other'] <- "Other Firmicutes"

df$Order[df$Phylum=="Proteobacteria" & df$Order=='Proteobacteria-Other'] <- "Other Proteobacteria"

df$group[df$group=='Unclassified-Other'] <- "Unclassified"

df$group[!df$Phylum %in% phylums] <- "Other Phyla"


df2 <- dplyr::select(df, axis, Phylum, Order, Abundance, group) %>%
  dplyr::mutate(Phylum=factor(Phylum, levels=c(phylums, "Others")),
         Order=fct_reorder(Order, 10*as.integer(Phylum) + grepl("Others", Order)))

sort(unique(df2$group))
colours <- c("green4","chartreuse1","chartreuse3","chartreuse4","green1", #Acidobacteriota
             "lightgoldenrod2","chocolate1","chocolate4","lightsalmon1", #Actinobacteriota
             "royalblue1","royalblue4", #Bacteroidota
             "seagreen1","seagreen4", #Dependentiae
             "hotpink1","hotpink3","hotpink4", #Firmicutes
             "blue","blue",#Myxococcota
             "khaki", #Other
             plasma(11), #Proteobacteria
             "black") #Unclassified

png("/Users/Max/Science/Amplicon/PHPP Application effects/Plots/prozentuale_abundanz_stackbar_16S_order1.png", width = 1600, height = 1200)
tiff("/Users/Max/Science/Amplicon/PHPP Application effects/Plots/three_trials_class.tiff", units="in", width=16, height=12, pointsize = 12, res=300)
ggplot(df2, aes(x=level_order1, y=Abundance, fill=group, order=group)) + 
  geom_bar(aes(fill=group), stat="identity", position="stack") + 
  mytheme + theme(axis.text.x = element_text(angle=75, vjust = 0.5)) +
  scale_fill_manual("", values=colours) +
  scale_x_discrete("") +
  theme(axis.text.x=element_text(angle=75, vjust=0.5, size=30)) +  # Vertical x-axis tick labels
  scale_y_continuous("Relative abundance [%]", breaks=seq(0,100,5), limits = c(0, 101)) +
  labs(y="Relative abundance") +
  guides(fill=guide_legend(nrow=length(unique(df2$group)))) + labs(fill = "Family") +
  geom_vline(xintercept = c(3.5), color = "black", size=2.5)
dev.off()


## 2 percent #####
remove(df)
remove(tmp1)
remove(tmp2)
tmp1 <- remainder
tmp2 <- tmp1[(mean <= 2), O := "Other"]

df <- tmp2
names(df)[names(df) == "P"] <- "Phylum"
names(df)[names(df) == "O"] <- "Order"
df$group <- paste0(df$Phylum, "-", df$Order, sep = "")
sort(unique(df$group))
cat(paste(shQuote(sort(unique(df$Phylum)), type="cmd"), collapse=", "))

#test
phylums <- c("Acidobacteriota","Actinobacteriota","Firmicutes","Proteobacteria")

df$Order[df$Phylum=="Acidobacteriota" & df$Order=='Acidobacteriota-Other'] <- "Other Acidobacteria"

df$Order[df$Phylum=="Actinobacteriota" & df$Order=='Actinobacteriota-Other'] <- "Other Actinobacteriota"

df$Order[df$Phylum=="Firmicutes" & df$Order=='Firmicutes-Other'] <- "Other Firmicutes"

df$Order[df$Phylum=="Proteobacteria" & df$Order=='Proteobacteria-Other'] <- "Other Proteobacteria"

df$group[df$group=='Unclassified-Other'] <- "Unclassified"

df$group[!df$Phylum %in% phylums] <- "Other Phyla"


df2 <- dplyr::select(df, axis, Phylum, Order, Abundance, group) %>%
  dplyr::mutate(Phylum=factor(Phylum, levels=c(phylums, "Others")),
                Order=fct_reorder(Order, 10*as.integer(Phylum) + grepl("Others", Order)))

sort(unique(df2$group))
colours <- c("chartreuse1","chartreuse2","chartreuse3","chartreuse4", #Acidobacteriota
             "chocolate1","chocolate3","chocolate4", #Actinobacteriota
             "hotpink1","hotpink3", #Firmicutes
             "khaki", #Other
             plasma(8)) #Proteobacteria

png("/Users/Max/OneDrive/Amplicon/PHPP Application effects/Plots/prozentuale_abundanz_stackbar_16S_order2.png", width = 1600, height = 1200)
tiff("/Users/Max/Science/Amplicon/PHPP Application effects/Plots/three_trials_class.tiff", units="in", width=16, height=12, pointsize = 12, res=300)
ggplot(df2, aes(x=level_order1, y=Abundance, fill=group, order=group)) + 
  geom_bar(aes(fill=group), stat="identity", position="stack") + 
  mytheme + theme(axis.text.x = element_text(angle=75, vjust = 0.5)) +
  scale_fill_manual("", values=colours) +
  scale_x_discrete("") + 
  theme(axis.text.x = element_text(angle=75, vjust=0.5, size=24), 
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +  
  scale_y_continuous("Relative abundance [%]", breaks=seq(0,100,5), limits = c(0, 101)) +
  guides(fill=guide_legend(nrow=length(unique(df2$group)))) + labs(fill = "Family") +
  geom_vline(xintercept = c(3.5,6.5), color = "black", size=2.5)
dev.off()


# Family #####
####################### Family -sorted ####################### 
DLR <- readRDS("/Users/Max/Science/Amplicon/DLR_temporal/phyloseq/ps_fp.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "F", NArm = FALSE)
K <- readRDS("/Users/Max/Science/Amplicon/V4a_temp/phyloseq/ps_fp.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "F", NArm = FALSE)
spat <- readRDS("/Users/Max/Science/Amplicon/V4a_spat/phyloseq/ps_fp.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "F", NArm = FALSE)

dat_taxa <- data.table(psmelt(K))
dat_taxa$F <- as.character(dat_taxa$F)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "F"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:19)]
dat_taxa$trial <- as.character("K")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder1 <- dat_taxa[(mean <= 2), F := "Other"]

dat_taxa <- data.table(psmelt(DLR))
dat_taxa$F <- as.character(dat_taxa$F)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "F"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:20)]
dat_taxa$trial <- as.character("DLR")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder2 <- dat_taxa[(mean <= 2), F := "Other"]
remainder2 <- remainder2[,-c(3:5)]

# dat_taxa <- data.table(psmelt(spat))
# dat_taxa$F <- as.character(dat_taxa$F)
# medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "F"]
# dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
# dat_taxa <- dat_taxa[,-c(1,4:13)]
# dat_taxa$trial <- as.character("Spat")
# dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
# remainder3 <- dat_taxa[(mean <= 2), F := "Other"]
dat_taxa <- data.table(psmelt(spat))
dat_taxa <- dat_taxa[,-c(1,4:16)]
colnames(dat_taxa) <- c("Sample","Abundance","P","C","O","F")
dat_taxa$F <- as.character(dat_taxa$F)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "F"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa$trial <- as.character("Spat")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder3 <- dat_taxa[(mean <= 2), F := "Other"]

remove(remainder)
remainder <- rbind(remainder1, remainder2, remainder3)
remainder <-remainder[!(remainder$axis=="Spat-B"),]
remainder <-remainder[!(remainder$axis=="DLR-B"),]

remainder$axis <- mapvalues(remainder$axis, from = c("Spat-L","K-L","DLR-L","Spat-T","K-T","DLR-T"),
                            to = c("Spatial - L","Temporal - L","ST - L","Spatial - T","Temporal - T","ST - T"))
level_order <- factor(remainder$axis, level = c("Spatial - L","Spatial - T","Temporal - L","Temporal - T","ST - L","ST - T"))

remove(df)
remove(tmp1)
remove(tmp2)
tmp1 <- remainder
tmp2 <- tmp1[(mean <= 2), F := "Other"]


ColourPalleteMulti <- function(df, group, subgroup){
  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
  category.start <- (scales::hue_pal(l = 100)(nrow(categories))) # Set the top of the colour pallete
  category.end  <- (scales::hue_pal(l = 20)(nrow(categories))) # set the bottom
  
  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
  return(colours)
}

df <- tmp2
names(df)[names(df) == "P"] <- "Phylum"
names(df)[names(df) == "F"] <- "Family"
df$group <- paste0(df$Phylum, "-", df$Family, sep = "")
sort(unique(df$group))
cat(paste(shQuote(sort(unique(df$Phylum)), type="cmd"), collapse=", "))

#test
phylums <- c("Acidobacteriota", "Actinobacteriota","Bacteroidota", "Dependentiae",
             "Firmicutes", "Proteobacteria", "Unclassified")


df$Family[df$Phylum=="Acidobacteriota" & df$Family=='Acidobacteriota-Other'] <- "Other Acidobacteria"

df$Family[df$Phylum=="Actinobacteriota" & df$Family=='Actinobacteriota-Other'] <- "Other Actinobacteriota"

df$Family[df$Phylum=="Bacteroidota" & df$Family=='Bacteroidota-Other'] <- "Other Bacteroidota"

df$Family[df$Phylum=="Dependentiae" & df$Family=='Dependentiae-Other'] <- "Other Dependentiae"

df$group[df$group=='Firmicutes-Other'] <- "Firmicutes"

df$Family[df$Phylum=="Proteobacteria" & df$Family=='Proteobacteria-Other'] <- "Other Proteobacteria"

df$group[df$group=='Unclassified-Other'] <- "Unclassified"

df$group[!df$Phylum %in% phylums] <- "Other Phyla"

# df$Family[df$Phylum=="Proteobacteria" & 
#            !df$Family %in% c('Acetobacterales','Betaproteobacteriales','Caulobacterales','Enterobacteriales',
#                             'Myxococcales','Pseudomonadales','Rhizobiales','Sphingomonadales','Xanthomonadales')] <- "Other Protobacteria"


df2 <- select(df, axis, Phylum, Family, Abundance, group) %>%
  mutate(Phylum=factor(Phylum, levels=c(phylums, "Others")),
         Family=fct_reorder(Family, 10*as.integer(Phylum) + grepl("Others", Family)))
df2$group[df2$group=="Proteobacteria-UnknownFamily"] <- "Proteobacteria-Unknown_Family"
df2$group[df2$group=="Proteobacteria-Unknown_Family"] <-"Proteobacteria-Unknown_Family"

sort(unique(df2$group))
colours <- c("chartreuse1","chartreuse3", #Acidobacteriota
             "chocolate1","chocolate4", #Actinobacteriota
             "royalblue1","royalblue3","royalblue4", #Bacteroidota
             "seagreen1","seagreen4", #Dependentiae
             "black", #Firmicutes
             "khaki", #Other
             plasma(13), #Proteobacteria
             "black") #Unclassified

png("/Users/Max/OneDrive/Amplicon/Spatio-temporal variation/Plots/prozentuale_abundanz_stackbar_16S_compartmentpsm_Family.png", width = 1600, height = 1000)
tiff("/Users/Max/OneDrive/Amplicon/Spatio-temporal variation/Plots/prozentuale_abundanz_stackbar_16S_compartmentpsm_Family.tiff", units="in", width=20, height=12, pointsize = 12, res=300)
ggplot(df2, aes(x=level_order, y=Abundance, fill=group, order=group)) + 
  geom_bar(aes(fill=group), stat="identity", position="stack") + 
  mytheme + theme(axis.text.x = element_text(angle=75, vjust = 0.5)) +
  scale_fill_manual("", values=colours) +
  scale_x_discrete("") +
  theme(
    axis.text.x = element_text(angle=55, vjust=0.5, size=30, colour = "black"), 
    axis.text.y = element_text(size = 22, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 26, face = "bold", colour = "black")
    ) +  # Vertical x-axis tick labels
  scale_y_continuous("Relative Abundance [%]", breaks=seq(0,100,5), limits = c(0, 101), expand = c(0,2)) +
  guides(fill=guide_legend(nrow=length(unique(df2$group)))) + labs(fill = "Family") +
  geom_vline(xintercept = c(2.5,4.5), color = "black", size=3)
dev.off()



# Genus  -----------------------------
temp <- readRDS("/Users/Max/Science/Amplicon/PHPP_temp/phyloseq/ps.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "G", NArm = FALSE)
conc <- readRDS("/Users/Max/Science/Amplicon/PHPP_conc/phyloseq/ps.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("Compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "G", NArm = FALSE)
straw <- readRDS("/Users/Max/Science/Amplicon/PHPP_strawberry/phyloseq/ps.RDS") %>%
  microbiome::transform("compositional") %>%
  merge_samples("Compartment") %>%
  transform_sample_counts(function(x) 100 * x/sum(x)) %>%
  tax_glom(taxrank = "G", NArm = FALSE)


dat_taxa <- data.table(psmelt(temp))
dat_taxa$G <- as.character(dat_taxa$G)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "G"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:17)]
dat_taxa$trial <- as.character("temp")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder1 <- dat_taxa[(mean <= 2), G := "Other"]
remainder1 <- remainder1[,-3]

dat_taxa <- data.table(psmelt(conc))
dat_taxa$G <- as.character(dat_taxa$G)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "G"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:12)]
dat_taxa$trial <- as.character("conc")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder2 <- dat_taxa[(mean <= 2), G := "Other"]

dat_taxa <- data.table(psmelt(straw))
dat_taxa$G <- as.character(dat_taxa$G)
medians <- dat_taxa[, mean := mean(Abundance, na.rm = TRUE), by = "G"]
dat_taxa[is.na(dat_taxa)] <- "Unclassified" 
dat_taxa <- dat_taxa[,-c(1,4:16)]
dat_taxa$trial <- as.character("straw")
dat_taxa$axis <- as.character(paste0(dat_taxa$trial,"-",dat_taxa$Sample))
remainder3 <- dat_taxa[(mean <= 2), G := "Other"]

#merge
remove(remainder)
remainder <- rbind(remainder1, remainder2, remainder3)
remainder$axis <- mapvalues(remainder$axis, from = c("temp-L", "temp-T", "temp-B", "conc-T", "conc-L", "conc-B", "straw-L", "straw-T", "straw-B"),
                            to = c("Temporal-L", "Temporal-T", "Temporal-B", "Concentration-T", "Concentration-L", "Concentration-B", "Strawberry-L", "Strawberry-T", "Strawberry-B"))
level_order1 <- factor(remainder$axis, level = c("Temporal-B","Concentration-B","Strawberry-B", "Temporal-L","Concentration-L","Strawberry-L", "Temporal-T","Concentration-T","Strawberry-T"))
level_order2 <- factor(remainder$axis, level = c("Temporal-B","Temporal-L","Temporal-T",
                                                 "Concentration-B","Concentration-L","Concentration-T",
                                                 "Strawberry-B","Strawberry-L","Strawberry-T"))

remove(df)
remove(tmp1)
remove(tmp2)
tmp1 <- remainder
tmp2 <- tmp1#[(mean <= 1), O := "Other"]

df <- tmp2
names(df)[names(df) == "P"] <- "Phylum"
names(df)[names(df) == "G"] <- "Genus"
df$group <- paste0(df$Phylum, "-", df$Genus, sep = "")
sort(unique(df$group))
cat(paste(shQuote(sort(unique(df$Phylum)), type="cmd"), collapse=", "))

#test
phylums <- c("Acidobacteriota","Actinobacteriota","Bacteroidota","Dependentiae","Firmicutes",
             "Myxococcota","Proteobacteria")

df$Genus[df$Phylum=="Acidobacteriota" & df$Genus=='Acidobacteriota-Other'] <- "Other Acidobacteria"

df$Genus[df$Phylum=="Actinobacteriota" & df$Genus=='Actinobacteriota-Other'] <- "Other Actinobacteriota"

df$Genus[df$Phylum=="Bacteroidota" & df$Genus=='Bacteroidota-Other'] <- "Other Bacteroidota"

df$Genus[df$Phylum=="Dependentiae" & df$Genus=='Dependentiae-Other'] <- "Other Dependentiae"

df$Genus[df$Phylum=="Firmicutes" & df$Genus=='Firmicutes-Other'] <- "Other Firmicutes"

df$Genus[df$Phylum=="Proteobacteria" & df$Genus=='Proteobacteria-Other'] <- "Other Proteobacteria"

df$group[df$group=='Unclassified-Other'] <- "Unclassified"

df$group[!df$Phylum %in% phylums] <- "Other Phyla"


df2 <- dplyr::select(df, axis, Phylum, Genus, Abundance, group) %>%
  dplyr::mutate(Phylum=factor(Phylum, levels=c(phylums, "Others")),
                Genus=fct_reorder(Genus, 10*as.integer(Phylum) + grepl("Others", Genus)))

sort(unique(df2$group))
colours <- c("green4","chartreuse1","chartreuse4",  #Acidobacteriota
             "chocolate1","chocolate4", #Actinobacteriota
             "royalblue4", #Bacteroidota
             "seagreen4", #Dependentiae
             "hotpink1","hotpink4", #Firmicutes
             "blue1",#Myxococcota
             "black", #Other
             plasma(10)) #Proteobacteria

png("/Users/Max/Science/Amplicon/PHPP Application effects/Plots/prozentuale_abundanz_stackbar_16S_genus.png", width = 1600, height = 1200)
tiff("/Users/Max/Science/Amplicon/PHPP Application effects/Plots/three_trials_genus.tiff", units="in", width=16, height=12, pointsize = 12, res=300)
ggplot(df2, aes(x=level_order1, y=Abundance, fill=group, order=group)) + 
  geom_bar(aes(fill=group), stat="identity", position="stack") + 
  mytheme + theme(axis.text.x = element_text(angle=75, vjust = 0.5)) +
  scale_fill_manual("", values=colours) +
  scale_x_discrete("") + 
  theme(axis.text.x = element_text(angle=75, vjust=0.5, size=24), 
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +  
  scale_y_continuous("Relative abundance [%]", breaks=seq(0,100,5), limits = c(0, 101)) +
  guides(fill=guide_legend(nrow=length(unique(df2$group)))) + labs(fill = "Family") +
  geom_vline(xintercept = c(3.5,6.5), color = "black", size=2.5)
dev.off()
