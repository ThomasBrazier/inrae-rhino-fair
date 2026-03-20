# New analyses for R1
# clear global environment: remove all variables
rm(list=ls(all=TRUE))
# library(rstudioapi)
# Get the directory of the file & set working directory
# filedir=dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(filedir)

#----------------------------------------------------------#
# Loading packages
# Check if packages are installed, install if necessary
# source("Sources/packages.R")
source("Sources/packages_R1.R")
library(tidyverse)


#==========================================================#
# Map of sampling locations ----
#==========================================================#

dyadsObsSelect_Pic = read.table("AssignationPic/outputs/dyadsObsSelect.txt", header=TRUE)
dyadsObsSelect_Thu = read.table("AssignationThu/outputs/dyadsObsSelect.txt", header=TRUE)


# Figure S1 ----

coordCol_Pic = read.table("Data/Pic/coordPic.txt",h=T)
# Pic_sampling_names = read_tsv("Data/Pic/Table_S1_EP.tsv")
# Pic_sampling_names$Colonie = Pic_sampling_names$sampling_site
# 
# coordCol_Pic %>%
#   left_join(Pic_sampling_names[,c("Colonie", "zarzoso_sampling_names")])

# coordCol_Pic = coordCol_Pic[-which(coordCol_Pic$Colonie=="M399" | coordCol_Pic$Colonie=="M1079" | coordCol_Pic$Colonie=="CXSGT"| coordCol_Pic$Colonie=="M1975"  | coordCol_Pic$Colonie=="M1979" | coordCol_Pic$Colonie=="P412"),]

colnames(coordCol_Pic)[1] = "Colony"
allel=read.table("Data/Pic/uniqueGenotypesWithInfo.txt",h=T)
allel$idcol = as.factor(allel$idcol)
colonies = as.character(allel$idcol)
coordCol_Pic = coordCol_Pic[coordCol_Pic$Colony %in% colonies,]

dyadsObsSelect_Pic = dyadsObsSelect_Pic[-which(dyadsObsSelect_Pic$offspring == dyadsObsSelect_Pic$father),]
dyadsObsSelect_Pic$xstart = NA
dyadsObsSelect_Pic$ystart = NA
dyadsObsSelect_Pic$xend = NA
dyadsObsSelect_Pic$yend = NA
for (i in 1:nrow(dyadsObsSelect_Pic)) {
  dyadsObsSelect_Pic$xstart[i] = coordCol_Pic$Long[which(coordCol_Pic$Colony == dyadsObsSelect_Pic$offspring[i])]
  dyadsObsSelect_Pic$ystart[i] = coordCol_Pic$Lat[which(coordCol_Pic$Colony == dyadsObsSelect_Pic$offspring[i])]
  dyadsObsSelect_Pic$xend[i] = coordCol_Pic$Long[which(coordCol_Pic$Colony == dyadsObsSelect_Pic$father[i])]
  dyadsObsSelect_Pic$yend[i] = coordCol_Pic$Lat[which(coordCol_Pic$Colony == dyadsObsSelect_Pic$father[i])]
}

plot(coordCol_Pic$Long, coordCol_Pic$Lat)
text(coordCol_Pic$Long, coordCol_Pic$Lat, coordCol_Pic$Colony)

summary(coordCol_Pic$Long)
summary(coordCol_Pic$Lat)

(map1 = ggplot(coordCol_Pic, aes(x = Long, y = Lat)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = Colony),
                   size = 3.5,
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.3, "lines")) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_curve(data = dyadsObsSelect_Pic, aes(x = xstart, y = ystart, xend = xend, yend = yend),
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "darkgrey",
             alpha = 0.5,
             size = 0.5,
             curvature = -0.3
  )  +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
          axis.title.x = element_text(color="black", size=14),
          axis.title.y = element_text(color="black", size=14),
          axis.text=element_text(size=14, colour="black"),
          legend.key = element_rect(fill = "white", colour = "white", size = 1),
          legend.key.height = unit(2,"line"),
          legend.key.width = unit(2,"line"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10, face = "bold"),
          legend.position = "none"))


# Use a map package for scale and north arrow
# Example
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# (sites <- data.frame(longitude = c(-80.144005, -80.109), latitude = c(26.479005, 
#                                                                       26.83)))

# arrow = data.frame(xstart = -80.144005,
#                    ystart = 26.479005,
#                    xend = -80.109,
#                    yend = 26.83)
# 
# ggplot(data = world) +
#   geom_sf() +
#   geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
#              shape = 23, fill = "darkred") +
#   coord_sf(xlim = c(-88, -78), ylim = c(24.5, 33), expand = FALSE) +
#   geom_curve(data = arrow, aes(x = xstart, y = ystart, xend = xend, yend = yend),
#              arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
#              color = "black",
#              alpha = 0.5,
#              size = 0.5,
#              curvature = -0.3
#   ) +
#   annotation_scale(location = "bl", width_hint = 0.5) +
#   annotation_north_arrow(location = "bl", which_north = "true", 
#                          pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
#                          style = north_arrow_fancy_orienteering)

max(coordCol_Pic$Long)
min(coordCol_Pic$Long)
max(coordCol_Pic$Lat)
min(coordCol_Pic$Lat)

(map1 = ggplot(data = world) +
  geom_sf() +
  geom_point(data = coordCol_Pic, aes(x = Long, y = Lat), size = 4, 
             shape = 23, fill = "darkred") +
  geom_text_repel(data = coordCol_Pic, aes(label = Colony, x = Long, y = Lat),
                  size = 5.5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_sf(xlim = c(2.7, 3.6), ylim = c(49.1, 49.7), expand = FALSE) +
  geom_curve(data = dyadsObsSelect_Pic, aes(x = xstart, y = ystart, xend = xend, yend = yend),
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black",
             alpha = 0.5,
             size = 0.5,
             curvature = -0.3
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering))



# THURINGIA MAP

coordCol_Thu = read.table("Data/Thu/coordThu.txt",h=T)
allel=read.table("Data/Thu/uniqueGenotypesWithInfo.txt",h=T)
allel$idcol = as.factor(allel$idcol)
colonies = as.character(allel$idcol)
coordCol_Thu = coordCol_Thu[coordCol_Thu$Colony %in% colonies,]
# coordCol_Thu = coordCol_Thu[coordCol_Thu$Colony != "Thu42",]

plot(coordCol_Thu$Long, coordCol_Thu$Lat)
text(coordCol_Thu$Long, coordCol_Thu$Lat, coordCol_Thu$Colony)

dyadsObsSelect_Thu = dyadsObsSelect_Thu[-which(dyadsObsSelect_Thu$offspring == dyadsObsSelect_Thu$father),]
dyadsObsSelect_Thu$xstart = NA
dyadsObsSelect_Thu$ystart = NA
dyadsObsSelect_Thu$xend = NA
dyadsObsSelect_Thu$yend = NA
for (i in 1:nrow(dyadsObsSelect_Thu)) {
  dyadsObsSelect_Thu$xstart[i] = coordCol_Thu$Long[which(coordCol_Thu$Colony == dyadsObsSelect_Thu$offspring[i])]
  dyadsObsSelect_Thu$ystart[i] = coordCol_Thu$Lat[which(coordCol_Thu$Colony == dyadsObsSelect_Thu$offspring[i])]
  dyadsObsSelect_Thu$xend[i] = coordCol_Thu$Long[which(coordCol_Thu$Colony == dyadsObsSelect_Thu$father[i])]
  dyadsObsSelect_Thu$yend[i] = coordCol_Thu$Lat[which(coordCol_Thu$Colony == dyadsObsSelect_Thu$father[i])]
}


(map2 = ggplot(coordCol_Thu, aes(x = Long, y = Lat)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = Colony),
                  size = 3.5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_curve(data = dyadsObsSelect_Thu, aes(x = xstart, y = ystart, xend = xend, yend = yend),
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "darkgrey",
             alpha = 0.5,
             size = 0.5,
             curvature = -0.3
  )  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", colour = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"),
        legend.position = "none"))


max(coordCol_Thu$Long)
min(coordCol_Thu$Long)
max(coordCol_Thu$Lat)
min(coordCol_Thu$Lat)

(map2 = ggplot(data = world) +
    geom_sf() +
    geom_point(data = coordCol_Thu, aes(x = Long, y = Lat), size = 4, 
               shape = 23, fill = "darkred") +
    geom_text_repel(data = coordCol_Thu, aes(label = Colony, x = Long, y = Lat),
                    size = 5.5,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines")) +
    xlab("Longitude") +
    ylab("Latitude") +
    coord_sf(xlim = c(10.1, 11.9), ylim = c(50.45, 51.6), expand = FALSE) +
    geom_curve(data = dyadsObsSelect_Thu, aes(x = xstart, y = ystart, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "black",
               alpha = 0.5,
               size = 0.5,
               curvature = -0.3
    ) +
    annotation_scale(location = "bl", width_hint = 0.5) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering))


ggarrange(map1, map2, ncol = 2, labels = c("A", "B"))

ggsave("Figures_R1/FigS1.jpeg",width = 18, height = 8)




#==========================================================#
# Summary statistics ----
#==========================================================#


allel=read.table("Data/Pic/uniqueGenotypesWithInfo.txt",h=T)
allel$idcol = as.factor(allel$idcol)

# Remove some colonies (too far away)
allel=allel[-which(allel$idcol=="M399" | allel$idcol=="M1079" | allel$idcol=="CXSGT"| allel$idcol=="M1975"  | allel$idcol=="M1979"),]
# allel=allel[allel$idcol != "M399",]

# Keep females only
# allel=allel[-which(allel$sexe=="M"),]
# Remove juveniles
# allel=allel[-which(allel$ageWhenFirstCaptur=="Juv"),]

# Remove individuals with uncomplete genotypes
# allel=allel[!apply(allel[,2:17],1,function (x) 0 %in% x),]

colonies = as.character(allel$idcol)

genotypes = data.frame(rha101 = paste(allel[,2], allel[,3], sep = "/"),
                       rha109 = paste(allel[,4], allel[,5], sep = "/"),
                       rha4 = paste(allel[,6], allel[,7], sep = "/"),
                       rha7 = paste(allel[,8], allel[,9], sep = "/"),
                       rhc108 = paste(allel[,10], allel[,11], sep = "/"),
                       rhc3 = paste(allel[,12], allel[,13], sep = "/"),
                       rhd102 = paste(allel[,14], allel[,15], sep = "/"),
                       rhd103 = paste(allel[,16], allel[,17], sep = "/"))


# matAlleles = df2genind(genotypes,
#                        sep="/",
#                        ploidy=2,
#                        type="codom",
#                        pop=colonies)

matAlleles = df2genind(genotypes,
                       sep="/",
                       ploidy=2,
                       NA.char = "0",
                       type="codom",
                       pop=colonies)
matAlleles
ploidy(matAlleles)
matAlleles@tab
  
basic_stats = hierfstat::basic.stats(matAlleles)
basic_stats$n.ind.samp

allrich = allelic.richness(matAlleles)

allrich = colMeans(as.data.frame(allrich$Ar))

df_Pic = data.frame(dataset = "Picardy",
                    sampling_site = colnames(basic_stats$n.ind.samp),
                mean_Ho = colMeans(basic_stats$Ho),
                mean_Hs = colMeans(basic_stats$Hs),
                mean_Fis = colMeans(basic_stats$Fis),
                mean_allelic_richness = allrich)

summary(df_Pic$mean_Ho)
summary(df_Pic$mean_allelic_richness)
summary(df_Pic$mean_Fis)

# DAPC
# dapc.clust =find.clusters(matAlleles,
#                           max.n.clust = length(unique(colonies)))
# dapc1 = dapc(matAlleles, dapc.clust$grp) 
# 
# scatter(dapc1)


allel=read.table("Data/Thu/uniqueGenotypesWithInfo.txt",
                 h=T)
# allel$idcol = as.factor(allel$idcol)

# enlever les males (femelles sont liees a la colonie beaucoup plus que les males)
# allel=allel[-which(allel$sexe=="M"),]
# enlever les juveniles
# allel=allel[-which(allel$ageWhenFirstCaptur=="Juv"),]
# Supprimer les individus avec un genotype incomplet
# i.e. individus avec un 0 parmi tous les alleles
#Enlever les individus avec missing data
# allel=allel[!apply(allel[,2:17],1,function (x) 0 %in% x),]


table(allel$idcol)
length(unique(allel$idcol))

# Remove Thu24 because only two samples
allel = allel[!(allel$idcol %in% c("Thu24", "Thu42")),]

table(allel$idcol)
length(unique(allel$idcol))



# Construction de l'objet genind qui sera transmis pour la Fst
colonies = as.character(allel$idcol)

genotypes = data.frame(rha101 = paste(allel[,2], allel[,3], sep = "/"),
                       rha109 = paste(allel[,4], allel[,5], sep = "/"),
                       rha4 = paste(allel[,6], allel[,7], sep = "/"),
                       rha7 = paste(allel[,8], allel[,9], sep = "/"),
                       rhc108 = paste(allel[,10], allel[,11], sep = "/"),
                       rhc3 = paste(allel[,12], allel[,13], sep = "/"),
                       rhd102 = paste(allel[,14], allel[,15], sep = "/"),
                       rhd103 = paste(allel[,16], allel[,17], sep = "/"))

matAlleles = df2genind(genotypes,
                       sep="/",
                       ploidy=2,
                       NA.char = "0",
                       type="codom",
                       pop=colonies)
matAlleles
ploidy(matAlleles)
d
basic_stats = hierfstat::basic.stats(matAlleles)
basic_stats$n.ind.samp

allrich = allelic.richness(matAlleles)

allrich = colMeans(as.data.frame(allrich$Ar))

df_Thu = data.frame(dataset = "Thuringia",
                    sampling_site = colnames(basic_stats$n.ind.samp),
                    mean_Ho = colMeans(basic_stats$Ho),
                    mean_Hs = colMeans(basic_stats$Hs),
                    mean_Fis = colMeans(basic_stats$Fis),
                mean_allelic_richness = allrich)

summary(df_Thu$mean_Ho)
summary(df_Thu$mean_allelic_richness)
summary(df_Thu$mean_Fis)

# DAPC
# dapc.clust =find.clusters(matAlleles,
#                           max.n.clust = length(unique(colonies)))
# dapc1 = dapc(matAlleles, dapc.clust$grp) 
# 
# scatter(dapc1)



df = bind_rows(df_Pic,
               df_Thu)


df %>%
  group_by(dataset) %>%
  summarise(avg_Ho = mean(mean_Ho),
            avg_Fis = mean(mean_Fis))

(df = df %>%
  select(-mean_Fis))

Pic_sampling_names = read_tsv("Data/Pic/Table_S1_EP.tsv")

df = df %>%
  left_join(Pic_sampling_names[,c("sampling_site", "zarzoso_sampling_names")])


df

# Table S1 Sum stats ----
write_tsv(df, "Tables_R1/Table_S1.tsv")




#==========================================================#
# Maternity and mating bonds ----
#==========================================================#

#-----------------------------------------------------------------------------#
# Distribution of offspring-father distances Picardy
#-----------------------------------------------------------------------------#
# The combined 8/8 and 6/8 dataset of inferred paternities AFTER filtering
dyadsObsSelect = read.table("AssignationPic/outputs/dyadsObsSelect.txt", header=TRUE)

# How many paternitites inferred AFTER filtering?
nrow(dyadsObsSelect)

# How many fathers?
table(dyadsObsSelect$fatherID)
length(unique(dyadsObsSelect$fatherID))

fatherInferred = read.table("AssignationPic/outputs/Paternities.txt",
                            sep = " ", header = T)
motherInferred = read.table("AssignationPic/outputs/Maternities.txt",
                            sep = " ", header = T)

allel=read.table("Data/Pic/uniqueGenotypesWithInfo.txt",h=T)
allel$idcol = as.factor(allel$idcol)

dyadsObsSelect$offspringID = fatherInferred$offspring
dyadsObsSelect$mother = NA
dyadsObsSelect$motherID = NA

for (i in 1:nrow(dyadsObsSelect)) {
  mom = motherInferred$mother[which(motherInferred$offspring == dyadsObsSelect$offspringID[i])]
  
  if (length(mom) > 0) {
    dyadsObsSelect$motherID[i] = mom
    dyadsObsSelect$mother[i] = as.character(allel$idcol[which(allel$idind == dyadsObsSelect$motherID[i])])
  }
}


pic = dyadsObsSelect

table(pic$motherID)
sum(!is.na(pic$motherID))
# How many couples? Recurrent mating bonds?
table(pic$motherID, pic$fatherID)
sort(table(pic$motherID, pic$fatherID), decreasing = T)

# Table S2 Pic inferred parents ----
write_tsv(dyadsObsSelect, "Tables_R1/Table_S2.tsv")


#-----------------------------------------------------------------------------#
# Distribution of offspring-father distances Thuringia ----
#-----------------------------------------------------------------------------#
# The combined 8/8 and 6/8 dataset of inferred paternities AFTER filtering
dyadsObsSelect = read.table("AssignationThu/outputs/dyadsObsSelect.txt", header=TRUE)

# How many paternitites inferred AFTER filtering?
nrow(dyadsObsSelect)

# How many fathers?
table(dyadsObsSelect$fatherID)
length(unique(dyadsObsSelect$fatherID))


fatherInferred = read.table("AssignationThu/outputs/Paternities.txt",
                            sep = " ", header = T)
motherInferred = read.table("AssignationThu/outputs/Maternities.txt",
                            sep = " ", header = T)

allel=read.table("Data/Thu/uniqueGenotypesWithInfo.txt",h=T)
allel$idcol = as.factor(allel$idcol)

dyadsObsSelect$offspringID = fatherInferred$offspring
dyadsObsSelect$mother = NA
dyadsObsSelect$motherID = NA

for (i in 1:nrow(dyadsObsSelect)) {
  mom = motherInferred$mother[which(motherInferred$offspring == dyadsObsSelect$offspringID[i])]
  
  if (length(mom) > 0) {
    dyadsObsSelect$motherID[i] = mom
    dyadsObsSelect$mother[i] = as.character(allel$idcol[which(allel$idind == dyadsObsSelect$motherID[i])])
  }
}

thu = dyadsObsSelect

table(thu$motherID)
sum(!is.na(thu$motherID))

# How many couples? Recurrent mating bonds?
table(thu$motherID, thu$fatherID)


# Table S3 Thu inferred parents ----
write_tsv(dyadsObsSelect, "Tables_R1/Table_S3.tsv")






#==========================================================#
# Pairwise Fst ----
#==========================================================#

#==========================================================#
# Picardy ----
#==========================================================#

#----------------------------------------------------------#
# MATRIX OF GENETIC DISTANCES

allel=read.table("Data/Pic/uniqueGenotypesWithInfo.txt",h=T)
allel$idcol = as.factor(allel$idcol)

# Remove some colonies (too far away)
allel=allel[-which(allel$idcol=="M399" | allel$idcol=="M1079" | allel$idcol=="CXSGT"| allel$idcol=="M1975"  | allel$idcol=="M1979"),]
# allel=allel[allel$idcol != "M399",]

# Keep females only
allel=allel[-which(allel$sexe=="M"),]
# Remove juveniles
allel=allel[-which(allel$ageWhenFirstCaptur=="Juv"),]

# Remove individuals with uncomplete genotypes
# allel=allel[!apply(allel[,2:17],1,function (x) 0 %in% x),]

colonies = as.character(allel$idcol)

genotypes = data.frame(rha101 = paste(allel[,2], allel[,3], sep = "/"),
                       rha109 = paste(allel[,4], allel[,5], sep = "/"),
                       rha4 = paste(allel[,6], allel[,7], sep = "/"),
                       rha7 = paste(allel[,8], allel[,9], sep = "/"),
                       rhc108 = paste(allel[,10], allel[,11], sep = "/"),
                       rhc3 = paste(allel[,12], allel[,13], sep = "/"),
                       rhd102 = paste(allel[,14], allel[,15], sep = "/"),
                       rhd103 = paste(allel[,16], allel[,17], sep = "/"))



matAlleles = df2genind(genotypes,
                       sep="/",
                       ploidy=2,
                       NA.char = "0",
                       type="codom",
                       pop=colonies)
matAlleles
ploidy(matAlleles)

# Pairwise Fst between colonies
# matGen=adegenet::pairwise.fst(matAlleles,res.type="matrix")
# Since a change in adegenet, the function moved to hierfstat
# matGen = genet.dist(matAlleles, method = "WC84")
matGen = genet.dist(genind2hierfstat(matAlleles), method = "WC84")

# matGen = pairwise.neifst(matAlleles)
matGen


# install_github("cam315/otuSummary")
library(otuSummary)

# tri_matrix is the lower triangular matrix
df_long = matrixConvert(matGen, colname = c("row", "col", "value"))

# Heatmap 
ggplot(df_long, aes(row, col, fill = value)) + 
  geom_tile() +
  geom_text(aes(label = round(value, 3))) +
  scale_fill_viridis_c() +
  theme_bw()

# Table S3 ----
# dir.create("Tables_R1")
# write_tsv(as.data.frame(as.matrix(matGen)),
#           "Tables_R1/Table_S3.tsv")

# install.packages('pheatmap') # if not installed already
library(pheatmap)

p = pheatmap(as.matrix(matGen),
         display_numbers = T,
         number_format = "%.3f",
         show_rownames = T,
         show_colnames = T)

gtPic = p$gtable

ggsave("Figures_R1/FstPic.jpeg", plot = gtPic,
       width=24, height=18,
       create.dir = T)

#==========================================================#
# Thuringia ----
#==========================================================#
# IBD for the adult females of the 19 colonies in Thuringia.

#----------------------------------------------------------#
# MATRIX OF GENETIC DISTANCES

allel=read.table("Data/Thu/uniqueGenotypesWithInfo.txt",
                 h=T)
# allel$idcol = as.factor(allel$idcol)

# enlever les males (femelles sont liees a la colonie beaucoup plus que les males)
allel=allel[-which(allel$sexe=="M"),]
# enlever les juveniles
allel=allel[-which(allel$ageWhenFirstCaptur=="Juv"),]
# Supprimer les individus avec un genotype incomplet
# i.e. individus avec un 0 parmi tous les alleles
#Enlever les individus avec missing data
# allel=allel[!apply(allel[,2:17],1,function (x) 0 %in% x),]


table(allel$idcol)
length(unique(allel$idcol))

# Remove Thu24 because only two samples
allel = allel[allel$idcol != "Thu24",]

table(allel$idcol)
length(unique(allel$idcol))



# Construction de l'objet genind qui sera transmis pour la Fst
colonies = as.character(allel$idcol)

genotypes = data.frame(rha101 = paste(allel[,2], allel[,3], sep = "/"),
                       rha109 = paste(allel[,4], allel[,5], sep = "/"),
                       rha4 = paste(allel[,6], allel[,7], sep = "/"),
                       rha7 = paste(allel[,8], allel[,9], sep = "/"),
                       rhc108 = paste(allel[,10], allel[,11], sep = "/"),
                       rhc3 = paste(allel[,12], allel[,13], sep = "/"),
                       rhd102 = paste(allel[,14], allel[,15], sep = "/"),
                       rhd103 = paste(allel[,16], allel[,17], sep = "/"))

matAlleles = df2genind(genotypes,
                       sep="/",
                       ploidy=2,
                       NA.char = "0",
                       type="codom",
                       pop=colonies)
matAlleles
table(colonies)
ploidy(matAlleles)

matGen = genet.dist(genind2hierfstat(matAlleles), method = "WC84")

# Distance genetique selon Rousset (1997)
matGen = matGen/(1-matGen)

matGen
max(matGen)
hist(matGen)


# Table S4 ----
# dir.create("Tables_R1")
# write_tsv(as.data.frame(as.matrix(matGen)),
#           "Tables_R1/Table_S4.tsv")


matGenDist = matGen / (1 - matGen)
pheatmap(as.matrix(matGenDist),
         display_numbers = T,
         number_format = "%.3f",
         show_rownames = T,
         show_colnames = T)


# Figure S7 ----
p = pheatmap(as.matrix(matGen),
             display_numbers = T,
             number_format = "%.3f",
             show_rownames = T,
             show_colnames = T)

gtThu = p$gtable

ggsave("Figures_R1/FstThu.jpeg",
       plot = gtThu,
       width=24, height=18,
       create.dir = T)


ggpubr::ggarrange(gtPic, gtThu, ncol = 2, labels = c("A", "B"))

ggsave("Figures_R1/FigS7.jpeg",
       width=16, height=8,
       create.dir = T)
