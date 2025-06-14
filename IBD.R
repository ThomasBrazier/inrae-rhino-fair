#==========================================================#
#      Lesser horseshoe
#     Isolation-by-distance
#==========================================================#


# clear global environment: remove all variables
rm(list=ls(all=TRUE))
# library(rstudioapi)
# Get the directory of the file & set working directory
# filedir=dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(filedir)

#----------------------------------------------------------#
# Loading packages
# Check if packages are installed, install if necessary
source("Sources/packages.R")


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

# Genetic distance following Rousset (1997)
matGen=matGen/(1-matGen)
matGen

#-------------------------------#
# MATRIX OF GEOGRAPHIC DISTANCES

# Distance geographique a recalculer
# Calcul d'une matrice de distances
coordCol=read.table("Data/Pic/coordPic.txt",h=T)
coordCol=coordCol[-which(coordCol$Colonie=="M399" | coordCol$Colonie=="M1079" | coordCol$Colonie=="CXSGT"| coordCol$Colonie=="M1975"  | coordCol$Colonie=="M1979" | coordCol$Colonie=="P412"),]
# coordCol=coordCol[-which(coordCol$Colonie=="M399" | coordCol$Colonie=="P412"),]
# coordCol=coordCol[coordCol$Colonie != "M399",]

lvl=unique(colonies)
# COlonies in the same order as in the genetic matrix
coord=data.frame(Colonie=rep(NA,length(lvl)),
                 Long=rep(NA,length(lvl)),
                 Lat=rep(NA,length(lvl)))
coord$Colonie=lvl
for (i in 1:nrow(coord)) {
  coord$Long[i]=coordCol$Long[which(coordCol$Colonie==coord$Colonie[i])]
  coord$Lat[i]=coordCol$Lat[which(coordCol$Colonie==coord$Colonie[i])]
}

# Calculating the matrix
matGeo=matrix(nrow=nrow(coord),ncol=nrow(coord))
colnames(matGeo)=coord$Colonie
rownames(matGeo)=coord$Colonie

# create distance matrix
matGeo=distm(cbind(coord$Long,coord$Lat), fun=distCosine)
# Matrix in km
matGeo=matGeo/1000
# log distance as in Rousset (1997)
matGeo=log(matGeo)
# Correcting negative values due to log
matGeo[which(matGeo<0)]=0


#-------------------------------#
# MANTEL TEST

# "For this purpose, we first quantified the linear relationship between genetic distance FST/(1-FST)
# and the logarithm of Euclidean geographic distance (Rousset 1997), i.e. isolation-by-distance (IBD)
# and tested its significance with a Mantel test using 10000 permutations" (Jan 2017)

matGen[matGen == 0] = NA
mantel.rtest(as.dist(matGeo,diag=FALSE,upper=FALSE),
             as.dist(matGen,diag=FALSE,upper=FALSE),
             nrepet=1000)

# PLOT IBD PATTERN
# df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen[lower.tri(matGen,diag=FALSE)]))
df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen))
#df=data.frame(Distance=as.vector(Dgeo),GeneDistance=as.vector(Dgen))
plot(df)
mod=lm(GeneDistance~Distance,data=df)
abline(lm(GeneDistance~Distance,data=df))
mod$coefficients[1]
mod$coefficients[2]

#-------------------------------#
# PLOT OF IBD PATTERN (13 colonies)

p1 = ggplot(data=df, aes(x=Distance, y=GeneDistance)) +
  geom_point(size=2)+
  geom_abline(intercept=mod$coefficients[1],slope=mod$coefficients[2],size=2) +
  xlab("Geographic distance (km)") + ylab("Genetic distance Fst/(1-Fst)") +
  scale_x_continuous(breaks=log(c(1, 2,5,10,25)), labels=c(0, 2,5,10,25), limits = log(c(1, 30))) +
  scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06), labels=c(0, 0.02, 0.04, 0.06), limits = c(0, 0.07)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"))
p1


#-------------------------------#
# Mean dispersal distance d'apres IBD pattern
#-------------------------------#

# Estimate of effective density
# 1000 sampled females ; 3382 sampled genotypes in 17 colonies, 2851 sampled gen in 13 colonies
# half of colonies were sampled (personal communication, Eric Petit)
# so 1000*2(half sampled)*2(for males)=4000 bats in the area (theoretical) or 5702/6764 (computed from data set)
# n effective = 10% * total = 400 or 570/675
# area of the sampling site = 500km2 (forêt) ou 2500km2 (conservatif, comprend tous les points de la zone)
# 
d=c(400/500,400/2500) # Effective density (ind.km-2)
# p=c(0.005826133,0.01095223) # slope of IBD regression (for the subset of 13 colonies or 17 colonies)
mod$coefficients

# d = 400/500
p = mod$coefficients[2]
# Formula the mean dispersal estimated from slope of IBD regression, given by Rousset (1997)
mean.dispersal = sqrt(1/(4*d*pi*p))
mean.dispersal


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

#-------------------------------#
# MATRICE DES DISTANCES GEOGRAPHIQUES

# Distance geographique a recalculer
# Calcul d'une matrice de distances
coordCol=read.table("Data/Thu/coordThu.txt",h=T)

plot(coordCol$Long, coordCol$Lat)
text(coordCol$Long, coordCol$Lat, coordCol$Colony)

coordCol = coordCol[coordCol$Colony %in% colonies,]

plot(coordCol$Long, coordCol$Lat)
text(coordCol$Long, coordCol$Lat, coordCol$Colony)

coordCol$Colonie

coordCol = coordCol[coordCol$Colony != "Thu42",]

lvl=unique(colonies)

# Creation du fichier contenant les coordonnees dans le bon ordre
# !!!! Meme ordre que dans la matrice des distances genetiques matGen
coord=data.frame(Colonie=rep(NA,length(lvl)),
                 Long=rep(NA,length(lvl)),
                 Lat=rep(NA,length(lvl)))
coord$Colonie=lvl

for (i in 1:nrow(coord)) {
  coord$Long[i]=coordCol$Long[which(coordCol$Colony==coord$Colonie[i])]
  coord$Lat[i]=coordCol$Lat[which(coordCol$Colony==coord$Colonie[i])]
}

matGeo=matrix(nrow=nrow(coord),ncol=nrow(coord))
colnames(matGeo)=coord$Colonie
rownames(matGeo)=coord$Colonie
# Calcul des distances paire a paire

# create distance matrix
matGeo=distm(cbind(coord$Long,coord$Lat), fun=distCosine)
# Matrix in km
matGeo=matGeo/1000
# Si on veut les distances en log selon Rousset (1997)
matGeo=log(matGeo)
# Correction des valeurs negatives apres le passage en log
matGeo[which(matGeo<0)]=0



#-------------------------------#
# MANTEL TEST

matGen[matGen == 0] = NA
mantel.rtest(as.dist(matGeo,diag=FALSE,upper=FALSE),
             as.dist(matGen,diag=FALSE,upper=FALSE),
             nrepet=1000)


# PLOT IBD PATTERN
# df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen[lower.tri(matGen,diag=FALSE)]))
df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen))

##### KEEPING ONLY DISTANCES INFERIOR TO 25KM (SHORT SCALE, same scale as IBD in Pic)
df = df[which(df$Distance < log(25)),]

#df=data.frame(Distance=as.vector(Dgeo),GeneDistance=as.vector(Dgen))
plot(df)
mod=lm(GeneDistance~Distance,data=df)
abline(lm(GeneDistance~Distance,data=df))
mod$coefficients[1]
mod$coefficients[2]

# save(df,file="Data/Pic/IBDPic.Rdata")
# save(mod,file="Data/Pic/modIBDPic.Rdata")

#-------------------------------#
# GGPLOT OF IBD PATTERN

p2 = ggplot(data=df, aes(x=Distance, y=GeneDistance)) +
  geom_point(size=2)+
  geom_abline(intercept=mod$coefficients[1],slope=mod$coefficients[2],size=2) +
  xlab("Geographic distance (km)") + ylab("") +
  scale_x_continuous(breaks=log(c(1, 2,5,10,25)), labels=c(0, 2,5,10,25), limits = log(c(1, 30))) +
  scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06), labels=c(0, 0.02, 0.04, 0.06), limits = c(0, 0.07)) +
  # ylim(0, 0.07) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"))
p2

#-------------------------------#
# Mean dispersal distance d'apres IBD pattern
#-------------------------------#

# Estimate of effective density
# 1000 sampled females ; 3382 sampled genotypes in 17 colonies, 2851 sampled gen in 13 colonies
# half of colonies were sampled (personal communication, Eric Petit)
# so 1000*2(half sampled)*2(for males)=4000 bats in the area (theoretical) or 5702/6764 (computed from data set)
# n effective = 10% * total = 400 or 570/675
# area of the sampling site = 500km2 (forêt) ou 2500km2 (conservatif, comprend tous les points de la zone)

mod$coefficients

d = 400/500
p = mod$coefficients[2]
# Formula the mean dispersal estimated from slope of IBD regression, given by Rousset (1997)
mean.dispersal=sqrt(1/(4*d*pi*p))
mean.dispersal




p = ggarrange(p1, p2, ncol = 2)
p

ggsave("Figures/FigS2.jpeg",
       dpi=320,units="cm",width=26,height=10,
       create.dir = T)
