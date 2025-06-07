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

# Supprimer les individus de la colonie M399 non prise en compte
allel=allel[-which(allel$idcol=="M399" | allel$idcol=="M1079" | allel$idcol=="CXSGT"| allel$idcol=="M1975"  | allel$idcol=="M1979"),]
#allel=allel[-which(allel$idcol=="M399"),]
# enlever les males (femelles sont liees a la colonie beaucoup plus que les males)
allel=allel[-which(allel$sexe=="M"),]
# enlever les juveniles
allel=allel[-which(allel$ageWhenFirstCaptur=="Juv"),]
# Supprimer les individus avec un genotype incomplet
# i.e. individus avec un 0 parmi tous les alleles
#Enlever les individus avec missing data
allel=allel[!apply(allel[,2:17],1,function (x) 0 %in% x),]
# Construction de l'objet genind qui sera transmis pour la Fst
colonies = as.character(allel$idcol)
matAlleles=df2genind(allel[,2:17],
                     sep=NULL,
                     ncode=3,
                     ploidy=2,
                     type="codom",
                     pop=colonies)

# Calcul de la Fst paire a paire entre colonies
# matGen=adegenet::pairwise.fst(matAlleles,res.type="matrix")
# Since a change in adegenet, the function moved to hierfstat
matGen = genet.dist(matAlleles, method = "Nei87")
matGen

# colnames(matGen)=levels(allel$idcol)[as.numeric(colnames(matGen))]
# rownames(matGen)=levels(allel$idcol)[as.numeric(rownames(matGen))]

# Distance genetique selon Rousset (1997)
matGen=matGen/(1-matGen)

#-------------------------------#
# MATRICE DES DISTANCES GEOGRAPHIQUES

# Distance geographique a recalculer
# Calcul d'une matrice de distances
coordCol=read.table("Data/Pic/coordPic.txt",h=T)
coordCol=coordCol[-which(coordCol$Colonie=="M399" | coordCol$Colonie=="M1079" | coordCol$Colonie=="CXSGT"| coordCol$Colonie=="M1975"  | coordCol$Colonie=="M1979" | coordCol$Colonie=="P412"),]
#coordCol=coordCol[-which(coordCol$Colonie=="M399" | coordCol$Colonie=="P412"),]
lvl=unique(colonies)
# Creation du fichier contenant les coordonnees dans le bon ordre
# !!!! Meme ordre que dans la matrice des distances genetiques matGen
coord=data.frame(Colonie=rep(NA,13),Long=rep(NA,13),Lat=rep(NA,13))
coord$Colonie=lvl
for (i in 1:nrow(coord)) {
  coord$Long[i]=coordCol$Long[which(coordCol$Colonie==coord$Colonie[i])]
  coord$Lat[i]=coordCol$Lat[which(coordCol$Colonie==coord$Colonie[i])]
}

# Calcul de la matrice de distance geographique
#   Passage des coordonnees en radians
#coord$Lat=coordCol$Lat*pi/180
#coord$Long=coordCol$Long*pi/180
# Matrice vide
matGeo=matrix(nrow=nrow(coord),ncol=nrow(coord)) # Matrix de 17*17 -> 17 colonies prises en compte
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


# library(geosphere)
# # create distance matrix
# matGeo=distm(cbind(coord$Long,coord$Lat), fun=distCosine)
# # Matrix in km
# matGeo=matGeo/1000

#-------------------------------#
# MANTEL TEST

# "For this purpose, we first quantified the linear relationship between genetic distance FST/(1-FST)
# and the logarithm of Euclidean geographic distance (Rousset 1997), i.e. isolation-by-distance (IBD)
# and tested its significance with a Mantel test using 10000 permutations" (Jan 2017)

# Pierre-Loup : femelles adultes only

# Test avec ade4
# Distances genetiques calculees par adegenet
# Distances EUCLIDIENNES !!!
# mGen=genind2genpop(matAlleles)
# Dgen=dist.genpop(mGen,method=2)
# Dgeo=dist(cbind(coord$Long,coord$Lat))
# Dgeo=log(Dgeo)
# ibd=mantel.randtest(Dgen,Dgeo)
# ibd

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

# save(df,file="Data/Pic/IBDPic.Rdata")
# save(mod,file="Data/Pic/modIBDPic.Rdata")

#-------------------------------#
# GRAPH GGPLOT OF IBD PATTERN (13 colonies)

# export 1800x1300
p1 = ggplot(data=df, aes(x=Distance, y=GeneDistance)) +
  geom_point(size=2)+
  geom_abline(intercept=mod$coefficients[1],slope=mod$coefficients[2],size=2) +
  xlab("Geographic distance (km)") + ylab("Genetic distance Fst/(1-Fst)") +
  scale_x_continuous(breaks=log(c(2,5,10,25)), labels=c(2,5,10,25)) +
  scale_y_continuous(breaks=c(0, 0.006, 0.012), labels=c(0, 0.006, 0.012)) +
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

# ggsave("IBD 13 colonies.ps",p,device="ps",path=graphdir,width=15,height=8,units="cm",dpi=150)
# ggsave("IBD 13 colonies.png",p,device="png",path=graphdir,width=15,height=8,units="cm",dpi=150)
# if (Sys.info()["sysname"]=="Windows"){
#   ggsave("IBD 13 colonies.wmf",p,device="wmf",path=graphdir,width=15,height=8,units="cm",dpi=150)
# }

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
# d=c(400/500,400/2500) # Effective density (ind.km-2)
# p=c(0.005826133,0.01095223) # slope of IBD regression (for the subset of 13 colonies or 17 colonies)

d = 400/500
p = mod$coefficients[2]
# Formula the mean dispersal estimated from slope of IBD regression, given by Rousset (1997)
mean.dispersal=sqrt(1/(4*d*pi*p))
mean.dispersal


#==========================================================#
# Thuringia ----
#==========================================================#

#----------------------------------------------------------#
# MATRIX OF GENETIC DISTANCES

allel=read.table("Data/Thu/uniqueGenotypesWithInfo.txt",h=T)
allel$idcol = as.factor(allel$idcol)


allel = allel[allel$idcol != "Thu42",]

# enlever les males (femelles sont liees a la colonie beaucoup plus que les males)
allel=allel[-which(allel$sexe=="M"),]
# enlever les juveniles
allel=allel[-which(allel$ageWhenFirstCaptur=="Juv"),]
# Supprimer les individus avec un genotype incomplet
# i.e. individus avec un 0 parmi tous les alleles
#Enlever les individus avec missing data
allel=allel[!apply(allel[,2:17],1,function (x) 0 %in% x),]
# Construction de l'objet genind qui sera transmis pour la Fst
colonies = as.character(allel$idcol)
matAlleles=df2genind(allel[,2:17],
                     sep=NULL,
                     ncode=3,
                     ploidy=2,
                     type="codom",
                     pop=colonies)

# Calcul de la Fst paire a paire entre colonies
# matGen=adegenet::pairwise.fst(matAlleles,res.type="matrix")
# Since a change in adegenet, the function moved to hierfstat
matGen = genet.dist(matAlleles, method = "Nei87")
matGen

# colnames(matGen)=levels(allel$idcol)[as.numeric(colnames(matGen))]
# rownames(matGen)=levels(allel$idcol)[as.numeric(rownames(matGen))]

# Distance genetique selon Rousset (1997)
matGen=matGen/(1-matGen)

#-------------------------------#
# MATRICE DES DISTANCES GEOGRAPHIQUES

# Distance geographique a recalculer
# Calcul d'une matrice de distances
coordCol=read.table("Data/Thu/coordThu.txt",h=T)

lvl=unique(colonies)
# Creation du fichier contenant les coordonnees dans le bon ordre
# !!!! Meme ordre que dans la matrice des distances genetiques matGen
coord=data.frame(Colonie=rep(NA,20),Long=rep(NA,20),Lat=rep(NA,20))
coord$Colonie=lvl
for (i in 1:nrow(coord)) {
  coord$Long[i]=coordCol$Long[which(coordCol$Colony==coord$Colonie[i])]
  coord$Lat[i]=coordCol$Lat[which(coordCol$Colony==coord$Colonie[i])]
}

# Calcul de la matrice de distance geographique
#   Passage des coordonnees en radians
#coord$Lat=coordCol$Lat*pi/180
#coord$Long=coordCol$Long*pi/180
# Matrice vide
matGeo=matrix(nrow=nrow(coord),ncol=nrow(coord)) # Matrix de 17*17 -> 17 colonies prises en compte
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


# library(geosphere)
# # create distance matrix
# matGeo=distm(cbind(coord$Long,coord$Lat), fun=distCosine)
# # Matrix in km
# matGeo=matGeo/1000

#-------------------------------#
# MANTEL TEST

# "For this purpose, we first quantified the linear relationship between genetic distance FST/(1-FST)
# and the logarithm of Euclidean geographic distance (Rousset 1997), i.e. isolation-by-distance (IBD)
# and tested its significance with a Mantel test using 10000 permutations" (Jan 2017)

# Pierre-Loup : femelles adultes only

# Test avec ade4
# Distances genetiques calculees par adegenet
# Distances EUCLIDIENNES !!!
# mGen=genind2genpop(matAlleles)
# Dgen=dist.genpop(mGen,method=2)
# Dgeo=dist(cbind(coord$Long,coord$Lat))
# Dgeo=log(Dgeo)
# ibd=mantel.randtest(Dgen,Dgeo)
# ibd

matGen[matGen == 0] = NA
mantel.rtest(as.dist(matGeo,diag=FALSE,upper=FALSE),
             as.dist(matGen,diag=FALSE,upper=FALSE),
             nrepet=1000)


# PLOT IBD PATTERN
# df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen[lower.tri(matGen,diag=FALSE)]))
df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=FALSE)]),GeneDistance=as.vector(matGen))

##### KEEPING ONLY DISTANCES INFERIOR TO 30KM (SHORT SCALE)
df = df[which(df$Distance < log(30)),]

#df=data.frame(Distance=as.vector(Dgeo),GeneDistance=as.vector(Dgen))
plot(df)
mod=lm(GeneDistance~Distance,data=df)
abline(lm(GeneDistance~Distance,data=df))
mod$coefficients[1]
mod$coefficients[2]

# save(df,file="Data/Pic/IBDPic.Rdata")
# save(mod,file="Data/Pic/modIBDPic.Rdata")

#-------------------------------#
# GRAPH GGPLOT OF IBD PATTERN

# export 1800x1300
p2 = ggplot(data=df, aes(x=Distance, y=GeneDistance)) +
  geom_point(size=2)+
  geom_abline(intercept=mod$coefficients[1],slope=mod$coefficients[2],size=2) +
  xlab("Geographic distance (km)") + ylab("") +
  scale_x_continuous(breaks=log(c(2,5,10,25)), labels=c(2,5,10,25)) +
  # scale_y_continuous(breaks=c(0, 0.02, 0.04, 0.06), labels=c(0, 0.02, 0.04, 0.06)) +
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

# ggsave("IBD 13 colonies.ps",p,device="ps",path=graphdir,width=15,height=8,units="cm",dpi=150)
# ggsave("IBD 13 colonies.png",p,device="png",path=graphdir,width=15,height=8,units="cm",dpi=150)
# if (Sys.info()["sysname"]=="Windows"){
#   ggsave("IBD 13 colonies.wmf",p,device="wmf",path=graphdir,width=15,height=8,units="cm",dpi=150)
# }

#-------------------------------#
# Mean dispersal distance d'apres IBD pattern
#-------------------------------#

# Estimate of effective density
# 1000 sampled females ; 3382 sampled genotypes in 17 colonies, 2851 sampled gen in 13 colonies
# half of colonies were sampled (personal communication, Eric Petit)
# so 1000*2(half sampled)*2(for males)=4000 bats in the area (theoretical) or 5702/6764 (computed from data set)
# n effective = 10% * total = 400 or 570/675
# area of the sampling site = 500km2 (forêt) ou 2500km2 (conservatif, comprend tous les points de la zone)

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
