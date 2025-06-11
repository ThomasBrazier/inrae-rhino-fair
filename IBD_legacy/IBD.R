#===============================#
#       Petit Rhinolophe
#     Isolation-by-distance
#===============================#

# clear global environment: remove all variables
rm(list=ls(all=TRUE))

#----------------------
# Loading packages
# SYSTEM
library(rstudioapi)
library(devtools)
library(MASS)
library(fitdistrplus)
# Libraries 'pscl' (zero-inflated functions) and 'extraDistr' could benefit to the model comparison
library(gamlss.dist)
library(ggplot2)
library(hierfstat)
library(adegenet)
library(geosphere)
library(ade4)
library(ggplot2)
library(adegenet)


#-------------------------------
# MATRICE DES DISTANCES GENETIQUES

# Distance genetique a calculer Fst/(1-Fst)
# Calculer matrice de Fst pour chaque colonie

# Pierre-Loup : femelles adultes only


allel=read.table("uniqueGenotypesWithInfo.txt",h=T)
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
allel=allel[!apply(allel[,2:13],1,function (x) 0 %in% x),]
# Construction de l'objet genind qui sera transmis pour la Fst
matAlleles=df2genind(allel[,2:13],sep=NULL,ncode=3,ploidy=2,type="codom",pop=as.numeric(allel$idcol))
# Calcul de la Fst paire a paire entre colonies
matGen=pairwise.fst(matAlleles,res.type="matrix")
colnames(matGen)=levels(allel$idcol)[as.numeric(colnames(matGen))]
rownames(matGen)=levels(allel$idcol)[as.numeric(rownames(matGen))]
# Distance genetique selon Rousset (1997)
matGen=matGen/(1-matGen)

#-------------------------------
# MATRICE DES DISTANCES GEOGRAPHIQUES

# Distance geographique a recalculer
# Calcul d'une matrice de distances
coordCol=read.table("coordpic.txt",h=T)
coordCol=coordCol[-which(coordCol$Colonie=="M399" | coordCol$Colonie=="M1079" | coordCol$Colonie=="CXSGT"| coordCol$Colonie=="M1975"  | coordCol$Colonie=="M1979" | coordCol$Colonie=="P412"),]
#coordCol=coordCol[-which(coordCol$Colonie=="M399" | coordCol$Colonie=="P412"),]
lvl=colnames(matGen)
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

#-------------------------------
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
###########
matGen[matGen == 0] = NA
mantel.rtest(as.dist(matGeo,diag=FALSE,upper=FALSE),as.dist(matGen,diag=FALSE,upper=FALSE),nrepet=10000)

# library(ecodist)
# mantel(as.dist(matGeo,diag=FALSE,upper=FALSE)~as.dist(matGen,diag=FALSE,upper=FALSE),nperm=10000)

#------------------
# CORRELOGRAMME
# correlogramme sur les distances genetiques et geographiques (km, non log)
# fournit IC, error bars et correlation de Mantel par permutation
# sur distances reelles, correlation genetique en fonction de chaque couple de colonie (distance)
# source("fonctionMantel.R")
# corr=mantel.correlogram(matGen,exp(matGeo),nclass = 4,nperm=10000)
# corr
# # PLOT CORRELOGRAMME MAISON
# library(ggplot2)
# df=data.frame(dist=corr[,2],mantel=corr[,1],upIC=corr[,3],lowIC=corr[,4],p=corr[,5],upH0=corr[,6],lowH0=corr[,7])
# # export 1450x1300
# ggplot(data=df, aes(x=dist, y=mantel)) +
#   geom_line()+
#   geom_errorbar(aes(ymin=lowIC, ymax=upIC), width=0,size=2) +
#   geom_line(aes(x=dist,y=upH0),linetype="dashed",size=1) +
#   geom_line(aes(x=dist,y=lowH0),linetype="dashed",size=1) +
#   geom_point(size=4,stroke=4,shape=21,fill=c("black","black","white","white","black","black"))+
#   geom_hline(yintercept=0,linetype="dashed",size=1.5) +
#   xlab("Geographic distance in km\n") + ylab("Mantel correlation") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=30,hjust = 0.5),
#         plot.margin = margin(.5, 2, .5, .5, "cm"),
#         axis.title.x = element_text(color="black", size=40),
#         axis.title.y = element_text(color="black", size=40),
#         axis.text=element_text(size=40, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(5,"line"),
#         legend.text=element_text(size=30),
#         legend.title=element_text(size=30))

# avec VEGAN, par classe
# library(vegan)
# corr=mantel.correlog(matGen,exp(matGeo),nperm=10000,n.class=6,r.type="spearman",cutoff=FALSE)$mantel.res
# plot(mantel.correlog(matGen,exp(matGeo),nperm=10000,n.class=6,r.type="spearman",cutoff=FALSE),
#      main="Autocorrelation between genetic distances and geographic distances (km)",
#      xlab="Distances",ylab="Mantel correlation",alpha=0.05)
# # autocorr jusqu'a 21m
# 
# # PLOT CORRELOGRAMME
# library(ggplot2)
# df=data.frame(dist=corr[,1],mantel=corr[,3])
# ggplot(data=df, aes(x=dist, y=mantel)) +
#   geom_line()+
#   geom_point(size=4,stroke=4,shape=21,fill=c("black","black","white","black","white","black"))+
#   geom_hline(yintercept=0,linetype="dashed",size=1) +
#   xlab("Geographic distance in km\n") + ylab("Mantel correlation") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=30,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=40),
#         axis.title.y = element_text(color="black", size=40),
#         axis.text=element_text(size=40, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(5,"line"),
#         legend.text=element_text(size=30),
#         legend.title=element_text(size=30))

# PLOT IBD PATTERN
df=data.frame(Distance=as.vector(matGeo[lower.tri(matGeo,diag=TRUE)]),GeneDistance=as.vector(matGen[lower.tri(matGen,diag=TRUE)]))
#df=data.frame(Distance=as.vector(Dgeo),GeneDistance=as.vector(Dgen))
plot(df)
mod=lm(GeneDistance~Distance,data=df)
abline(lm(GeneDistance~Distance,data=df))
mod$coefficients[1]
mod$coefficients[2]

save(df,file="IBDPic.Rdata")
save(mod,file="modIBDPic.Rdata")

#-------------------------------
# GRAPH GGPLOT OF IBD PATTERN (13 colonies)

# export 1800x1300
p=ggplot(data=df, aes(x=Distance, y=GeneDistance)) +
  geom_point(size=2)+
  geom_abline(intercept=mod$coefficients[1],slope=mod$coefficients[2],size=2) +
  xlab("Geographic distance (km)") + ylab("Genetic distance Fst/(1-Fst)") +
  scale_x_continuous(breaks=log(c(2,5,10,20)), labels=c(2,5,10,20)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=22, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=10,hjust = 0.5),
        axis.title.x = element_text(color="black", size=30),
        axis.title.y = element_text(color="black", size=30),
        axis.text=element_text(size=30, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=30),
        legend.title=element_text(size=30))
p

ggsave("IBD 13 colonies.ps",p,device="ps",path=graphdir,width=15,height=8,units="cm",dpi=150)
ggsave("IBD 13 colonies.png",p,device="png",path=graphdir,width=15,height=8,units="cm",dpi=150)
if (Sys.info()["sysname"]=="Windows"){
  ggsave("IBD 13 colonies.wmf",p,device="wmf",path=graphdir,width=15,height=8,units="cm",dpi=150)
}

# 2 graphs dans la meme image
# library(ggpubr)
# 
# ggarrange()

#-------------------------------#
# Mean dispersal distance d'apres IBD pattern
#-------------------------------

# Estimate of effective density
# 1000 sampled females ; 3382 sampled genotypes in 17 colonies, 2851 sampled gen in 13 colonies
# half of colonies were sampled (personal communication, Eric Petit)
# so 1000*2(half sampled)*2(for males)=4000 bats in the area (theoretical) or 5702/6764 (computed from data set)
# n effective = 10% * total = 400 or 570/675
# area of the sampling site = 500km2 (forêt) ou 2500km2 (conservatif, comprend tous les points de la zone)

d=c(400/500,400/2500) # Effective density (ind.km-2)
p=c(0.005826133,0.01095223) # slope of IBD regression (for the subset of 13 colonies or 17 colonies)

# Formula the mean dispersal estimated from slope of IBD regression, given by Rousset (1997)
mean.dispersal=sqrt(1/(4*d*pi*p))
mean.dispersal



#---------------#
# INBREEDING
#---------------


inbreeding(matAlleles,res.type="estimate")
