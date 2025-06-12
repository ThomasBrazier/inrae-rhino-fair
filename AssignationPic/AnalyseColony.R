#===============================#
#            COLONY
#         Colony on Data
#     Analyse des resultats
#===============================#

# DIRECTORIES
# Working directories
# clear global environment: remove all variables
rm(list=ls(all=TRUE))
# library(rstudioapi)
# Get the directory of the file & set working directory
# filedir=dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(filedir)

#----------------------------------------------------------#
# Loading packages
# Check if packages are installed, install if necessary
source("../Sources/packages.R")
source("fonctionsColony.R")

info = read.table("../Data/Pic/uniqueGenotypesWithInfo.txt", h = T)

#-------------------------------#
#   CLEAN DIRECTORIES & REMOVE USELESS FILES
#-------------------------------
# res=c("Colony99loci/","Colony79lociSansM399/","colAllCohorts50/","simFathersIncluded10/",
#       "simFathersIncluded20/","simFathersIncluded30/","simIncompleteGen/",
#       "simRandomPop/","simRandomPopIncomplete/","simreplicatedRuns2400/","simReplicatedRuns48/")
# for (j in 1:length(res)) {
#   cleanUp(res[j])
# }

#-------------------------------#
#   CONSTRUCTION DU DATASET
#   A PARTIR DE N RUNS DE COLONY
#-------------------------------

# Get all the juveniles assigned to a parent (father/mother)
# Build the juvenile files
# Then get paternities
# and build the paternity file
# excluding false paternitities

# Recuperation de tous les juveniles assignes a un parent
# construction du fichier des juveniles
# Recuperation des paternites
# Construction du fichier des paternites
# Exclusion des fausses assignations

# Arguments
# path = path to COLONY assignation output
# excl = the proportion of runs required to be kept as a father
# assignments below this threshold are removed

# Make it on the dataset with individuals with at least 7/9 loci, without colony M399
# and the dataset with only individuals with 9/9 loci
paternity=constructResults(path="Colony79lociSansM399/",excl=0.3)
paternity2=constructResults(path="Colony99loci/",excl=0.3)
# OR
# paternity=constructResults(path="Colony79lociSansM399/",excl=0.003)
# paternity2=constructResults(path="Colony99loci/",excl=0.003)

###### SELECTION de PATERNITY
###### RASSEMBLER LES DEUX JEUX DE DONNEES paternity ET paternity2
###### Merging paternity and paternity2 for a consensus dataset

# 1/ Conflits entre les deux resultats ?
# Combien de juveniles identiques retrouves
# How many identical juveniles founds ?
identique=sum((paternity[,1] %in% paternity2[,1]))
identique
totalpat=(nrow(paternity) + nrow(paternity2))-identique
identique/totalpat

# Conflits sur les juvenile identiques ?
# Compare le vecteur des peres 79 au vecteur des peres 99
father79=paternity[(paternity[,1] %in% paternity2[,1]),]
father99=paternity2[(paternity2[,1] %in% paternity[,1]),]
# Nombre de conflits sur un juvenile
count=c()
for (i in 1:nrow(father79)) {
  if(father79[i,2]!=father99[which(father99[,1]==father79[i,1]),2]){
    count=c(count,father79[i,1]) # compte des juveniles avec conflit
  }
}
length(count)

sum(sort(father79)!=sort(father99))

# 2/ Combiner les deux jeux de donnees
# 2/ Combine the two datasets
# Argument : "99" prend le jeu de donnees 9/9 loci complets
#           "79" prend le jeu de donnees 7/9 loci complets
#           "all" prend les deux jeux de donnees
Paternities=combineParentageMethods("all",sex="M")
write.table(Paternities,"Paternities.txt",col.name=TRUE,row.names=TRUE)

# MATERNITES
maternity=constructResultsMums(path="Colony79lociSansM399/",excl=0.3)
maternity2=constructResultsMums(path="Colony99loci/",excl=0.3)
descResults=describeResults(maternity)

Maternities=combineParentageMethods("all",sex='F')
write.table(Maternities,"Maternities.txt",col.name=TRUE,row.names=TRUE)

#-----------------------------#
# SIM RESULTS WITH INFOS - EXPORT DATA
#-----------------------------
# Save a file containing :
# - offspring ID
# - fathers and mothers ID
# - sex
# - correct/incorrect assignation
# - genotype

empiricalResultsWithInfo()

#-------------------------------#
#   DESCRIBE RESULTS
#   DESCRIPTION DES RESULTATS
#   A PARTIR DE N RUNS DE COLONY
#-------------------------------#
descResults=describeResults(Paternities)
write.table(descResults,paste(Rdir,"Assignation/resultsWithInfo.txt",sep=""),col.name=TRUE,row.names=FALSE)

length(unique(descResults$offspring)) # nombre de juveniles
length(unique(descResults$father)) # nombre de peres

# # Distribution du nombre d'occurence de chaque couple offsprings-fathers
# distBestF=distBestFather(path="Colony79lociSansM399/")
# # OR
# distBestF=distBestFather(path="Colony99loci/")
# 
# # PLOT : Distribution du nombre d'occurence de chaque couple offsprings-fathers
# # en fonction de la distribution attendue d'apres les simulations
# distSim=read.table(paste(Rdir,"Assignation/04_data/outputs/occurencesBestFather.txt",sep=""),header=TRUE)
# # OR
# distSim=read.table(paste(Rdir,"Assignation/04_data/outputs/occurencesBestFatherIncompleteGen.txt",sep=""),header=TRUE)
# 
# 
# # PLOT
# df=data.frame(occurence=distBestF$occurences/nrun)
# dfSim=data.frame(occurence=distSim$occurence[which(distSim$runs==29)])
# 
# hist(df$occurence,freq=FALSE,ylim=c(0,8))
# lines(density(as.numeric(unlist(df)),na.rm=TRUE))
# lines(density(as.numeric(unlist(dfSim))))
# 
# nrun=30
# 
# ggplot(data=df, aes(x=occurence)) +
#   geom_histogram(aes(y=..density..),binwidth=0.1,fill="grey",color="black") +
#   geom_line(aes(y=..density..),fill="grey",color="black",stat="density",size=1.3) +
#   geom_line(data=dfSim,aes(x=occurence,y=..density..),fill="grey",color="black",stat="density",size=2) +
#   
#   # geom_bar(aes(y=..count../457),fill="grey",color="black") +
#   # geom_line(data=dfSim,aes(x=occurence,y=..density../200),size=1.3,stat='density') +
#   xlim(0,1) +
#   scale_x_continuous(breaks=c(0.00,0.25,0.35,0.50,0.75,1.00)) +
#   geom_vline(xintercept=0.35,linetype="dashed",size=1) +
#   ggtitle(paste("Distribution of assignment frequencies\n(n=",nrow(df),")"))+
#   xlab("Assignment frequency") + ylab("Density") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=24),
#         axis.title.y = element_text(color="black", size=24),
#         axis.text=element_text(size=24, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(5,"line"),
#         legend.text=element_text(size=30),
#         legend.title=element_text(size=30))



#-------------------------------#
#   COMPUTE THE H0 NULL DISTRIBUTION
#   CONSTRUCTION DE LA DISTRIBUTION H0
#-------------------------------

# abscisse : distance offspring-father
# ordonnee : frequence cumulee

# DISTRIBUTION H0 AVEC INTERVALLE DE CONFIANCE
# 1/ 1000 tirages de distribution H0
# 2/ Frequence cumulee pour chaque classe de distance (classes de 1km)
# 3/ Pour chaque tirage et chaque classe de distance : calculer la freq cumulee (i.e. combien de valeurs <= distance / total)
# 4/ Prendre les quantiles a 2.5% et 97.5% parmi les 1000 freq cumulees echantillonnees pour chaque classe de distance

# DISTRIBUTION H0 WITH CONFIDENCE INTERVAL
# 1/ 1000 resampling
# 2/ Cumulative frequencies for each class of distance
# 3/ For each draw compute the cumulative frequency (i.e. how many values <= distance / total)
# 4/ Take the 2.5% et 97.5% quantiles
Maternities=read.table("Maternities.txt",header=TRUE)
Paternities=read.table("Paternities.txt",header=TRUE)

nit=1000 # nombre d'iterations
maxDist=60 # distance maximale de la distribution
distClasses=1 # intervalle des classes de distance
interval=seq(0,maxDist,distClasses)
# Resultats stockes dans un data frame
bootH0=as.data.frame(matrix(NA,ncol=nit,nrow=length(interval)))

# Dresser une liste de couples offspring-father
# Dresser une liste des vrais offsprings (tous les juveniles assignes a une mere ou un pere)

# 2 OPTIONS
# 1/ tous les vrais juveniles, i.e. juveniles pour lesquels il y a un pere ou mere retrouve
# offsp = c(maternity[, 1], Paternities[, 1]) # "79" loci
# offsp = c(maternity2[, 1], Paternities[, 1]) # "99" loci
offsp = c(Maternities[, 1], Paternities[, 1]) # "all" loci
offsp = as.data.frame(unique(offsp))

# 2/ tous les juveniles candidats
#offsp=read.table("Offsprings.txt",h=F)

# echanger les individus contre leurs colonies
info = read.table("uniqueGenotypesWithInfo.txt", h = T)
for (i in 1:nrow(offsp)) {
  offsp[i, 1] = as.character(info$idcol[which(info$idind == offsp[i, 1])])
}
# Dresser une liste des peres (tous les peres candidats)
fathers = read.table("Rhino_CMS.txt", h = F)
fathers = fathers$V1
for (i in 1:length(fathers)) {
  fathers[i] = as.character(info$idcol[which(info$idind == fathers[i])])
}
# enlever les individus (pere ou offspring) venant d'une des 5 colonies exclues
# dyads impossibles ,"M1979","M1975","CXSGT","M1079"
excludedCol = c("M399")
offsp = offsp[which(!(offsp[, 1] %in% excludedCol) &
                      !(offsp[, 1] %in% excludedCol)), 1]
fathers = fathers[which(!(fathers %in% excludedCol) &
                          !(fathers %in% excludedCol))]
dyadsH0 = data.frame(offspring = offsp)
colnames(dyadsH0) = "offspring"

pb = txtProgressBar(min=0, max=nit, style=3) 
for (it in 1:nit) {
  # calculer H0
  # Tirer pour chaque juvenile un pere parmi tous les peres, toutes colonies confondues
  dyadsH0$father = sample(fathers, size = nrow(dyadsH0), replace = TRUE)
  
  # pour chaque couple calculer la distance en km entre colonies
  # formule : ACOS(SIN(lat1)*SIN(lat2)+COS(lat1)*COS(lat2)*COS(lon2-lon1))*6371
  coordCol = read.table("coordpic.txt", h = T)
  # coordCol$Lat=coordCol$Lat*pi/180
  # coordCol$Long=coordCol$Long*pi/180
  library(geosphere)
  for (i in 1:nrow(dyadsH0)) {
    # dyadsH0$distance[i]=acos(sin(coordCol$Lat[which(coordCol$Colonie==as.character(dyadsH0$offspring[i]))])*sin(coordCol$Lat[which(coordCol$Colonie==as.character(dyadsH0$father[i]))])+cos(coordCol$Lat[which(coordCol$Colonie==as.character(dyadsH0$offspring[i]))])*cos(coordCol$Lat[which(coordCol$Colonie==as.character(dyadsH0$father[i]))])*cos(coordCol$Long[which(coordCol$Colonie==as.character(dyadsH0$father[i]))]-coordCol$Long[which(coordCol$Colonie==as.character(dyadsH0$offspring[i]))]))*6371
    # create distance matrix
    dyadsH0$distance[i] = distCosine(cbind(coordCol$Long[which(coordCol$Colonie ==
                                                                 as.character(dyadsH0$offspring[i]))], coordCol$Lat[which(coordCol$Colonie ==
                                                                                                                            as.character(dyadsH0$offspring[i]))]),
                                     cbind(coordCol$Long[which(coordCol$Colonie == as.character(dyadsH0$father[i]))], coordCol$Lat[which(coordCol$Colonie ==
                                                                                                                                           as.character(dyadsH0$father[i]))]),
                                     r = 6378137) / 1000
  }
  
  vecF=c() # vecteur des freq cumulees pour chaque classe de distance
  for (d in 1:length(interval)){
    vecF[d]=sum(dyadsH0$distance<=interval[d])/nrow(dyadsH0)
  }
  bootH0[,it]=vecF
  setTxtProgressBar(pb,it)
}
close(pb)
# Tableau de resultats contenant les bornes inf-sup de l'IC a 95% de la distribution H0
icH0=as.data.frame(matrix(NA,ncol=3,nrow=length(interval)))
icH0[,1]=interval
for (i in 1:nrow(icH0)) {
  icH0[i,2]=mean(unlist(bootH0[i,]))
  icH0[i,3]=quantile(unlist(bootH0[i,]),0.025)
  icH0[i,4]=quantile(unlist(bootH0[i,]),0.975)
}
colnames(icH0)=c("distance","meanFreq","lowerIC","upperIC")
write.table(icH0,"icH0.txt",col.name=TRUE,row.names=TRUE)


# BOOTSTRAP MEAN DISPERSAL DISTANCE UNDER H0
boots=numeric(10000)
for (i in 1:10000){
  boots[i]=mean(sample(dyadsH0$distance,replace=T))
}
hist(boots)
limBoot1=quantile(boots,c(.025,.975))
limBoot1
ICDispersalH0=limBoot1
mboot=mean(boots)
mboot
meanDispersalH0=mboot
etboot=sqrt(var(boots))
etboot
mboot+c(-1,1)*qnorm(1-0.025)*etboot
etDispersalH0=etboot


#----
#-------------------------------
#   MAKE THE OBSERVED DISTRIBUTION
#-------------------------------#
#----
# abscisse : distance offspring-father
# ordonnee : frequence cumulee

# Dresser une liste de couples offspring-father retrouves par Colony (parente validee selon critere a definir)
# offspm=read.table("04_data/outputs/colRun1/RhinoThomas2.Paternity",h=T,fill=TRUE)
# dyadsObs=data.frame(offspring=offspm$OffspringID,father=offspm$InferredDad1,distance=rep(NA,nrow(offspm)))
dyadsObs=read.table("Paternity.txt",h=F)
colnames(dyadsObs)=c("offspring","father")
# echanger les individus contre leurs colonies
info=read.table("uniqueGenotypesWithInfo.txt",h=T)
dyadsObs$fatherID=dyadsObs$father
for(i in 1:nrow(dyadsObs)){
  dyadsObs$offspring[i]=as.character(info$idcol[which(info$idind==dyadsObs$offspring[i])])
  dyadsObs$father[i]=as.character(info$idcol[which(info$idind==dyadsObs$father[i])])
}
# enlever les individus (pere ou offspring) venant d'une des colonies exclues
# dyads impossibles ,"M1979","M1975","CXSGT","M1079"
excludedCol=c("M399")
dyadsObs=dyadsObs[which(!(dyadsObs$offspring %in% excludedCol)),]
dyadsObs=dyadsObs[which(!(dyadsObs$father %in% excludedCol)),]
# pour chaque couple calculer la distance en km entre colonies
# formule : ACOS(SIN(lat1)*SIN(lat2)+COS(lat1)*COS(lat2)*COS(lon2-lon1))*6371
coordCol=read.table("coordpic.txt",h=T)
# coordCol$Lat=coordCol$Lat*pi/180
# coordCol$Long=coordCol$Long*pi/180
library(geosphere)
for (i in 1:nrow(dyadsObs)) {
  # dyadsH0$distance[i]=acos(sin(coordCol$Lat[which(coordCol$Colonie==as.character(dyadsH0$offspring[i]))])*sin(coordCol$Lat[which(coordCol$Colonie==as.character(dyadsH0$father[i]))])+cos(coordCol$Lat[which(coordCol$Colonie==as.character(dyadsH0$offspring[i]))])*cos(coordCol$Lat[which(coordCol$Colonie==as.character(dyadsH0$father[i]))])*cos(coordCol$Long[which(coordCol$Colonie==as.character(dyadsH0$father[i]))]-coordCol$Long[which(coordCol$Colonie==as.character(dyadsH0$offspring[i]))]))*6371
  # create distance matrix
  dyadsObs$distance[i]=distCosine(cbind(coordCol$Long[which(coordCol$Colonie==as.character(dyadsObs$offspring[i]))],coordCol$Lat[which(coordCol$Colonie==as.character(dyadsObs$offspring[i]))]),
                                 cbind(coordCol$Long[which(coordCol$Colonie==as.character(dyadsObs$father[i]))],coordCol$Lat[which(coordCol$Colonie==as.character(dyadsObs$father[i]))]),
                                 r=6378137)/1000
}

#-------------------------------#
#   EXPORT DATA TO GRAPH
#-------------------------------
write.table(dyadsH0,"dyadsH0.txt",col.names=TRUE,row.names=FALSE)
# Run observed distribution with a threshold of 0.3
write.table(dyadsObs,"dyadsObsSelect.txt",col.names=TRUE,row.names=FALSE)
# And then rerun it with a threshold of 0
write.table(dyadsObs,"dyadsObs.txt",col.names=TRUE,row.names=FALSE)

#-------------------------------#
#   MEAN DISPERSAL DISTANCE
#     BOOTSTRAP
#-------------------------------
# Intervalle de confiance de la moyenne : calcul par bootstrap
#########################

#creation d'un vecteur "boots" comprenant 1000 elements, pour le moment tous des "0"
boots=numeric(10000)
#petite boucle qui va reechantillonner 1000 fois avec remise les donnees comprises dans
#le vecteur dyadsObs$distance, et calculer a chaque fois la moyenne des valeurs reechantillonnees
for (i in 1:10000){
  boots[i]=mean(sample(dyadsObs$distance,replace=T))
}
#trace l'histogramme de boot
hist(boots)
#calcule l'intervalle de confiance comme les 2,5% et 97,5% percentiles de la distribution
#de boots
limBoot1=quantile(boots,c(.025,.975))
limBoot1
ICDispersal=limBoot1
abline(v=limBoot1,col="red")
#calcule la moyenne, l'ecart-type et l'erreur standard apres bootstrap et les range dans
#des variables
mboot=mean(boots)
mboot
meanDispersal=mboot
etboot=sqrt(var(boots))
etboot
mboot+c(-1,1)*qnorm(1-0.025)*etboot
etDisperal=etboot


#-------------------------------#
#   GRAPH
#-------------------------------

####### Import data
# ALL NECESSARY DATA IS INCLUDED IN THIS FILES WRITEN EACH TIME YOU RUNNED THE PREVIOUS SCRIPT.
# YOU DON'T HAVE TO RERUN PREVIOUS SCRIPT FOR COMPUTING THE PLOT
dyadsH0=read.table("dyadsH0.txt",header=TRUE)
dyadsObsSelect=read.table("dyadsObsSelect.txt",header=TRUE)
dyadsObs=read.table("dyadsObs.txt",header=TRUE)
meanDispersal=11.32188
icH0=read.table("icH0.txt",header=TRUE)

##### Fitted distribution are estimated in the Dispersal kernel directory
# Cumulative density function from fitdistrplus and exponential distribution
z=seq(0,100,by=0.01)
lambda=0.10013465
lambdaInf=0.07614281
lambdaSup=0.13780406
yfit=function( z ){1-exp(-1*lambda*z)}
yfitInf=function( z ){1-exp(-1*lambdaInf*z)}
yfitSup=function( z ){1-exp(-1*lambdaSup*z)}
# Cumulative density function from fitdistrplus and gamma distribution
a=1
b=9.970022
yfitGamma=function( z ){pgamma(z,shape=a,scale=b)}
# Cumulative density function from MasterBayes estimates
lambda2=0.4704699
lambda2Inf=0.8413
lambda2Sup=0.2828
ybayes=function( z ){1-exp(-1*lambda2*z)}
ybayesInf=function( z ){1-exp(-1*lambda2Inf*z)}
ybayesSup=function( z ){1-exp(-1*lambda2Sup*z)}

#####################
# DYADS RELATIVE FREQUENCY IN FUNCTION OF DISTANCE
#####################
df=data.frame(d=c(sort(dyadsH0$distance),sort(dyadsObs$distance),sort(dyadsObsSelect$distance)),
              freq=c(seq(1,length(dyadsH0$distance),1)/length(dyadsH0$distance),seq(1,length(dyadsObs$distance),1)/length(dyadsObs$distance),seq(1,length(dyadsObsSelect$distance),1)/length(dyadsObsSelect$distance)),
              dist=c(rep("H0",length(dyadsH0$distance)),rep("Obs",length(dyadsObs$distance)),rep("ObsSelect",length(dyadsObsSelect$distance)))
              )
write.table(df,"distancesDist.txt",row.names=FALSE,col.names=TRUE)
df1=df[which(df$dist=="H0"),]
df2=df[which(df$dist=="Obs"),]
df3=df[which(df$dist=="ObsSelect"),]

ggplot(data=df1, aes(x=d, y=freq)) +
  stat_function(fun = yfit)+
  stat_function(fun = yfitInf,lty=2)+
  stat_function(fun = yfitSup,lty=2)+
  # stat_function(fun = yfitGamma,size=1.5,lty=3)+
  stat_function(fun = ybayes,colour="grey")+
  stat_function(fun = ybayesInf,colour="grey",lty=2)+
  stat_function(fun = ybayesSup,colour="grey",lty=2)+
  # geom_point(aes(color="expected"),size=4)+
  geom_ribbon(data=icH0,aes(x=distance,y=meanFreq,ymin=lowerIC, ymax=upperIC),alpha=0.5) +
  # geom_point(data=df2,aes(color="observed"),size=4) +
  geom_point(data=df3,aes(color="observed with selection")) +
  # geom_line(data=icH0,aes(x=distance,y=lowerIC),color="grey",size=1.5) +
  # geom_line(data=icH0,aes(x=distance,y=upperIC),color="grey",size=1.5) +
  scale_colour_manual(values=c("black"))+
  ylim(0,1) +
  # scale_x_continuous(breaks=c(0,round(meanDispersal,digits=1),20,40,60),limits=c(0,60)) +
  # ggtitle(paste("Distribution of offspring-father distances\n(n=",nrow(df3),")"))+
  xlab("Geographic distance (km)") + ylab("Dyads relative frequency") +
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
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')

# ggsave(paste(graphdir,"/Distribution of offspring-father distances PicAll.tiff",sep=""),
#        device="tiff",dpi=320,units="cm",width=19,height=13)

################
# Violin plots
# df=rbind(df1[,c(1,3)],df2[,c(1,3)])
# ggplot(data=df, aes(x=dist, y=d,color=dist)) +
#   geom_violin(size=2)+
#   geom_boxplot(width=0.1) +
#   ggtitle(paste("Distribution of offspring-father distances\n(n=",length(dyadsObs$distance),")"))+
#   xlab("Distribution") + ylab("Distance (km)") +
#   scale_y_continuous(breaks=c(0,10,20,30,40,50,60),limits=c(0,60)) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=30),
#         axis.title.y = element_text(color="black", size=30),
#         axis.text=element_text(size=26, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(5,"line"),
#         legend.text=element_text(size=30),
#         legend.title=element_text(size=30),
#         legend.position='none')

# Droite de Henry en quantiles
# plot(x=cumsum(hist(dyadsObs$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density),
#      y=cumsum(hist(dyadsH0$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density),
#      xlim=c(0,1),ylim=c(0,1),
#      main="Droite de Henry entre distributions observées et attendues\ndes distances offspring-father",
#      xlab="Distribution observ?e",ylab="Distribution H0")
# abline(0,1)
# 
# df=data.frame(H0=cumsum(hist(dyadsObs$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density),
#               Obs=cumsum(hist(dyadsH0$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density))
# 
# ggplot(data=df, aes(x=H0, y=Obs)) +
#   geom_point(size=2)+
#   geom_abline(intercept=0,slope=1) +
#   ylim(0,1) + xlim(0,1) +
#   ggtitle(paste("Q-Q plot of distances offspring-father\n (n=",length(dyadsObs$distance),")"))+
#   xlab("Expected distribution under random mating") + ylab("Observed distribution") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=30),
#         axis.title.y = element_text(color="black", size=30),
#         axis.text=element_text(size=26, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(5,"line"),
#         legend.text=element_text(size=30),
#         legend.title=element_text(size=30))


# Test de Kolmogorov
# Test sur les frequences cumulees de deux distributions
?ks.test
#ks.test(x=cumsum(hist(dyadsObs$distance,plot=F,breaks=seq(0,130,10))$counts/sum(hist(dyadsObs$distance,plot=F,breaks=seq(0,130,10))$counts)),y=cumsum(hist(dyadsH0$distance,plot=F,breaks=seq(0,130,10))$counts/sum(hist(dyadsH0$distance,plot=F,breaks=seq(0,130,10))$counts)))

# Test sur les quantiles
# ks.test(x=cumsum(hist(dyadsObs$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density),
#         y=cumsum(hist(dyadsH0$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density))
ks.test(x=dyadsH0$distance,
        y=dyadsObsSelect$distance)

# EXCLUSION DES COLONIES
# excludedCol=c("M399","M1979","M1975","CXSGT","M1079")
# coordCol=coordCol[which(!(coordCol$Colonie %in% excludedCol) & !(coordCol$Colonie %in% excludedCol)),]



##############################################################
# Below is experimental - not required
##############################################################

#-------------------------------#
#   % OF PARENTS IN SAME COLONY AS OFFSPRING
#-------------------------------

#---------------#
# EXPECTED UNDER H0
#---------------
# Dresser une liste de couples offspring-father
# Dresser une liste des vrais offsprings (tous les juveniles assignes a une mere ou un pere)
# 2 OPTIONS
# 1/ tous les vrais juveniles, i.e. juveniles pour lesquels il y a un pere ou mere retrouve
offsp = c(Maternities[, 1], Paternities[, 1])
offsp = as.vector(t(unique(offsp)))

# 2/ tous les juveniles candidats
#offsp=read.table("Offsprings.txt",h=F)

# Faible nombre d'échantillons, donc estimation par reechantillonnage de la proportion de paternite intracoloniale sous H0
nperm=10000
somme=0 # somme totale des pater intra col pour nperm
listP=c() # liste des pater intra col pour chaque iteration

# echanger les individus contre leurs colonies
info=read.table("uniqueGenotypesWithInfo.txt",h=T)
for(i in 1:length(offsp)){
  offsp[i]=as.character(info$idcol[which(info$idind==offsp[i])])
}
# Dresser une liste des peres (tous les peres candidats)
fathers=paternity[,2]
fathers=as.vector(t(fathers))
for(i in 1:length(fathers)){
  fathers[i]=as.character(info$idcol[which(info$idind==fathers[i])])
}
# enlever les individus (pere ou offspring) venant d'une colonie exclue
# M399
excludedCol=c("M399")
offsp=offsp[which(!(offsp %in% excludedCol) & !(offsp %in% excludedCol))]
fathers=fathers[which(!(fathers %in% excludedCol) & !(fathers %in% excludedCol))]

# create progress bar
pb <- txtProgressBar(min = 0, max = nperm, style = 3)
for (n in 1:nperm) {
  # Distribution attendue pour chaque colonie
  # Pas de reechantillonnage sur les juveniles (lien avec la colonie)

  # Tirer pour chaque pere un juvenile parmi tous les juveniles, toutes colonies confondues
  dyadsH0=data.frame(father=fathers)
  colnames(dyadsH0)=c("father")
  dyadsH0$offspring=sample(offsp,size=nrow(dyadsH0),replace=TRUE)
  # condition : $intracol=TRUE si pere et offspring sont de la meme colonie 
  dyadsH0$intracol=(dyadsH0$offspring==dyadsH0$father)
  # Proportion of intracolonial paternity
  listP=c(listP,sum(dyadsH0$intracol==TRUE))
  somme=somme+sum(dyadsH0$intracol==TRUE)
  # update progress bar
  setTxtProgressBar(pb,n)
}
close(pb)

hist(listP)
# proportion de paternite intra-coloniale
propIntra=(somme/nperm)/length(fathers)
propIntra*100

mean(listP)/length(fathers) # proportion moy par bootstrap
#trace l'histogramme de boot
hist(listP)
#calcule l'intervalle de confiance comme les 2,5% et 97,5% percentiles de la distribution
limBoot1=quantile(listP,c(.025,.975))
limBoot1/length(fathers)
abline(v=limBoot1,col="red")


# Compter la proportion de TRUE parmi tous les juveniles de chaque colonie
# table(dyadsH0$intracol,
#       dyadsH0$offspring)
# H0ColonyPater=table(dyadsH0$intracol,dyadsH0$offspring)[2,]/colSums(table(dyadsH0$intracol,dyadsH0$offspring))

#---------------#
# OBSERVED
#---------------
# Dresser une liste de couples offspring-father retrouves par Colony (parente validee selon critere a definir)
# offspm=read.table("04_data/outputs/colRun1/RhinoThomas2.Paternity",h=T,fill=TRUE)
# dyadsObs=data.frame(offspring=offspm$OffspringID,father=offspm$InferredDad1,distance=rep(NA,nrow(offspm)))
dyadsObs=read.table("Paternity.txt",h=F)
colnames(dyadsObs)=c("offspring","father")
# echanger les individus contre leurs colonies
info=read.table("uniqueGenotypesWithInfo.txt",h=T)
for(i in 1:nrow(dyadsObs)){
  dyadsObs$offspring[i]=as.character(info$idcol[which(info$idind==dyadsObs$offspring[i])])
  dyadsObs$father[i]=as.character(info$idcol[which(info$idind==dyadsObs$father[i])])
}
# enlever les individus (pere ou offspring) venant d'une des 5 colonies exclues
# dyads impossibles ,"M1979","M1975","CXSGT","M1079"
excludedCol=c("M399")
dyadsObs=dyadsObs[which(!(dyadsObs$offspring %in% excludedCol)),]
dyadsObs=dyadsObs[which(!(dyadsObs$father %in% excludedCol)),]
# condition : $intracol=TRUE si pere et offspring sont de la meme colonie 
dyadsObs$intracol=(dyadsObs$offspring==dyadsObs$father)
# number of parents of same colony
sum(dyadsObs$intracol==TRUE)
# number of parents of a different colony
sum(dyadsObs$intracol==FALSE)
# Proportion of intracolonial paternity
sum(dyadsObs$intracol==TRUE)/length(dyadsObs$intracol)
mean(dyadsObs$intracol)

# ESTIMATION DU POURCENTAGE PAR BOOTSTRAP
nperm=10000
boots=numeric(nperm)
pb <- txtProgressBar(min = 0, max = nperm, style = 3)
for (i in 1:nperm){
  boots[i]=mean(sample(dyadsObs$intracol,replace=T)) # moyenne des TRUE equivalente a proportion
  setTxtProgressBar(pb, i)
}
close(pb)
hist(boots)
limBoot1=quantile(boots,c(.025,.975))
limBoot1
abline(v=limBoot1,col="red")
mboot=mean(boots)
mboot
etboot=sqrt(var(boots))
etboot

#---------------#
# TEST BINOMIAL
#---------------
propIntra
binom.test(40,82,propIntra)


#---------------#
# GRAPH
#---------------

##### GGPLOT2
# TOUTES COLONIES ENSEMBLE
require(ggplot2)
df=data.frame(distr=c("H0","Obs"),intracol=c(propIntra,sum(dyadsObs$intracol==TRUE)/length(dyadsObs$intracol)))
ggplot(data=df, aes(x=distr, y=intracol)) +
  geom_bar(stat="identity", width=0.5)+
  ylim(0,1)+ ggtitle("Intra-colonial paternity percentage\nobserved and expected under random mating",
                     subtitle=paste("n=",nrow(dyadsObs),sep=""))+
  xlab("Distribution") + ylab("Intra-colonial paternity percentage") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black", size=30,hjust = 0.5),
        axis.title.x = element_text(color="black", size=30),
        axis.title.y = element_text(color="black", size=30),
        axis.text=element_text(size=26, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=30),
        legend.title=element_text(size=30))


#-------------------------------#
#   ASSEMBLAGE DES RESULTATS
#   AVEC PERES INFERES ET GENOTYPES RECONSTRUITS PAR COLONY
#-------------------------------
# COLONY reconstruit le genotype du pere le plus probable quand il ne trouve pas de pere probable dans les peres candidats
# Cette information permet de faire certaines descriptions de la population
# Notamment le calcul du succes reproductif et de la variance du succes reproductif
# assignationsAll=AssignationWithInfos(dir="Colony79lociSansM399",runs=30,p=0.03)
# assignations=AssignationWithInfos(dir="Colony79lociSansM399",runs=30,p=0.3)
# 
# assignations=read.table("Assignations.txt",header=TRUE)


################# TRAVAIL A PARTIR DE PATERNITIES ET MATERNITIES
# DATA FRAME GENERES A PARTIR DES RESULTATS D'ASSIGNATION

assignations=rbind(Paternities,Maternities)
year=c()
# Annee de premiere capture ou naissance
for (i in 1:nrow(assignations)) {
  year[i]=info$yearWhenFirstCaptur[which(info$idind==assignations[i,1])]
}
# Classe d'age premiere capture
age=c()
for (i in 1:nrow(assignations)) {
  age[i]=as.character(info$ageWhenFirstCaptur[which(info$idind==assignations[i,2])])
}
# Classe d'age premiere capture
captur=c()
for (i in 1:nrow(assignations)) {
  captur[i]=as.character(info$yearWhenFirstCaptur[which(info$idind==assignations[i,2])])
}
assignations=cbind(assignations,c(rep("M",nrow(Paternities)),rep("F",nrow(Maternities))),year,age,captur)

colnames(assignations)=c("OffspringID","ParentID","sex","yearOffsp","ClassAgeParent","YearParent")
write.table(assignations,"Assignations.txt",row.names=FALSE,col.names=TRUE)

#-------------------------------#
#   NOMBRE D'OFFSPRING : MALES + FEMELLES
#-------------------------------
Rsuccess=Rsuccess()
write.table(Rsuccess,"Rsuccess.txt",row.names=FALSE,col.names=TRUE)

#-------------------------------#
#   MALE REPRODUCTIVE BIAS
#   TEST DU NOMBRE D'OFFSPRING ENTRE MALES/FEMELLES
#-------------------------------

# Comparaison de la moyenne du succes reproductif
# entre sampled breeding males (2013-2016) et sampled breeding females (2013-2016)
# Nombre de juveniles/an par pere s'etant reproduit
# un vecteur contenant les peres et le nombre de juveniles sur chaque annee (2013-2016)
ReprodSuccess=data.frame(FatherID=unique(assignations[which(assignations[,3]=="M"),2]),
                         NbOffspring2013=rep(NA,length(unique(assignations[which(assignations[,3]=="M"),2]))),
                         NbOffspring2014=rep(NA,length(unique(assignations[which(assignations[,3]=="M"),2]))),
                         NbOffspring2015=rep(NA,length(unique(assignations[which(assignations[,3]=="M"),2]))),
                         NbOffspring2016=rep(NA,length(unique(assignations[which(assignations[,3]=="M"),2]))))
for (i in 1:nrow(ReprodSuccess)) {
  ReprodSuccess[i,2]=length(which(assignations[,2]==ReprodSuccess[i,1] & assignations[,4]=="2013"))
  ReprodSuccess[i,3]=length(which(assignations[,2]==ReprodSuccess[i,1] & assignations[,4]=="2014"))
  ReprodSuccess[i,4]=length(which(assignations[,2]==ReprodSuccess[i,1] & assignations[,4]=="2015"))
  ReprodSuccess[i,5]=length(which(assignations[,2]==ReprodSuccess[i,1] & assignations[,4]=="2016"))
}

#-------------------------------
# FEMALES
ReprodSuccessF=data.frame(FatherID=unique(assignations[which(assignations[,3]=="F"),2]),
                          NbOffspring2013=rep(NA,length(unique(assignations[which(assignations[,3]=="F"),2]))),
                          NbOffspring2014=rep(NA,length(unique(assignations[which(assignations[,3]=="F"),2]))),
                          NbOffspring2015=rep(NA,length(unique(assignations[which(assignations[,3]=="F"),2]))),
                          NbOffspring2016=rep(NA,length(unique(assignations[which(assignations[,3]=="F"),2]))))
for (i in 1:nrow(ReprodSuccessF)) {
  ReprodSuccessF[i,2]=length(which(assignations[,2]==ReprodSuccessF[i,1] & assignations[,4]=="2013"))
  ReprodSuccessF[i,3]=length(which(assignations[,2]==ReprodSuccessF[i,1] & assignations[,4]=="2014"))
  ReprodSuccessF[i,4]=length(which(assignations[,2]==ReprodSuccessF[i,1] & assignations[,4]=="2015"))
  ReprodSuccessF[i,5]=length(which(assignations[,2]==ReprodSuccessF[i,1] & assignations[,4]=="2016"))
}


# TEST PAR PERMUTATION
#Test par permutation
# Difference des moyennes
#####################

# Distribution empirique de ces differences sous l'hypothese nulle qu'il n'y a pas de
# difference entre les deux groupes
# Difference de taille entre les deux echantillons
# Sample de deux echantillons males/femelles de meme taille pour etre equilibre ?
perm=numeric(10000)
for(i in 1:9999){
  size=length(as.vector(t(ReprodSuccess[,2:5])))
  v=sample(c(as.vector(t(ReprodSuccess[,2:5])),sample(as.vector(t(ReprodSuccessF[,2:5])),size=size,replace=FALSE)),replace=FALSE)
  perm[i]=mean(v[1:size])-mean(v[(size+1):(2*size)])
}
perm[10000]=mean(as.vector(t(ReprodSuccess[,2:5])))-mean(as.vector(t(ReprodSuccessF[,2:5])))
tap=table(abs(perm)>=abs(perm[10000]))
p=tap[2]/10000
p

limsum=quantile(perm,c(.025,.975))
hist(perm, main = "Histogramme de la difference des moyennes") # distribution de notre statistique de test sous H0 : p-valeur exacte. la différence des moyennes sous l'hypothèse H0
abline(v=limsum[1]) # quantile a 2.5%
abline(v=limsum[2]) # quantile a 97.5%
abline(v=perm[10000],lty=2) # valeur observee


#-------------------------------#
#   EFFECT OF VARIANCE IN REPRODUCTIVE SUCCESS
#   ON THE REPRODUCTIVE VARIANCE EFFECTIVE SIZE
#   HILL (1972) and BOUTEILLER PERRIN (2000)
#-------------------------------



#-------------------------------#
#   MATE SELECTION
#   MOTHERS MATE WITH THE SAME FATHERS OR NOT ?
#-------------------------------

# Femelles s'accouplent-elles plus souvent avec le meme male qu'au hasard ?
# Tracer liens males/femelles
# un data.frame :
# - ID offspring
# - ID father
# - ID mother
# A partir de assignations
mateSelect=data.frame(assignations[which(assignations[,3]=="M"),1:2],motherID=rep(NA,length(which(assignations[,3]=="M"))))
colnames(mateSelect)=c("offspID","fatherID","motherID")
for (i in 1:nrow(mateSelect)) {
  if(length(as.character(assignations[which(assignations[,3]=="F" & mateSelect$offspID[i]==assignations[,1]),2]))>0){
    mateSelect$motherID[i]=as.character(assignations[which(assignations[,3]=="F" & mateSelect$offspID[i]==assignations[,1]),2])
  }
}
mateSelect=mateSelect[!grepl("e[0-9]",mateSelect[,2]),]
mateSelect=mateSelect[!is.na(mateSelect$motherID),]
##### Maximal/observed/expected number of combinations
# max = number of males free for mating and assigned to a dyad by COLONY
length(which(!is.na(mateSelect$motherID)))
# observed = number of unique dyads
length(unique(paste(mateSelect[which(!is.na(mateSelect$motherID)),2],mateSelect[which(!is.na(mateSelect$motherID)),3])))
# expected = resampling between fathers and mothers for selected dyads ?


##### Number of offspring for each dyad
DyadNbOfOffspring=data.frame(dyad=unique(paste(mateSelect[which(!is.na(mateSelect$motherID)),2],mateSelect[which(!is.na(mateSelect$motherID)),3])),offspNb=rep(NA,length(unique(paste(mateSelect[which(!is.na(mateSelect$motherID)),2],mateSelect[which(!is.na(mateSelect$motherID)),3])))))
for (j in 1:nrow(DyadNbOfOffspring)) {
  DyadNbOfOffspring$offspNb[j]=length(which(as.character(DyadNbOfOffspring$dyad[j])==paste(mateSelect[which(!is.na(mateSelect$motherID)),2],mateSelect[which(!is.na(mateSelect$motherID)),3])))
}
# Nombre de couples n'ayant donne qu'un offsp
sum(DyadNbOfOffspring$offspNb==1)

# for (i in 1:nrow(mateSelect)) {
#   print(info[which(info$idind==mateSelect$motherID[i]),c(1,19,20)])
# }
# for (i in 1:nrow(mateSelect)) {
#   mateSelect[i,4]=info[which(info$idind==mateSelect$offspID[i]),21]
# }
# for (i in 1:nrow(mateSelect)) {
#   mateSelect[i,5]=paste(assignations[which(assignations$offspID==mateSelect$offspID[i] &
#                                              (assignations$parentID==mateSelect$motherID[i] & assignations$sex=="F")),4],collapse=" ")
# }

#-------------------------------#
#   NUMBER OF PARTNERS
#-------------------------------
# Number of partners for males
malesPartners=data.frame(FatherID=unique(mateSelect$fatherID),P=rep(NA,length(unique(mateSelect$fatherID))))
for (i in 1:nrow(malesPartners)) {
  malesPartners$P[i]=length(unique(mateSelect$motherID[which(mateSelect$fatherID==malesPartners$FatherID[i])]))
}

# Number of partners for females
femalesPartners=data.frame(MotherID=unique(mateSelect$motherID),P=rep(NA,length(unique(mateSelect$motherID))))
for (i in 1:nrow(femalesPartners)) {
  femalesPartners$P[i]=length(unique(mateSelect$fatherID[which(mateSelect$motherID==femalesPartners$MotherID[i])]))
}
table(malesPartners$P)
table(femalesPartners$P)

#-------------------------------#
#   DO OFFSPRING BREED ?
#-------------------------------
breedingOffsp=assignations[which(assignations[,5]=="Juv"),]
# Nb de juveniles s'etant reproduits
nrow(breedingOffsp)
# 2 cas de figures :
# - le juvenile est en fait un adulte non echantillonne en premiere session
# - le juvenile est vraiment un juvenile et donc une parente a pu lui etre assignee
# Trouver un parent du parent pour le verifier
breedingOffsp=cbind(breedingOffsp,rep(NA,nrow(breedingOffsp)))
colnames(breedingOffsp)[7]="trueOffsp"
for (i in 1:nrow(breedingOffsp)) {
  if(breedingOffsp[i,2] %in% assignations[,1]){
    breedingOffsp[i,7]=assignations[which(assignations[,1]==breedingOffsp[i,2]),2]
  }
}

