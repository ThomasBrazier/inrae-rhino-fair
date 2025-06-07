#==========================================================#
#            COLONY
#     Gametic dispersal
#     Gametic dispersal is the combination  of natal + breeding dispersal
#==========================================================#

#==========================================================#
# LOADING ENVIRONMENT
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
# Import data
#==========================================================#
dataPic = read.table("Data/Pic/uniqueGenotypesWithInfo.txt", header = TRUE)
dataThu = read.table("Data/Thu/uniqueGenotypesWithInfo.txt", header = TRUE)
# Migrant fathers
MigrantFathersPic = read.table("Tables/MigrantFathersPic.txt", header = TRUE)
MigrantFathersThu = read.table("Tables/MigrantFathersThu.txt", header = TRUE)


#==========================================================#
# Estimates of gametic dispersal
#==========================================================#

# Take the colony of origin of fathers instead of the sampling colony (residence colony)
# And re-compute all distances
# Once done, re-fit models on the distribution of distances and estimate distances


# Re-estimate Figure 2
# Figures_Gametic
# Figure 2; no need to re-compute H0 with actual fathers origins (birth colony of all fathers) -> random resampling
# Besides, re-do Dispersal kernel estimates and mean effective dispersal distance 

# Table 2
# fitDistR_model_comparison_Gametic
# Dispersal_Kernel_Gametic
# retropropagate dispersal kernel estimates to Figure 2

# Figure 3
# Figures_Gametic
# Dispersal kernel with new estimations of parameters



#==========================================================#
# New estimates of dispersal distribution
#==========================================================#
#----------------------------------------------------------#
# Picardy
#----------------------------------------------------------#

#----------------------------------------------------------#
#   CONSTRUCTION DE LA DISTRIBUTION OBSERVEE
#----------------------------------------------------------#

# abscisse : distance offspring-father
# ordonnee : frequence cumulee

# Dresser une liste de couples offspring-father retrouves par Colony (parente validee selon critere a definir)
# offspm=read.table("04_data/outputs/colRun1/RhinoThomas2.Paternity",h=T,fill=TRUE)
# dyadsObs=data.frame(offspring=offspm$OffspringID,father=offspm$InferredDad1,distance=rep(NA,nrow(offspm)))
dyadsObs = read.table("AssignationPic/resultsWithInfo.txt", header = TRUE)
dyadsObs = dyadsObs[,-c(5:16)]
# Gametic dispersal, hence exchange the sampling colony by the birth colony
# echanger les individus contre leurs colonies
info = read.table("AssignationPic/uniqueGenotypesWithInfo.txt",h=T)
MigrantFathersPic = read.table("Tables/MigrantFathersPic.txt", header = TRUE)
dyadsObs$birthColonyFather = NA
for(i in 1:nrow(dyadsObs)){
  dyadsObs$birthColonyFather[i] = as.character(MigrantFathersPic$ancestry[MigrantFathersPic$idind == dyadsObs$father[i]])
}
# enlever les individus (pere ou offspring) venant d'une des colonies exclues
# dyads impossibles ,"M1979","M1975","CXSGT","M1079"
excludedCol=c("M399")
dyadsObs = dyadsObs[which(!(dyadsObs$colonyOffsp %in% excludedCol)),]
dyadsObs = dyadsObs[which(!(dyadsObs$colonyFather %in% excludedCol)),]
dyadsObs = dyadsObs[which(!(dyadsObs$birthColonyFather %in% excludedCol)),]
# pour chaque couple calculer la distance en km entre colonies
# formule : ACOS(SIN(lat1)*SIN(lat2)+COS(lat1)*COS(lat2)*COS(lon2-lon1))*6371
coordCol=read.table("AssignationPic/coordpic.txt",h=T)
# coordCol$Lat=coordCol$Lat*pi/180
# coordCol$Long=coordCol$Long*pi/180
# Distance between birth colony of the juvenile and inferred birth colony of the father by STRUCTURE
dyadsObs$distance = NA
for (i in 1:nrow(dyadsObs)) {
  # create distance matrix
  dyadsObs$distance[i] = distCosine(as.matrix(cbind(coordCol$Long[which(coordCol$Colonie==as.character(dyadsObs$colonyOffsp[i]))],
                                          coordCol$Lat[which(coordCol$Colonie==as.character(dyadsObs$colonyOffsp[i]))])),
                                  as.matrix(cbind(coordCol$Long[which(coordCol$Colonie==as.character(dyadsObs$birthColonyFather[i]))],
                                        coordCol$Lat[which(coordCol$Colonie==as.character(dyadsObs$birthColonyFather[i]))])),
                                  r=6378137)/1000
}
dyadsObs$distance


#----------------------------------------------------------#
#   EXPORT NEW DATA TO GRAPH
#----------------------------------------------------------#
# Run observed distribution with a threshold of 0.3
# write.table(dyadsObs,"AssignationPic/dyadsObsSelect_gametic.txt",col.names=TRUE,row.names=FALSE)
dyadsObs = read.table("AssignationPic/dyadsObsSelect_gametic.txt", header = TRUE)

#----------------------------------------------------------#
#   MEAN DISPERSAL DISTANCE
#----------------------------------------------------------#
# Intervalle de confiance de la moyenne : calcul par bootstrap
#########################

#creation d'un vecteur "boots" comprenant 1000 elements, pour le moment tous des "0"
boots = numeric(10000)
#petite boucle qui va reechantillonner 1000 fois avec remise les donnees comprises dans
#le vecteur dyadsObs$distance, et calculer a chaque fois la moyenne des valeurs reechantillonnees
for (i in 1:10000){
  boots[i] = mean(sample(dyadsObs$distance,replace=T))
}
#trace l'histogramme de boot
hist(boots)
#calcule l'intervalle de confiance comme les 2,5% et 97,5% percentiles de la distribution
#de boots
mean(boots)
ICDispersal = quantile(boots,c(.025,.975))
ICDispersal
abline(v=ICDispersal,col="red")
#calcule la moyenne, l'ecart-type et l'erreur standard apres bootstrap et les range dans
#des variables
# meanDispersal = mean(boots)
# meanDispersal
# etDisperal=sqrt(var(boots))
# etDisperal
# meanDispersal+c(-1,1)*qnorm(1-0.025)*etDisperal


#----------------------------------------------------------#
#   FIGURE
#----------------------------------------------------------#

####### Import data
# ALL NECESSARY DATA IS INCLUDED IN THIS FILES WRITEN EACH TIME YOU RUNNED THE PREVIOUS SCRIPT.
# YOU DON'T HAVE TO RERUN PREVIOUS SCRIPT FOR COMPUTING THE PLOT
dyadsH0 = read.table("AssignationPic/outputs/dyadsH0.txt",header=TRUE)
dyadsObsSelect = read.table("AssignationPic/dyadsObsSelect_gametic.txt",header=TRUE)
meanDispersal = 17.56167
icH0 = read.table("AssignationPic/outputs/icH0.txt",header=TRUE)

        # ##### Fitted distribution are estimated in the Dispersal kernel directory
        # # Cumulative density function from fitdistrplus and exponential distribution
        # z=seq(0,100,by=0.01)
        # lambda=0.10013465
        # lambdaInf=0.07614281
        # lambdaSup=0.13780406
        # yfit=function( z ){1-exp(-1*lambda*z)}
        # yfitInf=function( z ){1-exp(-1*lambdaInf*z)}
        # yfitSup=function( z ){1-exp(-1*lambdaSup*z)}
        # # Cumulative density function from fitdistrplus and gamma distribution
        # a=1
        # b=9.970022
        # yfitGamma=function( z ){pgamma(z,shape=a,scale=b)}
        # # Cumulative density function from MasterBayes estimates
        # lambda2=0.4704699
        # lambda2Inf=0.8413
        # lambda2Sup=0.2828
        # ybayes=function( z ){1-exp(-1*lambda2*z)}
        # ybayesInf=function( z ){1-exp(-1*lambda2Inf*z)}
        # ybayesSup=function( z ){1-exp(-1*lambda2Sup*z)}
        # 
        # #####################
        # # DYADS RELATIVE FREQUENCY IN FUNCTION OF DISTANCE
        # #####################
        # df=data.frame(d=c(sort(dyadsH0$distance),sort(dyadsObsSelect$distance)),
        #               freq=c(seq(1,length(dyadsH0$distance),1)/length(dyadsH0$distance),seq(1,length(dyadsObsSelect$distance),1)/length(dyadsObsSelect$distance)),
        #               dist=c(rep("H0",length(dyadsH0$distance)),rep("ObsSelect",length(dyadsObsSelect$distance)))
        # )
        # write.table(df, "AssignationPic/outputs/distancesDistPic_gametic.txt", row.names=FALSE,col.names=TRUE)
        # df1 = df[which(df$dist=="H0"),]
        # # df2=df[which(df$dist=="Obs"),]
        # df3 = df[which(df$dist=="ObsSelect"),]
        # 
        # ggplot(data=df1, aes(x=d, y=freq)) +
        #   stat_function(fun = yfit)+
        #   stat_function(fun = yfitInf,lty=2)+
        #   stat_function(fun = yfitSup,lty=2)+
        #   # stat_function(fun = yfitGamma,size=1.5,lty=3)+
        #   stat_function(fun = ybayes,colour="grey")+
        #   stat_function(fun = ybayesInf,colour="grey",lty=2)+
        #   stat_function(fun = ybayesSup,colour="grey",lty=2)+
        #   # geom_point(aes(color="expected"),size=4)+
        #   geom_ribbon(data=icH0,aes(x=distance,y=meanFreq,ymin=lowerIC, ymax=upperIC),alpha=0.5) +
        #   # geom_point(data=df2,aes(color="observed"),size=4) +
        #   geom_point(data=df3,aes(color="observed with selection")) +
        #   # geom_line(data=icH0,aes(x=distance,y=lowerIC),color="grey",size=1.5) +
        #   # geom_line(data=icH0,aes(x=distance,y=upperIC),color="grey",size=1.5) +
        #   scale_colour_manual(values=c("black"))+
        #   ylim(0,1) +
        #   # scale_x_continuous(breaks=c(0,round(meanDispersal,digits=1),20,40,60),limits=c(0,60)) +
        #   # ggtitle(paste("Distribution of offspring-father distances\n(n=",nrow(df3),")"))+
        #   xlab("Geographic distance (km)") + ylab("Dyads relative frequency") +
        #   theme(axis.line = element_line(colour = "black"),
        #         panel.grid.major = element_blank(),
        #         panel.grid.minor = element_blank(),
        #         panel.border = element_blank(),
        #         panel.background = element_blank(),
        #         plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        #         axis.title.x = element_text(color="black", size=14),
        #         axis.title.y = element_text(color="black", size=14),
        #         axis.text=element_text(size=14, colour="black"),
        #         legend.key = element_rect(fill = "white", size = 1),
        #         legend.key.height = unit(2,"line"),
        #         legend.key.width = unit(5,"line"),
        #         legend.text=element_text(size=14),
        #         legend.title=element_text(size=14),
        #         legend.position='none')

# ggsave("Figures/Distribution of offspring-father distances PicAll GAMETIC.tiff",
#        device="tiff",dpi=320,units="cm",width=19,height=13)

# Test sur les quantiles
# ks.test(x=cumsum(hist(dyadsObs$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density),
#         y=cumsum(hist(dyadsH0$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density))
ks.test(x=dyadsH0$distance,
        y=dyadsObsSelect$distance)

#------------------------------------------------
# Thuringia
#------------------------------------------------
#----
#-------------------------------
#   CONSTRUCTION DE LA DISTRIBUTION OBSERVEE
#-------------------------------#
#----
# abscisse : distance offspring-father
# ordonnee : frequence cumulee

# Dresser une liste de couples offspring-father retrouves par Colony (parente validee selon critere a definir)
# offspm=read.table("04_data/outputs/colRun1/RhinoThomas2.Paternity",h=T,fill=TRUE)
# dyadsObs=data.frame(offspring=offspm$OffspringID,father=offspm$InferredDad1,distance=rep(NA,nrow(offspm)))
dyadsObs = read.table("AssignationThu/resultsWithInfo.txt", header = TRUE)
dyadsObs = dyadsObs[,-c(5:16)]
# Gametic dispersal, hence exchange the sampling colony by the birth colony
# echanger les individus contre leurs colonies
info = read.table("AssignationThu/uniqueGenotypesWithInfo.txt",h=T)
MigrantFathersThu = read.table("Tables/MigrantFathersThu.txt", header = TRUE)
dyadsObs$birthColonyFather = NA
for(i in 1:nrow(dyadsObs)){
  dyadsObs$birthColonyFather[i] = as.character(MigrantFathersThu$ancestry[MigrantFathersThu$idind == dyadsObs$father[i]])
}

# pour chaque couple calculer la distance en km entre colonies
# formule : ACOS(SIN(lat1)*SIN(lat2)+COS(lat1)*COS(lat2)*COS(lon2-lon1))*6371
coordCol = read.table("AssignationThu/coordThu.txt",h=T)
# coordCol$Lat=coordCol$Lat*pi/180
# coordCol$Long=coordCol$Long*pi/180
# Distance between birth colony of the juvenile and inferred birth colony of the father by STRUCTURE
dyadsObs$distance = NA
for (i in 1:nrow(dyadsObs)) {
  # create distance matrix
  # cat(i)
  dyadsObs$distance[i] = distCosine(as.matrix(cbind(coordCol$Long[which(coordCol$Colony==as.character(dyadsObs$colonyOffsp[i]))],
                                                    coordCol$Lat[which(coordCol$Colony==as.character(dyadsObs$colonyOffsp[i]))])),
                                    as.matrix(cbind(coordCol$Long[which(coordCol$Colony==as.character(dyadsObs$birthColonyFather[i]))],
                                                    coordCol$Lat[which(coordCol$Colony==as.character(dyadsObs$birthColonyFather[i]))])),
                                    r=6378137)/1000
}
dyadsObs$distance


#-------------------------------#
#   EXPORT NEW DATA TO GRAPH
#-------------------------------
# Run observed distribution with a threshold of 0.3
# write.table(dyadsObs,"AssignationThu/dyadsObsSelect_gametic.txt",
#             col.names=TRUE,row.names=FALSE)
dyadsObs = read.table("AssignationThu/dyadsObsSelect_gametic.txt", header = TRUE)

#-------------------------------#
#   MEAN DISPERSAL DISTANCE
#-------------------------------
# Intervalle de confiance de la moyenne : calcul par bootstrap
#########################

#creation d'un vecteur "boots" comprenant 1000 elements, pour le moment tous des "0"
boots = numeric(10000)
#petite boucle qui va reechantillonner 1000 fois avec remise les donnees comprises dans
#le vecteur dyadsObs$distance, et calculer a chaque fois la moyenne des valeurs reechantillonnees
for (i in 1:10000){
  boots[i] = mean(sample(dyadsObs$distance,replace=T))
}
#trace l'histogramme de boot
hist(boots)
#calcule l'intervalle de confiance comme les 2,5% et 97,5% percentiles de la distribution
#de boots
mean(boots)
ICDispersal = quantile(boots,c(.025,.975))
ICDispersal
abline(v=ICDispersal,col="red")
#calcule la moyenne, l'ecart-type et l'erreur standard apres bootstrap et les range dans
#des variables
# meanDispersal = mean(boots)
# meanDispersal
# etDisperal=sqrt(var(boots))
# etDisperal
# meanDispersal+c(-1,1)*qnorm(1-0.025)*etDisperal


#-------------------------------#
#   GRAPH
#-------------------------------

####### Import data
# ALL NECESSARY DATA IS INCLUDED IN THIS FILES WRITEN EACH TIME YOU RUNNED THE PREVIOUS SCRIPT.
# YOU DON'T HAVE TO RERUN PREVIOUS SCRIPT FOR COMPUTING THE PLOT
dyadsH0 = read.table("AssignationThu/outputs/dyadsH0.txt",header=TRUE)
dyadsObsSelect = read.table("AssignationThu/dyadsObsSelect_gametic.txt",header=TRUE)
meanDispersal = 22.24434
icH0 = read.table("AssignationThu/outputs/icH0.txt",header=TRUE)

      # ##### Fitted distribution are estimated in the Dispersal kernel directory
      # # Cumulative density function from fitdistrplus and exponential distribution
      # z=seq(0,100,by=0.01)
      # lambda=0.10013465
      # lambdaInf=0.07614281
      # lambdaSup=0.13780406
      # yfit=function( z ){1-exp(-1*lambda*z)}
      # yfitInf=function( z ){1-exp(-1*lambdaInf*z)}
      # yfitSup=function( z ){1-exp(-1*lambdaSup*z)}
      # # Cumulative density function from fitdistrplus and gamma distribution
      # a=1
      # b=9.970022
      # yfitGamma=function( z ){pgamma(z,shape=a,scale=b)}
      # # Cumulative density function from MasterBayes estimates
      # lambda2=0.4704699
      # lambda2Inf=0.8413
      # lambda2Sup=0.2828
      # ybayes=function( z ){1-exp(-1*lambda2*z)}
      # ybayesInf=function( z ){1-exp(-1*lambda2Inf*z)}
      # ybayesSup=function( z ){1-exp(-1*lambda2Sup*z)}
      # 
      # #####################
      # # DYADS RELATIVE FREQUENCY IN FUNCTION OF DISTANCE
      # #####################
      # df=data.frame(d=c(sort(dyadsH0$distance),sort(dyadsObsSelect$distance)),
      #               freq=c(seq(1,length(dyadsH0$distance),1)/length(dyadsH0$distance),seq(1,length(dyadsObsSelect$distance),1)/length(dyadsObsSelect$distance)),
      #               dist=c(rep("H0",length(dyadsH0$distance)),rep("ObsSelect",length(dyadsObsSelect$distance)))
      # )
      # write.table(df, "AssignationThu/outputs/distancesDistThu_gametic.txt", row.names=FALSE,col.names=TRUE)
      # df1 = df[which(df$dist=="H0"),]
      # # df2=df[which(df$dist=="Obs"),]
      # df3 = df[which(df$dist=="ObsSelect"),]
      # 
      # ggplot(data=df1, aes(x=d, y=freq)) +
      #   stat_function(fun = yfit)+
      #   stat_function(fun = yfitInf,lty=2)+
      #   stat_function(fun = yfitSup,lty=2)+
      #   # stat_function(fun = yfitGamma,size=1.5,lty=3)+
      #   stat_function(fun = ybayes,colour="grey")+
      #   stat_function(fun = ybayesInf,colour="grey",lty=2)+
      #   stat_function(fun = ybayesSup,colour="grey",lty=2)+
      #   # geom_point(aes(color="expected"),size=4)+
      #   geom_ribbon(data=icH0,aes(x=distance,y=meanFreq,ymin=lowerIC, ymax=upperIC),alpha=0.5) +
      #   # geom_point(data=df2,aes(color="observed"),size=4) +
      #   geom_point(data=df3,aes(color="observed with selection")) +
      #   # geom_line(data=icH0,aes(x=distance,y=lowerIC),color="grey",size=1.5) +
      #   # geom_line(data=icH0,aes(x=distance,y=upperIC),color="grey",size=1.5) +
      #   scale_colour_manual(values=c("black"))+
      #   ylim(0,1) +
      #   # scale_x_continuous(breaks=c(0,round(meanDispersal,digits=1),20,40,60),limits=c(0,60)) +
      #   # ggtitle(paste("Distribution of offspring-father distances\n(n=",nrow(df3),")"))+
      #   xlab("Geographic distance (km)") + ylab("Dyads relative frequency") +
      #   theme(axis.line = element_line(colour = "black"),
      #         panel.grid.major = element_blank(),
      #         panel.grid.minor = element_blank(),
      #         panel.border = element_blank(),
      #         panel.background = element_blank(),
      #         plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
      #         axis.title.x = element_text(color="black", size=14),
      #         axis.title.y = element_text(color="black", size=14),
      #         axis.text=element_text(size=14, colour="black"),
      #         legend.key = element_rect(fill = "white", size = 1),
      #         legend.key.height = unit(2,"line"),
      #         legend.key.width = unit(5,"line"),
      #         legend.text=element_text(size=14),
      #         legend.title=element_text(size=14),
      #         legend.position='none')

# ggsave("Figures/Distribution of offspring-father distances ThuAll GAMETIC.tiff",
#        device="tiff",dpi=320,units="cm",width=19,height=13)

# Test sur les quantiles
# ks.test(x=cumsum(hist(dyadsObs$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density),
#         y=cumsum(hist(dyadsH0$distance,plot=F,breaks=seq(0,ceiling(max(c(dyadsH0$distance,dyadsObs$distance))),1))$density))
ks.test(x=dyadsH0$distance,
        y=dyadsObsSelect$distance)


#------------------------------------------------#
# END
#------------------------------------------------#
