#===============================#
#            COLONY
#         Simulations
#    Analyse des resultats
#===============================#

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
source("fonctionsSim.R")



#-----------------------------#
# SENSITIVITY ANALYSES
# COUNT TABLES
#-----------------------------#
# Count individuals in each assignation category (true, type I or II error...)

# tableau de comptage des individus dans les differentes categories (true, type I or II error...)
# pour chaque simulation
# tableau de comparaison des nombres d'individus moyens et SD sur n runs
# saved dans le fichier simSensitivity.txt

# Take as arguments
# 1/ path = the path to the simulations output files
# 2/ run = number of runs used for estimations
tabSensitivity=sensitivityResults(path=c("simRandomPop",
                                         "simRandomPopIncomplete"),30)
View(tabSensitivity)

tabSensitivity=sensitivityResultsMothers(path=c("simRandomPop",
                                         "simRandomPopIncomplete"),30)
View(tabSensitivity)

#-----------------------------#
# CONSTRUCTION DE LA DISTRIBUTION DE LA STAT ~ NOMBRE DE RUNS
#-----------------------------#
# Compute statistics on COLONY output probabilities

# Pour les stats sur les probabilities de Colony
# en arguments :

# Take as arguments
# 1/ path, the path to the COLONY simulation outputs
# 2/ number of runs
# 3/ stats to compute : "mean","median","var"
distMean=statBestF("simRandomPopIncomplete",run=30,stat="mean")
write.table(distMean,file="meanBestFather.txt",quote=F,row.names = F,col.names = F)

distMedian=statBestF("simRandomPopIncomplete",run=30,stat="median")
write.table(distMedian,file="medianBestFather.txt",quote=F,row.names = F,col.names = F)

distVar=statBestF("simRandomPopIncomplete",run=30,stat="var")
write.table(distVar,file="varBestFather.txt",quote=F,row.names = F,col.names = F)

# To get the number of occurences of the best father

# Pour le nombre d'occurences du meilleur pere
# 1/ path, the path to the COLONY simulation outputs
# 2/ number of runs
distOcc=occurenceBestF("simRandomPop",30)
write.table(distOcc,file="occurencesBestFather.txt",quote=F,row.names = F,col.names = T)
distOcc=occurenceBestF("simRandomPopIncomplete",30)
write.table(distOcc,file="occurencesBestFatherIncompleteGen.txt",quote=F,row.names = F,col.names = T)

# distOcc=read.table("occurencesBestFather.txt",h=T)
# OR
distOcc=read.table("occurencesBestFatherIncompleteGen.txt",h=T)

colnames(distOcc)=c("runs","occurence","father")
dfT=distOcc[which(distOcc$father==TRUE),1:2]
dfF=distOcc[which(distOcc$father==FALSE),1:2]
colnames(dfT)=c("runs","occurence")
colnames(dfF)=c("runs","occurence")
pM=ggplot(data=dfT, aes(x=runs, y=occurence)) +
  geom_point(colour="Black",size=1)+
  geom_point(data=dfF,colour="darkgrey",size=1) +
  geom_hline(yintercept=0.3,linetype="dashed",size=1) +
  ylim(0.1,0.7) +
  # ggtitle("Assignation frequencies of the best father\n",
         # subtitle=paste("n=1200",sep = "")) +
  xlab("Number of runs") + ylab("Assignment frequency") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        #legend.key.height = unit(2,"line"),
        #legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
pM


# To get the number of occurences of the best mother

# Pour le nombre d'occurences du meilleur pere
# 1/ path, the path to the COLONY simulation outputs
# 2/ number of runs
distOcc=occurenceBestM("simRandomPop",30)
write.table(distOcc,file="occurencesBestMother.txt",quote=F,row.names = F,col.names = T)
distOcc=occurenceBestM("simRandomPopIncomplete",30)
write.table(distOcc,file="occurencesBestMotherIncompleteGen.txt",quote=F,row.names = F,col.names = T)

distOcc=read.table("occurencesBestMotherIncompleteGen.txt",h=T)

colnames(distOcc)=c("runs","occurence","father")
dfT=distOcc[which(distOcc$father==TRUE),1:2]
dfF=distOcc[which(distOcc$father==FALSE),1:2]
colnames(dfT)=c("runs","occurence")
colnames(dfF)=c("runs","occurence")
pF=ggplot(data=dfT, aes(x=runs, y=occurence)) +
  geom_point(colour="Black",size=1)+
  geom_point(data=dfF,colour="darkgrey",size=1) +
  geom_hline(yintercept=0.3,linetype="dashed",size=1) +
  ylim(0.1,0.7) +
  # ggtitle("Assignation frequencies of the best father\n",
          # subtitle=paste("n=1200",sep = "")) +
  xlab("Number of runs") + ylab("Assignment frequency") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        #legend.key.height = unit(2,"line"),
        #legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
pF

ggarrange(pM,pF,widths=1:1)
ggsave(paste(graphdir,"/Assignation frequencies of the best father Pic79.tiff",sep=""),
       device="tiff",dpi=320,units="cm",width=20,height=8)



#-----------------------------#
# SIM RESULTS WITH INFOS - EXPORT DATA
#-----------------------------#
# give a file containing :
# - offsrping ID
# - fathers and mothers ID
# - sex
# - correct/incorrect assignation
# - genotype

simResults=simResultsWithInfo(path="simRandomPop")
simResults=simResultsWithInfo(path="simRandomPopIncomplete")
















###################################
# EXPERIMENTAL                    #
###################################
# 
# #-------------------------------#
# # simReplicatedRuns48
# #-------------------------------#
# 
# #---------------------------------#
# # Tableau de comptage
# #---------------------------------#
# # fournir l'argument du chemin de la simulation
# # ex : simReplicatedRuns48
# tabAssigned48=countAssigned("simReplicatedRuns48")
# 
# # pour faire une moyenne sur plusieurs runs
# # SD = dispersion des donnees
# tabAssigned48_50runs=countAssignedRuns("simReplicatedRuns48",50)
# 
# #---------------------------------#
# # CALCUL DES STATS DE TESTS
# #---------------------------------
# # col1 : offspring ID
# # col2 : vrai pere
# # col3 : sampled (0/1)
# # col4 : meilleur pere retrouve
# # col5 : meilleur pere = vrai pere (0/1)
# # col6 : moyenne de probDad par Colony sur l'ensemble des runs pour le meilleur pere
# # col7 : mediane de probDad par Colony sur l'ensemble des runs pour le meilleur pere
# # col8 : variance de probDad par Colony sur l'ensemble des runs pour le meilleur pere
# 
# #---------------------------------#
# # FONCTION
# #---------------------------------
# # fournir en argument :
# # 1/ le nombre de runs (defaut = 50)
# # 2/ le chemin du dossier de la simulation
# # ex : simReplicatedRuns48
# tab1=statSim(50,"simReplicatedRuns48")
# tab2=statSimBestFather(50,"simReplicatedRuns48")
# 
# # A la recherche de la meilleure stat associee a la probabilite de trouver le bon pere
# # on compare les stats entre les bons peres et les mauvais peres
# falseDad=as.data.frame(split(tab1,tab1$trueDad)[1])
# colnames(falseDad)=c("offspring","realFather","sampled","meanProbObs","medianProbObs","varProbObs","bestFather","probBestFather","trueDad")
# trueDad=as.data.frame(split(tab1,tab1$trueDad)[2])
# colnames(trueDad)=c("offspring","realFather","sampled","meanProbObs","medianProbObs","varProbObs","bestFather","probBestFather","trueDad")
# 
# wilcox.test(meanProbObs~trueDad,data=tab1)
# wilcox.test(medianProbObs~trueDad,data=tab1)
# wilcox.test(varProbObs~trueDad,data=tab1)
# 
# #-------------------------------#
# # simReplicatedRuns2400
# #-------------------------------
# tabAssigned2400=countAssigned("simReplicatedRuns2400")
# tabAssigned2400_10runs=countAssignedRuns("simReplicatedRuns2400",50)
# 
# #----
# tab1_2400=statSim(50,"simReplicatedRuns2400")
# tab2_2400=statSimBestFather(50,"simReplicatedRuns2400")
# # A la recherche de la meilleure stat associee a la probabilite de trouver le bon pere
# # on compare les stats entre les bons peres et les mauvais peres
# falseDad=as.data.frame(split(tab2_2400,tab2_2400$trueDad)[1])
# colnames(falseDad)=c("offspring","realFather","sampled","bestFather","probBestFather","trueDad","meanProbObs","medianProbObs","varProbObs")
# trueDad=as.data.frame(split(tab2_2400,tab2_2400$trueDad)[2])
# colnames(trueDad)=c("offspring","realFather","sampled","bestFather","probBestFather","trueDad","meanProbObs","medianProbObs","varProbObs")
# 
# # trueDad$meanProbObs
# # falseDad$meanProbObs
# length(which(falseDad$meanProbObs<0.95))/nrow(falseDad)
# length(which(trueDad$meanProbObs<0.95))/nrow(trueDad)
# 
# # trueDad$medianProbObs
# # falseDad$medianProbObs
# length(which(falseDad$medianProbObs<0.95))/nrow(falseDad)
# length(which(trueDad$medianProbObs<0.95))/nrow(trueDad)
# 
# # trueDad$varProbObs
# # falseDad$varProbObs
# length(which(falseDad$varProbObs>0.05))/nrow(falseDad)
# length(which(trueDad$varProbObs>0.05))/nrow(trueDad)
# 
# wilcox.test(meanProbObs~trueDad,data=tab1_2400)
# wilcox.test(medianProbObs~trueDad,data=tab1_2400)
# wilcox.test(varProbObs~trueDad,data=tab1_2400)
# 
# #----
# # tests par permutation
# #----
# 
# 
# # Fournir un df avec une colonne probability et une colonne true dad
# 
# # Sur les valeurs de stat calculees pour chaque offspring
# # Offsrpings independants
# dataPerm=statSimBestFather(50,"simReplicatedRuns2400")
# # enlever les lignes ou le pere == NA
# dataPerm=dataPerm[-which(is.na(dataPerm$trueDad)),]
# 
# ###############################
# # Difference des moyennes des valeurs de la stat sur probability entre peres corrects et incorrects
# # Test par permutation sur la moyenne de probability
# 
# dal = numeric(10000) # table de la somme des valeurs avant/apres pour 10000 permutations
# for (n in 1:9999) {  # n repetitions
#   for (i in 1:nrow(dataPerm)) { # resample pere correct/incorrect pour chaque individu
#     dataPerm$perm = sample(dataPerm$trueDad, nrow(dataPerm), replace = FALSE)
#   }
#   dal[n]=mean(as.numeric(dataPerm$meanProbObs[which(dataPerm$perm==1)])) - mean(as.numeric(dataPerm$meanProbObs[which(dataPerm$perm==0)])) # somme des valeurs d'un groupe
# }
# # ajout de la valeur observee a la distribution
# dal[10000]=mean(dataPerm$meanProbObs[which(dataPerm$trueDad==1)]) - mean(dataPerm$meanProbObs[which(dataPerm$trueDad==0)])
# dal
# 
# # Table disjonctive qui compte le nombre de valeurs dans dal
# # qui sont superieures ou egales a la difference observee entre les moyennes
# tap=table(abs(dal)>=abs(mean(dataPerm$meanProbObs[which(dataPerm$trueDad==1)]) - mean(dataPerm$meanProbObs[which(dataPerm$trueDad==0)])))
# tap
# # calcul de la valeur de p
# p=tap[2]/10000
# p
# 
# #trouve les limites de l'intervalle de confiance pour la somme
# # des elements du premier groupe, sous l'hypothese nulle
# limsum=quantile(dal,c(.025,.975))
# limsum
# dal[10000]
# hist(dal, main = "Mean probability value differences between correct and incorrect fathers",
#      xlab="Differences of the mean probability values") # distribution de la difference des moyennes entre les deux groupes sous l'hypothèse H0
# abline(v=limsum[1], col = "red") # quantile a 2.5%
# abline(v=limsum[2], col = "red") # quantile a 97.5%
# abline(v=dal[10000],lty=2,lwd=2) # valeur observee
# 
# ###############################
# # Difference des medianes des valeurs de la stat sur probability entre peres corrects et incorrects
# # Test par permutation sur la moyenne de probability
# 
# dal = numeric(10000) # table de la somme des valeurs avant/apres pour 10000 permutations
# for (n in 1:9999) {  # n repetitions
#   for (i in 1:nrow(dataPerm)) { # resample pere correct/incorrect pour chaque individu
#     dataPerm$perm = sample(dataPerm$trueDad, nrow(dataPerm), replace = FALSE)
#   }
#   dal[n]=median(as.numeric(dataPerm$meanProbObs[which(dataPerm$perm==1)])) - median(as.numeric(dataPerm$meanProbObs[which(dataPerm$perm==0)])) # somme des valeurs d'un groupe
# }
# # ajout de la valeur observee a la distribution
# dal[10000]=median(dataPerm$meanProbObs[which(dataPerm$trueDad==1)]) - median(dataPerm$meanProbObs[which(dataPerm$trueDad==0)])
# dal
# 
# # Table disjonctive qui compte le nombre de valeurs dans dal
# # qui sont superieures ou egales a la difference observee entre les moyennes
# tap=table(abs(dal)>=abs(median(dataPerm$meanProbObs[which(dataPerm$trueDad==1)]) - median(dataPerm$meanProbObs[which(dataPerm$trueDad==0)])))
# tap
# # calcul de la valeur de p
# p=tap[2]/10000
# p
# 
# #trouve les limites de l'intervalle de confiance pour la somme
# # des elements du premier groupe, sous l'hypothese nulle
# limsum=quantile(dal,c(.025,.975))
# limsum
# dal[10000]
# hist(dal, main = "Median probability value differences between correct and incorrect fathers",
#      xlab="Differences of the median probability values") # distribution de la difference des moyennes entre les deux groupes sous l'hypothèse H0
# abline(v=limsum[1], col = "red") # quantile a 2.5%
# abline(v=limsum[2], col = "red") # quantile a 97.5%
# abline(v=dal[10000],lty=2,lwd=2) # valeur observee
# 
# 
# 
# ###############################
# # Difference des variances des valeurs de la stat sur probability entre peres corrects et incorrects
# # Test par permutation sur la moyenne de probability
# 
# dal = numeric(10000) # table de la somme des valeurs avant/apres pour 10000 permutations
# for (n in 1:9999) {  # n repetitions
#   for (i in 1:nrow(dataPerm)) { # resample pere correct/incorrect pour chaque individu
#     dataPerm$perm = sample(dataPerm$trueDad, nrow(dataPerm), replace = FALSE)
#   }
#   dal[n]=var(as.numeric(dataPerm$meanProbObs[which(dataPerm$perm==1)])) - var(as.numeric(dataPerm$meanProbObs[which(dataPerm$perm==0)])) # somme des valeurs d'un groupe
# }
# # ajout de la valeur observee a la distribution
# dal[10000]=var(dataPerm$meanProbObs[which(dataPerm$trueDad==1)]) - var(dataPerm$meanProbObs[which(dataPerm$trueDad==0)])
# dal
# 
# # Table disjonctive qui compte le nombre de valeurs dans dal
# # qui sont superieures ou egales a la difference observee entre les moyennes
# tap=table(abs(dal)>=abs(var(dataPerm$meanProbObs[which(dataPerm$trueDad==1)]) - var(dataPerm$meanProbObs[which(dataPerm$trueDad==0)])))
# tap
# # calcul de la valeur de p
# p=tap[2]/10000
# p
# 
# #trouve les limites de l'intervalle de confiance pour la somme
# # des elements du premier groupe, sous l'hypothese nulle
# limsum=quantile(dal,c(.025,.975))
# limsum
# dal[10000]
# hist(dal, main = "Variance probability value differences between correct and incorrect fathers",
#      xlab="Differences of the var probability values") # distribution de la difference des moyennes entre les deux groupes sous l'hypothèse H0
# abline(v=limsum[1], col = "red") # quantile a 2.5%
# abline(v=limsum[2], col = "red") # quantile a 97.5%
# abline(v=dal[10000],lty=2,lwd=2) # valeur observee
# 
# 
# 



