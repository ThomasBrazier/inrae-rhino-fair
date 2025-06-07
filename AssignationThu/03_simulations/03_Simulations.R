#===============================#
#            COLONY
#         Simulations
#===============================#

# DIRECTORIES
if (Sys.info()["sysname"]=="Darwin") {
  setwd("/Volumes/Samsung_T5/INRA RHINO/R/AssignationThu/03_simulations/")
  wd="/Volumes/Samsung_T5/INRA RHINO/R/AssignationThu/03_simulations/"
  graphdir="/Volumes/Samsung_T5/INRA RHINO/R/Graph/"
  Rdir="/Volumes/Samsung_T5/INRA RHINO/R/"
  source("fonctionsPrepaSim.R")
}else {
  if (Sys.info()["sysname"]=="Windows"){
    setwd("E:/INRA RHINO/R/AssignationThu/03_simulations/")
    wd="E:/INRA RHINO/R/AssignationThu/03_simulations/"
    graphdir="E:/INRA RHINO/R/Graph/"
    Rdir="E:/INRA RHINO/R/"
    source("fonctionsPrepaSim.R")
  }
}


#-------------------------------#
# Replicated Runs
#-------------------------------
# Creation des fichiers de config de la simulation
# Chaque simulation doit etre lancee dans COLONY Windows Only

# Replication du jeu de donnees .Dat simule en n exemplaires identiques avec une graine aleatoire differente
replicateSim("Pop2400rep3",30)
simShellCmd("Pop2400rep3",30)


#----------------------------------#
# FREQUENCES ALLELIQUES ALLEMANDES
#----------------------------------
library(related)
info=read.table("uniqueGenotypesWithInfo.txt",h=T)
write.table(info[,1:17],"genotypeData.txt",sep=" ",col.names = FALSE,row.names = FALSE)
allelFreq=readgenotypedata("genotypeData.txt")$freqs
f=c()
for (i in 1:8) {
  f=c(f,paste(unlist(t(round(allelFreq[[i]],digits=4))[2,]),collapse=" "))
}
write.table(f,"FreqAllelique.txt",col.names=F,row.names=F)

#-------------------------------#
# CREATION D'UNE POP THEORIQUE ALEATOIRE
#-------------------------------

##### MATRICE D'ACCOUPLEMENT ALEATOIRE
# Procedure de force brute
# Tirage de matrices jusqu'a obtenir matrice conforme
# n = nombre de sous-structure independantes dans la matrice
# dont depend le nombre d'individus adultes = n*24 parents
# t = nombre de generations
Mating=MatingMatrix()
# OR
Mating=as.matrix(read.table("MatriceAccouplement.txt",header=FALSE,sep=" "))

##### CONSTRUCTION DE L'INPUT POUR SIMULATION
# Met en forme le fichier Input.Par
# avec les parametres choisis et la matrice d'accouplement
# 1/ n transmis a la fonction MatingMatrix()
# 2/ t transmis a la fonction MatingMatrix()
# 3/ pop, taille de la popualtion simulee, qui sera arrondie a l'inferieur pour etre un multiple de la taille de la mating matrix
ConstructSimInput(n=10,t=4,pop=4000,name="Simulation")

##### LANCEMENT DE LA SIMULATION DU MODULE COLONY
# LaunchSim()
# Lance le module simulation de COLONY
# nettoie le dossier en enlevant tous les fichiers inutiles




##### CREATION D'UN FICHIER COLONY.DAT SUIVANT UN SCENARIO
# CreateScenario()
# Scenario d'echantillonnage et de genotypage
# 1/ Proba de detection
# Les individus sont retires aleatoirement du fichier COLONY2.DAT selon un tirage sans remise
# Proportion d'individus retires de chaque categorie fixee en argument
# - Unsampled, sous la forme d'un vecteur c(males, femelles, juv)
# Si les valeurs sont superieures a 1, elles sont prises comme des effectifs et non des proportions
# 2/ Erreurs de genotypages
# Sur les individus restants, une proportion d'entre eux est aleatoirement selectionnee pour avoir un genotype incomplet
# Les alleles sont retires selon une probabilite fixee en argument
# - propIncomplete (si proportion > 1, alors on considere que c'est un effectif a retirer)
# - probMissingLocus
CreateScenario(name="simRandomPopThu",unsampled=c(0.44,0.1,0.1),propIncomplete=0,probMissingLocus=0)
CreateScenario(name="simRandomPopThuIncomplete",unsampled=c(0.44,0.1,0.1),propIncomplete=0.2,probMissingLocus=2/16)
# Replication du jeu de donnees .Dat simule en n exemplaires identiques avec une graine aleatoire differente
replicateSim("simRandomPopThu",30)
simShellCmd("simRandomPopThu",30)
replicateSim("simRandomPopThuIncomplete",30)
simShellCmd("simRandomPopThuIncomplete",30)







