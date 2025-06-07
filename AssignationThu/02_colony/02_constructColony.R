#===============================#
#            Assignation
#   Preparation Inputs de COLONY
#
#     Construction du fichier COLONY.dat
#===============================#


# D'apres un script fourni par Pierre-Loup JAN

# DIRECTORIES

if (Sys.info()["sysname"]=="Darwin") {
  wd="/Volumes/Samsung_T5/INRA RHINO/R/AssignationThu/02_colony/"
  graphdir="/Volumes/Samsung_T5/INRA RHINO/R/Graph/"
  Rdir="/Volumes/Samsung_T5/INRA RHINO/R/"
  setwd(wd)
  source(paste(wd,"fonctionPrepaColony.R",sep=""))
  source(paste(wd,"fonctionExcluded.R",sep=""))
}else {
  if (Sys.info()["sysname"]=="Windows"){
    wd="E:/INRA RHINO/R/AssignationThu/02_colony/"
    graphdir="E:/INRA RHINO/R/Graph/"
    Rdir="E:/INRA RHINO/R/"
    setwd(wd)
    source(paste(wd,"fonctionPrepaColony.R",sep=""))
    source(paste(wd,"fonctionExcluded.R",sep=""))
  }
}




#-----------------------#
#       A executer
#-----------------------#
# Effacer les anciens fichiers
nettoyage=c("Rhino_OFS.txt","Rhino_CMS.txt","Rhino_CFS.txt","Rhino_ExcludedMothers.txt","Rhino_ExcludedFathers.txt","Rhino_ExcludedMaternalSibs.txt")
file.remove(paste(wd,nettoyage,sep=""))

prepaColony() # prepare les fichiers OFS, CMS et CFS
excluded() # prepare fichiers d'exclusion

# Creation des .Dat
constructColony(30) # pour obtenir le fichier Colony.dat en n replicats (n graines aleatoires differentes)
shellCommand(30) # pour obtenir le fichier qsColThomas avec n lignes de commandes








