#===============================#
#            Assignation
#   Preparation Inputs de COLONY
#
#     Construction du fichier COLONY.dat
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
source("../../Sources/packages.R")
source("fonctionPrepaColony.R")
source("fonctionExcluded.R")


unique=read.table("../uniqueGenotypesWithInfo.txt",h=T)

#-----------------------#
#       A executer
#-----------------------#
# Effacer les anciens fichiers
# nettoyage=c("Rhino_OFS.txt","Rhino_CMS.txt","Rhino_CFS.txt","Rhino_ExcludedMothers.txt","Rhino_ExcludedFathers.txt","Rhino_ExcludedMaternalSibs.txt")
# file.remove(paste(nettoyage,sep=""))

prepaColony() # prepare les fichiers OFS, CMS et CFS
excluded() # prepare fichiers d'exclusion

# Creation des .Dat
constructColony(30) # pour obtenir le fichier Colony.dat en n replicats (n graines aleatoires differentes)
shellCommand(30) # pour obtenir le fichier qsColThomas avec n lignes de commandes








