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
source("fonctionExcluded.R")
source("fonctionPrepaColony.R")


# Copy the file "uniqueGenotypesWithInfo.txt"
file.copy("../01_prepadata/uniqueGenotypesWithInfo.txt", "uniqueGenotypesWithInfo.txt")


# Effacer les anciens fichiers
cleanup=c("Rhino_OFS.txt","Rhino_CMS.txt","Rhino_CFS.txt","Rhino_ExcludedMothers.txt","Rhino_ExcludedFathers.txt","Rhino_ExcludedMaternalSibs.txt")
file.remove(cleanup)

prepaColony() # prepare les fichiers OFS, CMS et CFS
excluded() # prepare fichiers d'exclusion

# Creation des .Dat
constructColony(30) # pour obtenir le fichier Colony.dat en n replicats (n graines aleatoires differentes)
# Pensez a changer le repertoire /home sur le serveur dans qsCol (1ere ligne bash)
shellCommand(30) # pour obtenir le fichier qsCol avec n lignes de commandes








