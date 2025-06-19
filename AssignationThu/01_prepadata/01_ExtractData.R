#===============================#
#            Assignation
#   Extraction des donnees de la BDD
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
source("fonctionsDataThu.R")


#-------------------------------#
# Lancer la fonction
#-------------------------------#
# Convertit le fichier All_Genotypes_FINAL_man.csv en uniqueGenotypesWithInfos.txt
unique=convertToUnique(falseloci=2,sexmark=FALSE)

#-------------------------------#
# Construction du tableau des effectifs
#-------------------------------#
# effectifs=effectifsTable(unique)
# View(effectifs)

#-------------------------------#
# Analyse du jeu de donnees uniqueGenotypesWithInfo
#-------------------------------#
totalF=sum(unique$sexe=="F") # nombre total de femelles
totalF
totalM=sum(unique$sexe=="M") # nombre total de males
totalM

# n total
totalF+totalM

# n juveniles
sum(unique$ageWhenFirstCaptur=="Juv")

# sexe ratio (% de males)
totalM/(totalF+totalM)*100

#-------------------------------#
# Routines de test
#-------------------------------
#source(paste(wd,"01fonction_tests.R",sep=""))
# Renseigner l'objet "unique", le nombre de tirages et si vous voulez un rapport plus detaille
# tests(unique,10,verbose=F)



