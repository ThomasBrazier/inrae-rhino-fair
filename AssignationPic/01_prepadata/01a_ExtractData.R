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
source("fonctionsPrepaData.R")


#-------------------------------#
# Get the data
#-------------------------------#

BDD=read.table("allInfosBDD.txt",h=T,sep="\t")
# Construit le jeu de donnees
# argument 1 : BDD contient toutes les donnees
# argument 2 : falseloci= le nombre max de loci incorrects acceptes avant d'exclure l'individu
# e.g. falseloci=2 signifie que tous les ind avec 0,1 ou 2 loci incomplets sont conserves
# argument 3 : withcolony=TRUE pour avoir les donnees d'echantillonnage de chaque session
# argument 4 : sexmark=TRUE pour garder le marqueur sexuel dby8
unique=UniqueGenotypes(BDD,falseloci=2,withcolony=T,sexmark=F)
# le data set pret est dans uniqueGenotypesWithInfo.txt

#-------------------------------#
# Sample sizes
#-------------------------------#
effectifs=effectifsTable(unique)
View(effectifs)



