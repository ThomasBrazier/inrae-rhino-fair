#===============================#
#            COLONY
#         Simulations
#===============================#
# clear global environment: remove all variables
rm(list=ls(all=TRUE))

#----------------------------------------------------------#
# Loading packages
# Check if packages are installed, install if necessary
source("../../Sources/packages.R")
source("fonctionsPrepaSim.R")


#-------------------------------#
# CREATION D'UNE POP THEORIQUE ALEATOIRE
#-------------------------------

##### MATRICE D'ACCOUPLEMENT ALEATOIRE
# Procedure de force brute
# Tirage de matrices jusqu'a obtenir matrice conforme
# n = nombre de sous-structure independantes dans la matrice
# dont depend le nombre d'individus adultes = n*24 parents
# t = nombre de generations
###### Vous n'avez pas a utiliser cette fonction pour construire votre pop, elle est deja integree
# dans ConstructSimInput()
# Les commandes sont presentees a titre informatif pour visualiser a quoi ressemble la matrice d'accouplement produite
# Mating=MatingMatrix()
# OR
# Mating=as.matrix(read.table("MatriceAccouplement.txt",header=FALSE,sep=" "))

library(rstudioapi)
wd=dirname(rstudioapi::getSourceEditorContext()$path)

##### CONSTRUCTION DE L'INPUT POUR SIMULATION
# Met en forme le fichier Input.Par
# avec les parametres choisis et la matrice d'accouplement
# 1/ n transmis a la fonction MatingMatrix()
# 2/ t transmis a la fonction MatingMatrix()
# 3/ pop, taille de la popualtion simulee, qui sera arrondie a l'inferieur pour etre un multiple de la taille de la mating matrix
ConstructSimInput(n=10,t=4,pop=4000,name="Simulation")

##### LANCEMENT DE LA SIMULATION DU MODULE COLONY
# Lancer le module simulation de COLONY manuellement


# conserver uniquement les fichiers :

##### CREATION D'UN FICHIER COLONY.DAT SUIVANT UN SCENARIO
# CreateScenario()
# Scenario d'echantillonnage et de genotypage
# 1/ Proba de detection
# Les individus sont retires aleatoirement du fichier COLONY.DAT selon un tirage sans remise
# Proportion d'individus retires de chaque categorie fixee en argument
# - Unsampled, sous la forme d'un vecteur c(males, femelles, juv)
# Si les valeurs sont superieures a 1, elles sont prises comme des effectifs et non des proportions
# 2/ Erreurs de genotypages
# Sur les individus restants, une proportion d'entre eux est aleatoirement selectionnee pour avoir un genotype incomplet
# Les alleles sont retires selon une probabilite fixee en argument
# - propIncomplete (si proportion > 1, alors on considere que c'est un effectif a retirer)
# - probMissingLocus
CreateScenario(name="Simulation",unsampled=c(0.44,0.1,0.1),propIncomplete=0,probMissingLocus=0)
# Parfois un message d'erreur est produit
# Warning message:
#   In .Internal(grepl(as.character(pattern), x, ignore.case, FALSE,  :
#                        closing unused connection 3 (COLONY2.DAT)
# Ne pas en tenir compte et relancer la commande une seconde fois

replicateSim("Simulation",30)
simShellCmd("Simulation",30)

# Exemple de la Picardie
CreateScenario(name="simRandomPop",unsampled=c(0.44,0.1,0),propIncomplete=0,probMissingLocus=0)
CreateScenario(name="simRandomPopIncomplete",unsampled=c(0.44,0.1,0),propIncomplete=0.2,probMissingLocus=2/16)
CreateScenario(name="simRandomPopFullySampled",unsampled=c(0,0,0),propIncomplete=0,probMissingLocus=0)
# Replication du jeu de donnees .Dat simule en n exemplaires identiques avec une graine aleatoire differente
replicateSim("simRandomPop",30)
simShellCmd("simRandomPop",30)
replicateSim("simRandomPopIncomplete",30)
simShellCmd("simRandomPopIncomplete",30)
replicateSim("simRandomPopFullySampled",30)
simShellCmd("simRandomPopFullySampled",30)








