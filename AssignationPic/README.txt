#===============================#
#            Colony
#   Assignation de parente
#===============================#


Dans ce dossier sont stockes tous les scripts pour parametrer Colony, preparer le jeu de donnees, ainsi que toutes les donnees d'inputs et d'outputs.

#-------------------------------#
# Dossiers
#-------------------------------#

### Dossier 01_prepadata

Ce dossier contient les scripts R necessaires pour manipuler les donnees de la BDD.

Le script 01_ExtractData contient toutes les commandes necessaires.

    * La fonction dbExtract() renvoie un objet BDD qui contient toutes les donnees extraites de la base de donnees PostgreSQL

    * La fonction UniqueGenotypes(BDD,falseloci,withcolony=F) permet d'ecrire le fichier contenant les genotypes uniques
    
    * Il contient quelques instructions pour decrire le jeu de donnees (effectifs, sexe-ratio)

    * Il fournit le fichier uniqueGenotypesWithInfos.txt necessaire pour la deuxieme partie.

### Dossier 02_colony

Ce dossier contient les scripts R necessaires pour configurer Colony et preparer les inputs.

Il permet de preparer :
- CMS.txt (liste des peres avec genotype)
- CFS.txt (liste des meres avec genotype)
- OFS.txt (liste des descendants avec genotype)
- ExcludedMothers.txt (matrice des meres exclues)
- ExcludedFathers.txt (matrice des peres exclus)
- ExcludedOffsprings.txt (matrice des offsprings exclus)
- Colony.Dat
- shellColony.txt

Les scripts 01_prepaColony, 02_ExcludedMothers, 03_ExcludedFathers et 04_ExlcudedOffsprings sont faits pour etre executes dans l'ordre, les fichiers R de fonctions ne sont pas a executer.

#### 03_simulations

Ce dossier contient les scripts R necessaires pour configurer et preparer les inputs de datasets simules.

Chaque jeu de donnees simule est place dans /inputs/"nom de la simulation"

### 04_Data

Ce dossier contient les donnees d'outputs de Colony, pour les simulations et les donnees observees.

Il doit contenir :
- /outputs     # contient les resultats de Colony sur le serveur
      - /Colony # Resultats des 50 runs de Colony sur donnees observees
      - /simReplicatedRuns48      # simulation 50 runs sur petite population (48 parents)
      - /simReplicatedRuns48rep1  # simulation replicat 1 avec ech aleatoire
      - /simReplicatedRuns48rep2  # simulation replicat 2 avec ech aleatoire
      - /simReplicatedRuns48rep3  # simulation replicat 3 avec ech aleatoire
      - /simReplicatedRuns2400    # simulation 50 runs sur grande population (2400 parents)
      - /simIncompleteLoci        # simulation quand 2 loci incomplets
      - /simFatherIncluded10      # simulation de la probabilite 0.1 que le pere ne soit pas echantillonne
      - /simFatherIncluded20      # simulation de la probabilite 0.2 que le pere ne soit pas echantillonne
      - /simFatherIncluded30      # simulation de la probabilite 0.3 que le pere ne soit pas echantillonne
      
      
      Chacun de ces dossiers contient :
      - Colony.Dat
      - qsColThomas (qui lance la tache de calcul)
      - CommColPutty (pour nettoyer le repertoire de travail et afficher les resultats)
      - MLcol.R
      - parentscol.R
      - preppedfiles.R
      - uniqueGenotypesWithInfo.txt
      - tous les fichiers de resultat


#-------------------------------#
# Serveur de calcul
#-------------------------------#

hostname : 147.99.162.220
login : tbrazier
mdp : thobra

Apres le transfert d'un fichier windows sur le serveur :
dos2unix file

Creer des fichiers avec le script de Gilles :
colonyAutoDat.py -i MMCohorts_Run1.Dat -m ?? -min ??? -max ?????

-o /home/users/petit/ouTuVeux/
et tu travailles dans ton home ;-)

m : nombre de replicats
min et max : bornes entre lesquelles sont tires les graines pour le tirage aleatoire pour COLONY


Commandes pour lancer colony (sur PuTTY) :
bash qsColThomas pour travailler directement sur le serveur (attention, session est occupee tout le temps du calcul)
OU
qsub qscolThomas pour utiliser le cluster de calcul

Pour savoir ou en est le serveur :
qstat


#-------------------------------#
# Gestion de la base de donnees
#-------------------------------#

### Creation d'une base de donnees en local

Demander a Gilles Maubert

Installer PostgreSQL et PGadmin

Fichier : petit_rhinolophe.sql, extrait du serveur PVE 05.

Creer un role "pljan"


### Connexion en local

Base de donnee en PostGreSQL

server : localhost
name : petit_rhinolophe
username : postgres
password : postgres

Code R pour la connexion :
# Nouvelle methode de connexion SQL
library(RPostgreSQL)
# driver PostGreSQL
drv <- dbDriver("PostgreSQL")
# cree une connexion avec la BDD
conn <- dbConnect(drv, host="localhost", user="postgres", password="postgres", dbname="petit_rhinolophe")
# pour se deconnecter
dbDisconnect(conn)


### Connexion sur le serveur PVE 05  ## Ne pas utiliser pour l'instant

Base de donnee en PostGreSQL

server : 
name : petit_rhinolophe
username : 
password : 

Code R pour la connexion :
# Nouvelle methode de connexion SQL
library(RPostgreSQL)
# driver PostGreSQL
drv <- dbDriver("PostgreSQL")
# cree une connexion avec la BDD
conn <- dbConnect(drv, host="", user="", password="", dbname="")
# pour se deconnecter
dbDisconnect(conn)