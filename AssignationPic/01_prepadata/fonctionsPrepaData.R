#===============================#
# FUNCTIONS FOR PREPARING DATA
#===============================#


#-------------------------------#
# dbExtract : 
# Recuperation des donnees de la base de donnees en postGreSQL
#-------------------------------#

dbExtract=function() {
  
  #-------------------------------#
  # Fonction pour l'importation des donnees de la base de donnees
  #-------------------------------
  timeIn=Sys.time()

  # Nouvelle methode de connexion SQL
  # installer packages RPostgreSQL et DBI
  library(RPostgreSQL)
  # driver PostGreSQL
  drv <- dbDriver("PostgreSQL")
  # cree une connexion avec la BDD
  conn <- dbConnect(drv, host="localhost", user="postgres", password="postgres", dbname="petit_rhinolophe")
  
  #-------------------------------#
  # Creation d'un fichier qui contient toutes les infos du jeu de donnees
  #-------------------------------#
  
  #-------------------------------#
  # Recuperation de toutes les donnees
  #-------------------------------#
  
  # jointure SQL de genotype, corres_col et corres_ind en utilisant les criteres idgen, year_sample et idind
  allInfos=dbGetQuery(conn,"select distinct idind,dby8,dby8b,rha101,rha101b,rha109,rha109b,rha4,rha4b,rha7,rha7b,rhc108,rhc108b,rhc3,rhc3b,
                      rhd102,rhd102b,rhd103,rhd103b,idcol,year_sample,col2013_1,repet2013_1,col2013_2,
                      repet2013_2,col2014_1,repet2014_1,col2014_2,repet2014_2,col2015_1,repet2015_1,col2015_2,repet2015_2,
                      col2016_1,repet2016_1,col2016_2,repet2016_2 from genotype inner join corres_col using(idgen,year_sample) inner join corres_ind using(idgen) inner join ind using(idind)")
  write.table(allInfos,paste(wd,"allInfosBDD.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  
  
  # Deconnexion de la BDD
  dbDisconnect(conn)
  
  timeOut=Sys.time()
  cat("Recuperation des donnees en",difftime(timeOut,timeIn),"secondes.")
  # renvoie les tableaux
  return(allInfos)
  
}


#-------------------------------#
# UniqueGenotypes : 
# Fonction pour construire le jeu de donnees avec les individus (genotypes uniques)
#-------------------------------#
# Construit le jeu de donnees
# argument 1 : BDD contient toutes les donnees
# argument 2 : falseloci= le nombre max de loci incorrects acceptes avant d'exclure l'individu
# e.g. falseloci=2 signifie que tous les ind avec 0,1 ou 2 loci incomplets sont conserves
# argument 3 : withcolony=TRUE pour avoir les donnees d'echantillonnage de chaque session
# argument 4 : sexmark=TRUE pour garder le marqueur sexuel dby8

UniqueGenotypes=function(BDD,falseloci=1,withcolony=F,sexmark=F) {
  timeIn=Sys.time()
  
  #-------------------------------#
  # Creer un fichier de tous les genotypes uniques
  #-------------------------------#
  
  # Structure : Names (idind), Alleles, colonie (idcol), sexe, age premiere capture (juv,adult), annee de premiere capture
  uniqueGenotypesWithInfo=matrix(nrow=nrow(BDD),ncol=39)
  colnames(uniqueGenotypesWithInfo)=c("idind","dby8","dby8b","rha101","rha101b","rha109","rha109b","rha4","rha4b","rha7",
                                      "rha7b","rhc108","rhc108b","rhc3","rhc3b","rhd102","rhd102b","rhd103","rhd103b",
                                      "idcol","sexe","ageWhenFirstCaptur","yearWhenFirstCaptur","col2013_1","repet2013_1",
                                      "col2013_2","repet2013_2","col2014_1","repet2014_1","col2014_2","repet2014_2",
                                      "col2015_1","repet2015_1","col2015_2","repet2015_2","col2016_1","repet2016_1",
                                      "col2016_2","repet2016_2")
  uniqueGenotypesWithInfo=as.data.frame(uniqueGenotypesWithInfo)
  
  
  # Remplir le tableau avec les infos connues d'apres BDD
  uniqueGenotypesWithInfo[,1:20]=BDD[,1:20]
  uniqueGenotypesWithInfo[,23:39]=BDD[,21:37]
  
  # Meme principe pour la colonie et l'annee de capture -> on ne conserve que la premiere colonie ou il a ete capture
  # ATTENTION: comme seule la colonie de premiere capture est notee,
  # il faut aller chercher dans les colonnes 24:39 pour voir les mouvements entre colonies
  
  # Ne garder que les colonies etudiees
  uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[(uniqueGenotypesWithInfo$idcol %in%
                                                     c('CXSGT','M244','M809','M811','M813','M815','M831','M1079','M1102','M1106',
                                                       'M1112','M1113','M1551','M1975','M1979','PP1','PP3')),]
  
  
  #-------------------------------#
  # Supprimer les doublons d'idind
  #-------------------------------#
  n2013i=length(unique(uniqueGenotypesWithInfo$idind[which(uniqueGenotypesWithInfo$yearWhenFirstCaptur==2013)]))
  # on ne souhaite conserver que la premiere apparition de l'individu (premiere annee de capture)
  # Ajouter les lignes pour chaque individu qui repondent a :
  # / 1 idind a la fois
  # / min des annees de sampling
  
  allInfosSansDb=uniqueGenotypesWithInfo[-(1:nrow(uniqueGenotypesWithInfo)),]
  for (i in unique(uniqueGenotypesWithInfo$idind)) {
    doublons=uniqueGenotypesWithInfo[which(uniqueGenotypesWithInfo$idind==i),]
    allInfosSansDb=rbind(allInfosSansDb,doublons[which.min(doublons$yearWhenFirstCaptur),])
    
    #   allInfosSansDb=rbind(allInfosSansDb,
    #                        uniqueGenotypesWithInfo[which(uniqueGenotypesWithInfo$idind==i
    #                        & uniqueGenotypesWithInfo$yearWhenFirstCaptur==min(uniqueGenotypesWithInfo$yearWhenFirstCaptur[which(uniqueGenotypesWithInfo$idind==i)],na.rm=T)),])
  }
  
  # TEST de l'algorithme pour enlever les doublons
  # compter les effectifs uniques totaux et effectifs uniques 2013
  # conditions de validite
  # 1 / nombre de lignes correspond au nombre d'individus
  # 1 / n avant = n apres
  if (nrow(allInfosSansDb)!=length(unique(uniqueGenotypesWithInfo$idind))) {
    warning("Pas le bon nombre d'individus apres la procedure de suppression des doublons.")
  } else{uniqueGenotypesWithInfo=allInfosSansDb}
  
  
  # 2 / n2013 avant = n2013 apres
  if (length(unique(uniqueGenotypesWithInfo$idind[which(uniqueGenotypesWithInfo$yearWhenFirstCaptur==2013)]))!=n2013i) {
    warning("Probleme dans la suppression des doublons, des individus de 2013 ont disparu.")
  }
  
  #-------------------------------#
  # Enlever les donnees manquantes
  #-------------------------------#
  
  # Option 1
  # On garde l'individu quand (9-"falseloci") loci ont ete correctement genotypes : 7/9 d'apres Jan (2017)
  # grouper les alleles deux a deux (produit des deux alleles)
  # si 1 des alleles=0 -> FALSE (0), si  les 2 alleles != 0 -> TRUE (1)
  # si la somme des loci faux pour une ligne > "falseloci" alors le genotype est rejete
  # P.S. : plus tard, les ind avec loci manquants sur le marqueur sexuel seront retires egalement
  loci=uniqueGenotypesWithInfo[,seq(2,18,2)]*uniqueGenotypesWithInfo[,seq(3,19,2)]
  uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[(rowSums(apply(loci,c(1,2),function(x) 0 %in% x),na.rm=TRUE))<=falseloci,]
  
  # Option 2
  # Methode Pierre-Loup
  # On enleve la ligne quand il y a au moins un locus a 0
  #uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[-which(0 %in% uniqueGenotypesWithInfo[2:19,]),]
  # uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[!apply(uniqueGenotypesWithInfo[2:19],1,function (x) 0 %in% x),]
  
  #-------------------------------#
  # Identifier les femelles
  #-------------------------------#
  
  # D'abord nettoyer le jeu de donnees en enlevant les ind sans marqueur sexuel (dbY8 ou dbY8b == 0)
  # en prenant la precaution de verifier si il existe des valeurs nulles dans les colonnes concernees
  if(prod(uniqueGenotypesWithInfo$dby8)==0 | prod(uniqueGenotypesWithInfo$dby8b)==0){
    uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[-which(uniqueGenotypesWithInfo$dby8==0 | uniqueGenotypesWithInfo$dby8b==0),]
  }
  
  # Option 1
  # Si l'individu n'a aucun de ces trois alleles, c'est une femelle (Zarsozo-Lacoste 2017)
  #colf=col[col$dby8!=128&col$dby8!=131&col$dby8!=133,]
  uniqueGenotypesWithInfo$sexe[which(uniqueGenotypesWithInfo$dby8!=128&uniqueGenotypesWithInfo$dby8!=131&uniqueGenotypesWithInfo$dby8!=133)]="F"
  
  # Option 2
  # Les femelles ont une difference absolue entre leurs deux alleles inferieure a 10
  #uniqueGenotypesWithInfo$sexe[which(abs(uniqueGenotypesWithInfo$dby8-uniqueGenotypesWithInfo$dby8b)<10)]="F"
  
  #-------------------------------#
  # Identifier les males
  #-------------------------------#  
  ###pour ne prendre que les males
  
  # Option 1
  # Si l'individu possede un de ces trois alleles, c'est un male (Zarsozo-Lacoste 2017)
  #colm=col[col$dby8==128|col$dby8==131|col$dby8==133,]
  uniqueGenotypesWithInfo$sexe[which(uniqueGenotypesWithInfo$dby8==128|uniqueGenotypesWithInfo$dby8==131|uniqueGenotypesWithInfo$dby8==133)]="M"
  
  # Option 2
  #uniqueGenotypesWithInfo$sexe[which(abs(uniqueGenotypesWithInfo$dby8-uniqueGenotypesWithInfo$dby8b)>=10)]="M"
  
  #-------------------------------#
  # Identifier les juveniles
  #-------------------------------#
  for (i in 1:nrow(uniqueGenotypesWithInfo)) {
    year=uniqueGenotypesWithInfo$yearWhenFirstCaptur[i]
    if (is.na(uniqueGenotypesWithInfo[i,which(colnames(uniqueGenotypesWithInfo)==paste("col",year,"_1",sep=""))])) {
      if (is.na(uniqueGenotypesWithInfo[i,which(colnames(uniqueGenotypesWithInfo)==paste("col",year,"_2",sep=""))])){
        # c'est une erreur, annee de premier echantillonnage ne correspond pas
        warning("Attention : annee de premier echantillonnage ne correspond pas")
      } else {
        # echantillonne pour la premier fois en deuxieme session donc considere comme juvenile
        uniqueGenotypesWithInfo$ageWhenFirstCaptur[i]="Juv"
      }
    } else {
      if (!(is.na(uniqueGenotypesWithInfo[i,which(colnames(uniqueGenotypesWithInfo)==paste("col",year,"_1",sep=""))]))) {
        # echantillonne pour la premiere fois en premiere session donc considere comme adulte
        uniqueGenotypesWithInfo$ageWhenFirstCaptur[i]="Adult"
      }
    }
  }
  
  #-------------------------------#
  # Enlever les dates d'echantillonnage
  #------------------------------- 
  uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[,-seq(25,39,2)]
  
  if (withcolony==F) {
    uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[,-c(24:31)]
  }
  
  #-------------------------------#
  # Enlever le marqueur sexuel dbY8
  #-------------------------------
  if (sexmark==F) {
    uniqueGenotypesWithInfo=uniqueGenotypesWithInfo[,-c(2,3)]
  }
  

  
  #-------------------------------#
  # Ecriture du fichier
  #-------------------------------
  write.table(uniqueGenotypesWithInfo,"uniqueGenotypesWithInfo.txt",row.names=F,col.names=T,quote=F,sep="\t")
  
  timeOut=Sys.time()
  cat("Fichier des genotypes uniques construit en",difftime(timeOut,timeIn),"secondes.")
  
  return(uniqueGenotypesWithInfo)
}


#-------------------------------#
# Tableau des effectifs
#-------------------------------#
# Dans cette version les individus ne sont comptes qu'une seule fois : la premiere annee de capture

effectifsTable=function(unique) {
  
  JuvF2013=0
  JuvF2014=0
  JuvF2015=0
  JuvF2016=0
  JuvM2013=0
  JuvM2014=0
  JuvM2015=0
  JuvM2016=0
  AdultF2013=0
  AdultF2014=0
  AdultF2015=0
  AdultF2016=0
  AdultM2013=0
  AdultM2014=0
  AdultM2015=0
  AdultM2016=0
  
  for (y in 2013:2016) {
    l=which(unique$yearWhenFirstCaptur==y)
    for (i in l) { # pour chaque individu (ligne) de l'annee
      assign(paste(unique$ageWhenFirstCaptur[i],unique$sexe[i],y,sep = ""),
             get(paste(unique$ageWhenFirstCaptur[i],unique$sexe[i],y,sep = ""))+1)
    }
  }
  JuvF=c(JuvF2013,JuvF2014,JuvF2015,JuvF2016) # Nombre de nouveaux juveniles femelles par an, 2013:2016
  JuvM=c(JuvM2013,JuvM2014,JuvM2015,JuvM2016) # Nombre de nouveaux juveniles males par an
  AdultF=c(AdultF2013,AdultF2014,AdultF2015,AdultF2016) # Nombre de nouveaux adultes femelles par an
  AdultM=c(AdultM2013,AdultM2014,AdultM2015,AdultM2016) # Nombre de nouveaux adultes males par an
  
  years=c(2013,2014,2015,2016)
  sumPerYear=JuvF+JuvM+AdultF+AdultM
  effectifs=cbind(years,JuvF,JuvM,AdultF,AdultM,sumPerYear)
  sumPerClass=c(NA,sum(JuvF),sum(JuvM),sum(AdultF),sum(AdultM),sum(effectifs[,6],na.rm=T))
  effectifs=rbind(effectifs,sumPerClass)
  write.table(effectifs,"effectifs.txt",row.names=F,col.names=T,quote=F,sep="\t")
  return(effectifs)
}


tests=function(unique,n,verbose=F) {
  
  #-------------------------------#
  # Fonction lanÃ§ant des tests automatiques
  # pour comparer le jeu de donnees et la BDD
  #-------------------------------#
  
  library(RPostgreSQL)
  
  
  #-------------------------------#
  # Routine en boucle
  #-------------------------------#
  for (i in 1:n) {
    # driver PostGreSQL
    drv <- dbDriver("PostgreSQL")
    # cree une connexion avec la BDD
    conn <- dbConnect(drv, host="localhost", user="postgres", password="postgres", dbname="petit_rhinolophe")
    
    # tirage d'un idind au hasard
    id=sample(unique$idind,1)
    # si le flag passe a faux durant la procedure, une erreur est notee, avec l'idind enregistre
    flag=TRUE
    err=c() # liste des ind incorrects
    
    #-------------------------------#
    # Verification des genotypes
    #-------------------------------#
    
    # compare genotype (fichier sortie) et genotype (BDD) ~ idind
    # recupere le genotype dans unique
    gen=unique[which(unique$idind==id),2:17]
    genBDD=dbGetQuery(conn,paste("select rha101,rha101b,rha109,rha109b,rha4,rha4b,rha7,rha7b,rhc108,rhc108b,
                                 rhc3,rhc3b,rhd102,rhd102b,rhd103,rhd103b
                                 FROM genotype inner join corres_col using(idgen,year_sample) inner join corres_ind using(idgen) inner join ind using(idind)
                                 WHERE idind=",id,";",sep=""))
    
    # Deconnexion de la BDD
    dbDisconnect(conn)
    
    
    if (verbose==T) {
      print(gen)
      print(genBDD)
    }
    
    if (sum(gen!=genBDD[1,])>0) {
      flag=FALSE
      err=c(err,id)
      cat(i,": Erreur de genotype sur l'individu",id,"! Different de la BDD.\n")
    }
    
    
    
  }
  
  #-------------------------------#
  # Verification des colonies et donnees d'echantillonnage
  #-------------------------------#
  
  # verification manuelle sur un sous echantillonnage
  cat("Verification des colonies :\n")
  for (i in 1:n) {
    # driver PostGreSQL
    drv <- dbDriver("PostgreSQL")
    # cree une connexion avec la BDD
    conn <- dbConnect(drv, host="localhost", user="postgres", password="postgres", dbname="petit_rhinolophe")
    # tirage d'un idind au hasard
    id=sample(unique$idind,1)
    
    cat(unique$idind[which(unique$idind==id)]," : ",as.character(unique$idcol[which(unique$idind==id)])," : ",unique$ageWhenFirstCaptur[which(unique$idind==id)]," : ",unique$yearWhenFirstCaptur[which(unique$idind==id)],"\n")
    print(unique[which(unique$idind==id),22:29])
    samplingBDD=dbGetQuery(conn,paste("select distinct col2013_1,col2013_2,
                                      col2014_1,col2014_2,col2015_1,col2015_2,
                                      col2016_1,col2016_2 from genotype inner join corres_col using(idgen,year_sample)
                                      inner join corres_ind using(idgen) inner join ind using(idind)
                                      WHERE idind=",id,sep=""))
    print(samplingBDD)
    cat("\n")
    
    # Deconnexion de la BDD
    dbDisconnect(conn)
  }
  
  
  if (flag) {
    cat("Pas d'erreur detectee sur l'ensemble des tests.")
  }
  
}


