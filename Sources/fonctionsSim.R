#===============================#
#            COLONY
#         Simulations
#    Analyse des resultats
#===============================#



#---------------------------------#
# COUNTING ASSIGNATIONS
#---------------------------------#

# fournir l'argument du chemin de la simulation
# ex : simReplicatedRuns48
countAssigned=function(path=""){
  # sur run 1
  
  # comptage en 5 categories des juveniles
  # question : pour un jeune, est-ce qu'on retrouve le bon pere ?
  # 1/ pere est echantillone
  #   1.a/ bon pere (P)
  #   1.b/ pas de pere retrouve
  #   1.c/ mauvais pere assigne (P)
  # 2/ pere non echantillone
  #   2.a/ on trouve un pere (P)
  #   2.b/ pas de pere retrouve
  
  # A partir de ces categories, on veut maximiser le contraste en  1.a et 2.b
  # Pour trois categories, on a des valeurs de probability fournies par Colony (P)
  
  
  require(data.table)
  # Recuperation de la liste des offsprings
  # generer un fichier texte avec les offsprings
  offsprings=read.table(paste(path,"/offsprings.txt",sep=""),sep=",",h=F)
  dyads=as.data.frame(offsprings[,1])
  colnames(dyads)="offspring"
  # Recuperation de la liste des paternites
  if (max(count.fields(paste(path,"/",path,"C_1.Paternity",sep="")))==3) {
    paternity=read.table(paste(path,"/",path,"C_1.Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
  } else {
    paternity=fread(paste(path,"/",path,"C_1.Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
    paternity=paternity[,1:3]
  }
  # generer un fichier txt avec la liste des peres echantillonnes
  fathers=read.table(paste(path,"/fathers.txt",sep=""),sep=",",h=F)
  fathers=fathers[,1]
  
  
  # reconstitution des dyads
  # pour chaque juvenile, le meilleur pere trouve est assigne (selon la valeur de probability fournie par Colony)
  # C'est le premier pere indique par Colony
  for (i in 1:nrow(dyads)) {
    # vrai pere
    dyads$realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",dyads$offspring[i])
    # on affecte les meilleurs peres aux juveniles quand ils existent
    if (dyads$offspring[i] %in% paternity$OffspringID) {
      dyads$bestFather[i]=paternity$InferredDad1[which(paternity$OffspringID==dyads$offspring[i])]
      dyads$probDad[i]=paternity$ProbDad1[which(paternity$OffspringID==dyads$offspring[i])]
    } else {
      dyads$bestFather[i]=NA
      dyads$probDad[i]=NA
    }
    # ajout de l'indicateur echantillonne ou pas
    if (dyads$realFather[i] %in% fathers) {
      dyads$sampled[i]=1
    } else {
      dyads$sampled[i]=0
    }
  }
  
  # tableau de comptage
  results=data.frame(offsprings=c(sum(dyads$sampled==1 & dyads$realFather==dyads$bestFather,na.rm=T),
                                  sum(dyads$sampled==1 & is.na(dyads$bestFather),na.rm=T),
                                  sum(dyads$sampled==1 & dyads$realFather!=dyads$bestFather & !is.na(dyads$bestFather),na.rm=T),
                                  sum(dyads$sampled==0 & !is.na(dyads$bestFather),na.rm=T),
                                  sum(dyads$sampled==0 & is.na(dyads$bestFather),na.rm=T)
  ))
  rownames(results)=c("MsampledMtrue","MsampledMnotfound","MsampledMfalse","MnotsampledMfound","MnotsampledMnotfound")
  # Proportion de chaque categorie
  # Nombre d'offsprings divise par nombre total d'offsprings ech dans la pop
  results$propObserved=results$offsprings/nrow(offsprings)
  # Proportion attendue dans la population
  # Dans notre cas, il s'agit du nombre d'offsprings attendu dans une categorie
  # divise par le nombre total d'offsprings dans la pop
  # par exemple tous les offsprings dont le pere a ete echantillonne et qui sont bien assignes
  # divise par le nombre total d'offsprings -> proportion attendue de resultats vrais dans la population totale
  # en consequence, la somme des proportions fait 1
  results$propExpected=c(sum(dyads$sampled==1)/nrow(offsprings),0,0,0,sum(dyads$sampled==1)/nrow(offsprings))
  # Probabilite fournie par Colony
  results$probDadColony=c(mean(dyads$probDad[which(dyads$sampled==1 & dyads$realFather==dyads$bestFather)],na.rm=T),
                          NA,
                          mean(dyads$probDad[which(dyads$sampled==1 & dyads$realFather!=dyads$bestFather & !is.na(dyads$bestFather))],na.rm=T),
                          mean(dyads$probDad[which(dyads$sampled==0 & !is.na(dyads$bestFather))],na.rm=T),
                          NA
  )
  return(results)
}

# pour faire une moyenne sur plusieurs runs
# SD = dispersion des donnees
countAssignedRuns=function(path="",run=50){
  # sur n runs
  
  # comptage en 5 categories des juveniles
  # question : pour un jeune, est-ce qu'on retrouve le bon pere ?
  # 1/ pere est echantillone
  #   1.a/ bon pere (P)
  #   1.b/ pas de pere retrouve
  #   1.c/ mauvais pere assigne (P)
  # 2/ pere non echantillone
  #   2.a/ on trouve un pere (P)
  #   2.b/ pas de pere retrouve
  
  # A partir de ces categories, on veut maximiser le contraste en  1.a et 2.b
  # Pour trois categories, on a des valeurs de probability fournies par Colony (P)
  
  require(data.table)
  # Recuperation de la liste des offsprings
  # generer un fichier texte avec les offsprings
  offsprings=read.table(paste(path,"/offsprings.txt",sep=""),sep=",",h=F)
  dyads=as.data.frame(offsprings[,1])
  colnames(dyads)="offspring"
  # generer un fichier txt avec la liste des peres echantillonnes
  fathers=read.table(paste(path,"/fathers.txt",sep=""),sep=",",h=F)
  fathers=fathers[,1]
  
  # tableaux vides
  effectifs=data.frame(c(0,0,0,0,0))
  propObserved=data.frame(c(0,0,0,0,0))
  probDadColony=data.frame(c(0,0,0,0,0))
  
  for (r in 1:run) {
    # Recuperation de la liste des paternites
    if (max(count.fields(paste(path,"/",path,"C_",r,".Paternity",sep="")))==3) {
      paternity=read.table(paste(path,"/",path,"C_",r,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      paternity=fread(paste(path,"/",path,"C_",r,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      paternity=paternity[,1:3]
    }
    
    
    # reconstitution des dyads
    # pour chaque juvenile, le meilleur pere trouve est assigne (selon la valeur de probability fournie par Colony)
    # C'est le premier pere indique par Colony
    for (i in 1:nrow(dyads)) {
      # vrai pere
      dyads$realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",dyads$offspring[i])
      # on affecte les meilleurs peres aux juveniles quand ils existent
      if (dyads$offspring[i] %in% paternity$OffspringID) {
        dyads$bestFather[i]=paternity$InferredDad1[which(paternity$OffspringID==dyads$offspring[i])]
        dyads$probDad[i]=paternity$ProbDad1[which(paternity$OffspringID==dyads$offspring[i])]
      } else {
        dyads$bestFather[i]=NA
        dyads$probDad[i]=NA
      }
      # ajout de l'indicateur echantillonne ou pas
      if (dyads$realFather[i] %in% fathers) {
        dyads$sampled[i]=1
      } else {
        dyads$sampled[i]=0
      }
    }
    
    # tableau de comptage
    effectifs[,r]=c(sum(dyads$sampled==1 & dyads$realFather==dyads$bestFather,na.rm=T),
                    sum(dyads$sampled==1 & is.na(dyads$bestFather),na.rm=T),
                    sum(dyads$sampled==1 & dyads$realFather!=dyads$bestFather & !is.na(dyads$bestFather),na.rm=T),
                    sum(dyads$sampled==0 & !is.na(dyads$bestFather),na.rm=T),
                    sum(dyads$sampled==0 & is.na(dyads$bestFather),na.rm=T)
    )
    # Proportion de chaque categorie
    # Nombre d'offsprings divise par nombre total d'offsprings ech dans la pop
    propObserved[,r]=effectifs[,r]/sum(effectifs[,r])
    
    # Probabilite fournie par Colony
    probDadColony[,r]=c(mean(dyads$probDad[which(dyads$sampled==1 & dyads$realFather==dyads$bestFather)],na.rm=T),
                        NA,
                        mean(dyads$probDad[which(dyads$sampled==1 & dyads$realFather!=dyads$bestFather & !is.na(dyads$bestFather))],na.rm=T),
                        mean(dyads$probDad[which(dyads$sampled==0 & !is.na(dyads$bestFather))],na.rm=T),
                        NA
    )
  }
  
  
  # Proportion attendue dans la population
  # Dans notre cas, il s'agit du nombre d'offsprings attendu dans une categorie
  # divise par le nombre total d'offsprings dans la pop
  # par exemple tous les offsprings dont le pere a ete echantillonne et qui sont bien assignes
  # divise par le nombre total d'offsprings -> proportion attendue de resultats vrais dans la population totale
  # en consequence, la somme des proportions fait 1
  results=data.frame(offspring=rowSums(effectifs),
                     meanPropObserved=apply(propObserved,1,mean),
                     sdPropObserved=apply(propObserved,1,sd),
                     propExpected=c(sum(dyads$sampled==1)/nrow(offsprings),0,0,0,sum(dyads$sampled==1)/nrow(offsprings)),
                     meanProbDadColony=apply(probDadColony,1,mean))
  rownames(results)=c("MsampledMtrue","MsampledMnotfound","MsampledMfalse","MnotsampledMfound","MnotsampledMnotfound")
  
  return(results)
}

#----
#---------------------------------#
# CALCUL DES STATS DE TESTS
#---------------------------------#
#----

# col1 : offspring ID
# col2 : vrai pere
# col3 : sampled (0/1)
# col4 : meilleur pere retrouve
# col5 : meilleur pere = vrai pere (0/1)
# col6 : moyenne de probDad par Colony sur l'ensemble des runs pour le meilleur pere
# col7 : mediane de probDad par Colony sur l'ensemble des runs pour le meilleur pere
# col8 : variance de probDad par Colony sur l'ensemble des runs pour le meilleur pere

#---------------------------------#
# FONCTION
#---------------------------------#
# fournir en argument :
# 1/ le nombre de runs (defaut = 50)
# 2/ le chemin du dossier de la simulation
# ex : simReplicatedRuns48

statSim=function(run=50,path=""){
  library(data.table)
  run=run+1
  offsprings=read.table(paste(path,"/offsprings.txt",sep=""),sep=",",h=F)
  dyads=as.data.frame(offsprings[,1])
  colnames(dyads)="offspring"
  # Recuperation de la liste des paternites
  if (max(count.fields(paste(path,"/",path,"C_1.Paternity",sep="")))==3) {
    paternity=read.table(paste(path,"/",path,"C_1.Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
  } else {
    paternity=fread(paste(path,"/",path,"C_1.Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
    paternity=paternity[,1:3]
  }
  # generer un fichier txt avec la liste des peres echantillonnes
  fathers=read.table(paste(path,"/fathers.txt",sep=""),sep=",",h=F)
  fathers=fathers[,1]
  # reconstitution des dyads
  # pour chaque juvenile, le meilleur pere trouve est assigne (selon la valeur de probability fournie par Colony)
  # C'est le premier pere indique par Colony
  for (i in 1:nrow(dyads)) {
    # vrai pere
    dyads$realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",dyads$offspring[i])
    # ajout de l'indicateur echantillonne ou pas
    if (dyads$realFather[i] %in% fathers) {
      dyads$sampled[i]=1
    } else {
      dyads$sampled[i]=0
    }
  }
  dyads
  
  # ProbDad : liste des probability de Colony a chaque run 1-50
  probDad=as.data.frame(dyads$offspring)
  colnames(probDad)="OffspringID"
  for (i in 1:(run-1)) {
    # proportions d'assignations pour paternite
    if (max(count.fields(paste(path,"/",path,"C_",i,".Paternity",sep="")))==3) {
      pat=read.table(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      pat=fread(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      pat=pat[,1:3]
    }
    
    for (j in 1:nrow(probDad)) {
      if(probDad$OffspringID[j] %in% pat$OffspringID){
        probDad[j,i+1]=pat$ProbDad1[which(pat$OffspringID==as.character(probDad$OffspringID[j]))]
      } else {
        probDad[j,i+1]=NA
      }
    }
  }
  
  # ListFathers : liste des meilleurs peres assignes a chaque run 1-50
  ListFathers=as.data.frame(dyads$offspring)
  colnames(ListFathers)="OffspringID"
  for (i in 1:(run-1)) {
    # proportions d'assignations pour paternite
    if (max(count.fields(paste(path,"/",path,"C_",i,".Paternity",sep="")))==3) {
      pat=read.table(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      pat=fread(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      pat=pat[,1:3]
    }
    
    for (j in 1:nrow(ListFathers)) {
      if(ListFathers$OffspringID[j] %in% pat$OffspringID){
        ListFathers[j,i+1]=pat$InferredDad1[which(pat$OffspringID==as.character(probDad$OffspringID[j]))]
      } else {
        ListFathers[j,i+1]=NA
      }
    }
  }
  # nombre d'occurences du meilleur pere
  # identifier le meilleur pere : celui qui revient le plus souvent dans les assignations
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre dyads et ListFathers
  for (i in 1:nrow(dyads)){
    uniq=unique(na.omit(as.character(ListFathers[i,2:run])))
    
    if(length(uniq)==0){
      dyads$bestFather[i]=NA
    }else{
      best=c(NA,0)
      for (j in 1:length(uniq)) {
        occur=length(which(ListFathers[i,2:run]==uniq[j]))
        if(occur>best[2]){
          best[1]=uniq[j]
          best[2]=occur
        }
      }
      dyads$bestFather[i]=best[1]
    }
  }
  # compter le nombre d'occurences du meilleur pere/nombre de runs
  for(i in 1:nrow(dyads)){
    dyads$probBestFather[i]=length(which(ListFathers[i,2:run]==as.character(dyads$bestFather[i])))/(run-1)
  }
  
  # true dad ?
  # vrai pere est trouve si le meilleur pere (le plus frequent) correspond au vrai pere
  dyads$trueDad=as.numeric(dyads$realFather==dyads$bestFather)
  
  # colonne meanProbObs
  dyads$meanProbObs=apply(probDad[2:run],1,function(x) mean(x,na.rm=T))
  # colonne varProbObs
  dyads$medianProbObs=apply(probDad[2:run],1,function(x) median(x,na.rm=T))
  # colonne varProbObs
  dyads$varProbObs=apply(probDad[2:run],1,function(x) var(x,na.rm=T))
  
  
  
  return(dyads)
}

statSimBestFather=function(run=50,path=""){
  library(data.table)
  run=run+1
  offsprings=read.table(paste(path,"/offsprings.txt",sep=""),sep=",",h=F)
  dyads=as.data.frame(offsprings[,1])
  colnames(dyads)="offspring"
  # Recuperation de la liste des paternites
  if (max(count.fields(paste(path,"/",path,"C_1.Paternity",sep="")))==3) {
    paternity=read.table(paste(path,"/",path,"C_1.Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
  } else {
    paternity=fread(paste(path,"/",path,"C_1.Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
    paternity=paternity[,1:3]
  }
  # generer un fichier txt avec la liste des peres echantillonnes
  fathers=read.table(paste(path,"/fathers.txt",sep=""),sep=",",h=F)
  fathers=fathers[,1]
  # reconstitution des dyads
  # pour chaque juvenile, le meilleur pere trouve est assigne (selon la valeur de probability fournie par Colony)
  # C'est le premier pere indique par Colony
  for (i in 1:nrow(dyads)) {
    # vrai pere
    dyads$realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",dyads$offspring[i])
    # ajout de l'indicateur echantillonne ou pas
    if (dyads$realFather[i] %in% fathers) {
      dyads$sampled[i]=1
    } else {
      dyads$sampled[i]=0
    }
  }
  dyads
  
  # ProbDad : liste des probability de Colony a chaque run 1-50
  probDad=as.data.frame(dyads$offspring)
  colnames(probDad)="OffspringID"
  for (i in 1:(run-1)) {
    # proportions d'assignations pour paternite
    if (max(count.fields(paste(path,"/",path,"C_",i,".Paternity",sep="")))==3) {
      pat=read.table(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      pat=fread(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      pat=pat[,1:3]
    }
    
    for (j in 1:nrow(probDad)) {
      if(probDad$OffspringID[j] %in% pat$OffspringID){
        probDad[j,i+1]=pat$ProbDad1[which(pat$OffspringID==as.character(probDad$OffspringID[j]))]
      } else {
        probDad[j,i+1]=NA
      }
    }
  }
  
  # ListFathers : liste des meilleurs peres assignes a chaque run 1-50
  ListFathers=as.data.frame(dyads$offspring)
  colnames(ListFathers)="OffspringID"
  for (i in 1:(run-1)) {
    # proportions d'assignations pour paternite
    if (max(count.fields(paste(path,"/",path,"C_",i,".Paternity",sep="")))==3) {
      pat=read.table(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      pat=fread(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      pat=pat[,1:3]
    }
    
    for (j in 1:nrow(ListFathers)) {
      if(ListFathers$OffspringID[j] %in% pat$OffspringID){
        ListFathers[j,i+1]=pat$InferredDad1[which(pat$OffspringID==as.character(probDad$OffspringID[j]))]
      } else {
        ListFathers[j,i+1]=NA
      }
    }
  }
  
  # nombre d'occurences du meilleur pere
  # identifier le meilleur pere : celui qui revient le plus souvent dans les assignations
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre dyads et ListFathers
  for (i in 1:nrow(dyads)){
    uniq=unique(na.omit(as.character(ListFathers[i,2:run])))
    
    if(length(uniq)==0){
      dyads$bestFather[i]=NA
    }else{
      best=c(NA,0)
      for (j in 1:length(uniq)) {
        occur=length(which(ListFathers[i,2:run]==uniq[j]))
        if(occur>best[2]){
          best[1]=uniq[j]
          best[2]=occur
        }
      }
      dyads$bestFather[i]=best[1]
    }
  }
  # compter le nombre d'occurences du meilleur pere/nombre de runs
  for(i in 1:nrow(dyads)){
    dyads$probBestFather[i]=length(which(ListFathers[i,2:run]==as.character(dyads$bestFather[i])))/(run-1)
  }
  
  # true dad ?
  # vrai pere est trouve si le meilleur pere (le plus frequent) correspond au vrai pere
  dyads$trueDad=as.numeric(dyads$realFather==dyads$bestFather)
  
  
  # colonne meanProbObs
  # moyenne des probability de COLONY pour le MEILLEUR PERE
  for (i in 1:nrow(dyads)) {
    dyads$meanProbObs[i]=mean(as.numeric(probDad[i,which(ListFathers[i,2:run]==dyads$bestFather[i])+1]),na.rm=T)
  }
  # colonne varProbObs
  # mediane des probability de COLONY pour le MEILLEUR PERE
  for (i in 1:nrow(dyads)) {
    dyads$medianProbObs[i]=median(as.numeric(probDad[i,which(ListFathers[i,2:run]==dyads$bestFather[i])+1]),na.rm=T)
  }  # colonne varProbObs
  # variance des probability de COLONY pour le MEILLEUR PERE
  for (i in 1:nrow(dyads)) {
    dyads$varProbObs[i]=var(as.numeric(probDad[i,which(ListFathers[i,2:run]==dyads$bestFather[i])+1]),na.rm=T)
  } 
  
  # FIN
  return(dyads)
}

#----
# proportions d'assignations pour paternite
#----
# paternity=read.table("simReplicatedRuns48/simReplicatedRuns48C_1.Paternity",h=T,colClasses = c("character","character","numeric"))
# summary(paternity)
# str(paternity)
# 
# for (i in 1:nrow(paternity)) {
#   # pere correspond au juvenile & proba > 0.8 (recommandation Walling 2010)
#   paternity$Match[i]=(grepl(paste("^",paternity$InferredDad1[i],sep=""),paternity$OffspringID[i]) & paternity$ProbDad1[i] > 0.8)
# }
# 
# # liste de peres correctement assignes/nombre total de peres
# match=unique(paternity$InferredDad1[which(paternity$Match==TRUE)])
# 
# # liste de peres incorrectement assignes/nombre total de peres
# mismatch=unique(paternity$InferredDad1[which(paternity$Match==FALSE)])
# 
# # prop of peres correctement assignes pour toutes les inferences
# # (confiance la plus conservative, seuil le plus strict et proportions justes i.e. = 1)
# # peres correctement assignes qui ne sont pas dans la liste des inferences incorrectes
# match=match[!(match %in% mismatch)]
# length(match)/numberOfDads
# 
# # prop of peres incorrectement assignes pour au moins une inference
# length(mismatch)/numberOfDads
# 
# # prop of not assigned
# # nombre de peres non assignes/nombre total de peres
# # i.e. total-(match + mismatch)
# (numberOfDads-length(match)-length(mismatch))/numberOfDads
# 
# 
# # population level confidence
# # selon la methode decrite par Walling 2010
# # le seuil de confiance au niveau de la population est defini comme
# # la moyenne des confiances individuelles pour les individus assignes (assignation acceptee par Colony)
# mean(paternity$ProbDad1[which(paternity$ProbDad1 > 0.8)],na.rm=T)
# 
# 
# # moyenne sur run 1-50
# #----
# library(data.table)
# # moyenne des proportions et la population-level confidence sur l'ensemble des 50 runs
# result50=data.frame(match=rep(NA,50),mismatch=rep(NA,50),notassigned=rep(NA,50),mean=rep(NA,50))
# for (j in 1:50) {
#   print(j)
#   # proportions d'assignations pour paternite
#   if (max(count.fields(paste("simReplicatedRuns48C_",j,".Paternity",sep="")))==3) {
#     paternity=read.table(paste("simReplicatedRuns48C_",j,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
#   } else {
#     paternity=fread(paste("simReplicatedRuns48C_",j,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
#     paternity=paternity[,1:3]
#   }
#   for (i in 1:nrow(paternity)) {
#     # pere correspond au juvenile & proba > 0.8 (recommandation Walling 2010)
#     paternity$Match[i]=(grepl(paste("^",paternity$InferredDad1[i],sep=""),paternity$OffspringID[i]) & paternity$ProbDad1[i] > 0.8)
#   }
#   
#   # liste de peres correctement assignes/nombre total de peres
#   match=unique(paternity$InferredDad1[which(paternity$Match==TRUE)])
#   
#   # liste de peres incorrectement assignes/nombre total de peres
#   mismatch=unique(paternity$InferredDad1[which(paternity$Match==FALSE)])
#   
#   # prop of peres correctement assignes pour toutes les inferences
#   # (confiance la plus conservative, seuil le plus strict et proportions justes i.e. = 1)
#   # peres correctement assignes qui ne sont pas dans la liste des inferences incorrectes
#   match=match[!(match %in% mismatch)]
#   result50$match[j]=length(match)/numberOfDads
#   
#   # prop of peres incorrectement assignes pour au moins une inference
#   result50$mismatch[j]=length(mismatch)/numberOfDads
#   
#   # prop of not assigned
#   # nombre de peres non assignes/nombre total de peres
#   # i.e. total-(match + mismatch)
#   result50$notassigned[j]=(numberOfDads-length(match)-length(mismatch))/numberOfDads
#   
#   
#   # population level confidence
#   # selon la methode decrite par Walling 2010
#   # le seuil de confiance au niveau de la population est defini comme
#   # la moyenne des confiances individuelles pour les individus assignes (assignation acceptee par Colony)
#   result50$mean[j]=mean(paternity$ProbDad1[which(paternity$ProbDad1 > 0.8)],na.rm=T)
# }
# result50
# apply(result50,2,mean)
# 
# 
# # 50 runs pooles ensemble
# #----
# listOfFathers=c()
# confidence=c()
# nbmatch=c()
# for (j in 1:50){
#   if (max(count.fields(paste("simReplicatedRuns48C_",j,".Paternity",sep="")))==3) {
#     pat=read.table(paste("simReplicatedRuns48C_",j,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
#   } else {
#     pat=fread(paste("simReplicatedRuns48C_",j,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
#     pat=pat[,1:3]
#   }
#   for (i in 1:nrow(pat)) {
#     # pere correspond au juvenile & proba > 0.8 (recommandation Walling 2010)
#     pat$Match[i]=(grepl(paste("^",pat$InferredDad1[i],sep=""),pat$OffspringID[i]) & pat$ProbDad1[i] > 0.8)
#   }
#   match=unique(pat$InferredDad1[which(pat$Match==TRUE)])
#   mismatch=unique(pat$InferredDad1[which(pat$Match==FALSE)])
#   match=match[!(match %in% mismatch)]
#   print(match)
#   listOfFathers=c(listOfFathers,match)
#   confidence=c(confidence,pat$ProbDad1[which(pat$ProbDad1 > 0.8)])
#   nbmatch=c(nbmatch,length(unique(match)))
# }
# nbmatch
# mean(nbmatch)
# sd(nbmatch)
# listOfFathers
# unique(listOfFathers)
# mean(confidence)
# 
# 
# 



#----
#----
# CONSTRUCTION DE LA DISTRIBUTION DE LA STAT ~ NOMBRE DE RUNS
#----

# en arguments :
# 1/ path, le chemin du dossier d'outputs
# 2/ nombre de runs
# 3/ stat : "mean","median","var","bestFather"
# 4/ method : "exhaustive" (tous les points) ou "sample" (tirage aleatoire de 50 points par run)
# 5/ bestF : TRUE -> calcul de la stat sur le meilleur pere pour un offspring parmi tous les runs (freq la plus grande)
#            FALSE -> calcul de la stat sur tous les peres pour un offspring

makeDistributionStat=function(path="",run=50,stat="mean",method="exhaustive",bestF=TRUE,TD=TRUE){
  library(data.table)
  offsprings=read.table(paste(path,"/offsprings.txt",sep=""),sep=",",h=F)
  # ProbDad : liste des probability de Colony a chaque run 1-50
  probDad=as.data.frame(offsprings[,1])
  colnames(probDad)="OffspringID"
  # listDad : liste des peres en miroir a probDad
  listDad=as.data.frame(offsprings[,1])
  colnames(listDad)="OffspringID"
  for (i in 2:51) {
    probDad[,i]=NA
  }
  for (i in 1:run) {
    # proportions d'assignations pour paternite
    if (max(count.fields(paste(path,"/",path,"C_",i,".Paternity",sep="")))==3) {
      pat=read.table(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      pat=fread(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      pat=pat[,1:3]
    }
    
    for (j in 1:nrow(probDad)) {
      if(probDad$OffspringID[j] %in% pat$OffspringID){
        probDad[j,i+1]=pat$ProbDad1[which(pat$OffspringID==as.character(probDad$OffspringID[j]))]
        listDad[j,i+1]=pat$InferredDad1[which(pat$OffspringID==as.character(listDad$OffspringID[j]))]
      } else {
        probDad[j,i+1]=NA
        listDad[j,i+1]=NA
      }
    }
  }
  
  
  # ListFathers : liste des meilleurs peres assignes a chaque run 1-50
  ListFathers=as.data.frame(offsprings[,1])
  colnames(ListFathers)="OffspringID"
  for (i in 1:(run-1)) {
    # proportions d'assignations pour paternite
    if (max(count.fields(paste(path,"/",path,"C_",i,".Paternity",sep="")))==3) {
      pat=read.table(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      pat=fread(paste(path,"/",path,"C_",i,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      pat=pat[,1:3]
    }
    
    for (j in 1:nrow(ListFathers)) {
      if(ListFathers$OffspringID[j] %in% pat$OffspringID){
        ListFathers[j,i+1]=pat$InferredDad1[which(pat$OffspringID==as.character(probDad$OffspringID[j]))]
      } else {
        ListFathers[j,i+1]=NA
      }
    }
  }
  
  # nombre d'occurences du meilleur pere
  # identifier le meilleur pere : celui qui revient le plus souvent dans les assignations
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre dyads et ListFathers
  bestFather=c()
  for (i in 1:nrow(ListFathers)){
    uniq=unique(na.omit(as.character(ListFathers[i,2:run])))
    
    if(length(uniq)==0){
      bestFather[i]=NA
    }else{
      best=c(NA,0)
      for (j in 1:length(uniq)) {
        occur=length(which(ListFathers[i,2:run]==uniq[j]))
        if(occur>best[2]){
          best[1]=uniq[j]
          best[2]=occur
        }
      }
      bestFather[i]=best[1]
    }
  }
  # compter le nombre d'occurences du meilleur pere/nombre de runs
  probBestFather=c()
  for(i in 1:length(bestFather)){
    probBestFather[i]=length(which(ListFathers[i,2:run]==as.character(bestFather[i])))/(run-1)
  }
  
  # true dad ?
  realFather=c()
  for (i in 1:nrow(ListFathers)) {
    # vrai pere
    realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",ListFathers$OffspringID[i])
  }
  # vrai pere est trouve si le meilleur pere (le plus frequent) correspond au vrai pere
  trueDad=(realFather==bestFather)
  
  
  # enlever les valeurs de probability qui ne correspondent pas au meilleur pere
  if(bestF==TRUE){
    # exclure les peres qui ne sont pas le meilleur pere -> valeur NA
    # on compare la liste bestFather a chaque run de chaque ligne de listDad
    # listDad est un data.frame miroir de probDad
    for (i in nrow(probDad)) {
      probDad[i,(which(listDad[i,2:run+1]!=bestFather[i]))+1]=NA
    }
  }else{
    if(bestF!=FALSE){
      warning("Wrong choice for bestFather parameter !")
    }
  }
  
  if(method=="exhaustive"){
    # Individus independants, calcul de la valeur pour chaque individu
    
    # creation d'un tableau de n*m lignes (n individus * m runs) et 2 colonnes
    dist=matrix(nrow=nrow(probDad),ncol=run)
    colnames(dist)=seq(1,run,1)
    if(stat=="mean"){
      # pour chaque individu i, la colonne j de probDad est la valeur de la stat associee a probability
      # sur un tirage sans remise de m runs parmi la ligne i
      for (i in 1:nrow(dist)) {
        for (j in 1:run) {
          if(is.na(sum(probDad[i,2:(run+1)]))){
            dist[i,j]=NA
          }else{
            dist[i,j]=mean(as.numeric(sample(na.omit(probDad[i,2:(run+1)]),size=j,replace=FALSE)),na.rm=TRUE)
          }
        }
      }
    }else{
      if(stat=="median"){
        # pour chaque individu i, la colonne j de probDad est la valeur de la stat associee a probability
        # sur un tirage sans remise de m runs parmi la ligne i
        for (i in 1:nrow(dist)) {
          for (j in 1:run) {
            if(is.na(sum(probDad[i,2:(run+1)]))){
              dist[i,j]=NA
            }else{
              dist[i,j]=median(as.numeric(sample(na.omit(probDad[i,2:(run+1)]),size=j,replace=FALSE)),na.rm=TRUE)
            }
          }
        }
      }else{
        if(stat=="var"){
          # pour chaque individu i, la colonne j de probDad est la valeur de la stat associee a probability
          # sur un tirage sans remise de m runs parmi la ligne i
          for (i in 1:nrow(dist)) {
            for (j in 1:run) {
              if(is.na(sum(probDad[i,2:(run+1)]))){
                dist[i,j]=NA
              }else{
                dist[i,j]=var(as.numeric(sample(na.omit(probDad[i,2:(run+1)]),size=j,replace=FALSE)),na.rm=TRUE)
              }
            }
          }
        }else{
          warning("Wrong choice of stat !")
        }
      }
    }
    
    # SI PARAM TD==TRUE
    # Affichage separe de la distribution de la stat pour vrai pere et mauvais pere
    if(TD==TRUE){
      # la liste des vrais peres existe deja dans realFather
      # comparer avec liste des meilleurs peres et former une liste de tri : TRUE/FALSE
      # c'est la liste trueDad
      vecTrue=rep(trueDad,run)
      # afficher TRUE en rouge (vrai pere) et FALSE en bleu (mauvais pere)
      vecDist=c()
      vecRuns=c()
      for (l in 1:run) {
        vecDist=c(vecDist,dist[,l])
        vecRuns=c(vecRuns,rep(l,nrow(dist)))
      }
      plot(x=vecRuns[which(vecTrue==TRUE)],y=vecDist[which(vecTrue==TRUE)],xlab="Nombre de runs",ylab="Metrique (exhaustive)",ylim=c(0,1),xlim=c(1,run),
           main=paste("Distribution de la statistique ",stat,"\nen fonction du nombre de runs\npour les vrais peres",sep=""),
           col="red")
      plot(x=vecRuns[which(vecTrue==FALSE)],y=vecDist[which(vecTrue==FALSE)],xlab="Nombre de runs",ylab="Metrique (exhaustive)",ylim=c(0,1),xlim=c(1,run),
           main=paste("Distribution de la statistique ",stat,"\nen fonction du nombre de runs\npour les mauvais peres",sep=""),
           col="blue")
      if(stat=="mean" | stat=="median"){
        abline(h=0.95,col="red")
      }else{
        abline(h=0.05,col="red")
      }
    }else{
      if(TD==FALSE){
        vecDist=c()
        vecRuns=c()
        for (l in 1:run) {
          vecDist=c(vecDist,dist[,l])
          vecRuns=c(vecRuns,rep(l,nrow(dist)))
        }
        plot(x=vecRuns,y=vecDist,xlab="Nombre de runs",ylab="Metrique (exhaustive)",ylim=c(0,1),xlim=c(1,run),
             main=paste("Distribution de la statistique ",stat,"\nen fonction du nombre de runs ",sep=""))
      }else{
        warning("Wrong choice for parameter TD (affichage separe des distribution vrai pere/mauvais pere) : TRUE/FALSE only !")
      }
    }
    
  }else{
    if(method=="sample"){
      # Les individus ne sont pas independants dans cette methode !!!
      # Tirage avec remise de j juveniles (j = nombre de juveniles)
      #     Parmi cette sous pop, tirage de r valeurs de probability par juvenile
      #     et calcul de la valeur de stat sur ces r valeurs (r =  nb de runs)
      if(TD==TRUE){
        # separer le jeu de donnees en vrai pere/faux pere
        juvTrueDad=probDad[which(realFather==bestFather),]
        juvFalseDad=probDad[which(realFather!=bestFather & !is.na(bestFather)),]
        
        j1=nrow(juvTrueDad)
        j2=nrow(juvFalseDad)
        
        runseq=c(1,seq(5,run,5))
        dist1=data.frame(runs=NA,stats=NA)
        dist2=data.frame(runs=NA,stats=NA)
        
        # Construction de la distribution pour les vrais peres
        # pour chaque nombre de runs 1~50
        for (r in runseq) {
          # bootstrap de 100 valeurs pour chaque nombre de runs
          for (i in 1:100) {
            # tirage de j juveniles
            juv=juvTrueDad[sample(seq(1,j1,1),size=j1,replace=TRUE),]
            for (g in 1:j1) {
              # tirage de r valeurs par juvenile et calcul de la stat associee
              if(stat=="mean"){
                statValues=mean(na.omit(as.vector(t(juv[,2:(run+1)]))))
                lim=c(0.8,1)
              }else{
                if(stat=="median"){
                  statValues=median(na.omit(as.vector(t(juv[,2:(run+1)]))))
                  lim=c(0.8,1)
                }else{
                  if(stat=="var"){
                    statValues=var(na.omit(as.vector(t(juv[,2:(run+1)]))))
                    lim=c(0,0.2)
                  }else{
                    warning("Wrong choice of stat !")
                  }
                }
              }
              
            }
            dist1=rbind(dist1,c(r,statValues))
          }
        }
        plot(dist1$runs,dist1$stat,xlab="Nombre de runs",ylab="Metrique (100 valeurs par reechantillonnage)",
             ylim=lim,col="red")
        
        
        # Construction de la distribution pour les faux peres
        # pour chaque nombre de runs 1~50
        for (r in runseq) {
          # bootstrap de 100 valeurs pour chaque nombre de runs
          for (i in 1:100) {
            # tirage de j juveniles
            juv=juvFalseDad[sample(seq(1,j2,1),size=j2,replace=TRUE),]
            for (g in 1:j2) {
              # tirage de r valeurs par juvenile et calcul de la stat associee
              if(stat=="mean"){
                statValues=mean(na.omit(as.vector(t(juv[,2:(run+1)]))))
                lim=c(0.8,1)
              }else{
                if(stat=="median"){
                  statValues=median(na.omit(as.vector(t(juv[,2:(run+1)]))))
                  lim=c(0.8,1)
                }else{
                  if(stat=="var"){
                    statValues=var(na.omit(as.vector(t(juv[,2:(run+1)]))))
                    lim=c(0,0.2)
                  }else{
                    warning("Wrong choice of stat !")
                  }
                }
              }
              
            }
            dist2=rbind(dist2,c(r,statValues))
          }
        }
        plot(dist2$runs,dist2$stat,xlab="Nombre de runs",ylab="Metrique (100 valeurs par reechantillonnage)",
             ylim=lim,col="blue")
        
        # Distribution a exporter
        dist=rbind(dist1,dist2)
        dist$trueDad=c(rep(TRUE,nrow(dist1)),rep(FALSE,nrow(dist2)))
        
      } else{
        if(TD==FALSE){
          j=nrow(probDad)
          runseq=c(1,seq(5,run,5))
          dist=data.frame(runs=NA,stats=NA)
          # pour chaque nombre de runs 1~50
          for (r in runseq) {
            # bootstrap de 100 valeurs pour chaque nombre de runs
            for (i in 1:100) {
              # tirage de j juveniles
              juv=probDad[sample(seq(1,j,1),size=j,replace=TRUE),]
              for (g in 1:j) {
                # tirage de r valeurs par juvenile et calcul de la stat associee
                if(stat=="mean"){
                  statValues=mean(na.omit(as.vector(t(juv[,2:(run+1)]))))
                  lim=c(0.8,1)
                }else{
                  if(stat=="median"){
                    statValues=median(na.omit(as.vector(t(juv[,2:(run+1)]))))
                    lim=c(0.8,1)
                  }else{
                    if(stat=="var"){
                      statValues=var(na.omit(as.vector(t(juv[,2:(run+1)]))))
                      lim=c(0,0.2)
                    }else{
                      warning("Wrong choice of stat !")
                    }
                  }
                }
                
              }
              dist=rbind(dist,c(r,statValues))
            }
          }
          plot(dist$runs,dist$stat,xlab="Nombre de runs",ylab="Metrique (100 valeurs par reechantillonnage)",
               ylim=lim)
        }else{
          warning("Wrong choice for parameter TD (affichage separe des distribution vrai pere/mauvais pere) : TRUE/FALSE only !")
        }
      }
      
      
      
      
      
      ############# OBSOLETE ###################
      # # un tableau renvoyé par 
      # # selection de 50 echantillons pour chaque nombre de runs
      # # abscisse : nombre de runs par pas de 5
      # # ordonnee : distribution de la stat de test (nuage de points)
      # # pour le meilleur pere selon Colony (best father le plus fréquent sur n runs)
      # 
      # # creation d'un tableau de 50 lignes (50 individus) et n/5 colonnes (n runs par pas de 5 runs)
      # dist=matrix(nrow=50,ncol=(run/5)+1)
      # colnames(dist)=c(1,seq(5,run,5))
      # # pour chaque valeur de nombre de runs (1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50) (chaque colonne)
      # #     tirage aleatoire de 50 individus (chaque ligne) dans probDad, la liste de toutes les probability
      # #     (condition 1 : tirage parmi les individus qui ne sont pas NA)
      # #         parmi ces 50 individus tires, on calcule la stat de test sur n valeurs (non NA) tirees au hasard
      # nbRuns=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
      # if(stat=="mean"){
      #   for (r in 1:(round(run/5)+1)) {
      #     index=sample(x=which(rowSums(probDad[,2:51],na.rm=T)>0),size=50,replace=T) # les indices des individus tires
      #     for (i in 1:length(index)) {
      #       dist[i,r]=mean(as.numeric(probDad[index[i],sample(x=(which(probDad[index[i],2:51]>0)+1),size=nbRuns[r],replace=T)]),na.rm=T)
      #     }
      #   }
      # }else{
      #   if(stat=="median"){
      #     for (r in 1:(round(run/5)+1)) {
      #       index=sample(x=which(rowSums(probDad[,2:51],na.rm=T)>0),size=50,replace=T) # les indices des individus tires
      #       for (i in 1:length(index)) {
      #         dist[i,r]=median(as.numeric(probDad[index[i],sample(x=(which(probDad[index[i],2:51]>0)+1),size=nbRuns[r],replace=T)]),na.rm=T)
      #       }
      #     }
      #   }else{
      #     if(stat=="var"){
      #       for (r in 1:(round(run/5)+1)) {
      #         index=sample(x=which(rowSums(probDad[,2:51],na.rm=T)>0),size=50,replace=T) # les indices des individus tires
      #         for (i in 1:length(index)) {
      #           dist[i,r]=var(as.numeric(probDad[index[i],sample(x=(which(probDad[index[i],2:51]>0)+1),size=nbRuns[r],replace=T)]),na.rm=T)
      #         }
      #       }
      #     }else{
      #       
      #     }
      #   }
      # }
      # graphDist=data.frame(x=c(rep(1,50),rep(5,50),rep(10,50),rep(15,50),rep(20,50),rep(25,50),rep(30,50),rep(35,50),
      #                          rep(40,50),rep(45,50),rep(50,50)),
      #                      y=c(dist[,1],dist[,2],dist[,3],dist[,4],dist[,5],dist[,6],dist[,7],dist[,8],dist[,9],
      #                          dist[,10],dist[,11]))
      # plot(graphDist,xlab="Nombre de runs (par pas de 5)",ylab="Metrique (50 valeurs)",ylim=c(0,1),xlim=c(0,50),
      #      main="Distribution de la statistique\nen fonction du nombre de runs ")
    }else{
      warning("Wrong choice of method !")
    }
  }
  
  
  return(dist)
}


#----
#----
# OCCURENCE BEST FATHER ~ NOMBRE DE RUNS
#----
# Compter le nombre d'occurences du meilleur pere, en fonction du nombre de runs
# Les runs sont independants, ce sont les replicats
# Produit un tableau de r colonnes (nb de runs) et n lignes (nb de replicats)
# Produire un graph
# x : nombre de runs
# y : nombre de runs ou le meiller pere est retrouve
# deux series de mesures : vrai pere/mauvais pere
# abline : expected x=y
# Les NA sont consideres comme des 0 (pere incorrect est equivalent a pas de pere trouve)

occurenceBestF=function(path="",run=50){
  timeIn=Sys.time()
  library(data.table)
  cat("Assemblage des donnees.\n")
  offsprings=read.table(paste(path,"/offsprings.txt",sep=""),sep=" ",h=F)
  # ProbDad : liste des probability de Colony a chaque run 1-50
  probDad=as.data.frame(offsprings[,1])
  colnames(probDad)="OffspringID"
  # listDad : liste des peres en miroir a probDad
  listDad=as.data.frame(offsprings[,1])
  colnames(listDad)="OffspringID"
  for (i in 2:51) {
    probDad[,i]=NA
  }
  for (i in 1:run) {
    # proportions d'assignations pour paternite
    if (max(count.fields(paste(path,"/",path,"_",i,".Paternity",sep="")))==3) {
      pat=read.table(paste(path,"/",path,"_",i,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      pat=fread(paste(path,"/",path,"_",i,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      pat=pat[,1:3]
    }
    
    for (j in 1:nrow(probDad)) {
      if(probDad$OffspringID[j] %in% pat$OffspringID){
        probDad[j,i+1]=pat$ProbDad1[which(pat$OffspringID==as.character(probDad$OffspringID[j]))]
        listDad[j,i+1]=pat$InferredDad1[which(pat$OffspringID==as.character(listDad$OffspringID[j]))]
      } else {
        probDad[j,i+1]=NA
        listDad[j,i+1]=NA
      }
    }
  }
  
  # identifier le meilleur pere : celui qui revient le plus souvent dans les assignations
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre bestFather et listDad
  cat("Identification du meilleur pere.\n")
  bestFather=c()
  for (i in 1:nrow(listDad)){
    uniq=unique(na.omit(as.character(listDad[i,2:run])))
    
    if(length(uniq)==0){
      bestFather[i]=NA
    }else{
      best=c(NA,0)
      for (j in 1:length(uniq)) {
        occur=length(which(listDad[i,2:run]==uniq[j]))
        if(occur>best[2]){
          best[1]=uniq[j]
          best[2]=occur
        }
      }
      bestFather[i]=best[1]
    }
  }
  
  # true dad ?
  # Liste des lignes qui contiennent le vrai pere (TRUE) ou un mauvais pere (FALSE) pour trier listDad
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre realFather et listDad
  realFather=c()
  for (i in 1:nrow(listDad)) {
    # vrai pere
    realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",listDad$OffspringID[i])
  }
  
  cat("Contruction du tableau de comptage.\n")
  # Pour chaque nombre de runs, on fait un tableau de presence/absence du meilleur pere trouve
  countDad=data.frame()
  for (n in 1:run){
    for (i in 1:nrow(listDad)){
      if (is.na(listDad[i,(n+1)]) | is.na(bestFather[i])){
        countDad[i,n]=0
      }else{
        if (bestFather[i]==listDad[i,(n+1)]) {
          countDad[i,n]=1
        }else{
          countDad[i,n]=0
        }
      }
    }
  }
  # vrai pere est trouve si le meilleur pere (le plus frequent) correspond au vrai pere
  # partage du tableau de comptage en deux, selon categorie
  trueDad=countDad[which(realFather==bestFather),]
  falseDad=countDad[which(realFather!=bestFather),]
  
  
  # STRUCTURE : construction du tableau
  # fournir un tableau pour vrai pere trouve et un autre pour mauvais pere trouve
  # pour chaque nombre de runs (r), on construit un bootstrap de 100 valeurs
  #     on sample n runs avec remise, 1OO fois
  #           sur chaque sample (donc 1:100), on fait un bootstrap de la valeur pour
  #           chaque condition : vrai pere trouve/mauvais pere trouve
  cat("Reechantillonnage du nombre d'occurences du meilleur pere.\n")
  require("boot")
  meanfun=function(data, j){
    d=data[j]
    return(mean(d,na.rm=TRUE))   
  }
  # Sur les vrais peres
  cat("Peres corrects :\n")
  boT=data.frame()
  progress=txtProgressBar(min = 1, max = run, style = 3)
  for (r in 1:run) {
    for (i in 1:100) {
      index=sample(seq(1,run,1),size=r,replace=TRUE) # reechantillonnage avec remise des indices des runs
      # bootstrapping : nombre moyen d'occurence
      # moyenne des 0/1 sur toutes les valeurs des runs tires (index)
      # stocke dans un tableau ou x=nb de runs=r et y=1:100=i (replicats)
      boT[i,r]=boot(as.vector(t(trueDad[,index])), statistic=meanfun, R=1000)$t0
    }
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(progress,r)
  }
  close(progress)
  # ET sur les peres incorrects
  cat("Peres incorrects :\n")
  boF=data.frame()
  progress=txtProgressBar(min = 1, max = run, style = 3)
  for (r in 1:run) {
    for (i in 1:100) {
      index=sample(seq(1,run,1),size=r,replace=TRUE) # reechantillonnage avec remise des indices des runs
      # bootstrapping : nombre moyen d'occurence
      # moyenne des 0/1 sur toutes les valeurs des runs tires (index)
      # stocke dans un tableau ou x=nb de runs=r et y=1:100=i (replicats)
      boF[i,r]=boot(as.vector(t(falseDad[,index])), statistic=meanfun, R=1000)$t0
    }
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(progress,r)
  }
  close(progress)
  
  cat("Fonction graphique.\n")
  # produit 2 df avec le nombre d'occurence des vrais peres et des mauvais peres
  # avec x=nombre de runs et y=proportion de runs ou meilleur pere est retrouve
  dfT=data.frame(runs=NA,occurence=NA)
  for (i in 1:nrow(boT)) {
    for (n in 1:run) {
      dfT=rbind(dfT,c(n,boT[i,n]))
    }
  }
  dfT=dfT[-1,]
  dfF=data.frame(runs=NA,occurence=NA)
  for (i in 1:nrow(boF)) {
    for (n in 1:run) {
      dfF=rbind(dfF,c(n,boF[i,n]))
    }
  }
  dfF=dfF[-1,]
  
  # graph avec ggplot2
  library(ggplot2)
  print(ggplot(data=dfT, aes(x=runs, y=occurence)) +
          geom_point(aes(color="Correct father"))+
          geom_point(data=dfF,aes(color="Incorrect father")) +
          labs(color="Distribution") +
          scale_colour_manual(values=c("black","darkgrey"))+
          geom_hline(yintercept=0.3,linetype="dashed",size=1) +
          ggtitle("Assignation frequencies of the best father\nfor correct or incorrect fathers",
                  subtitle=paste("n=",nrow(offsprings),sep = "")) +
          xlab("Number of runs") + ylab("Assignation frequency of the best father") +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
                plot.subtitle = element_text(color="black",size=30,hjust = 0.5),
                axis.title.x = element_text(color="black", size=30),
                axis.title.y = element_text(color="black", size=30),
                axis.text=element_text(size=26, colour="black"),
                legend.key = element_rect(fill = "white", size = 1),
                legend.key.height = unit(2,"line"),
                legend.key.width = unit(5,"line"),
                legend.text=element_text(size=30),
                legend.title=element_text(size=30)))
  
  # renvoie les deux distributions concatenees en un seul df
  dist=rbind(dfT,dfF)
  dist$father=c(rep(TRUE,nrow(dfT)),rep(FALSE,nrow(dfF)))
  
  timeOut=Sys.time()
  cat("Execution en",difftime(timeOut,timeIn,units="mins"),"minutes.\n")
  return(dist)
}


#----
# OCCURENCE BEST MOTHER ~ NOMBRE DE RUNS
#----
# Compter le nombre d'occurences de la meilleure mere, en fonction du nombre de runs
# Les runs sont independants, ce sont les replicats
# Produit un tableau de r colonnes (nb de runs) et n lignes (nb de replicats)
# Produire un graph
# x : nombre de runs
# y : nombre de runs ou la meilleure mere est retrouvee
# deux series de mesures : vraie mere/mauvais mere
# abline : expected x=y
# Les NA sont consideres comme des 0 (mere incorrecte est equivalent a pas de mere trouvee)

##### Quick adaptation : variables names are still masculine.

occurenceBestM=function(path="",run=50){
  timeIn=Sys.time()
  library(data.table)
  cat("Assemblage des donnees.\n")
  offsprings=read.table(paste(path,"/offsprings.txt",sep=""),sep=" ",h=F)
  # ProbDad : liste des probability de Colony a chaque run 1-50
  probDad=as.data.frame(offsprings[,1])
  colnames(probDad)="OffspringID"
  # listDad : liste des peres en miroir a probDad
  listDad=as.data.frame(offsprings[,1])
  colnames(listDad)="OffspringID"
  for (i in 2:51) {
    probDad[,i]=NA
  }
  for (i in 1:run) {
    # proportions d'assignations pour paternite
    if (max(count.fields(paste(path,"/",path,"_",i,".Maternity",sep="")))==3) {
      pat=read.table(paste(path,"/",path,"_",i,".Maternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      pat=fread(paste(path,"/",path,"_",i,".Maternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      pat=pat[,1:3]
    }
    
    for (j in 1:nrow(probDad)) {
      if(probDad$OffspringID[j] %in% pat$OffspringID){
        probDad[j,i+1]=pat$ProbMum1[which(pat$OffspringID==as.character(probDad$OffspringID[j]))]
        listDad[j,i+1]=pat$InferredMum1[which(pat$OffspringID==as.character(listDad$OffspringID[j]))]
      } else {
        probDad[j,i+1]=NA
        listDad[j,i+1]=NA
      }
    }
  }
  
  # identifier le meilleur pere : celui qui revient le plus souvent dans les assignations
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre bestFather et listDad
  cat("Identification de la meilleure mere.\n")
  bestFather=c()
  for (i in 1:nrow(listDad)){
    uniq=unique(na.omit(as.character(listDad[i,2:run])))
    
    if(length(uniq)==0){
      bestFather[i]=NA
    }else{
      best=c(NA,0)
      for (j in 1:length(uniq)) {
        occur=length(which(listDad[i,2:run]==uniq[j]))
        if(occur>best[2]){
          best[1]=uniq[j]
          best[2]=occur
        }
      }
      bestFather[i]=best[1]
    }
  }
  
  # true dad ?
  # Liste des lignes qui contiennent le vrai pere (TRUE) ou un mauvais pere (FALSE) pour trier listDad
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre realFather et listDad
  realFather=c()
  for (i in 1:nrow(listDad)) {
    # vraie mere
    tmp=sub(pattern = "M[0-9]*","",listDad$OffspringID[i])
    realFather[i]=sub(pattern="C+[0-9]+","",tmp)
  }
  
  cat("Contruction du tableau de comptage.\n")
  # Pour chaque nombre de runs, on fait un tableau de presence/absence du meilleur pere trouve
  countDad=data.frame()
  for (n in 1:run){
    for (i in 1:nrow(listDad)){
      if (is.na(listDad[i,(n+1)]) | is.na(bestFather[i])){
        countDad[i,n]=0
      }else{
        if (bestFather[i]==listDad[i,(n+1)]) {
          countDad[i,n]=1
        }else{
          countDad[i,n]=0
        }
      }
    }
  }
  # vrai pere est trouve si le meilleur pere (le plus frequent) correspond au vrai pere
  # partage du tableau de comptage en deux, selon categorie
  trueDad=countDad[which(realFather==bestFather),]
  falseDad=countDad[which(realFather!=bestFather),]
  
  
  # STRUCTURE : construction du tableau
  # fournir un tableau pour vrai pere trouve et un autre pour mauvais pere trouve
  # pour chaque nombre de runs (r), on construit un bootstrap de 100 valeurs
  #     on sample n runs avec remise, 1OO fois
  #           sur chaque sample (donc 1:100), on fait un bootstrap de la valeur pour
  #           chaque condition : vrai pere trouve/mauvais pere trouve
  cat("Reechantillonnage du nombre d'occurences du meilleur pere.\n")
  require("boot")
  meanfun=function(data, j){
    d=data[j]
    return(mean(d,na.rm=TRUE))   
  }
  # Sur les vrais peres
  cat("Meres correctes :\n")
  boT=data.frame()
  progress=txtProgressBar(min = 1, max = run, style = 3)
  for (r in 1:run) {
    for (i in 1:100) {
      index=sample(seq(1,run,1),size=r,replace=TRUE) # reechantillonnage avec remise des indices des runs
      # bootstrapping : nombre moyen d'occurence
      # moyenne des 0/1 sur toutes les valeurs des runs tires (index)
      # stocke dans un tableau ou x=nb de runs=r et y=1:100=i (replicats)
      boT[i,r]=boot(as.vector(t(trueDad[,index])), statistic=meanfun, R=1000)$t0
    }
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(progress,r)
  }
  close(progress)
  # ET sur les peres incorrects
  cat("Meres incorrectes :\n")
  boF=data.frame()
  progress=txtProgressBar(min = 1, max = run, style = 3)
  for (r in 1:run) {
    for (i in 1:100) {
      index=sample(seq(1,run,1),size=r,replace=TRUE) # reechantillonnage avec remise des indices des runs
      # bootstrapping : nombre moyen d'occurence
      # moyenne des 0/1 sur toutes les valeurs des runs tires (index)
      # stocke dans un tableau ou x=nb de runs=r et y=1:100=i (replicats)
      boF[i,r]=boot(as.vector(t(falseDad[,index])), statistic=meanfun, R=1000)$t0
    }
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(progress,r)
  }
  close(progress)
  
  cat("Fonction graphique.\n")
  # produit 2 df avec le nombre d'occurence des vrais peres et des mauvais peres
  # avec x=nombre de runs et y=proportion de runs ou meilleur pere est retrouve
  dfT=data.frame(runs=NA,occurence=NA)
  for (i in 1:nrow(boT)) {
    for (n in 1:run) {
      dfT=rbind(dfT,c(n,boT[i,n]))
    }
  }
  dfT=dfT[-1,]
  dfF=data.frame(runs=NA,occurence=NA)
  for (i in 1:nrow(boF)) {
    for (n in 1:run) {
      dfF=rbind(dfF,c(n,boF[i,n]))
    }
  }
  dfF=dfF[-1,]
  
  # graph avec ggplot2
  library(ggplot2)
  print(ggplot(data=dfT, aes(x=runs, y=occurence)) +
          geom_point(aes(color="Correct mother"))+
          geom_point(data=dfF,aes(color="Incorrect mother")) +
          labs(color="Distribution") +
          scale_colour_manual(values=c("black","darkgrey"))+
          geom_hline(yintercept=0.3,linetype="dashed",size=1) +
          ggtitle("Assignation frequencies of the best mother\nfor correct or incorrect fathers",
                  subtitle=paste("n=",nrow(offsprings),sep = "")) +
          xlab("Number of runs") + ylab("Assignation frequency of the best mother") +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
                plot.subtitle = element_text(color="black",size=30,hjust = 0.5),
                axis.title.x = element_text(color="black", size=30),
                axis.title.y = element_text(color="black", size=30),
                axis.text=element_text(size=26, colour="black"),
                legend.key = element_rect(fill = "white", size = 1),
                legend.key.height = unit(2,"line"),
                legend.key.width = unit(5,"line"),
                legend.text=element_text(size=30),
                legend.title=element_text(size=30)))
  
  # renvoie les deux distributions concatenees en un seul df
  dist=rbind(dfT,dfF)
  dist$mother=c(rep(TRUE,nrow(dfT)),rep(FALSE,nrow(dfF)))
  
  timeOut=Sys.time()
  cat("Execution en",difftime(timeOut,timeIn,units="mins"),"minutes.\n")
  return(dist)
}

#----------------------------#
# STAT BEST FATHER ~ NUMBER OF RUNS
#----------------------------#
# calculer la valeur de la stat choisie pour le meilleur pere, en fonction du nombre de runs
# Les runs sont independants, ce sont les replicats
# Produit un tableau de r colonnes (nb de runs) et n lignes (nb de replicats)
# Produire un graph
# x : nombre de runs
# y : nombre de runs ou le meiller pere est retrouve
# deux series de mesures : vrai pere/mauvais pere
# abline : expected x=y
# Les NA ne sont pas pris en compte

# Trois stats : mean, median, var

statBestF=function(path="",run=50,stat="mean"){
  timeIn=Sys.time()
  library(data.table)
  cat("Assemblage des donnees.\n")
  offsprings=read.table(paste(path,"/offsprings.txt",sep=""),sep=" ",h=F)
  # ProbDad : liste des probability de Colony a chaque run 1-50
  probDad=as.data.frame(offsprings[,1])
  colnames(probDad)="OffspringID"
  # listDad : liste des peres en miroir a probDad
  listDad=as.data.frame(offsprings[,1])
  colnames(listDad)="OffspringID"
  for (i in 2:51) {
    probDad[,i]=NA
  }
  for (i in 1:run) {
    # proportions d'assignations pour paternite
    if (max(count.fields(paste(path,"/",path,"_",i,".Paternity",sep="")))==3) {
      pat=read.table(paste(path,"/",path,"_",i,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
    } else {
      pat=fread(paste(path,"/",path,"_",i,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
      pat=pat[,1:3]
    }
    
    for (j in 1:nrow(probDad)) {
      if(probDad$OffspringID[j] %in% pat$OffspringID){
        probDad[j,i+1]=pat$ProbDad1[which(pat$OffspringID==as.character(probDad$OffspringID[j]))]
        listDad[j,i+1]=pat$InferredDad1[which(pat$OffspringID==as.character(listDad$OffspringID[j]))]
      } else {
        probDad[j,i+1]=NA
        listDad[j,i+1]=NA
      }
    }
  }
  
  # identifier le meilleur pere : celui qui revient le plus souvent dans les assignations
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre bestFather et listDad
  cat("Identification du meilleur pere.\n")
  bestFather=c()
  for (i in 1:nrow(listDad)){
    uniq=unique(na.omit(as.character(listDad[i,2:run])))
    
    if(length(uniq)==0){
      bestFather[i]=NA
    }else{
      best=c(NA,0)
      for (j in 1:length(uniq)) {
        occur=length(which(listDad[i,2:run]==uniq[j]))
        if(occur>best[2]){
          best[1]=uniq[j]
          best[2]=occur
        }
      }
      bestFather[i]=best[1]
    }
  }
  
  # true dad ?
  # Liste des lignes qui contiennent le vrai pere (TRUE) ou un mauvais pere (FALSE) pour trier listDad
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre realFather et listDad
  realFather=c()
  for (i in 1:nrow(listDad)) {
    # vrai pere
    realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",listDad$OffspringID[i])
  }
  
  cat("Contruction du tableau de comptage.\n")
  # Pour chaque nombre de runs, on fait un tableau de presence/absence du meilleur pere trouve
  countDad=data.frame()
  for (n in 1:run){
    for (i in 1:nrow(listDad)){
      if (is.na(listDad[i,(n+1)]) | is.na(bestFather[i])){
        countDad[i,n]=0
      }else{
        if (bestFather[i]==listDad[i,(n+1)]) {
          countDad[i,n]=1
        }else{
          countDad[i,n]=0
        }
      }
    }
  }
  # vrai pere est trouve si le meilleur pere (le plus frequent) correspond au vrai pere
  # partage du tableau de comptage en deux, selon categorie
  trueDad=probDad[which(realFather==bestFather),2:(run+1)]
  falseDad=probDad[which(realFather!=bestFather),2:(run+1)]
  
  
  # STRUCTURE : construction du tableau des stats
  # fournir un tableau pour vrai pere trouve et un autre pour mauvais pere trouve
  # pour chaque nombre de runs (r), on construit un bootstrap de 100 valeurs
  #     on sample n runs avec remise, 1OO fois
  #           sur chaque sample (donc 1:100), on fait un bootstrap de la valeur pour
  #           chaque condition : vrai pere trouve/mauvais pere trouve
  cat("Calcul de la stat du meilleur pere par bootstrap.\n")
  require("boot")
  if(stat=="mean"){
    fun=function(data, j){
      d=data[j]
      return(mean(d,na.rm=TRUE))   
    }
  }else{
    if(stat=="median"){
      fun=function(data, j){
        d=data[j]
        return(median(d,na.rm=TRUE))   
      }
    }else{
      if(stat=="var"){
        fun=function(data, j){
          d=data[j]
          return(var(d,na.rm=TRUE))   
        }
      }else{
        warning("Wrong choice for 'stat' parameter.")
      }
    }
  }
  # Sur les vrais peres
  cat("Peres corrects :\n")
  boT=data.frame()
  progress=txtProgressBar(min = 1, max = run, style = 3)
  for (r in 1:run) {
    for (i in 1:100) {
      index=sample(seq(1,run,1),size=r,replace=TRUE) # reechantillonnage avec remise des indices des runs
      # bootstrapping : nombre moyen d'occurence
      # moyenne des 0/1 sur toutes les valeurs des runs tires (index)
      # stocke dans un tableau ou x=nb de runs=r et y=1:100=i (replicats)
      boT[i,r]=boot(as.vector(t(trueDad[,index])), statistic=fun, R=1000)$t0
    }
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(progress,r)
  }
  close(progress)
  # ET sur les peres incorrects
  cat("Peres incorrects :\n")
  boF=data.frame()
  progress=txtProgressBar(min = 1, max = run, style = 3)
  for (r in 1:run) {
    for (i in 1:100) {
      index=sample(seq(1,run,1),size=r,replace=TRUE) # reechantillonnage avec remise des indices des runs
      # bootstrapping : nombre moyen d'occurence
      # moyenne des 0/1 sur toutes les valeurs des runs tires (index)
      # stocke dans un tableau ou x=nb de runs=r et y=1:100=i (replicats)
      boF[i,r]=boot(as.vector(t(falseDad[,index])), statistic=fun, R=1000)$t0
    }
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(progress,r)
  }
  close(progress)
  
  cat("Fonction graphique.\n")
  # produit 2 df avec le nombre d'occurence des vrais peres et des mauvais peres
  # avec x=nombre de runs et y=proportion de runs ou meilleur pere est retrouve
  dfT=data.frame(runs=NA,occurence=NA)
  for (i in 1:nrow(boT)) {
    for (n in 1:run) {
      dfT=rbind(dfT,c(n,boT[i,n]))
    }
  }
  dfT=dfT[-1,]
  dfF=data.frame(runs=NA,occurence=NA)
  for (i in 1:nrow(boF)) {
    for (n in 1:run) {
      dfF=rbind(dfF,c(n,boF[i,n]))
    }
  }
  dfF=dfF[-1,]
  
  # graph avec ggplot2
  library(ggplot2)
  print(ggplot(data=dfT, aes(x=runs, y=occurence)) +
          geom_point(aes(color="Correct father"))+
          geom_point(data=dfF,aes(color="Incorrect father")) +
          labs(color="Distribution") +
          scale_colour_manual(values=c("black","darkgrey"))+
          ggtitle(paste("Bootstrapped value of the",stat,"statistic\nfor the best father\nfor a number of colony runs"),
                  subtitle=paste("n=",nrow(listDad),sep = "")) +
          xlab("Number of runs") + ylab("Bootstrapped value of the statistic") +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(color="black", size=34, face="bold.italic",hjust = 0.5),
                plot.subtitle = element_text(color="black",size=30,hjust = 0.5),
                axis.title.x = element_text(color="black", size=30),
                axis.title.y = element_text(color="black", size=30),
                axis.text=element_text(size=26, colour="black"),
                legend.key = element_rect(fill = "white", size = 1),
                legend.key.height = unit(2,"line"),
                legend.key.width = unit(5,"line"),
                legend.text=element_text(size=30),
                legend.title=element_text(size=30)))
  
  # renvoie les deux distributions concatenees en un seul df
  dist=rbind(dfT,dfF)
  dist$father=c(rep(TRUE,nrow(dfT)),rep(FALSE,nrow(dfF)))
  
  timeOut=Sys.time()
  cat("Execution en",difftime(timeOut,timeIn,units="mins"),"minutes.\n")
  return(dist)
}

#----------------------------#
# SENSITIVITY RESULTS
#----------------------------#
# produit un tableau qui contient les effectifs dans chaque categorie d'assignation
# (true, type I error, type II error, no dad)
# pour chaque simulation, une colonne mean number et une colonne SD
# en argument :
# 1/ path=c(), la liste de tous les dossiers de simulation a integrer au tableau
# 2/ nombre de runs

# comptage en 5 categories des juveniles
# question : pour un jeune, est-ce qu'on retrouve le bon pere ?
# 1/ pere est echantillone
#   1.a/ bon pere (P)
#   1.b/ pas de pere retrouve
#   1.c/ mauvais pere assigne (P)
# 2/ pere non echantillone
#   2.a/ on trouve un pere (P)
#   2.b/ pas de pere retrouve

sensitivityResults=function(path=c(""),run=50){
  # creation du tableau
  # 5 lignes : categories de resultats attendus
  # 2n+2 colonnes ou n est le nombre de paths differents
  tab=matrix(nrow=5,ncol=(2*length(path)+2))
  tab[,1]=c("MsampledMtrue","MsampledMnotfound","MsampledMfalse","MnotsampledMfound","MnotsampledMnotfound")
  tab[,2]=c(0.5,0,0,0,0.5)
  colnames(tab)=c("Result","ExpectedProportion",paste(unlist(lapply(path,function(x) rep(x,2))),c("Mean","SD"),sep=""))
  # parcours de chaque dossier de simulation, 1 par 1
  for (p in 1:length(path)) {
    cat("Processing",path[p],"\n")
    require(data.table)
    # Recuperation de la liste des offsprings
    # generer independamment un fichier texte avec les offsprings
    # d'apres .DAT
    offsprings=read.table(paste(path[p],"/offsprings.txt",sep=""),sep=" ",h=F)
    dyads=as.data.frame(offsprings[,1])
    colnames(dyads)="offspring"
    # generer un fichier txt avec la liste des peres echantillonnes
    # d'apres .DAT
    fathers=read.table(paste(path[p],"/fathers.txt",sep=""),sep=" ",h=F)
    fathers=as.data.frame(fathers[,1])
    
    # tableaux vides
    effectifs=data.frame(c(0,0,0,0,0))
    propObserved=data.frame(c(0,0,0,0,0))
    
    for (r in 1:run) {
      # Recuperation de la liste des paternites
      if (max(count.fields(paste(path[p],"/",path[p],"_",r,".Paternity",sep="")))==3) {
        paternity=read.table(paste(path[p],"/",path[p],"_",r,".Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
      } else {
        paternity=fread(paste(path[p],"/",path[p],"_",r,".Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
        paternity=paternity[,1:3]
      }
      
      # reconstitution des dyads
      # pour chaque juvenile, le meilleur pere trouve est assigne (selon la valeur de probability fournie par Colony)
      # C'est le premier pere indique par Colony
      for (i in 1:nrow(dyads)) {
        # vrai pere
        dyads$realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",dyads$offspring[i])
        # on affecte les meilleurs peres aux juveniles quand ils existent
        if (dyads$offspring[i] %in% paternity$OffspringID) {
          dyads$bestFather[i]=paternity$InferredDad1[which(paternity$OffspringID==dyads$offspring[i])]
        } else {
          dyads$bestFather[i]=NA
        }
        # ajout de l'indicateur echantillonne ou pas
        if (dyads$realFather[i] %in% fathers[,1]) {
          dyads$sampled[i]=1
        } else {
          dyads$sampled[i]=0
        }
      }
      
      # tableau de comptage
      effectifs[,r]=c(sum(dyads$sampled==1 & dyads$realFather==dyads$bestFather,na.rm=T),
                      sum(dyads$sampled==1 & is.na(dyads$bestFather),na.rm=T),
                      sum(dyads$sampled==1 & dyads$realFather!=dyads$bestFather & !is.na(dyads$bestFather),na.rm=T),
                      sum(dyads$sampled==0 & !is.na(dyads$bestFather),na.rm=T),
                      sum(dyads$sampled==0 & is.na(dyads$bestFather),na.rm=T)
      )
      # Proportion de chaque categorie
      # Nombre d'offsprings divise par nombre total d'offsprings ech dans la pop
      propObserved[,r]=effectifs[,r]/sum(effectifs[,r])
      cat(r,"runs :",propObserved[,r],"\n")
    }
    tab[,(p*2+1)]=apply(propObserved,1,mean)
    tab[,(p*2+2)]=apply(propObserved,1,sd)
    cat("Mean =",tab[,(p*2+1)],"and SD =",tab[,(p*2+2)],"\n")
  }
  
  # ecriture du tableau
  write.table(tab,"simSensitivity.txt",row.names=FALSE,col.names=TRUE)
  return(tab)
}

#----------------------------#
# SENSITIVITY RESULTS FOR MOTHERS
#----------------------------
# produit un tableau qui contient les effectifs dans chaque categorie d'assignation
# (true, type I error, type II error, no dad)
# pour chaque simulation, une colonne mean number et une colonne SD
# en argument :
# 1/ path=c(), la liste de tous les dossiers de simulation a integrer au tableau
# 2/ nombre de runs

# comptage en 5 categories des juveniles
# question : pour un jeune, est-ce qu'on retrouve le bon pere ?
# 1/ pere est echantillone
#   1.a/ bon pere (P)
#   1.b/ pas de pere retrouve
#   1.c/ mauvais pere assigne (P)
# 2/ pere non echantillone
#   2.a/ on trouve un pere (P)
#   2.b/ pas de pere retrouve

sensitivityResultsMothers=function(path=c(""),run=50){
  # creation du tableau
  # 5 lignes : categories de resultats attendus
  # 2n+2 colonnes ou n est le nombre de paths differents
  tab=matrix(nrow=5,ncol=(2*length(path)+2))
  tab[,1]=c("MsampledMtrue","MsampledMnotfound","MsampledMfalse","MnotsampledMfound","MnotsampledMnotfound")
  tab[,2]=c(0.9,0,0,0,0.1)
  colnames(tab)=c("Result","ExpectedProportion",paste(unlist(lapply(path,function(x) rep(x,2))),c("Mean","SD"),sep=""))
  # parcours de chaque dossier de simulation, 1 par 1
  for (p in 1:length(path)) {
    cat("Processing",path[p],"\n")
    require(data.table)
    # Recuperation de la liste des offsprings
    # generer independamment un fichier texte avec les offsprings
    # d'apres .DAT
    offsprings=read.table(paste(path[p],"/offsprings.txt",sep=""),sep=" ",h=F)
    dyads=as.data.frame(offsprings[,1])
    colnames(dyads)="offspring"
    # generer un fichier txt avec la liste des peres echantillonnes
    # d'apres .DAT
    fathers=read.table(paste(path[p],"/mothers.txt",sep=""),sep=" ",h=F)
    fathers=as.data.frame(fathers[,1])
    
    # tableaux vides
    effectifs=data.frame(c(0,0,0,0,0))
    propObserved=data.frame(c(0,0,0,0,0))
    
    for (r in 1:run) {
      # Recuperation de la liste des paternites
      if (max(count.fields(paste(path[p],"/",path[p],"_",r,".Maternity",sep="")))==3) {
        paternity=read.table(paste(path[p],"/",path[p],"_",r,".Maternity",sep=""),h=T,colClasses = c("character","character","numeric"))
      } else {
        paternity=fread(paste(path[p],"/",path[p],"_",r,".Maternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
        paternity=paternity[,1:3]
      }
      
      # reconstitution des dyads
      # pour chaque juvenile, le meilleur pere trouve est assigne (selon la valeur de probability fournie par Colony)
      # C'est le premier pere indique par Colony
      for (i in 1:nrow(dyads)) {
        # vrai pere
        # dyads$realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",dyads$offspring[i])
        # vraie mere
        tmp=sub(pattern = "M[0-9]*","",dyads$offspring[i])
        dyads$realFather[i]=sub(pattern="C+[0-9]+","",tmp)
        
        # on affecte les meilleurs peres aux juveniles quand ils existent
        if (dyads$offspring[i] %in% paternity$OffspringID) {
          dyads$bestFather[i]=paternity$InferredMum1[which(paternity$OffspringID==dyads$offspring[i])]
        } else {
          dyads$bestFather[i]=NA
        }
        # ajout de l'indicateur echantillonne ou pas
        if (dyads$realFather[i] %in% fathers[,1]) {
          dyads$sampled[i]=1
        } else {
          dyads$sampled[i]=0
        }
      }
      
      # tableau de comptage
      effectifs[,r]=c(sum(dyads$sampled==1 & dyads$realFather==dyads$bestFather,na.rm=T),
                      sum(dyads$sampled==1 & is.na(dyads$bestFather),na.rm=T),
                      sum(dyads$sampled==1 & dyads$realFather!=dyads$bestFather & !is.na(dyads$bestFather),na.rm=T),
                      sum(dyads$sampled==0 & !is.na(dyads$bestFather),na.rm=T),
                      sum(dyads$sampled==0 & is.na(dyads$bestFather),na.rm=T)
      )
      # Proportion de chaque categorie
      # Nombre d'offsprings divise par nombre total d'offsprings ech dans la pop
      propObserved[,r]=effectifs[,r]/sum(effectifs[,r])
      cat(r,"runs :",propObserved[,r],"\n")
    }
    tab[,(p*2+1)]=apply(propObserved,1,mean)
    tab[,(p*2+2)]=apply(propObserved,1,sd)
    cat("Mean =",tab[,(p*2+1)],"and SD =",tab[,(p*2+2)],"\n")
  }
  
  # ecriture du tableau
  write.table(tab,"simSensitivityMothers.txt",row.names=FALSE,col.names=TRUE)
  return(tab)
}


#-----------------------------#
# SIM RESULTS WITH INFOS - EXPORT DATA
#-----------------------------
# give a file containing :
# - offsrping ID
# - fathers and mothers ID
# - sex
# - correct/incorrect assignation
# - genotype

simResultsWithInfo=function(path=""){
  require(data.table)
  # Recuperation de la liste des offsprings
  # generer un fichier texte avec les offsprings
  offsprings=read.table(paste(path,"/offsprings.txt",sep=""),sep=" ",h=F)
  dyads=as.data.frame(offsprings[,1])
  colnames(dyads)="offspring"
  # Recuperation de la liste des paternites
  if (max(count.fields(paste(path,"/",path,"_1.Paternity",sep="")))==3) {
    paternity=read.table(paste(path,"/",path,"_1.Paternity",sep=""),h=T,colClasses = c("character","character","numeric"))
  } else {
    paternity=fread(paste(path,"/",path,"_1.Paternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
    paternity=paternity[,1:3]
  }
  # Recuperation de la liste des maternites
  if (max(count.fields(paste(path,"/",path,"_1.Maternity",sep="")))==3) {
    maternity=read.table(paste(path,"/",path,"_1.Maternity",sep=""),h=T,colClasses = c("character","character","numeric"))
  } else {
    maternity=fread(paste(path,"/",path,"_1.Maternity",sep=""),h=T,select=c(1:5),fill=TRUE,colClasses = c("character","character","numeric","character","numeric"))
    maternity=maternity[,1:3]
  }
  # For each offspring, mother and father are searched and true/false assignment is checked at the same time
  dyads$father=rep(NA,nrow(dyads))
  dyads$mother=rep(NA,nrow(dyads))
  dyads$assignmentFather=rep(NA,nrow(dyads)) # 0 if false assignment, 1 if true
  dyads$assignmentMother=rep(NA,nrow(dyads)) # 0 if false assignment, 1 if true
  
  # Liste des lignes qui contiennent le vrai pere (TRUE) ou un mauvais pere (FALSE) pour trier listDad
  # !!!! on utilise les numeros de lignes : ils doivent concorder entre realFather et listDad
  realFather=c()
  for (i in 1:nrow(dyads)) {
    # vrai pere
    realFather[i]=sub(pattern = "F[0-9]*C*[0-9]*","",dyads$offspring[i])
  }
  # Same operation for the true mother
  realMother=c()
  for (i in 1:nrow(dyads)) {
    # vraie mere
    tmp=sub(pattern = "M[0-9]*","",dyads$offspring[i])
    realMother[i]=sub(pattern="C+[0-9]+","",tmp)
  }
  for (i in 1:nrow(dyads)){
    if (sum(paternity$OffspringID==dyads$offspring[i])){
      dyads$father[i]=paternity$InferredDad1[which(paternity$OffspringID==dyads$offspring[i])]
    }
    if (sum(maternity$OffspringID==dyads$offspring[i])){
      dyads$mother[i]=maternity$InferredMum1[which(maternity$OffspringID==dyads$offspring[i])]
    }
    if (!is.na(dyads$father[i])){
      if(dyads$father[i]==realFather[i]){
        dyads$assignmentFather[i]=1
      } else {
        dyads$assignmentFather[i]=0
      }
    }
    if (!is.na(dyads$mother[i])){
      if(dyads$mother[i]==realMother[i]){
        dyads$assignmentMother[i]=1
      } else {
        dyads$assignmentMother[i]=0
      }
    }
  }
  
  # Add the genotype of the father
  genoFather=read.table(paste(datadir,path,"/fathers.txt",sep=""),header=FALSE,sep=" ")
  genoMother=read.table(paste(datadir,path,"/mothers.txt",sep=""),header=FALSE,sep=" ")
  # columns 6 to 21  are for father genotype
  # columns 22 to 37 are for mother genotype
  dyads[,6:37]=rep(NA,nrow(dyads))
  for (i in 1:nrow(dyads)) {
    if (!is.na(dyads$father[i])){
      dyads[i,6:21]=genoFather[which(genoFather[,1]==dyads$father[i]),2:17]
    }
    if (!is.na(dyads$mother[i])){
      dyads[i,22:37]=genoMother[which(genoMother[,1]==dyads$mother[i]),2:17]
    }
  }
  
  # Export
  cat(length(which(dyads$assignmentFather==1)),"over",length(which(dyads$assignmentFather==0 | dyads$assignmentFather==1)),"fathers assigned :",length(which(dyads$assignmentFather==1))/length(which(dyads$assignmentFather==0 | dyads$assignmentFather==1))*100,"%\n")
  cat(length(which(dyads$assignmentMother==1)),"over",length(which(dyads$assignmentMother==0 | dyads$assignmentMother==1)),"mothers assigned :",length(which(dyads$assignmentMother==1))/length(which(dyads$assignmentMother==0 | dyads$assignmentMother==1))*100,"%\n")
  write.table(dyads,paste(path,"_resultsWithInfo.txt",sep=""),col.names=TRUE,row.names=FALSE)
  return(dyads)
}