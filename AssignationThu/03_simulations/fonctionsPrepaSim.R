#===============================#
#            COLONY
#         Simulations
#===============================#


#-------------------------------
# Fonction ReplicateSim
#-------------------------------#

# Le run de Colony en simulation prend un fichier Input.Par en entree
# Il construit un fichier .Dat a partir de Input.Par
# Cette fonction permet de repliquer ce fichier .Dat
# en n exemplaires identiques, sauf pour la graine aleatoire

# Arguments
# 1/ le chemin d'acces du fichier dans le repertoire R
# 2/ le nombre de replicats

replicateSim=function(path="",r=2) {
  #wd="F:/R/Assignation/03_simulations/"
  
  if (r>=2) {
    # pour ajouter des fichiers numerotes sans toucher au fichier 1
    pb=txtProgressBar(min = 2, max = r, style = 3)
    for (i in 2:r) {
      # supprime version precedente
      if (file.exists(paste(wd,path,"/COLONY2_",(i),".DAT",sep=""))){
        file.remove(paste(wd,path,"/COLONY2_",(i),".DAT",sep=""))
      }
      # lit le fichier
      dat=readLines(paste(wd,path,"/COLONY2_",(i-1),".DAT",sep=""))
      # modifie dataset name et output file name
      dat[1]=paste(path,"_",i,"       ! Dataset name, Length<51",sep="")
      dat[2]=paste(path,"_",i,"       ! Main output file name, Length<21",sep="")
      # modifie la graine (ligne 5)
      dat[5]=paste(round(runif(1,1000,9999),0),"      ! Seed for random number generator",sep="")
      # ecrit le nouveau fichier
      writeLines(dat,paste(wd,path,"/COLONY2_",(i),".DAT",sep=""))
      Sys.sleep(0.01)
      # update progress bar
      setTxtProgressBar(pb,i)
    } 
  }
  close(pb)
}

#-------------------------------
# Fonction simShellCmd
#-------------------------------#
#Pour lancer Colony sur tous les fichiers DAT de data simulees

#Commandes utilisateur
simShellCmd=function(direct="",n=1){
  Commandes=data.frame(a="#!/bin/bash
                       #Script Torque
                       #PBS -l ncpus=2
                       #PBS -l mem=10000mb
                       #PBS -l cput=600:00:00
                       #PBS -V
                       #PBS -m abe
                       #PBS -M thomas.brazier@inra.fr
                       #PBS -l walltime=200:00:00
                       #PBS -q batch
                       # cri_job_type = serial
                       # cri_initialdir =	/usr/local/torque/SERIAL/IN
                       # cri_finaldir =	/usr/local/torque/SERIAL/OUT",
                       stringsAsFactors = F)
  Commandes=rbind(Commandes,paste("cd /home/users/tbrazier/",direct,sep=""))
  for (r in 1:n) {
    Commandes=rbind(Commandes,paste("colony2s.gnu.out IFN:COLONY2_",r,".DAT",sep=""))
  }
  write.table(Commandes,paste(wd,direct,"/qsColThomas",sep=""),quote=F,row.names = F,col.names = F)
}


#-------------------------------#
# CREATION D'UNE POP THEORIQUE ALEATOIRE
#-------------------------------

##### MATRICE D'ACCOUPLEMENT ALEATOIRE
# Procedure de force brute
# Tirage de matrices jusqu'a obtenir matrice conforme
# n = nombre de sous-structure independantes dans la matrice
# dont depend le nombre d'individus adultes = n*24 parents
# t = nombre de generations
MatingMatrix=function(n=10,t=4){
  timeIn=Sys.time()
  
  # Creation de la matrice totale
  #n=10 # nombre de sous-structure independantes dans la matrice
  # indices de 12 en 12
  index=seq(1,n*12,12)
  Mat=matrix(0,nrow=n*12,ncol=n*12)
  
  for (sub in index) {
    success=FALSE
    c=0
    # Calcul d'une matrice 12*12
    while (success==FALSE) {
      # Tirage de matrices carrees aléatoires selon regles de tirage :
      # Matrice representant donc un accouplement au hasard entre individus
      # t = nombre de generations
      # - 3*t valeurs a tirer parmi {0, 1, 2, 3} pour les males
      # - t valeurs a tirer parmi {0, 1, 2, 3} pour les femelles
      # matrice composee de sous-structure de 12*12 (nb de valeurs pour les males)
      # Matrice composee de n repetitions de ces sous-structures 12*12
      #print(c)
      M=matrix(0,nrow=12,ncol=12)
      for (j in 1:6) { # 1:8 parce qu'une fraction des pères ne se reproduit pas (valeur arbitraire)
        M[j,]=rbinom(n=12,size=4,0.15)
      }
      
      # Matrice soumise a une verification d'un ensemble de regles pour representer une population
      # biologiquement coherente
      # 1/ total(somme des colonnes) = total(somme des lignes)
      # 2/ somme ligne <= 3*t (pas plus de 3 juveniles par an par male)
      # 3/ somme colonne <=t (pas plus de 1 juvenile par an par femelle)
      if(sum(rowSums(M))!=sum(colSums(M))){
      }else{
        if(sum(colSums(M)>=t)>0){
        }else{
          if(sum(rowSums(M)>=(3*t))>0){
          }else{
            success=TRUE
          }
        }
      }
      c=c+1
    }
    write.table(x=M,file=paste("Matrice",sub,".txt",sep=""),row.names = FALSE,col.names = FALSE)
    
    # Ajout de la sous matrice dans la matrice totale
    Mat[sub:(sub+11),sub:(sub+11)]=M
  }
  write.table(x=Mat,file=paste("MatriceAccouplement.txt",sep=""),row.names = FALSE,col.names = FALSE)
  
  return(Mat)
  timeOut=Sys.time()
  cat("Generation aleatoire en",difftime(timeOut,timeIn,units="mins"),"minutes.\n")
}

##### CONSTRUCTION DE L'INPUT POUR SIMULATION
# Met en forme le fichier Input.Par
# avec les parametres choisis et la matrice d'accouplement
# 1/ n transmis a la fonction MatingMatrix()
# 2/ t transmis a la fonction MatingMatrix()
# 3/ pop, taille de la popualtion simulee, qui sera arrondie a l'inferieur pour etre un multiple de la taille de la mating matrix
ConstructSimInput=function(n=10,t=4,pop=4000,name="",mating=NA){
  # Directory
  cat("Deplacement dans le repertoire.\n")
  setwd(paste(wd,name,"/",sep=""))
  
  cat("Chargement d'une matrice d'accouplement.\n")
  if (is.na(mating)){
    # Construit d'abord la matrice d'accouplement
    M=MatingMatrix(n,t)
  }else{
    M=mating
  }
  
  cat("Calcul de la taille de population.\n")
  # Number of mating matrices
  number=floor(4000/(nrow(M)*2))
  # Recalculer la population adulte arrondie pour etre multiple des matrices d'accouplement
  pop=nrow(M)*2*number
  cat(pop," parents.\n")
  cat(sum(M)*number," juveniles.\n")
  
  #-----------------------#
  # Construction du fichier Input.Par
  #-----------------------
  library(dplyr)
  
  cat("Construction du fichier Input.Par\n")
  
  input=data.frame()
  
  # Parametres
  input=rbind(input,data.frame(a=paste(name," ! Project (output file) name",sep=""),b=""))
  
  input=rbind(input,data.frame(a=c(paste(1, "         ! Number of replicates",sep=""),
                                   "1         ! 0/1/2=Pair-Likelihood-Score(PLS)/Full-Likelihood(FL)/FL-PLS-combined(FPLS) method",
                                   "1         ! 0/1/2/3=Low/Medium/High/VeryHigh precision",
                                   "2 0       !2/1=Dioecious/Monoecious,Selfing rate for monoecious",
                                   paste(number,"         ! Number of mating matrices"),
                                   paste(nrow(M)," ",ncol(M),"         ! Dads and mums in mating matrices",sep="")),
                               b="")
  )
  
  # Mating structure
  for (i in 1:nrow(M)) {
    input=rbind(input,data.frame(a=paste(M[i,],collapse=" "),b=""))
  }
  
  input=rbind(input,data.frame(a="",b=""))
  # known dad and known mum
  for (i in 1:nrow(M)) {
    input=rbind(input,data.frame(a=paste(rep(0,ncol(M)),collapse=" "),b=""))
  }
  
  input=rbind(input,data.frame(a="",b=""))
  # known dad and unknown mum
  for (i in 1:nrow(M)) {
    input=rbind(input,data.frame(a=paste(rep(0,ncol(M)),collapse=" "),b=""))
  }
  
  input=rbind(input,data.frame(a="",b=""))
  # unknown dad and known mum
  for (i in 1:nrow(M)) {
    input=rbind(input,data.frame(a=paste(rep(0,ncol(M)),collapse=" "),b=""))
  }
  
  input=rbind(input,data.frame(a="",b=""))
  # Dads and mums included
  input=rbind(input,data.frame(a=paste(rep(1,ncol(M)),collapse=" "),b=""))
  input=rbind(input,data.frame(a=paste(rep(1,ncol(M)),collapse=" "),b=""))
  
  input=rbind(input,data.frame(a="",b=""))
  input=rbind(input,data.frame(a="",b=""))
  input=rbind(input,data.frame(a=c("0 0         ! Unrelated candidate males and females",
                                   "0.5 0.94         ! Assumed probabilities of father and mother included",
                                   "8          ! Number of loci",
                                   "0          ! Probability of missing genotypes",
                                   "0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05                !Dropout rates",
                                   "0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01                !Other error rates",
                                   "0 0 0 0 0 0 0 0                !Codominant/Dominant (0/1) markers",
                                   # pensez a adapter le nombre d'allele a chaque population en fonction des freq alleliques
                                   "9 8 15 10 12 8 11 8         !#alleles/locus, must be=2 for dominant marker",
                                   "3                              !0/1/2/3=Uniform/Equal/Triangular/other allele freq. distr."),
                               b=""))
  
  # Frequences alleliques
  # Doivent etre comprises entre 0.001 et 1
  # freq en dessous sont eliminees de la simulation
  options(scipen=50,digits=5)
  fAllel=read.table("FreqAllelique.txt",fill=TRUE)
  for (i in 1:nrow(fAllel)) {
    input=rbind(input,data.frame(a=paste(fAllel[i,which(!is.na(fAllel[i,]))],collapse=" "),
                                 b=""))
  }
  
  input=rbind(input,data.frame(a=c("2                              !1/n=HaploDiploid/n-ploid species",
                                   "0  1                           !1/0 =Mono/Polygamy for males & females",
                                   paste(round(runif(1,1000,9999),0),"      ! Seed for random number generator",sep=""),
                                   "1 1 1                          !I,R,R : I=0,1,2,3 for No,weak,medium,strong sibship size prior R,R=mean paternal & maternal sibship size",
                                   "0                              !B, 0/1=Clone inference =No/Yes",
                                   "1                              !B, 0/1=Scale full sibship=No/Yes",
                                   "1                              !1/0 (Y/N) for known allele frequency",
                                   "0                              !1/0 for updating allele freq. or not",
                                   "1                              !#replicate runs",
                                   "2                              !0/1/2/3/4=VeryShort/Short/Medium/Long/VeryLong run",
                                   "-1                             !Map length in Morgans. <0 for infinite length",
                                   "0                              !Inbreeding coefficient of parents",
                                   "0                              !0/1=N/Y for allowing inbreeding in Colony",
                                   "1                              !0/1 for Windows DOS/GUI run, parameter in Colony",
                                   "1                              !#iterates/#second for Windows DOS/GUI run, parameter in Colony"),
                               b="")
  )
  
  cat("Ecriture du fichier...\n")
  write.table(input,paste(wd,name,"/input3.Par",sep=""),quote=F,row.names = F,col.names = F)
  
  cat("Construction du fichier de projet.\n")
  project=data.frame()
  project=rbind(project,data.frame(a=c("Output file path & name : C:/ZSL/Colony/Simulation",
                                       "Number of loci : 8",
                                       paste(" Number of offspring in the sample : ",sum(M)*number,sep=""),
                                       "Outbreeding (0) or inbreeding (1) model : 0",
                                       paste("Number of male candidates : ",pop,sep=""),
                                       paste("Number of female candidates : ",pop,sep=""),
                                       "Number of known paternal sibships : 0",
                                       "Number of known maternal sibships : 0",
                                       "Number of offspring with excluded fathers : 0",
                                       "Number of offspring with excluded mothers : 0",
                                       "Male mating system : Polygamy",
                                       "Female mating system : Monogamy",
                                       "Number of threads : 1",
                                       "Number of Excluded Paternal Sibships : 0",
                                       "Number of Excluded Maternal Sibships : 0",
                                       "Dioecious (2) or monoecious (1) : 2",
                                       "Seed for random number generator : 1234",
                                       "Allele frequency : No updating by accounting for the inferred relationship",
                                       "Species : Diploid",
                                       "Sibship size prior : Yes",
                                       "paternal & maternal sibship sizes : 1  1",
                                       "Known population allele frequency : Yes",
                                       "Number of run : 1",
                                       "Length of run : Medium",
                                       "Monitor intermiediate results by : Every 1 second",
                                       "Prob. a dad is included in the male candidates : 0.50",
                                       "Prob. a mum is included in the female candidates : 0.94",
                                       "Project data input produced : 20/06/2018"),
                                   b=""))
  
  cat("Ecriture du fichier...\n")
  write.table(project,paste(wd,name,"/ProjectInformation.txt",sep=""),quote=F,row.names = F,col.names = F)
  
  cat("Fin.\n")
  setwd(wd)
}


##### LANCEMENT DE LA SIMULATION DU MODULE COLONY
# LaunchSim()
# Lance le module simulation de COLONY
# nettoie le dossier en enlevant tous les fichiers inutiles




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
CreateScenario=function(name="",unsampled=c(0,0,0),propIncomplete=0,probMissingLocus=0) {
  # Directory
  cat("Deplacement dans le repertoire.\n")
  setwd(paste(wd,name,"/",sep=""))
  
  cat("Ouverture du fichier.\n")
  f=file("COLONY2.DAT")
  colony=readLines(f)
  close(f)
  # Parcourt une fois toutes les lignes pour identifier les lignes d'interet, stockees dans trois vecteurs
  # - males (ligne commence par 'M[0-9]*,' et n'est pas suivie de F)
  # - femelles (ligne commence par 'F[0-9]*,')
  # - juveniles (Ligne commence par 'M[0-9]*F[0-9]*C1,')
  IDmales=c()
  IDfemales=c()
  IDoffspring=c()
  for (i in 1:length(colony)){
    if(grepl('F[0-9]*,',colony[i])){
      IDfemales[length(IDfemales)+1]=i
    }else{
      if(grepl('M[0-9]*,',colony[i])){
        IDmales[length(IDmales)+1]=i
      }else{
        if(grepl('M[0-9]*F[0-9]*C1,',colony[i])){
          IDoffspring[length(IDoffspring)+1]=i
        }
      }
    }
  }
  # Verification des indices
  # Juveniles, males et femelles se suivent, avec une ligne vide entre chacun
  # ET
  # length(males) == length(females)
  if((max(IDoffspring)+5 != min(IDmales)) | (max(IDmales)+2 != min(IDfemales))){
    warning("Individus manquants : les indices ne se suivent pas")
  }
  if(length(IDmales) != length(IDfemales)){
    warning("Individus manquants : pas le meme nombre de males et femelles")
  }
  
  # Echantillonnage aleatoire des individus qui vont disparaitre
  if (prod(unsampled<=1)==1){
    IDmales=sample(IDmales,size=round(unsampled[1]*length(IDmales)),replace=FALSE)
    IDfemales=sample(IDfemales,size=round(unsampled[2]*length(IDfemales)),replace=FALSE)
    IDoffspring=sample(IDoffspring,size=round(unsampled[3]*length(IDoffspring)),replace=FALSE)
  }else{
    IDmales=sample(IDmales,size=unsampled[1],replace=FALSE)
    IDfemales=sample(IDfemales,size=unsampled[2],replace=FALSE)
    IDoffspring=sample(IDoffspring,size=unsampled[3],replace=FALSE)
  }
  # SUPPRESSION
  colony=colony[-c(IDmales,IDfemales,IDoffspring)]
  
  # MISSING GENOTYPES
  # Reprend tous les indices une fois les individus enleves
  
  # Parcourt une fois toutes les lignes pour identifier les lignes d'interet, stockees dans trois vecteurs
  # - males (ligne commence par 'M[0-9]*,' et n'est pas suivie de F)
  # - femelles (ligne commence par 'F[0-9]*,')
  # - juveniles (Ligne commence par 'M[0-9]*F[0-9]*C1,')
  IDmales=c()
  IDfemales=c()
  IDoffspring=c()
  for (i in 1:length(colony)){
    if(grepl('F[0-9]*,',colony[i])){
      IDfemales[length(IDfemales)+1]=i
    }else{
      if(grepl('M[0-9]*,',colony[i])){
        IDmales[length(IDmales)+1]=i
      }else{
        if(grepl('M[0-9]*F[0-9]*C[0-9],',colony[i])){
          IDoffspring[length(IDoffspring)+1]=i
        }
      }
    }
  }
  # Remplacer virgules par des espaces
  for (x in c(IDmales,IDfemales,IDoffspring)) {
    colony[x]=gsub(","," ",x=colony[x])
  }
  
  # Mettre a jour les infos sur le nombre d'individus
  # Offspring
  colony[3]=sub("[0-9]+",as.character(length(IDoffspring)),x=colony[3])
  # Parents
  colony[max(IDoffspring)+1]=" "
  colony[max(IDoffspring)+2]=paste(1-unsampled[1],"  ",1-unsampled[2],"  !Prob that the dad & mum of an offspring included in candidates")
  colony[max(IDoffspring)+3]=paste(length(IDmales),"  ",length(IDfemales),"  !Number of candidate males & females")
  
  # Echantillonnage aleatoire des lignes avec genotype incomplet
  if (propIncomplete<=1){
    IDmales=sample(IDmales,size=round(propIncomplete*length(IDmales)),replace=FALSE)
    IDfemales=sample(IDfemales,size=round(propIncomplete*length(IDfemales)),replace=FALSE)
    IDoffspring=sample(IDoffspring,size=round(propIncomplete*length(IDoffspring)),replace=FALSE)
  }else{
    IDmales=sample(IDmales,size=propIncomplete,replace=FALSE)
    IDfemales=sample(IDfemales,size=propIncomplete,replace=FALSE)
    IDoffspring=sample(IDoffspring,size=propIncomplete,replace=FALSE)
  }
  missingGen=c(IDmales,IDfemales,IDoffspring)
  # Echantillonnage aleatoire des loci manquant (1:16)
  # pour chaque ligne dans le vecteur des genotypes manquants, retire certains loci aleatoirement
  for (j in missingGen) {
    v=colony[j]
    values=gregexpr(' [0-9]{1}',v)[[1]][1:16]+1
    suppr=sample(values,size=round(probMissingLocus*16),replace=FALSE)
    for (i in 1:length(suppr)) {
      substr(v,start=suppr[i],stop=suppr[i])="0"
    }
    colony[j]=v
  }
  
  # Fin du programme et retour du dataset
  cat("Ecriture du fichier...")
  write.table(colony,"COLONY2_1.DAT",quote=F,row.names = F,col.names = F)
  
  #return(colony)
  
  cat("Fin.\n")
  setwd(wd)
}