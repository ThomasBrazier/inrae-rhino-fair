#===============================#
#            Assignation
#   Preparation Inputs de COLONY
#===============================#


# D'apres un script fourni par Eric Petit

unique=read.table("uniqueGenotypesWithInfo.txt",h=T)

prepaColony=function(){
  #######################################################################
  # preparing the files for colony(see MM_pedigrees.docx) #
  #######################################################################
  # OFS
  # Offsprings potentiels
  index=which(unique$ageWhenFirstCaptur=="Juv")
  length(index)
  Rhino_OFS1=unique[index,c(1:17)]
  write.table(Rhino_OFS1,"Rhino_OFS.txt",quote=F,row.names=F,col.names=F)
  
  # CMS
  # Males candidats
  # Tous les males + tous les juveniles males nes avant 2017
  index=c(
    which(unique$ageWhenFirstCaptur=="Juv" & unique$yearWhenFirstCaptur=="2015" & unique$sexe=="M")
    ,
    which(unique$ageWhenFirstCaptur=="Juv" & unique$yearWhenFirstCaptur=="2016" & unique$sexe=="M")
    ,
    which(unique$ageWhenFirstCaptur=="Adult" & unique$sexe=="M")
  )
  length(index)
  Rhino_CMS=unique[index,c(1:17)]
  write.table(Rhino_CMS,"Rhino_CMS.txt",quote=F,row.names=F,col.names=F)
  
  # CFS
  # Femelles candidates
  # Toutes les femelles + femelles juv nees avant 2017
  index=c(
    which(unique$ageWhenFirstCaptur=="Juv" & unique$yearWhenFirstCaptur=="2015" & unique$sexe=="F")
    ,
    which(unique$ageWhenFirstCaptur=="Juv" & unique$yearWhenFirstCaptur=="2016" & unique$sexe=="F")
    ,
    which(unique$ageWhenFirstCaptur=="Adult" & unique$sexe=="F")
  )
  length(index)
  Rhino_CFS=unique[index,c(1:17)]
  write.table(Rhino_CFS,"Rhino_CFS.txt",quote=F,row.names=F,col.names=F)
  
}

#-----------------------#
#       ConstructColony
#-----------------------

constructColony=function(n=1){
  
  timeIn=Sys.time()
  
  library(dplyr)
  
  for (r in 1:n) {
    inputcolony=data.frame()
    
    # Parametres
    inputcolony=rbind(inputcolony,data.frame(a=paste("RhinoThomas",r," ! Project name",sep=""),b=""))
    inputcolony=rbind(inputcolony,data.frame(a=paste("RhinoThomas",r," ! Output file name",sep=""),b=""))
    
    inputcolony=rbind(inputcolony,
                      data.frame(a=c(paste(nrow(read.table(paste(wd,"Rhino_OFS.txt",sep=""))), "         ! Number of offspring in the sample",sep=""),
                                     "8         ! Number of loci",
                                     paste(round(runif(1,1000,9999),0),"      ! Seed for random number generator",sep=""),
                                     "0         ! 0/1=Not updating/updating allele frequency",
                                     "2         ! 2/1=Dioecious/Monoecious species",
                                     "0         ! 0/1=Inbreeding absent/present",
                                     "0         ! 0/1=Diploid species/HaploDiploid species",
                                     "0 0       ! 0/1=Polygamy/Monogamy for males & females",
                                     "0         ! 0/1 = Clone inference = No/Yes",
                                     "1         ! 0/1=Scale full sibship=No/Yes",
                                     "1 1 1     ! 0/1/2/3/4=No/Weak/Medium/Strong sibship prior; 4=Optimal sibship prior for Ne",
                                     "0         ! 0/1=Unknown/Known population allele frequency",
                                     "1         ! Number of runs",
                                     "2         ! 1/2/3/4 = Short/Medium/Long/VeryLong run",
                                     "0         ! 0/1=Monitor method by Iterate",
                                     "10000     ! Monitor interval in Iterate",
                                     "0         ! 0/1=DOS/Windows version",
                                     "1         ! 0/1/2=Pair-Likelihood-Score(PLS)/Full-Likelihood(FL)/FL-PLS-combined(FPLS) method",
                                     "1         ! 0/1/2/3=Low/Medium/High/VeryHigh precision",
                                     "rha101 rha109  rha4  rha7  rhc108  rhc3  rhd102  rhd103",
                                     "0  0 0 0 0 0 0 0",
                                     "0.07 0.05  0.07  0.07  0.09  0.10  0.06  0.05",
                                     "0.01 0.02  0.03  0.03  0.02  0.05  0.02  0.03",
                                     ""),b="")
    )
    
    # #Offsprings
    inputcolony=rbind(inputcolony,
                      data.frame(a=apply(read.table(paste(wd,"Rhino_OFS.txt",sep="")),1,function (x) paste(x,collapse=" ")),b="!Offspring ID and genotypes")
    )
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    
    # Proba de trouver des parents
    inputcolony=rbind(inputcolony,
                      data.frame(a=c("0.5  1"),b="! Proba de trouver des parents"))
    #   # Nombre de parents
    inputcolony=rbind(inputcolony,
                      data.frame(a=paste(nrow(read.table(paste(wd,"Rhino_CMS.txt",sep=""))),"  ",
                                         nrow(read.table(paste(wd,"Rhino_CFS.txt",sep=""))),sep=""),b="!Numbers of candidate males and females" )
    )
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    

    #   candidate males
    inputcolony=rbind(inputcolony,
                      data.frame(a=apply(read.table(paste(wd,"Rhino_CMS.txt",sep="")),1,function (x) paste(x,collapse=" ")),b="!Candidate males ID and genotypes")
    )
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    # Candidate females
    inputcolony=rbind(inputcolony,
                      data.frame(a=apply(read.table(paste(wd,"Rhino_CFS.txt",sep="")),1,function (x) paste(x,collapse=" ")),b="!Candidate females ID and genotypes")
    )
    
    
    # Zero
    # maj de Colony
    
    # Structure des donnees
    # !Number of offspring with known paternity, exclusion threshold
    # !IDs of known offspring-father dyad
    inputcolony=rbind(inputcolony,
                      data.frame(a="",b=""))
    inputcolony=rbind(inputcolony,
                      data.frame(a="0  0",b="!Number of offspring with known paternity, exclusion threshold"))    
    # !Number of offspring with known maternity, exclusion threshold
    inputcolony=rbind(inputcolony,
                      data.frame(a="",b=""))
    inputcolony=rbind(inputcolony,
                      data.frame(a="0  0",b="!Number of offspring with known maternity, exclusion threshold"))
    # !Number of known paternal sibship
    # !Size of known paternal sibship, and IDs of offspring in the sibship
    inputcolony=rbind(inputcolony,
                      data.frame(a="",b=""))
    inputcolony=rbind(inputcolony,
                      data.frame(a="0   ",b="!Number of known paternal sibship"))
    # !Number of known maternal sibship
    # !Size of known maternal sibship, and IDs of offspring in the sibship
    inputcolony=rbind(inputcolony,
                      data.frame(a="",b=""))
    inputcolony=rbind(inputcolony,
                      data.frame(a="0   ",b="!Number of known maternal sibship"))
    
    # !Number of offspring with known excluded paternity
    # !Offspring ID, number of excluded males, the IDs of excluded males
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    # Nombre d'excluded
    
    inputcolony=rbind(inputcolony,
                      data.frame(a=paste(nrow(read.table(paste(wd,"Rhino_ExcludedFathers.txt",sep=""),h=F,fill=TRUE))),b="!Number of offspring with known excluded paternity")
                      )
    
    #Excludes
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    
    #Excludes
    Exclude=apply(read.table(paste(wd,"Rhino_ExcludedFathers.txt",sep=""),sep="\n"),1,function (x) paste(x,collapse=" "))
    
    inputcolony=rbind(inputcolony,
                      data.frame(a=Exclude,b="")
    )

    # !Number of offspring with known excluded maternity
    # !Offspring ID, number of excluded females, the IDs of excluded females
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    #inputcolony=rbind(inputcolony,data.frame(a="0",b="!Number of offspring with known excluded maternity"))
    
    # Nombre d'excluded
    
    inputcolony=rbind(inputcolony,
                      data.frame(a=paste(length(apply(read.table(paste(wd,"Rhino_ExcludedMothers.txt",sep=""),sep="\n"),1,function (x) paste(x,collapse=" ")))),b="!Number of offspring with known excluded maternity")
    )

    #Excludes
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))

    #Excludes
    Exclude=apply(read.table(paste(wd,"Rhino_ExcludedMothers.txt",sep=""),sep="\n"),1,function (x) paste(x,collapse=" "))

    inputcolony=rbind(inputcolony,
                      data.frame(a=Exclude,b="")
    )
    
    
    # !Number of offspring with known excluded paternal sibships
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    inputcolony=rbind(inputcolony,data.frame(a="0  ",b="!Number of offspring with known excluded paternal sibships"))
    
    # !Number of offspring with known excluded maternal sibships
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    #♥inputcolony=rbind(inputcolony,data.frame(a="0",b="!Number of offspring with known excluded maternal sibships"))
    
    # Nombre d'excluded
    inputcolony=rbind(inputcolony,
                      data.frame(a=paste(length(apply(read.table(paste(wd,"Rhino_ExcludedMaternalSibs.txt",sep=""),sep="\n"),1,function (x) paste(x,collapse=" ")))),b="!Number of offspring with known excluded maternal sibships"))

    #Excludes
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))

   #Excludes
    Exclude=apply(read.table(paste(wd,"Rhino_ExcludedMaternalSibs.txt",sep=""),sep="\n"),1,function (x) paste(x,collapse=" "))

    inputcolony=rbind(inputcolony,
                      data.frame(a=Exclude,b="")
    )
    
    
    # #Fin
    # inputcolony=rbind(inputcolony,
    #                   data.frame(a=c("!_________________NOTE_________________!",
    #                                  "Input file generated by an R script",
    #                                  "itself generated by Pierre-Loup JAN",
    #                                  "and modified by Thomas BRAZIER"),b=""))
    # 
    write.table(inputcolony,paste(wd,"inputs/Colony",r,".Dat",sep=""),quote=F,row.names = F,col.names = F)
    
  }
  
  timeOut=Sys.time()
  cat("Construction des Colony.Dat en",difftime(timeOut,timeIn,units="secs"),"secondes.\n")
}


#--------------------------#
# Ecriture des commandes Shell
#--------------------------#
#!/bin/bash
#Script Torque
#PBS -l ncpus=1
#PBS -l mem=10000mb
#PBS -l cput=200:00:00
#PBS -V
#PBS -m abe
#PBS -M thomas.brazier@inra.fr
#PBS -l walltime=200:00:00
#PBS -q batch
# cri_job_type = serial
# cri_initialdir =	/usr/local/torque/SERIAL/IN
# cri_finaldir =	/usr/local/torque/SERIAL/OUT

#Commandes utilisateur
shellCommand=function(n=1){
  Commandes=data.frame(a="#!/bin/bash
#Script Torque
#PBS -l ncpus=2
#PBS -l mem=10000mb
#PBS -l cput=200:00:00
#PBS -V
#PBS -m abe
#PBS -M thomas.brazier@inra.fr
#PBS -l walltime=200:00:00
#PBS -q batch
# cri_job_type = serial
# cri_initialdir =	/usr/local/torque/SERIAL/IN
# cri_finaldir =	/usr/local/torque/SERIAL/OUT
cd /home/users/tbrazier/Colony\n",
                       stringsAsFactors = F)
  for (r in 1:n) {
    Commandes=rbind(Commandes,paste("colony2s.gnu.out IFN:",paste("Colony",r,".Dat",sep=""),sep=""))
  }
  write.table(Commandes,"inputs/qsColThomas",quote=F,row.names = F,col.names = F)
}
