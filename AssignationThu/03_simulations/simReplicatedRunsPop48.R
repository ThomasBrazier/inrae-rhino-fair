#===============================#
#            COLONY
#  fonction simulation de replicats de runs
#===============================#

# working dir INRA
#wd="W:/R/Assignation/03_simulations/"
#wd="C:/Users/tbrazier/Desktop/R/Assignation/03_simulations/"
# working dir HOME
#wd="~/Dropbox/M1 EFCE/INRA RHINO/R/Assignation/03_simulations/"
#setwd(wd)


# Cette fonction cree un jeu de donnees simulees pour tester l'influence du nombre de runs sur l'efficacite de l'assignation

# Le parametre d'interet est n (nombre de replicats, nombre de runs avec une graine differente)
# On fait n runs avec le MEME dataset simule une fois et copie (n - 1) fois
# (changement de graine a chaque copie de dataset)

simReplicatedRunsPop48=function() {
  timeIn=Sys.time()
    inputcolony=data.frame()
    
    # Parametres
    inputcolony=rbind(inputcolony,
                      data.frame(a="simReplicatedRunsPop50   ! Project name",b="")
    )
    inputcolony=rbind(inputcolony,
                      data.frame(a="1   ! Number of replicates",b="")
    )
    
    # Number of replicates of mating matrices
    # Total number of offsprings = # of mating matrices * sum of offspring in a single mating matrix (12)
    # total number of mums = #mating matrices * #mums in a matrix (12)
    # total number of dads = #mating matrices * #dads in a matrix (12)
    # Pour 48 parents
    inputcolony=rbind(inputcolony,
                           data.frame(a=c("1    !0/1/2=Pair-likelihood Score(PLS)/Full likelihood(FL)/FL-PLS combined (FPLS) method",
                                          "1    !0/1/2/3 for low/medium/high/very high precision",
                                          "2 0  !2/1=Dioecious/Monoecious Selfing rate for monoecious",
                                          "2    !Number of mating matrices",
                                          "12 12 !dads & mums in a mating structure"),b="")
    )

    # Matrix of family structure : number of offsprings per dyad male-female
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    
    # inputcolony=rbind(inputcolony,
    #                   data.frame(a=c(paste(matingstr)),b="")
    # ) 
    inputcolony=rbind(inputcolony,
                      data.frame(a=c("1 0 0 0 0 0 0 0 0 0 0 0 ! matrix of family structure",
                                     "0 1 1 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 1 1 1 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 1 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 1 1 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 1 1 1",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0"),b="")
    ) 

    
    #full siblings of known dad i (row, 1-12) and known mum j (col, 1-12)
    # a priori 0
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    
    inputcolony=rbind(inputcolony,
                      data.frame(a=c("0 0 0 0 0 0 0 0 0 0 0 0 ! full siblings of known dad i (row, 1-12) and known mum j (col, 1-12)",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0"),b="")
    )      
    
    #full siblings of known dad i (row, 1-12) and unknown mum j (col, 1-12)
    # a priori 0
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    
    inputcolony=rbind(inputcolony,
                      data.frame(a=c("0 0 0 0 0 0 0 0 0 0 0 0 ! full siblings of known dad i (row, 1-12) and unknown mum j (col, 1-12)",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0"),b="")
    )  
    
    #full siblings of unknown dad (row, 1-12) and known mum j (col, 1-12)
    # a priori 0
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    
    inputcolony=rbind(inputcolony,
                      data.frame(a=c("0 0 0 0 0 0 0 0 0 0 0 0 ! full siblings of unknown dad i (row, 1-12) and known mum j (col, 1-12)",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0",
                                     "0 0 0 0 0 0 0 0 0 0 0 0"),b="")
    )   
    
    #Dad 1-9 included (1) in or excluded (0) from candidate list
    #mum 1-18 included (1) in or excluded (0) from candidate list
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    
    inputcolony=rbind(inputcolony,
                      data.frame(a=c("1 1 1 1 0 0 0 0 0 0 0 0 !Dad 1-12 included (1) in or excluded (0) from candidate list 33%",
                                     "1 1 1 1 1 1 1 1 1 1 1 0 !mum 1-12 included (1) in or excluded (0) from candidate list 91%"),b="")
    ) 
    
    
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    
    # !# unrelated candidate males & females -> nul car on se considere en metapop fermee
    # !Assumed prbs of fathers & mothers included in candidates -> liee a la proba de detection
    # !Prob. of missing genotypes -> valeur par defaut 0.05
    inputcolony=rbind(inputcolony, 
                      data.frame(a=c("0 0       !# unrelated candidate males & females",
                                     "0.5  1    !Assumed prbs of fathers & mothers included in candidates",
                                     "8         !Number of Loci",
                                     "0.05      !Prob. of missing genotypes"),b="")
    )    
    
    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    # 8 loci
    # markers = rha101 rha109  rha4  rha7  rhc108  rhc3  rhd102  rhd103
    # alleles a chaque locus calcule d'apres le jeu de donnees
    inputcolony=rbind(inputcolony,
                      data.frame(a=c("0.05 0.05  0.05  0.05  0.05  0.05  0.05  0.05 !Drop rate for each locus",
                                     "0.01 0.01  0.01  0.01  0.01  0.01  0.01  0.01 !OtherErrorRate for each locus",
                                     "0 0 0 0 0 0 0 0                               !Marker types, 0/1=codominant/dominant",
                                     "14 10 21 14 14 6 11 10                        !# alleles at each locus"),b="")
    )

    inputcolony=rbind(inputcolony,data.frame(a="",b=""))
    
    inputcolony=rbind(inputcolony,
                      data.frame(a=c("0       !0/1/2/3=Uniform/Equal/Triangular/other allele freq. distr.",
                                     "2       !1/n=HaploDiploid/n-ploid species",
                                     "0 1     ! 0/1=Polygamy/Monogamy for males & females",
                                     paste(round(runif(1,1000,9999),0),"      ! Seed for random number generator",sep=""),
                                     "1  1  1 !0/1/2/3=No/Weak/Medium/Strong sibship prior",
                                     "0       !B, 0/1=Clone inference =No/Yes",
                                     "1       !B, 0/1=Scale full sibship=No/Yes",
                                     "0       !1/0 (Y/N) for known allele frequency",
                                     "0       !1/0 for updating allele freq. or not",
                                     "1       ! number of replicate runs",
                                     "0       !0/1/2/3/4=VeryShort/Short/Medium/Long/VeryLong run",
                                     "-1      !Map length in Morgans. <0 for infinite length",
                                     "0.0     !Inbreeding coefficient of parents",
                                     "0       !0/1=N/Y for allowing inbreeding in Colony",
                                     "0       !0/1 for Windows DOS/GUI run, parameter in Colony",
                                     "10000   !#iterates/#second for Windows DOS/GUI run, parameter in Colony"),b="")
    )
    
    # #Fin
    write.table(inputcolony,"simReplicatedRunsPop48/Colony.Par",quote=F,row.names = F,col.names = F)

  timeOut=Sys.time()
  cat("Construction des fichiers de simulation en",difftime(timeOut,timeIn),"secondes.\n")
}


#--------------------------#
# Ecriture des commandes Shell
# pour lancer les .Dat dans Colony
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
shellCommandSim1=function(n=1){
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
cd/home/users/tbrazier/Colony\n",
                       stringsAsFactors = F)
  for (r in 1:n) {
    Commandes=rbind(Commandes,paste("qsub colony2s.gnu.out IFN:",paste("Colony",r,".Dat",sep=""),sep=""))
  }
  write.table(Commandes,"simReplicatedRuns/qsColThomas",quote=F,row.names = F,col.names = F)
}