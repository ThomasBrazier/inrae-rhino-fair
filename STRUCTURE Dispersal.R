#===============================#
#            COLONY
#     Natal dispersal
#     Migration of fathers estimated with STRUCTURE
#===============================#

#==========================================================
# LOADING ENVIRONMENT
#==========================================================

# clear global environment: remove all variables
rm(list=ls(all=TRUE))
# library(rstudioapi)
# Get the directory of the file & set working directory
# filedir=dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(filedir)

#----------------------------------------------------------#
# Loading packages
# Check if packages are installed, install if necessary
source("Sources/packages.R")


#=================================================
# Import data
#=================================================
dataPic = read.table("Data/Pic/uniqueGenotypesWithInfo.txt", header = TRUE,
                     stringsAsFactors = T)
dataThu = read.table("Data/Thu/uniqueGenotypesWithInfo.txt", header = TRUE,
                     stringsAsFactors = T)

resPic = read.table("Data/Pic/resultsWithInfo.txt", header = TRUE,
                    stringsAsFactors = T)
resThu = read.table("Data/Thu/resultsWithInfo.txt", header = TRUE,
                    stringsAsFactors = T)


dir.create("STRUCTURE/input")
dir.create("STRUCTURE/output")

dir.create("STRUCTURE/output/Pic")
dir.create("STRUCTURE/output/Thu")

# Use only male individuals, as genetic structure can be different among sexes in sex-biased dispersal
# See more in Dussex et al. 2016
# To gain more information about populations, females juveniles, not dispersing, were added
# dataPic = dataPic[-which(dataPic$sexe == "F" & dataPic$ageWhenFirstCaptur == "Adult"),]
# dataThu = dataThu[-which(dataThu$sexe == "F" & dataThu$ageWhenFirstCaptur == "Adult"),]


# Dispatch in pre- and post-dispersal subsets
# Pre-dispersal are female adults
# Post-dispersal are male adults

#=================================================
# Convert data to STRUCTURE input format
#=================================================
# Data file format
# Individuals in lines
# in one row, where each locus is in two consecutive columns.
# Columns are
# 1/ Labels of individuals (string)
# 2/ PopData (integer), id of the colony
# 3/ PopFlag (Optional; 0 or 1) A Boolean flag which indicates whether to use the PopData when using learning samples
# n+16 / Genotype data

# Create the individuals infos
# Assign to pre- or post-dispersal
# pre- is for learning (adult females)
# post- is for assignment (adult males)
samplePic = dataPic[which(dataPic$ageWhenFirstCaptur == "Adult"),]
dfPic = data.frame(indivID = samplePic$idind,
                   PopID = as.numeric(samplePic$idcol),
                   PopFlag = as.numeric(samplePic$sexe == "F")
)
# Colonies are encoded in integers
# Translate from integers to names
colonynamesPic = levels(samplePic$idcol)[as.numeric(samplePic$idcol)]
# Add genotypes
dfPic = cbind(dfPic, samplePic[,2:17])
# Save the input file
write.table(dfPic, file = "STRUCTURE/input/inputPic.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

#---------------------------------------
# And for Thuringia
#---------------------------------------
sampleThu = dataThu[which(dataThu$ageWhenFirstCaptur == "Adult"),]
dfThu = data.frame(indivID = sampleThu$idind,
                   PopID = as.numeric(sampleThu$idcol),
                   PopFlag = as.numeric(sampleThu$sexe == "F")
)
# Colonies are encoded in integers
# Translate from integers to names
colonynamesThu = levels(sampleThu$idcol)[as.numeric(sampleThu$idcol)]
# Add genotypes
dfThu = cbind(dfThu, sampleThu[,2:17])
# Save the input file
write.table(dfThu, file = "STRUCTURE/input/inputThu.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)


#=================================================
# STRUCTURE multiple runs
#=================================================
# K was known, as the number of colonies.
# "Then, the second step of the analysis involved assigning individuals
# to each of the K subpopulations. Samples were
# placed into the respective subpopulation based upon the
# highest percentage of membership (q). When comparingresults to other assignment test approaches, a threshold value of q≥0.90 was chosen"
# @cegelskiAssessingPopulationStructure2003

# USEPOPINFO model: pre-specify the population of origin of some individuals to
# assist ancestry estimation for individuals of unknown origin. A second way to use
# the USEPOPINFO model is to define “learning samples” that are pre-defined as coming from
# particular clusters. structure is then used to cluster the remaining individuals.
# Learning samples are implemented using the PopFlag column in the data file.
# The predefined population is used for those individuals for whom PopFlag=1 (and whose PopData is
# in (1...K)). The PopData value is ignored for individuals for whom PopFlag=0.


# 1000000 iterations with a 500000 burnin and a sampling interval of 100

# Command line options:
#   
# -m mainparams
# -e extraparams
# -s stratparams
# -K MAXPOPS 
# -L NUMLOCI
# -N NUMINDS
# -i input file
# -o output file
# -D SEED

# Settings different from default values:
#define ONEROWPERIND 1    // (B) store data for individuals in a single line
#define MISSING     0    // (int) value given to missing genotype data

#define POPFLAG   1     // (B) Input file contains a flag which says 
# whether to use popinfo when USEPOPINFO==1

#define GENSBACK    0  //(int) For use when inferring whether an indiv-
# idual is an immigrant, or has an immigrant an-
#   cestor in the past GENSBACK generations.  eg, if 
# GENSBACK==2, it tests for immigrant ancestry 
# back to grandparents.

# Define a pool of pre-dispersal individuals (juveniles) and a pool of post-dispersal individuals to assign to existing clusters (adults)
# define PFROMPOPFLAGONLY 1 // (B) only use individuals with POPFLAG=1 to update	P.
# This is to enable use of a reference set of 
# individuals for clustering additional "test" 
# individuals.

#define ANCESTDIST   1  // (B) collect data about the distribution of an-
# cestry coefficients (Q) for each individual

#define PRINTQHAT    1  // (B) Q-hat printed to a separate file.  Turn this 
# on before using STRAT.

# Analysis performed on all individuals.
# Settings passed in command line
# -m mainparams
# -e extraparams
# -K MAXPOPS 
# -L NUMLOCI
# -N NUMINDS
# -i input file
# -o output file
# -D SEED
niter = 1000000
burnin = 1000000
sampling = 100 # Though STRUCTURE does not subsample the MCMC chain, it can be useful for subsequent analyses
# Run 10 independent MCMC chains for Picardy
for (i in 1:10) {
  cat("============================\n")
  cat("Run", i, "\n")
  cat("============================\n\n\n")
  system(paste("rm STRUCTURE/output/Pic/trace", i, sep = "")) # Remove past trace files
  # Trace files are the output of STRUCTURE written in a file
  system("chmod +x ./STRUCTURE/structure_linux")
  system(paste("./STRUCTURE/structure_linux -m STRUCTURE/mainparams -e STRUCTURE/extraparams -K ", length(unique(colonynamesPic)),
               " -L 8 -N ", nrow(dfPic),
               " -i STRUCTURE/input/inputPic.txt -o STRUCTURE/output/Pic/run", i, " -D ", sample(1:1000, 1), " >> STRUCTURE/output/Pic/trace", i, sep = ""))
}

# Run for Thuringia as well
for (i in 1:10) {
  cat("============================\n")
  cat("Run", i, "\n")
  cat("============================\n\n\n")
  system(paste("./STRUCTURE/structure_linux -m STRUCTURE/mainparams -e STRUCTURE/extraparams -K ", length(unique(colonynamesThu)),
               " -L 8 -N ", nrow(dfThu),
               " -i STRUCTURE/input/inputThu.txt -o STRUCTURE/output/Thu/run", i, " -D ", sample(1:1000, 1), " >> STRUCTURE/output/Thu/trace", i, sep = ""))
}



#=================================================
# Diagnostic - Assess convergence within independant MCMC chains
#=================================================
# Extract sampling chains from trace files
# Use the 'structure_MCMC_samplingchains.sh' scripts in output directories

# Now sampling chains were extracted and we could retrieve the trace properly

# Load traces from the 10 independant MCMC chains
i = 1
trace = read.table(paste("STRUCTURE/output/Pic/sampling_chain_", i, sep = ""), header = TRUE)
for (i in 2:10) {
  cat(i, "\n")
  df = read.table(paste("STRUCTURE/output/Pic/sampling_chain_", i, sep = ""), header = TRUE)
  trace = rbind(trace, df)
}
write.table(trace, "STRUCTURE/output/Pic/trace.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)
tracePic = read.table("STRUCTURE/output/Pic/trace.txt", header = TRUE)

# Load traces from the 10 independant MCMC chains
i = 1
trace = read.table(paste("STRUCTURE/output/Thu/sampling_chain_", i, sep = ""), header = TRUE)
for (i in 2:10) {
  cat(i, "\n")
  df = read.table(paste("STRUCTURE/output/Thu/sampling_chain_", i, sep = ""), header = TRUE)
  trace = rbind(trace, df)
}
write.table(trace, "STRUCTURE/output/Thu/trace.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)
traceThu = read.table("STRUCTURE/output/Thu/trace.txt", header = TRUE)



##### DIAGNOSTICS FOR CONVERGENCE OF SINGLE CHAINS
# see Hamra et al. 2013
# Trace + autocorr + density plots
# Visualize traces, with and without discarding burnin of 5,000,000
tracePic$Run = as.factor(tracePic$Run)
ggplot(data = tracePic, aes(x = Rep, y = Ln.Like)) + geom_line(aes(colour = Run))
# Remove NA from burnin
ggplot(data = tracePic[which(!is.na(tracePic$Ln.Like)),], aes(x = Rep, y = Ln.Like)) + geom_line(aes(colour = Run))
# Overall convergence attained for burnin of 1,000,000 and MCMC sampling of 1,000,000
# Yet a problem of local minima with departure from the chain before retunring to a stable chain.
# Besides, check convergence of alpha parameter
ggplot(data = tracePic[which(!is.na(tracePic$Alpha)),], aes(x = Rep, y = Alpha)) + geom_line(aes(colour = Run))
# Convergence for parameters seems globally achievable, yet with local minima


#========================================================================
# TESTS OF CONVERGENCE
#========================================================================
# Subsequently, we ran a battery of diagnostic tests of convergence on the
# sampling chains of Alpha, r and Ln Likelihood of K
#========================================================================
# Geweke Diagnostic of MC Chain Stability
#------------------------------------------------------------------------
# if the whole chain is stationary, the means of the values early and late in
# the sequence should be similar
# convergence diagnostic 'Z' is the difference between the 2 means divided by
# the asymptotic standard error of their difference
# values of 'Z' near the extreme tails of the N(0,1) indicates lack of
# convergence
# can also estimate p-value of 'Z' from the normal distribution
# yields the probability that the divided chain means are different
?geweke.diag
#------------------------------------------------------------------------
# Convergence tests in Picardy
i = 1
tracePic[tracePic$Run == i,]
mcmcPic = as.mcmc(tracePic[(tracePic$Run == 1 & !is.na(tracePic$Ln.Like)),c(3, 21)], start = 100, stop = 2000000, thin = 100)
traceplot(mcmcPic, lty = 1)
geweke.diag(mcmcPic)$z[1]
geweke.diag(mcmcPic)$z[2]
# Alpha Ln.Like  for run 1
# 0.3434 -0.9362 
# Convergence for Alpha, not for Ln.Lik

convergencePic = data.frame(Run = c(1:10), Alpha = NA, Ln.Lik = NA)
for (i in convergencePic$Run) {
  mcmcPic = as.mcmc(tracePic[(tracePic$Run == i & !is.na(tracePic$Ln.Like)),c(3, 21)], start = 100, stop = 2000000, thin = 100)
  convergencePic$Alpha[i] = geweke.diag(mcmcPic)$z[1]
  convergencePic$Ln.Lik[i] = geweke.diag(mcmcPic)$z[2]
}
# Which is the best under Geweke diagnostic?
convergencePic
which(abs(convergencePic$Alpha) < 0.6 & abs(convergencePic$Ln.Lik) < 0.6)

# Run 2
i = 2
mcmcPic = as.mcmc(tracePic[(tracePic$Run == i & !is.na(tracePic$Ln.Like)),c(3, 21)], start = 100, stop = 2000000, thin = 100)
geweke.diag(mcmcPic)
?raftery.diag
raftery.diag(mcmcPic)
traceplot(mcmcPic, lty = 1)



#========================================================================
# THURINGIA
#========================================================================

##### DIAGNOSTICS FOR CONVERGENCE OF SINGLE CHAINS
# see Hamra et al. 2013
# Trace + autocorr + density plots
# Visualize traces, with and without discarding burnin of 5,000,000
traceThu$Run = as.factor(traceThu$Run)
# ggplot(data = tracePic, aes(x = Rep, y = Ln.Like)) + geom_line(aes(colour = Run))
# Remove NA from burnin
ggplot(data = traceThu[which(!is.na(traceThu$Ln.Like)),], aes(x = Rep, y = Ln.Like)) + geom_line(aes(colour = Run))
# Overall convergence attained for burnin of 1,000,000 and MCMC sampling of 1,000,000
# Yet a problem of local minima with departure from the chain before retunring to a stable chain.
# Besides, check convergence of alpha parameter
ggplot(data = traceThu[which(!is.na(traceThu$Alpha)),], aes(x = Rep, y = Alpha)) + geom_line(aes(colour = Run))
# Convergence for parameters seems globally achievable


#========================================================================
# TESTS OF CONVERGENCE
#========================================================================
# Subsequently, we ran a battery of diagnostic tests of convergence on the
# sampling chains of Alpha, r and Ln Likelihood of K
#========================================================================
# Geweke Diagnostic of MC Chain Stability
#------------------------------------------------------------------------
# if the whole chain is stationary, the means of the values early and late in
# the sequence should be similar
# convergence diagnostic 'Z' is the difference between the 2 means divided by
# the asymptotic standard error of their difference
# values of 'Z' near the extreme tails of the N(0,1) indicates lack of
# convergence
# can also estimate p-value of 'Z' from the normal distribution
# yields the probability that the divided chain means are different
?geweke.diag
#------------------------------------------------------------------------
# Convergence tests in Picardy
i = 1
traceThu[traceThu$Run == i,]
mcmcThu = as.mcmc(traceThu[(traceThu$Run == 1 & !is.na(traceThu$Ln.Like)),c(3, 24)], start = 100, stop = 2000000, thin = 100)
traceplot(mcmcPic, lty = 1)
geweke.diag(mcmcThu)$z[1]
geweke.diag(mcmcThu)$z[2]
# Alpha Ln.Like  for run 1
# 0.3434 -0.9362 
# Convergence for Alpha, not for Ln.Lik

convergenceThu = data.frame(Run = c(1:10), Alpha = NA, Ln.Lik = NA)
for (i in convergenceThu$Run) {
  mcmcThu = as.mcmc(traceThu[(traceThu$Run == i & !is.na(traceThu$Ln.Like)),c(3, 24)], start = 100, stop = 2000000, thin = 100)
  convergenceThu$Alpha[i] = geweke.diag(mcmcThu)$z[1]
  convergenceThu$Ln.Lik[i] = geweke.diag(mcmcThu)$z[2]
}
# Which is the best under Geweke diagnostic?
convergenceThu
which(abs(convergenceThu$Alpha) < 0.6 & abs(convergenceThu$Ln.Lik) < 0.9)
# Run 2
i = 7
mcmcThu = as.mcmc(traceThu[(tracePic$Run == i & !is.na(traceThu$Ln.Like)),c(3, 24)], start = 100, stop = 2000000, thin = 100)
geweke.diag(mcmcThu)
raftery.diag(mcmcThu)
traceplot(mcmcThu, lty = 1)




# So far, run 2 was the best for Picardy and run 7 for Thuringia
# Minimizing the Geweke diagnostic criterion for both Alpha and Ln Likelihood



#=================================================
# Diagnostic - Assess convergence among different independant MCMC chains
#=================================================
# Not relevant because we choose the best sampling chain and do not average them

#=================================================
# Output natal dispersal distances of males involved in paternities
#=================================================

# Dussex et al. 2016
# "An individual was flagged as a possible migrant
# when the most probable village of membership was not its village of capture." Supplementary

#------------------------------------------------------------------------
# PICARDY
#------------------------------------------------------------------------
# Retrieve lines in the output file concerning only males involved in paternities
# idind = 2960 # Single  individual for testing
i = 2 # Select the run to estimate ancestry probabilities
individuals = unique(resPic$father) # List of individuals to retrieve (inferred fathers)

# Run a for loop on the list of individuals
# Re-init result file
system("rm STRUCTURE/output/Pic/results.txt")
for (idind in individuals) {
  # Copy the lines for individuals of interest in a result file
  system(paste("grep ' ", idind, " ' STRUCTURE/output/Pic/run", i, "_f >> STRUCTURE/output/Pic/results.txt", sep = ""))
}
system("cat STRUCTURE/output/Pic/results.txt")

PopAssignPic = read.table("STRUCTURE/output/Pic/results.txt", header = FALSE)
PopAssignPic = PopAssignPic[,-c(1, 3, 5)] # Remove unwanted columns
colnames(PopAssignPic) = c("idind", "idcol", levels(as.factor(colonynamesPic))) # Rename columns
PopAssignPic$idcol = (levels(as.factor(colonynamesPic))[PopAssignPic$idcol])
PopAssignPic = PopAssignPic[,1:19]

# Now which is the inferred ancestry and at which probability
PopAssignPic$ancestry = NA
PopAssignPic$probancestry = NA
for (i in 1:nrow(PopAssignPic)) {
  PopAssignPic$ancestry[i] = as.character(colnames(PopAssignPic[,3:19])[which.max(PopAssignPic[i,3:19])])
  PopAssignPic$probancestry[i] = max(PopAssignPic[i,3:19]) 
}
# Assess the validity of the inference
# Is the probability of the second probable ancestry much lower than the first one?

# Is there natal dispersal?
# i.e. ancestry and colony are the same
PopAssignPic$nataldisp = (PopAssignPic$idcol == PopAssignPic$ancestry)
table(PopAssignPic$nataldisp)
# 50/50 chances of migrants/non-migrants
mean(PopAssignPic$probancestry)
# Are probabilities of ancestry for migrant or non-migrant the same?
mean(PopAssignPic$probancestry[PopAssignPic$nataldisp == TRUE])
mean(PopAssignPic$probancestry[PopAssignPic$nataldisp == FALSE])
hist(PopAssignPic$probancestry[PopAssignPic$nataldisp == TRUE], breaks = 20, xlim = c(0, 0.7))
hist(PopAssignPic$probancestry[PopAssignPic$nataldisp == FALSE], breaks = 20, xlim = c(0, 0.7))
# Test it with Mann-Whitney
wilcox.test(PopAssignPic$probancestry[PopAssignPic$nataldisp == TRUE], PopAssignPic$probancestry[PopAssignPic$nataldisp == FALSE])
# Slightly significant for individuals of interest
# Astonishingly, fathers that are migrants have higher probabilities of assignment than fathers assigned to their residency colony

# Export results
dir.create("Tables")
write.table(PopAssignPic, "Tables/MigrantFathersPic.txt", col.names = TRUE, row.names = FALSE,
            quote = FALSE)
PopAssignPic = read.table("Tables/MigrantFathersPic.txt", header = TRUE)




# Yang et al. (2005): At least 15-20 markers to achieve a confident assignement. Huge loss of accuracy under 15 markers.
# If the most informative markers are selected, then we can reach >99% assignment with only 4-6 markers.
# Besides, STRUCTURE need a sufficient sample size to detect informations
# Broquet et al. 2009
# "The bias and precision of the results from Structure and BayesAss is
# nearly indistinguishable from our method (Fig. 8). Our
# method gave slightly more accurate results than Structure, while BayesAss returned the least-biased estimates
# (especially for philopatry). However, for BayesAss, the
# confidence intervals for the estimates with the largest error
# were often too small, which meant that they did not
# include the true value (this occurred in 17% of immigration
# estimates and 58% of philopatry estimates, as compared to
# 13% and 14%, respectively, with our approach)."
# STRUCTURE. Underestimation of philopatry.

# # Explore results for all individuals assigned (PopFlag == 0)
# i = 1 # Select the run to estimate ancestry probabilities
# individuals = unique(dfPic$indivID[dfPic$PopFlag == 0]) # List of individuals to retrieve (inferred fathers)
# 
# # Run a for loop on the list of individuals
# # Re-init result file
# system("rm STRUCTURE/output/Pic/resultsAll.txt")
# for (idind in individuals) {
#   # Copy the lines for individuals of interest in a result file
#   system(paste("grep ' ", idind, "    (' STRUCTURE/output/Pic/run", i, "_f >> STRUCTURE/output/Pic/resultsAll.txt", sep = ""))
# }
# # system("cat STRUCTURE/output/Pic/resultsAll.txt")
# 
# PopAssign = read.table("STRUCTURE/output/Pic/resultsAll.txt", header = FALSE)
# PopAssign = PopAssign[,-c(1, 3, 5)] # Remove unwanted columns
# colnames(PopAssign) = c("idind", "idcol", levels(as.factor(colonynamesPic))) # Rename columns
# PopAssign$idcol = (levels(as.factor(colonynamesPic))[PopAssign$idcol])
# PopAssign = PopAssign[,1:19]
# 
# # Now which is the inferred ancestry and at which probability
# PopAssign$ancestry = NA
# PopAssign$probancestry = NA
# for (i in 1:nrow(PopAssign)) {
#   PopAssign$ancestry[i] = as.character(colnames(PopAssign[,3:19])[which.max(PopAssign[i,3:19])])
#   PopAssign$probancestry[i] = max(PopAssign[i,3:19]) 
# }
# # Assess the validity of the inference
# # Is the probability of the second probable ancestry much lower than the first one?
# 
# # Is there natal dispersal?
# # i.e. ancestry and colony are the same
# PopAssign$nataldisp = (PopAssign$idcol != PopAssign$ancestry)
# table(PopAssign$nataldisp)
# # Proportion of migrants
# sum(PopAssign$nataldisp == TRUE)/nrow(PopAssign)
# # 68% of migrants seems a high percentage
# 
# # Are probabilities of ancestry for migrant or non-migrant the same?
# mean(PopAssign$probancestry[PopAssign$nataldisp == TRUE])
# mean(PopAssign$probancestry[PopAssign$nataldisp == FALSE])
# hist(PopAssign$probancestry[PopAssign$nataldisp == TRUE], breaks = 20, xlim = c(0, 0.7))
# hist(PopAssign$probancestry[PopAssign$nataldisp == FALSE], breaks = 20, xlim = c(0, 0.7))
# # Test it with Mann-Whitney
# wilcox.test(PopAssign$probancestry[PopAssign$nataldisp == TRUE], PopAssign$probancestry[PopAssign$nataldisp == FALSE])
# # Not significant for individuals of interest




#------------------------------------------------------------------------
# PICARDY
#------------------------------------------------------------------------
# Retrieve lines in the output file concerning only males involved in paternities
# idind = 2960 # Single  individual for testing
i = 7 # Select the run to estimate ancestry probabilities
individuals = unique(resThu$father) # List of individuals to retrieve (inferred fathers)

# Run a for loop on the list of individuals
# Re-init result file
system("rm STRUCTURE/output/Thu/results.txt")
for (idind in individuals) {
  # Copy the lines for individuals of interest in a result file
  system(paste("grep ' ", idind, " ' STRUCTURE/output/Thu/run", i, "_f >> STRUCTURE/output/Thu/results.txt", sep = ""))
}
system("cat STRUCTURE/output/Thu/results.txt")

PopAssignThu = read.table("STRUCTURE/output/Thu/results.txt", header = FALSE)
PopAssignThu = PopAssignThu[,-c(1, 3, 5)] # Remove unwanted columns
colnames(PopAssignThu) = c("idind", "idcol", levels(as.factor(colonynamesThu))) # Rename columns
PopAssignThu$idcol = (levels(as.factor(colonynamesThu))[PopAssignThu$idcol])
PopAssignThu = PopAssignThu[,1:19]

# Now which is the inferred ancestry and at which probability
PopAssignThu$ancestry = NA
PopAssignThu$probancestry = NA
for (i in 1:nrow(PopAssignThu)) {
  PopAssignThu$ancestry[i] = as.character(colnames(PopAssignThu[,3:19])[which.max(PopAssignThu[i,3:19])])
  PopAssignThu$probancestry[i] = max(PopAssignThu[i,3:19]) 
}
# Assess the validity of the inference
# Is the probability of the second probable ancestry much lower than the first one?

# Is there natal dispersal?
# i.e. ancestry and colony are the same
PopAssignThu$nataldisp = (PopAssignThu$idcol == PopAssignThu$ancestry)
table(PopAssignThu$nataldisp)
# 7 non-migrants for 11 migrants
mean(PopAssignThu$probancestry)
# Are probabilities of ancestry for migrant or non-migrant the same?
mean(PopAssignThu$probancestry[PopAssignThu$nataldisp == TRUE])
mean(PopAssignThu$probancestry[PopAssignThu$nataldisp == FALSE])
hist(PopAssignThu$probancestry[PopAssignThu$nataldisp == TRUE], breaks = 20, xlim = c(0, 0.7))
hist(PopAssignThu$probancestry[PopAssignThu$nataldisp == FALSE], breaks = 20, xlim = c(0, 0.7))
# Test it with Mann-Whitney
wilcox.test(PopAssignThu$probancestry[PopAssignThu$nataldisp == TRUE], PopAssignThu$probancestry[PopAssignThu$nataldisp == FALSE])
# Not significant for individuals of interest
# Besides probabilities associated to migrnats and n on-migrants are high
# especially compared to those of Picardy
# Maybe sampling is better and we have more polymorphism/structure in genetic markers
# More power to detect migrants

# Export results
write.table(PopAssignThu, "Tables/MigrantFathersThu.txt", col.names = TRUE, row.names = FALSE,
            quote = FALSE)
PopAssignThu = read.table("Tables/MigrantFathersThu.txt", header = TRUE)



