#=============================================================================#
#     Automatically generates all the results necessary for the manuscript
#=============================================================================#


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



#=============================================================================#
# Dataset ----
#=============================================================================#

#-----------------------------------------------------------------------------#
# Picardy
#-----------------------------------------------------------------------------#
# Of the 3,706 unique genotypes in Picardy (France), 1,602 individuals were labelled as offspring (43%). 72% of
# individuals considered as adults were females (58% to 92% females depending on colony). 21% of individuals
# had an incomplete genotype (1 or 2 missing loci). The mean sample size per colony was 218±36 (s.e.m.)
# individuals for cumulated years 2013-2016, ranging from a minimal sample size of 53 bats to a maximal
# sample size of 746 bats.

genotypes = read.table("Data/Pic/uniqueGenotypesWithInfo.txt", header = TRUE)

# number of genotypes
nrow(genotypes)
# Numbers of juveniles and adults
table(genotypes$ageWhenFirstCaptur)
table(genotypes$ageWhenFirstCaptur) / nrow(genotypes)

# Number of adults females and males
table(genotypes$ageWhenFirstCaptur, genotypes$sexe)
table(genotypes$sexe[which(genotypes$ageWhenFirstCaptur == "Adult")])
table(genotypes$sexe[which(genotypes$ageWhenFirstCaptur == "Adult")]) / length(genotypes$sexe[which(genotypes$ageWhenFirstCaptur == "Adult")])

# Number of incomplete genotypes
n_incomplete = apply(genotypes[,2:17], 1, function(x) {sum(x == 0) / 2})
table(n_incomplete)

# Mean sample size of a colony
table(genotypes$idcol)
mean(table(genotypes$idcol))
min(table(genotypes$idcol))
max(table(genotypes$idcol))



#-----------------------------------------------------------------------------#
# Thuringia
#-----------------------------------------------------------------------------#
# In Thuringia (Germany), there were 3,913 unique genotypes with 2,128 potential
# offspring (54%). 62% of individuals considered as adults were females. 19% of individuals had an incomplete
# genotype (1 or 2 missing loci). The mean sample size per colony was 126±22 (s.e.m.) for years 2015-2017
# (range: 10 to 385).

genotypes = read.table("Data/Thu/uniqueGenotypesWithInfo.txt", header = TRUE)

# number of genotypes
nrow(genotypes)
# Numbers of juveniles and adults
table(genotypes$ageWhenFirstCaptur)
table(genotypes$ageWhenFirstCaptur) / nrow(genotypes)

# Number of adults females and males
table(genotypes$ageWhenFirstCaptur, genotypes$sexe)
table(genotypes$sexe[which(genotypes$ageWhenFirstCaptur == "Adult")])
table(genotypes$sexe[which(genotypes$ageWhenFirstCaptur == "Adult")]) / length(genotypes$sexe[which(genotypes$ageWhenFirstCaptur == "Adult")])

# Number of incomplete genotypes
n_incomplete = apply(genotypes[,2:17], 1, function(x) {sum(x == 0) / 2})
table(n_incomplete)

# Mean sample size of a colony
table(genotypes$idcol)
# Thu42 was excluded (1 individual)
mean(table(genotypes$idcol[which(genotypes$idcol != "Thu42")]))
min(table(genotypes$idcol[which(genotypes$idcol != "Thu42")]))
max(table(genotypes$idcol[which(genotypes$idcol != "Thu42")]))



#=============================================================================#
# Simulations ----
#=============================================================================#
# The proportion of successful father assignments was low with a high prevalence of type II errors (no father
# found for a given offspring, Table S1). In simulations, a conservative frequency value of 0.3 in Picardy and
# 0.25 in Thuringia accurately excluded wrong fathers when at least 20 runs of COLONY were used, even for
# the least powerful analysis with 6/8 complete loci (Fig. 1).

# Table S1: Mean assignment frequency and standard deviation (SD) for each category of result in power
# simulations (simulated datasets). True fathers were known. (a) and (e) were correct results. (b) and (d) were
# type I errors, (c) was type II error. Expected frequencies were calculated under the hypothesis of 100% of
# correct results. Frequencies were calculated for 3,892 offspring and 30 runs.


#-----------------------------------------------------------------------------#
# Picardy
#-----------------------------------------------------------------------------#
simulations_info = read.table("AssignationPic/simRandomPop_resultsWithInfo.txt",
                              header = T)
simulations_incomplete_info = read.table("AssignationPic/simRandomPopIncomplete_resultsWithInfo.txt",
                              header = T)
# - offspring ID
# - fathers and mothers ID
# - sex
# - correct/incorrect assignation
# - genotype

# How many offsprings?
nrow(simulations_info)



simulation_results = read.table("AssignationPic/outputs/simSensitivity.txt",
                                    header = T)
results_Pic = simulation_results[c(1,3,2,4,5),]


#-----------------------------------------------------------------------------#
# Thuringia
#-----------------------------------------------------------------------------#

simulation_results = read.table("AssignationThu/outputs/simSensitivity.txt",
                                header = T)


results_Thu = simulation_results[c(1,3,2,4,5),]


#-----------------------------------------------------------------------------#
# Comestic for output
#-----------------------------------------------------------------------------#
results_sim = cbind(results_Pic,
                    results_Thu[,3:6])

results_sim$ExpectedProportion = c(0.56, 0, 0, 0, 0.44)

colnames(results_sim) = c("Results",
                          "Expected proportion",
                          "Mean Picardy 8/8 loci",
                          "SD Picardy 8/8 loci",
                          "Mean Picardy 6/8 loci",
                          "SD Picardy 6/8 loci",
                          "Mean Thuringia 8/8 loci",
                          "SD Thuringia 8/8 loci",
                          "Mean Thuringia 6/8 loci",
                          "SD Thuringia 6/8 loci")

write.table(results_sim, "Tables/TableS1.tsv", sep = "\t",
            col.names = T, row.names = F)


#=============================================================================#
# Parentage assignment ----
#=============================================================================#

# We pooled inferences from genotypes with no more than 2 missing loci (≥ 6/8 loci) with the inferences from
# the complete genotypes (8/8 loci). Ambiguous assignments (different fathers for the same offspring) were
# resolved in favour of the most powerful dataset (8/8 complete loci). COLONY found 530 unique paternities
# (33% of putative offspring) in the French dataset. Based on simulation results (Fig. 1), we removed inferred
# paternities found in less than 30% of 30 runs, resulting in 82 paternities being considered as true paternities.
# There were 17 fathers responsible for these 82 paternities over four years. From the German dataset with
# the same procedure, 62 paternities were retained with only 18 fathers over three years.

#-----------------------------------------------------------------------------#
# Distribution of offspring-father distances Picardy
#-----------------------------------------------------------------------------#
# The combined 8/8 and 6/8 dataset of inferred paternities AFTER filtering
dyadsObsSelect = read.table("AssignationPic/outputs/dyadsObsSelect.txt", header=TRUE)

# How many paternitites inferred AFTER filtering?
nrow(dyadsObsSelect)

# How many fathers?
table(dyadsObsSelect$fatherID)
length(unique(dyadsObsSelect$fatherID))


#-----------------------------------------------------------------------------#
# Distribution of offspring-father distances Thuringia
#-----------------------------------------------------------------------------#
# The combined 8/8 and 6/8 dataset of inferred paternities AFTER filtering
dyadsObsSelect = read.table("AssignationThu/outputs/dyadsObsSelect.txt", header=TRUE)

# How many paternitites inferred AFTER filtering?
nrow(dyadsObsSelect)

# How many fathers?
table(dyadsObsSelect$fatherID)
length(unique(dyadsObsSelect$fatherID))




#=============================================================================#
# Mating dispersal ----
#=============================================================================#

# Based on the empirical distribution of mating dispersal distances, the Kolmogorov-Smirnov test between the
# two distributions of offspring-father distances (expected and observed) showed a significant deviation (D =
# 0.40 in France and 0.53 in Germany, p < 0.001 in both cases, n = 82 and 62) with more zero distances than
# expected by chance (Fig. 2). The mean mating dispersal distance was 11.3 km (IC95 = [7.4; 15.7], 10,000
# bootstraps, n = 82) in France and 11.6 km (IC95 = [7.3; 16.6], 10,000 bootstraps, n = 62) in Germany.
# However, the median distance were much lower in these two populations (7.7 km, IC95 = [5.5, 8.6] and 4.7
# km, IC95 = [0, 7.8], respectively). Both empirical distributions were globally congruent (Fig. 2) and the
# MasterBayes Bayesian approach yielded similar results (Supplementary Information).


#-----------------------------------------------------------------------------#
# Distribution of offspring-father distances Picardy
#-----------------------------------------------------------------------------#
# The expected random distribution
dyadsH0 = read.table("AssignationPic/outputs/dyadsH0.txt", header=TRUE)

# The combined 8/8 and 6/8 dataset of inferred paternities AFTER filtering
dyadsObsSelect = read.table("AssignationPic/outputs/dyadsObsSelect.txt", header=TRUE)

# The Kolmogorov-Smirnov test between the
# two distributions of offspring-father distances (expected and observed)
ks.test(dyadsH0$distance, dyadsObsSelect$distance)


# The mean mating dispersal distance
# Bootstrap implies some randomness in the results
nboot = 10000
boot = numeric(nboot)

for (i in 1:nboot) {
  boot[i] = mean(sample(dyadsObsSelect$distance, replace = T))
}

mean(boot)
quantile(boot, c(0.025, 0.975))


# The median mating dispersal distance
for (i in 1:nboot) {
  boot[i] = median(sample(dyadsObsSelect$distance, replace = T))
}

mean(boot)
quantile(boot, c(0.025, 0.975))




#-----------------------------------------------------------------------------#
# Distribution of offspring-father distances Thuringia
#-----------------------------------------------------------------------------#
# The expected random distribution
dyadsH0 = read.table("AssignationThu/outputs/dyadsH0.txt", header=TRUE)

# The combined 8/8 and 6/8 dataset of inferred paternities AFTER filtering
dyadsObsSelect = read.table("AssignationThu/outputs/dyadsObsSelect.txt", header=TRUE)

# The Kolmogorov-Smirnov test between the
# two distributions of offspring-father distances (expected and observed)
ks.test(dyadsH0$distance, dyadsObsSelect$distance)

# The mean mating dispersal distance
# Bootstrap implies some randomness in the results
nboot = 10000
boot = numeric(nboot)

for (i in 1:nboot) {
  boot[i] = mean(sample(dyadsObsSelect$distance, replace = T))
}

mean(boot)
quantile(boot, c(0.025, 0.975))


# The median mating dispersal distance
for (i in 1:nboot) {
  boot[i] = median(sample(dyadsObsSelect$distance, replace = T))
}

mean(boot)
quantile(boot, c(0.025, 0.975))



#=============================================================================#
# Natal dispersal ----
#=============================================================================#

# Assignment to natal colonies with STRUCTURE indicated a high proportion of natal dispersers. Of 17
# fathers in Picardy, eight were inferred as natal dispersers, yet their assignment probabilities were low with
# a mean probability of 0.13. There were 11 dispersing fathers in Thuringia (of a total of 18 fathers) with a
# mean assignment probability of 0.52. All sampling chains of STRUCTURE globally converged to stationarity
# and Geweke Z-scores of the selected replicate for Picardy were -0.45 (α) and 0.44 (Log-Likelihood). Geweke
# Z-scores for Thuringia were 0.06 (α) and 0.84 (Log-Likelihood).


#-----------------------------------------------------------------------------#
# Natal dispersal in Picardy
#-----------------------------------------------------------------------------#
PopAssignPic = read.table("Tables/MigrantFathersPic.txt", header = TRUE)

# Number of natal dispersers
nrow(PopAssignPic)
table(PopAssignPic$nataldisp)

# Mean assignment probability
mean(PopAssignPic$probancestry)

##### DIAGNOSTICS FOR CONVERGENCE OF SINGLE CHAINS
tracePic = read.table("STRUCTURE/output/Pic/trace.txt", header = TRUE)

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

#------------------------------------------------------------------------#
# Convergence tests in Picardy
convergencePic = data.frame(Run = c(1:10), Alpha = NA, Ln.Lik = NA)

for (i in convergencePic$Run) {
  mcmcPic = as.mcmc(tracePic[(tracePic$Run == i & !is.na(tracePic$Ln.Like)),c(3, 21)], start = 100, stop = 2000000, thin = 100)
  convergencePic$Alpha[i] = geweke.diag(mcmcPic)$z[1]
  convergencePic$Ln.Lik[i] = geweke.diag(mcmcPic)$z[2]
}
# Which is the best under Geweke diagnostic?
convergencePic
# which(abs(convergencePic$Alpha) < 0.6 & abs(convergencePic$Ln.Lik) < 0.6)

# Selected Run = 2
i = 2
convergencePic[2,]
mcmcPic = as.mcmc(tracePic[(tracePic$Run == i & !is.na(tracePic$Ln.Like)),c(3, 21)], start = 100, stop = 2000000, thin = 100)
geweke.diag(mcmcPic)
raftery.diag(mcmcPic)
traceplot(mcmcPic, lty = 1)






#-----------------------------------------------------------------------------#
# Natal dispersal in Picardy
#-----------------------------------------------------------------------------#
PopAssignPic = read.table("Tables/MigrantFathersThu.txt", header = TRUE)

# Number of natal dispersers
nrow(PopAssignPic)
table(PopAssignPic$nataldisp)

# Mean assignment probability
mean(PopAssignPic$probancestry)

##### DIAGNOSTICS FOR CONVERGENCE OF SINGLE CHAINS
traceThu = read.table("STRUCTURE/output/Thu/trace.txt", header = TRUE)

# see Hamra et al. 2013
# Trace + autocorr + density plots
# Visualize traces, with and without discarding burnin of 5,000,000
traceThu$Run = as.factor(traceThu$Run)
ggplot(data = traceThu, aes(x = Rep, y = Ln.Like)) + geom_line(aes(colour = Run))
# Remove NA from burnin
ggplot(data = traceThu[which(!is.na(traceThu$Ln.Like)),], aes(x = Rep, y = Ln.Like)) + geom_line(aes(colour = Run))
# Overall convergence attained for burnin of 1,000,000 and MCMC sampling of 1,000,000
# Yet a problem of local minima with departure from the chain before retunring to a stable chain.
# Besides, check convergence of alpha parameter
ggplot(data = traceThu[which(!is.na(traceThu$Alpha)),], aes(x = Rep, y = Alpha)) + geom_line(aes(colour = Run))
# Convergence for parameters seems globally achievable, yet with local minima

#------------------------------------------------------------------------#
# Convergence tests in Picardy
convergenceThu = data.frame(Run = c(1:10), Alpha = NA, Ln.Lik = NA)

for (i in convergenceThu$Run) {
  mcmcPic = as.mcmc(traceThu[(traceThu$Run == i & !is.na(traceThu$Ln.Like)),c(3, 21)], start = 100, stop = 2000000, thin = 100)
  convergenceThu$Alpha[i] = geweke.diag(mcmcPic)$z[1]
  convergenceThu$Ln.Lik[i] = geweke.diag(mcmcPic)$z[2]
}
# Which is the best under Geweke diagnostic?
convergenceThu
# which(abs(convergenceThu$Alpha) < 0.6 & abs(convergenceThu$Ln.Lik) < 0.6)

# Selected Run = 7
i = 7
convergenceThu[7,]
mcmcPic = as.mcmc(traceThu[(traceThu$Run == i & !is.na(traceThu$Ln.Like)),c(3, 21)], start = 100, stop = 2000000, thin = 100)
geweke.diag(mcmcPic)
raftery.diag(mcmcPic)
traceplot(mcmcPic, lty = 1)






#=============================================================================#
# Gametic dispersal ----
#=============================================================================#

# The mean gametic dispersal distance (i.e. pairwise distance between the natal colony and the colony where
# the offspring was born resulting from both natal and mating dispersal) was slightly higher than for mating
# dispersal alone: 17.6 km (IC95 = [14.4; 20.7], 10,000 bootstraps, n = 82) in France and 22.2 km (IC95 = [15.2;
# 29.6], 10,000 bootstraps, n = 62) in Germany. However, as for mating dispersal, the median gametic dispersal
# distances were lower than the mean (16.2 km, IC95 = [10.7, 18.9] and 6.9 km, IC95 = [4.1, 8.7], respectively),
# revealing a strong skewness of the distribution of distances. 


#-----------------------------------------------------------------------------#
# Distribution of offspring-father distances Picardy
#-----------------------------------------------------------------------------#
# The combined 8/8 and 6/8 dataset of inferred paternities AFTER filtering
dyadsObsSelect = read.table("AssignationPic/dyadsObsSelect_gametic.txt", header=TRUE)

# The mean mating dispersal distance
# Bootstrap implies some randomness in the results
nboot = 10000
boot = numeric(nboot)

for (i in 1:nboot) {
  boot[i] = mean(sample(dyadsObsSelect$distance, replace = T))
}

mean(boot)
quantile(boot, c(0.025, 0.975))


# The median mating dispersal distance
for (i in 1:nboot) {
  boot[i] = median(sample(dyadsObsSelect$distance, replace = T))
}

mean(boot)
quantile(boot, c(0.025, 0.975))




#-----------------------------------------------------------------------------#
# Distribution of offspring-father distances Thuringia
#-----------------------------------------------------------------------------#
# The distances natal colony of the father - offspring
dyadsObsSelect = read.table("AssignationThu/dyadsObsSelect_gametic.txt", header=TRUE)

# The mean mating dispersal distance
# Bootstrap implies some randomness in the results
nboot = 10000
boot = numeric(nboot)

for (i in 1:nboot) {
  boot[i] = mean(sample(dyadsObsSelect$distance, replace = T))
}

mean(boot)
quantile(boot, c(0.025, 0.975))


# The median mating dispersal distance
for (i in 1:nboot) {
  boot[i] = median(sample(dyadsObsSelect$distance, replace = T))
}

mean(boot)
quantile(boot, c(0.025, 0.975))






# Upon the individual distances of gametic disper-
# sal, model comparison based on AICc/BIC showed that hurdle models consistently outperformed all single
# continuous distributions (Table 1). When the conditional probability of mating inside/outside the colony and
# dispersal distances outside the colony were modelled separately in the hurdle model, the best fitting distribu-
# tion in Picardy was the binomial + Weibull (Fig. S1). The proportion of gametic dispersal outside the natal
# colony, estimated by the binomial parameter, was 0.78 (IC95 = [0.68; 0.87], 1,000 bootstraps). The shape
# and scale parameters of the fitted Weibull were 1.90 (IC95 = [1.54; 2.37], 1,000 bootstraps) and 25.3 (IC95
# = [20.9; 29.4], 1,000 bootstraps). In Thuringia, hurdle models outperformed all other ones and Log-Normal
# was the best fit, yet with very similar AICc/BIC values for different fat-tailed families (∆ AICc/BIC ≤ 3;
# Table 1; Fig.S1). We selected the Weibull model for comparison between Picardy and Thuringia (Fig. 3)
# and compared it to the negative exponential dispersal kernel inferred by Lehnen et al. (2021) on a large
# European meta-population (mean distance = 16.77 km). The binomial parameter in Thuringia was slightly
# lower than Picardy, at 0.70 (IC95 = [0.56; 0.82], 1,000 bootstraps). The fitted Weibull had a shape of 0.99
# (IC95 = [0.87; 1.16], 1,000 bootstraps) and a scale of 31.8 (IC95 = [22.4; 42.5], 1,000 bootstraps).

# Need the distancesDist.txt output file as dataset
distPic = read.table("AssignationPic/outputs/distancesDistPic_gametic.txt", header = TRUE)
distThu = read.table("AssignationThu/outputs/distancesDistThu_gametic.txt", header = TRUE)

dataPic = distPic[which(distPic$dist == "ObsSelect"), 1]
plotdist(distPic[which(distPic$dist == "ObsSelect"), 1], histo = TRUE, demp = TRUE)
descdist(distPic[which(distPic$dist == "ObsSelect"), 1], discrete=FALSE, boot=1000)



#-----------------------------------------------------------------------------#
# Table 1 ----
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# Picardy
#-----------------------------------------------------------------------------#
Table1a = data.frame(model = c("Gaussian",
                              "Negative exponential",
                              "Exponential power",
                              "Weibull",
                              "Gamma"),
                    n_parameters = c(2, 1, 3, 2, 2),
                    Thuringia_AIC = NA,
                    Thuringia_BIC = NA,
                    Picardy_AIC = NA,
                    Picardy_BIC = NA)


# Eligible functions are continuous, with a random variable X defined on 0 to +Inf
# Gaussian and exponential, usually poorer fits, are used as references
fit_e = fitdist(dataPic, distr = "exp")
summary(fit_e)
fit_gauss = fitdist(dataPic, distr = "norm")
summary(fit_gauss)
# Based on previous graphical representation, we try to fit weibull and gamma distributions
fit_w  = fitdist(dataPic, distr = "weibull", lower = c(0,0), start = list(shape = 1, scale = 1))
summary(fit_w)
fit_g  = fitdist(dataPic, distr = "gamma", lower = c(0,0), start = list(shape = 1, rate = 1))
summary(fit_g)
# Weibull and Gamma not fitting !
# They are defined on values > 0

# Exponential Power can be loaded from the package 'gnorm' described here: https://cran.r-project.org/web/packages/gnorm/vignettes/gnormUse.html
# Generalized Normal Distribution is the other name of the exponential power
fit_expower = fitdist(dataPic, distr = "gnorm", lower = c(0,0,0), start = list(mu = 0.5, alpha = 1, beta = 1))
summary(fit_expower)



# LogLikelihood
fit_gauss$loglik
fit_e$loglik
fit_expower$loglik
fit_w$loglik
fit_g$loglik
# AIC values compared for model choice
fit_gauss$aic
fit_e$aic
fit_expower$aic
fit_w$aic
fit_g$aic
# BIC values
fit_gauss$bic
fit_e$bic
fit_expower$bic
fit_w$bic
fit_g$bic
# Based on AIC values, Log normal, Weibull and Gamma fit the best the data (AIC ~ -301, -303 and -330)
# The Gaussian, Exponential and Exponential power are not adequate (AIC = 667, 543 and 568) and poorly fit graphically

Table1a$Picardy_AIC = c(fit_gauss$aic,
                         fit_e$aic,
                         fit_expower$aic,
                         fit_w$aic,
                         fit_g$aic)

Table1a$Picardy_BIC = c(fit_gauss$bic,
                         fit_e$bic,
                         fit_expower$bic,
                         fit_w$bic,
                         fit_g$bic)



#-----------------------------------------------------------------------------#
# Thuringia
#-----------------------------------------------------------------------------#
dataThu = distThu[which(distThu$dist == "ObsSelect"), 1]
plotdist(distThu[which(distThu$dist == "ObsSelect"), 1], histo = TRUE, demp = TRUE)
descdist(distThu[which(distThu$dist == "ObsSelect"), 1], discrete=FALSE, boot=1000)


# Eligible functions are continuous, with a random variable X defined on 0 to +Inf
# Gaussian and exponential, usually poorer fits, are used as references
fit_e = fitdist(dataThu, distr = "exp")
summary(fit_e)
fit_gauss = fitdist(dataThu, distr = "norm")
summary(fit_gauss)
# Based on previous graphical representation, we try to fit weibull and gamma distributions
fit_w  = fitdist(dataThu, distr = "weibull", lower = c(0,0), start = list(shape = 1, scale = 1))
summary(fit_w)
fit_g  = fitdist(dataThu, distr = "gamma", lower = c(0,0), start = list(shape = 1, rate = 1))
summary(fit_g)
# Weibull and Gamma not fitting !
# They are defined on values > 0

# Exponential Power can be loaded from the package 'gnorm' described here: https://cran.r-project.org/web/packages/gnorm/vignettes/gnormUse.html
# Generalized Normal Distribution is the other name of the exponential power
fit_expower = fitdist(dataThu, distr = "gnorm", lower = c(0,0,0), start = list(mu = 0.5, alpha = 1, beta = 1))
summary(fit_expower)



# LogLikelihood
fit_gauss$loglik
fit_e$loglik
fit_expower$loglik
fit_w$loglik
fit_g$loglik
# AIC values compared for model choice
fit_gauss$aic
fit_e$aic
fit_expower$aic
fit_w$aic
fit_g$aic
# BIC values
fit_gauss$bic
fit_e$bic
fit_expower$bic
fit_w$bic
fit_g$bic
# Based on AIC values, Log normal, Weibull and Gamma fit the best the data (AIC ~ -301, -303 and -330)
# The Gaussian, Exponential and Exponential power are not adequate (AIC = 667, 543 and 568) and poorly fit graphically

Table1a$Thuringia_AIC = c(fit_gauss$aic,
                         fit_e$aic,
                         fit_expower$aic,
                         fit_w$aic,
                         fit_g$aic)

Table1a$Thuringia_BIC = c(fit_gauss$bic,
                        fit_e$bic,
                        fit_expower$bic,
                        fit_w$bic,
                        fit_g$bic)










#-----------------------------------------------------------------------------#
# HURDLE MODELS: Modelling two different processes... ----
#-----------------------------------------------------------------------------#
# This is a two-parts model, with a binomial process of intra/extra colonial mating (i.e. probability of a 0)
# and a continuous distribution of distances
# Some informations:
# https://seananderson.ca/2014/05/18/gamma-hurdle/
# https://stats.stackexchange.com/questions/320924/is-a-hurdle-model-really-one-model-or-just-two-separate-sequential-models
# https://stackoverflow.com/questions/36312771/making-zero-inflated-or-hurdle-model-with-r


Table1b = data.frame(model = c("Binomial + Gaussian",
                               "Binomial + Negative exponential",
                               "Binomial + Exponential power",
                               "Binomial + Weibull",
                               "Binomial + Gamma",
                               "Binomial + Log-normal"),
                     n_parameters = c(3, 2, 4, 3, 3, 3),
                     Thuringia_AIC = NA,
                     Thuringia_BIC = NA,
                     Picardy_AIC = NA,
                     Picardy_BIC = NA)


#----------------------------------------------------------------------------#
# Picardy
#----------------------------------------------------------------------------#
# Fitting the intercept, the probability of intra-colonial paternity (i.e. null value)
# with a binomial law
non_zero = ifelse(dataPic > 0, 0, 1)
(fit_bin = glm(non_zero ~ 1, family = binomial(link = logit))) # binomial: probability of 0 value
(bin_coef = plogis(coef(fit_bin)[[1]])) # prob = 0.51
summary(fit_bin)
logLik(fit_bin)
AIC(fit_bin)
BIC(fit_bin)
# The intercept is the y value from which the fitted distribution of non-null values must start

# Estimate CI for binomial coefficient with bootstrap
paramPic = data.frame(parameter = character(), mean = numeric(), lower = numeric(), upper = numeric())
k = 3 # Number of parameters estimated
threshold = 0.025
boot = numeric()
nboot = 1000
for (i in 1:nboot) {
  print(i)
  fit_bin_boot = glm(sample(non_zero, replace = TRUE) ~ 1, family = binomial(link = logit))
  boot[i] = (bin_coef = plogis(coef(fit_bin_boot)[[1]]))
}
mean(boot)
quantile(boot, threshold)
quantile(boot, 1 - threshold)

# The probability of mating outside the colony:
1 - mean(boot)
1 - quantile(boot, threshold)
1 - quantile(boot, 1 - threshold)


# THEN... Fitting only the non-null values
# Making a zero-truncated dataset
non_zero = ifelse(dataPic > 0, 1, 0)
dataPic_pos = subset(dataPic, non_zero == 1)

# Gaussian and exponential as references
fit_e = fitdist(dataPic_pos, distr = "exp")
summary(fit_e)
fit_gauss = fitdist(dataPic_pos, distr = "norm")
summary(fit_gauss)
# Based on previous graphical representation, we try to fit weibull and gamma distributions
fit_w  = fitdist(dataPic_pos, distr = "weibull")
summary(fit_w)
fit_g  = fitdist(dataPic_pos, distr = "gamma")
summary(fit_g)
# A lognormal law is defined by x>0
fit_ln  = fitdist(dataPic_pos, distr = "lnorm")
summary(fit_ln)
fit_expower = fitdist(dataPic_pos, distr = "gnorm", lower = c(0,0,0), start = list(mu = 0.5, alpha = 1, beta = 1))
summary(fit_expower)
# Exp power fails to fit


# Bootstrapping CI for parameters of the chosen model, here Weibull
# Shape and scale
boot_shape = numeric()
boot_scale = numeric()
nboot = 1000
for (i in 1:nboot) {
  print(i)
  fit_w_boot  = fitdist(sample(dataPic_pos, replace = TRUE), distr = "weibull")
  boot_shape[i] = fit_w_boot$estimate[1]
  boot_scale[i] = fit_w_boot$estimate[2]
}
mean(boot_shape)
quantile(boot_shape, 0.025)
quantile(boot_shape, 0.975)
mean(boot_scale)
quantile(boot_scale, 0.025)
quantile(boot_scale, 0.975)




# AIC values compared for model choice
fit_e$aic
fit_gauss$aic
fit_w$aic
fit_g$aic
fit_ln$aic
fit_expower$aic


# Full model bin + continuous dist
# Recompute logLik, AIC and BIC

# LogLik is the sum of the two LogLik estimated separately
(LLik_gauss = logLik(fit_bin) + logLik(fit_gauss))
(LLik_e = logLik(fit_bin) + logLik(fit_e))
(LLik_expower = logLik(fit_bin) + logLik(fit_expower))
(LLik_w = logLik(fit_bin) + logLik(fit_w))
(LLik_g = logLik(fit_bin) + logLik(fit_g))
(LLik_ln = logLik(fit_bin) + logLik(fit_ln))

# AIC : choose the model which minimises 2p – 2ln(L)
# where p is number of parameters
# Function AIC
# AIC = 2p – 2ln(L)
# (AIC_gauss = AIC(LLik_gauss))
(AIC_gauss = 2*3-2*LLik_gauss)
# (AIC_e = AIC(LLik_e))
(AIC_e = 2*2-2*LLik_e)
# (AIC_expower = AIC(LLik_expower))
(AIC_expower = 2*3-2*LLik_expower)
# (AIC_w = AIC(LLik_w))
(AIC_w = 2*3-2*LLik_w)
# (AIC_g = AIC(LLik_g))
(AIC_g = 2*3-2*LLik_g)
# (AIC_ln = AIC(LLik_ln))
(AIC_ln = 2*4-2*LLik_ln)


# BIC : choose the model which minimises p.ln(N) - 2ln(L)
# where p is number of parameters and N number of observations
# BIC = p.ln(N) - 2ln(L)
# (BIC_gauss = BIC(LLik_gauss))
(BIC_gauss = 3*log(length(fit_bin$y))-2*LLik_gauss)
# (BIC_e = BIC(LLik_e))
(BIC_e = 2*log(length(fit_bin$y))-2*LLik_e)
# (BIC_expower = BIC(LLik_expower))
(BIC_expower = 4*log(length(fit_bin$y))-2*LLik_expower)
# (BIC_w = BIC(LLik_w))
(BIC_w = 3*log(length(fit_bin$y))-2*LLik_w)
# (BIC_g = BIC(LLik_g))
(BIC_g = 3*log(length(fit_bin$y))-2*LLik_g)
# (BIC_ln = BIC(LLik_ln))
(BIC_ln = 3*log(length(fit_bin$y))-2*LLik_ln)


Table1b$Picardy_AIC = c(AIC_gauss,
                        AIC_e,
                        AIC_expower,
                        AIC_w,
                        AIC_g,
                        AIC_ln)
Table1b$Picardy_BIC = c(BIC_gauss,
                        BIC_e,
                        BIC_expower,
                        BIC_w,
                        BIC_g,
                        BIC_ln)



#----------------------------------------------------------------------------#
# Picardy
#----------------------------------------------------------------------------#
# Fitting the intercept, the probability of intra-colonial paternity (i.e. null value)
# with a binomial law
non_zero = ifelse(dataThu > 0, 0, 1)
(fit_bin = glm(non_zero ~ 1, family = binomial(link = logit))) # binomial: probability of 0 value
(bin_coef = plogis(coef(fit_bin)[[1]])) # prob = 0.51
summary(fit_bin)
logLik(fit_bin)
AIC(fit_bin)
BIC(fit_bin)
# The intercept is the y value from which the fitted distribution of non-null values must start

# Estimate CI for binomial coefficient with bootstrap
paramThu = data.frame(parameter = character(), mean = numeric(), lower = numeric(), upper = numeric())
k = 3 # Number of parameters estimated
threshold = 0.025
boot = numeric()
nboot = 1000
for (i in 1:nboot) {
  print(i)
  fit_bin_boot = glm(sample(non_zero, replace = TRUE) ~ 1, family = binomial(link = logit))
  boot[i] = (bin_coef = plogis(coef(fit_bin_boot)[[1]]))
}
mean(boot)
quantile(boot, threshold)
quantile(boot, 1 - threshold)

# The probability of mating outside the colony:
1 - mean(boot)
1 - quantile(boot, threshold)
1 - quantile(boot, 1 - threshold)


# THEN... Fitting only the non-null values
# Making a zero-truncated dataset
non_zero = ifelse(dataThu > 0, 1, 0)
dataThu_pos = subset(dataThu, non_zero == 1)

# Gaussian and exponential as references
fit_e = fitdist(dataThu_pos, distr = "exp")
summary(fit_e)
fit_gauss = fitdist(dataThu_pos, distr = "norm")
summary(fit_gauss)
# Based on previous graphical representation, we try to fit weibull and gamma distributions
fit_w  = fitdist(dataThu_pos, distr = "weibull")
summary(fit_w)
fit_g  = fitdist(dataThu_pos, distr = "gamma")
summary(fit_g)
# A lognormal law is defined by x>0
fit_ln  = fitdist(dataThu_pos, distr = "lnorm")
summary(fit_ln)
fit_expower = fitdist(dataThu_pos, distr = "gnorm", lower = c(0,0,0), start = list(mu = 0.5, alpha = 1, beta = 1))
summary(fit_expower)
# Exp power fails to fit


# Bootstrapping CI for parameters of the chosen model, here Weibull
# Shape and scale
boot_shape = numeric()
boot_scale = numeric()
nboot = 1000
for (i in 1:nboot) {
  print(i)
  fit_w_boot  = fitdist(sample(dataThu_pos, replace = TRUE), distr = "weibull")
  boot_shape[i] = fit_w_boot$estimate[1]
  boot_scale[i] = fit_w_boot$estimate[2]
}
mean(boot_shape)
quantile(boot_shape, 0.025)
quantile(boot_shape, 0.975)
mean(boot_scale)
quantile(boot_scale, 0.025)
quantile(boot_scale, 0.975)




# AIC values compared for model choice
fit_e$aic
fit_gauss$aic
fit_w$aic
fit_g$aic
fit_ln$aic
fit_expower$aic


# Full model bin + continuous dist
# Recompute logLik, AIC and BIC

# LogLik is the sum of the two LogLik estimated separately
(LLik_gauss = logLik(fit_bin) + logLik(fit_gauss))
(LLik_e = logLik(fit_bin) + logLik(fit_e))
(LLik_expower = logLik(fit_bin) + logLik(fit_expower))
(LLik_w = logLik(fit_bin) + logLik(fit_w))
(LLik_g = logLik(fit_bin) + logLik(fit_g))
(LLik_ln = logLik(fit_bin) + logLik(fit_ln))

# AIC : choose the model which minimises 2p – 2ln(L)
# where p is number of parameters
# Function AIC
# AIC = 2p – 2ln(L)
# (AIC_gauss = AIC(LLik_gauss))
(AIC_gauss = 2*3-2*LLik_gauss)
# (AIC_e = AIC(LLik_e))
(AIC_e = 2*2-2*LLik_e)
# (AIC_expower = AIC(LLik_expower))
(AIC_expower = 2*3-2*LLik_expower)
# (AIC_w = AIC(LLik_w))
(AIC_w = 2*3-2*LLik_w)
# (AIC_g = AIC(LLik_g))
(AIC_g = 2*3-2*LLik_g)
# (AIC_ln = AIC(LLik_ln))
(AIC_ln = 2*4-2*LLik_ln)


# BIC : choose the model which minimises p.ln(N) - 2ln(L)
# where p is number of parameters and N number of observations
# BIC = p.ln(N) - 2ln(L)
# (BIC_gauss = BIC(LLik_gauss))
(BIC_gauss = 3*log(length(fit_bin$y))-2*LLik_gauss)
# (BIC_e = BIC(LLik_e))
(BIC_e = 2*log(length(fit_bin$y))-2*LLik_e)
# (BIC_expower = BIC(LLik_expower))
(BIC_expower = 4*log(length(fit_bin$y))-2*LLik_expower)
# (BIC_w = BIC(LLik_w))
(BIC_w = 3*log(length(fit_bin$y))-2*LLik_w)
# (BIC_g = BIC(LLik_g))
(BIC_g = 3*log(length(fit_bin$y))-2*LLik_g)
# (BIC_ln = BIC(LLik_ln))
(BIC_ln = 3*log(length(fit_bin$y))-2*LLik_ln)



Table1b$Thuringia_AIC = c(AIC_gauss,
                        AIC_e,
                        AIC_expower,
                        AIC_w,
                        AIC_g,
                        AIC_ln)
Table1b$Thuringia_BIC = c(BIC_gauss,
                        BIC_e,
                        BIC_expower,
                        BIC_w,
                        BIC_g,
                        BIC_ln)



Table1 = rbind(Table1a, Table1b)

write.table(Table1, "Tables/Table1.tsv", sep = "\t",
            col.names = , row.names = F)


#=============================================================================#
# END
#=============================================================================#