#===============================#
#            COLONY
#     Effective dispersal kernel
#     Estimate parameters
#     on a model-based approach
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


#---------------------------------
# Import data
#---------------------------------
# Need the distancesDist.txt output file as dataset
distPic = read.table("AssignationPic/outputs/distancesDistPic_gametic.txt", header = TRUE)
distThu = read.table("AssignationThu/outputs/distancesDistThu_gametic.txt", header = TRUE)

#=================================#
# DISTRIBUTION CHOICE
#=================================
# Continuous distribution, with a lot of zeros
# Properties of the random variable (e.g. range, expected values) suggest:
# - that it is not a beta, which is defined in [0,1]
# - it is not count data, hence not Poisson or binomial
# - it could be a gamma or a Weibull
# - it is closer to a lognormal, but all values must be positive with a lognormal

# Classical dispersal probability functions were Gaussian and (Negative) Exponential
# 2 models that should be used as reference against more realistic models

# General mixture dispersal kernel
# Well suited when 2 patterns suspected, for short and long distances
# Clobert et al; 2012

# Although, we can characterize the empirical distribution by its high proportion of null values.
# Maybe we should test some zero inflated distributions
# Zero-inflated Discrete Weibull (survival analysis, in times)
# Zero-inflated Gamma

#---------------------------------#
# Descriptive statistics
#---------------------------------
dataPic = distPic[which(distPic$dist == "ObsSelect"), 1]
plotdist(distPic[which(distPic$dist == "ObsSelect"), 1], histo = TRUE, demp = TRUE)
descdist(distPic[which(distPic$dist == "ObsSelect"), 1], discrete=FALSE, boot=1000)

# Skewness and kurtosis on the Cullen Frey graph
# For kurtosis, it is close to logistic, beta, lognormal, gamma and weibull
# For skewness, it is close to beta, lognormal, gamma and weibull but far from logistic
# Bootstrap gives the best expectation as beta, and a parallel tendency to logn, gamma and weibull
png("Figures/Pic Density GAMETIC.png",width=1000,height=500)
plotdist(distPic[which(distPic$dist == "ObsSelect"), 1], histo = TRUE, demp = TRUE)
dev.off()
png("Figures/Pic Theoretical expectations GAMETIC.png",width=1000,height=1000)
descdist(distPic[which(distPic$dist == "ObsSelect"), 1], discrete=FALSE, boot=1000)
dev.off()
# summary statistics
# ------
#   min:  0   max:  55.59229 
# median:  0.6885013 
# mean:  9.96432 
# estimated sd:  13.92085 
# estimated skewness:  1.473765 
# estimated kurtosis:  4.427601 

#---------------------------------#
# Fitting - First try with continuous functions
#---------------------------------
?fitdist
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

# Exponential power well suited for kernel shape comparisons
# Power-law (fat tail) and exponential power are also important distribution functions in dispersal

# Exponential Power can be loaded from the package 'gnorm' described here: https://cran.r-project.org/web/packages/gnorm/vignettes/gnormUse.html
# Generalized Normal Distribution is the other name of the exponential power
fit_expower = fitdist(dataPic, distr = "gnorm", lower = c(0,0,0), start = list(mu = 0.5, alpha = 1, beta = 1))
summary(fit_expower)

# Half normal loaded from 'fdrtool'
library(fdrtool)
fit_halfnorm = fitdist(dataPic, distr = "halfnorm", start = list(theta=sqrt(pi/2)))
summary(fit_halfnorm)


# An error is catched 'Not defined for negative values of alpha and/or beta.'
# Hence, no estimates of Std Error and Correlation matrix

# Power-law
# with custom power-law distribution
# source("Sources/powerlaw.R")
# # fit_power = fitdist(dataPic, distr = "powerlaw", lower = c(0,1), start = list(xmin = 0.1, alpha = 2))
# fit_power = fitdist(dataPic, distr = "powerlaw", method = "qme", probs = c(1/3,2/3), lower = c(0,0), start = list(xmin = 0.5, alpha = 0.5))
# summary(fit_power)
# Not working wit the 'mle' method, so use the quantile method 'qme' which does not provide log likelood
# And very poor fit to the data
# Hence, discarding the power law from our approach

# A lognormal law is not possible with full distribution of values, since log-normal is defined by X > 0
# fit_ln  = fitdist(dataPic, distr = "lnorm")
# summary(fit_ln)

# Our first approach is to consider that we have a continuous distribution,
# Hence, null-valuies could be considered as very short dispersal distances,
# allowing us to perform a slight transformation of data
# for performance of optimisation algorithms only
# Then try with transformed data:
# So, trying to fit to data with a small noise (10^-6)
    # dataPicTransformed = dataPic+0.000001
    # # Gaussian and exponential as references
    # fit_e = fitdist(dataPicTransformed, distr = "exp")
    # summary(fit_e)
    # fit_gauss = fitdist(dataPicTransformed, distr = "norm")
    # summary(fit_gauss)
    # # Based on previous graphical representation, we try to fit weibull and gamma distributions
    # fit_w  = fitdist(dataPicTransformed, distr = "weibull")
    # summary(fit_w)
    # fit_g  = fitdist(dataPicTransformed, distr = "gamma")
    # summary(fit_g)
    # # A lognormal law is defined by x>0
    # fit_ln  = fitdist(dataPicTransformed, distr = "lnorm")
    # summary(fit_ln)
    # library(gnorm)
    # fit_expower = fitdist(dataPicTransformed, distr = "gnorm", lower = c(0,0,0), start = list(mu = 0.5, alpha = 1, beta = 1))
    # summary(fit_expower)

#-------------------------
# Goodness of fit plots
# Check which distribution law fits best the data
# The Q-Q plot emphasizes the lack-of-fit at the distribution tails
# while the P-P plot emphasizes the lack-of-fit at the distribution center.
# par(mfrow=c(2,2))
# plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Exponential Power")
plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Exponential Power")
denscomp(list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend, ylim = c(0, 0.2))
cdfcomp (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
# par(mfrow=c(1,1))

png("Figures/Pic Fitdist AllValues GAMETIC.png",width=1000,height=1000)
par(mfrow=c(2,2))
plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Exponential Power")
denscomp(list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend, ylim = c(0, 0.2))
cdfcomp (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
par(mfrow=c(1,1))
dev.off()

#------------------------
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



# It emphasizes the need to deal with zero values
# - zero-inflated or zero-adjusted distributions
# - fitting only effective dispersal (i.e. distances > 0), meaning a new distribution
# - fitting a hurdle model treating zeros as a different distribution than the non-zero values (i.e. model zeros separately)

# https://stats.stackexchange.com/questions/279273/zero-inflated-distributions-what-are-they-really
# https://seananderson.ca/2014/05/18/gamma-hurdle/




          # #=================================#
          # # Test of a zero-adjusted distribution in fitdistR
          # #=================================
          # # One advantage of zero-adjusted dist is to compare them to non zero-adjusted on the same dataset (AIC can be compared)
          # 
          # # Zero adjusted Gamma distribution
          # # mu, location parameter
          # # sigma, scale parameter
          # # nu, the probability of zero, must be between 0 and 1
          # fit_ZAGA = fitdist(dataPic, distr = "ZAGA", start = list(mu = 2, sigma = 2, nu = 0.5))
          # # fit_gamlss = fitdist(dataPic, distr = dZAGA(x = freqPic), lower = c(0,0,0), upper = c(100,100,1), start = list(mu = 2, sigma = 2, nu = 0.5))
          # summary(fit_ZAGA)
          # plot.legend = c("ZAGA")
          # denscomp(list(fit_ZAGA), legendtext = plot.legend, ylim = c(0, 1))
          # cdfcomp (list(fit_ZAGA), legendtext = plot.legend)
          # qqcomp  (list(fit_ZAGA), legendtext = plot.legend)
          # ppcomp  (list(fit_ZAGA), legendtext = plot.legend)
          # 
          # # The zero-adjusted Gamma distribution works very well with the data
          # # Comparing with other distributions...
          # 
          # # dZALG, the zero-adjusted logarithmic distribution
          # # mu, positive means, between 0 and 1
          # # sigma, probabilities at 0, between 0 and 1
          # fit_ZALG = fitdist(dataPic, distr = "ZALG", upper = c(1,1), lower = c(0,0), start = list(mu = 0.5, sigma = 0.5))
          # summary(fit_ZALG)
          # 
          # # dZAIG, the zero adjusted Inverse Gaussian distribution
          # # mu, location parameter
          # # sigma, scale parameter
          # # nu, the probability of zero, must be between 0 and 1
          # fit_ZAIG = fitdist(dataPic, distr = "ZAIG", start = list(mu = 1, sigma = 1, nu = 0.5))
          # summary(fit_ZAIG)
          # 
          # #-------------------------
          # # Goodness of fit plots
          # # Check which distribution law fits best the data
          # # The Q-Q plot emphasizes the lack-of-fit at the distribution tails
          # # while the P-P plot emphasizes the lack-of-fit at the distribution center.
          # # par(mfrow=c(2,2))
          # plot.legend = c("ZAGA", "ZALG", "ZAIG")
          # denscomp(list(fit_ZAGA, fit_ZALG, fit_ZAIG), legendtext = plot.legend, ylim = c(0, 0.2))
          # cdfcomp (list(fit_ZAGA, fit_ZALG, fit_ZAIG), legendtext = plot.legend)
          # qqcomp  (list(fit_ZAGA, fit_ZALG, fit_ZAIG), legendtext = plot.legend)
          # ppcomp  (list(fit_ZAGA, fit_ZALG, fit_ZAIG), legendtext = plot.legend)
          # # par(mfrow=c(1,1))
          # 
          # 
          # png("Figures/Pic Fitdist ZeroAdjusted.png",width=1000,height=1000)
          # par(mfrow=c(2,2))
          # plot.legend = c("ZAGA", "ZALG", "ZAIG")
          # denscomp(list(fit_ZAGA, fit_ZALG, fit_ZAIG), legendtext = plot.legend, ylim = c(0, 0.2))
          # cdfcomp (list(fit_ZAGA, fit_ZALG, fit_ZAIG), legendtext = plot.legend)
          # qqcomp  (list(fit_ZAGA, fit_ZALG, fit_ZAIG), legendtext = plot.legend)
          # ppcomp  (list(fit_ZAGA, fit_ZALG, fit_ZAIG), legendtext = plot.legend)
          # par(mfrow=c(1,1))
          # dev.off()
          # 
          # #------------------------
          # # AIC values compared for model choice
          # fit_ZAGA$aic
          # fit_ZALG$aic
          # fit_ZAIG$aic
          # # ZAGA is the best fit of all distributions (AIC = 447),
          # # and 2nd best fits are ZAIG (AIC = 471) and ZALG (AIC = 483)
          # # while the best fit for non-adjusted dist are exp (AIC = 543) and Weibull (AIC = 545)
          # 
          # 
          # 
          # 
          # 
          # #=================================#
          # # Constructing Zero-adjusted functions (Weibull, exponential, Gaussian)
          # #=================================
          # # # (4) defining your own distribution functions, here for the Gumbel distribution
          # # # for other distributions, see the CRAN task view 
          # # # dedicated to probability distributions
          # # #
          # # 
          # # dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
          # # pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
          # # qgumbel <- function(p, a, b) a-b*log(-log(p))
          # # 
          # # fitgumbel <- fitdist(serving, "gumbel", start=list(a=10, b=10))
          # # summary(fitgumbel)
          # # plot(fitgumbel)
          # 



#=================================#
# HURDLE MODELS: Modelling two different processes...
#=================================
# This is a two-parts model, with a binomial process of intra/extra colonial mating (i.e. probability of a 0)
# and a continuous distribution of distances
# Some informations:
# https://seananderson.ca/2014/05/18/gamma-hurdle/
# https://stats.stackexchange.com/questions/320924/is-a-hurdle-model-really-one-model-or-just-two-separate-sequential-models
# https://stackoverflow.com/questions/36312771/making-zero-inflated-or-hurdle-model-with-r

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

paramPic$parameter = as.character(paramPic$parameter)
paramPic[1,] = c(as.character("Binomial"), mean(boot), quantile(boot, 0.025), quantile(boot, 0.975))
paramPic

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
# Half normal loaded from 'fdrtool'
library(fdrtool)
fit_halfnorm = fitdist(dataPic_pos, distr = "halfnorm", start = list(theta=sqrt(pi/2)))
summary(fit_halfnorm)


# Parameters saved for figures
paramPicKernel = data.frame(dist = c("exp", "norm", "weibull", "gamma", "lnorm", "expower"),
                            param1 = unlist(c(fit_e$estimate, fit_gauss$estimate[1], fit_w$estimate[1], fit_g$estimate[1], fit_ln$estimate[1], fit_expower$estimate[1])),
                            param2 = unlist(c(NA, fit_gauss$estimate[2], fit_w$estimate[2], fit_g$estimate[2], fit_ln$estimate[2], fit_expower$estimate[2])),
                            param3 = unlist(c(NA, NA, NA, NA, NA, fit_expower$estimate[3])))
write.table(paramPicKernel, file = "Tables/paramPicKernel GAMETIC.txt", quote = F, col.names = T, row.names = F)

# LogLik and AIC can be the sum of the LogLik/AIC of the binomial model and the conditional model (zero-truncated) estimated separately. (McDowell, 2003)

# Bootstrapping CI for parameters of the chosen model, here Weibull
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

paramPic[2,] = c("Shape", mean(boot_shape), quantile(boot_shape, threshold), quantile(boot_shape, 1 - threshold))
paramPic[3,] = c("Scale", mean(boot_scale), quantile(boot_scale, threshold), quantile(boot_scale, 1 -threshold))
paramPic
write.table(paramPic, file = "Tables/ParamPic GAMETIC.txt", quote = F, col.names = T, row.names = F)

#-------------------------
# Goodness of fit plots
# Check which distribution law fits best the data
# The Q-Q plot emphasizes the lack-of-fit at the distribution tails
# while the P-P plot emphasizes the lack-of-fit at the distribution center.
# par(mfrow=c(2,2))
# plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Exponential Power")
plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Log Normal", "Exponential Power")
denscomp(list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend, ylim = c(0, 0.2))
cdfcomp (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
# par(mfrow=c(1,1))

png("Figures/Pic Fitdist PosValues GAMETIC.png",width=1000,height=1000)
par(mfrow=c(2,2))
plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Log Normal", "Exponential Power")
denscomp(list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend, ylim = c(0, 0.2))
cdfcomp (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
par(mfrow=c(1,1))
dev.off()

#------------------------
# AIC values compared for model choice
fit_e$aic
fit_gauss$aic
fit_w$aic
fit_g$aic
fit_ln$aic
fit_expower$aic

#------------------------
# Full model bin + continuous dist
# Recompute logLik, AIC and BIC

# LogLik is the sum of the two LogLik estimated separately
(LLik_gauss = logLik(fit_bin) + logLik(fit_gauss))
(LLik_e = logLik(fit_bin) + logLik(fit_e))
(LLik_expower = logLik(fit_bin) + logLik(fit_expower))
(LLik_w = logLik(fit_bin) + logLik(fit_w))
(LLik_g = logLik(fit_bin) + logLik(fit_g))
(LLik_ln = logLik(fit_bin) + logLik(fit_ln))

(LLik_hnorm = logLik(fit_bin) + logLik(fit_halfnorm))

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

(AIC_hnorm = 2*2-2*LLik_hnorm)

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

(BIC_hnorm = 4*log(length(fit_bin$y))-2*LLik_hnorm)



#=================================#
# MODEL VALIDITY CHECKS
#=================================

# Normality or residuals


# Variance



#=================================#
# PARAMETER ESTIMATION
#=================================

# Parameter estimation of the chosen model(s)



# Bootstrap procedure for confidence intervals



#=================================#
# FITTING QUALITY ASSESSMENT
#=================================



#=================================#
# MAKING TABLES
#=================================






########################################################################################
#
#   Thuringia
#
########################################################################################


#---------------------------------#
# Descriptive statistics
#---------------------------------
dataThu = distThu[which(distThu$dist == "ObsSelect"), 1]
plotdist(distThu[which(distThu$dist == "ObsSelect"), 1], histo = TRUE, demp = TRUE)
descdist(distThu[which(distThu$dist == "ObsSelect"), 1], discrete=FALSE, boot=1000)

#---------------------------------#
# Fitting - First try with continuous functions
#---------------------------------
fit_e = fitdist(dataThu, distr = "exp")
summary(fit_e)
fit_gauss = fitdist(dataThu, distr = "norm")
summary(fit_gauss)
fit_w  = fitdist(dataThu, distr = "weibull", lower = c(0,0), start = list(shape = 1, scale = 1))
summary(fit_w)
fit_g  = fitdist(dataThu, distr = "gamma", lower = c(0,0), start = list(shape = 1, rate = 1))
summary(fit_g)
fit_expower = fitdist(dataThu, distr = "gnorm", lower = c(0,0,0), start = list(mu = 0.5, alpha = 1, beta = 1))
summary(fit_expower)
# fit_power = fitdist(dataThu, distr = "powerlaw", method = "qme", probs = c(1/3,2/3), lower = c(0,0), start = list(xmin = 0.5, alpha = 0.5))
# summary(fit_power)

plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Exponential Power")
denscomp(list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend, ylim = c(0, 0.2))
cdfcomp (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
# par(mfrow=c(1,1))

png("Figures/Thu Fitdist AllValues GAMETIC.png",width=1000,height=1000)
par(mfrow=c(2,2))
plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Exponential Power")
denscomp(list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend, ylim = c(0, 0.2))
cdfcomp (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_expower), legendtext = plot.legend)
par(mfrow=c(1,1))
dev.off()

#------------------------
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

#=================================#
# HURDLE MODELS: Modelling two different processes...
#=================================
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

paramThu = data.frame(parameter = character(), mean = numeric(), lower = numeric(), upper = numeric())
paramThu$parameter = as.character(paramThu$parameter)
paramThu[1,] = c(as.character("Binomial"), mean(boot), quantile(boot, threshold), quantile(boot, 1 - threshold))

# THEN... Fitting only the non-null values
dataThu_pos = subset(dataThu, non_zero == 0)
hist(dataThu_pos)
fit_e = fitdist(dataThu_pos, distr = "exp")
summary(fit_e)
fit_gauss = fitdist(dataThu_pos, distr = "norm")
summary(fit_gauss)
fit_w  = fitdist(dataThu_pos, distr = "weibull")
summary(fit_w)
fit_g  = fitdist(dataThu_pos, distr = "gamma")
summary(fit_g)
fit_ln  = fitdist(dataThu_pos, distr = "lnorm")
summary(fit_ln)
fit_expower = fitdist(dataThu_pos, distr = "gnorm", start = list(mu = 0.05, alpha = 0.1, beta = 0.1))
summary(fit_expower)

# Parameters saved for figures
paramThuKernel = data.frame(dist = c("exp", "norm", "weibull", "gamma", "lnorm", "expower"),
                            param1 = unlist(c(fit_e$estimate, fit_gauss$estimate[1], fit_w$estimate[1], fit_g$estimate[1], fit_ln$estimate[1], fit_expower$estimate[1])),
                            param2 = unlist(c(NA, fit_gauss$estimate[2], fit_w$estimate[2], fit_g$estimate[2], fit_ln$estimate[2], fit_expower$estimate[2])),
                            param3 = unlist(c(NA, NA, NA, NA, NA, fit_expower$estimate[3])))
write.table(paramThuKernel, file = "Tables/ParamThuKernel GAMETIC.txt", quote = F, col.names = T, row.names = F)

# LogLik and AIC can be the sum of the LogLik/AIC of the binomial model and the conditional model (zero-truncated) estimated separately. (McDowell, 2003)

#-------------------------
# Goodness of fit plots
plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Log Normal", "Exponential Power")
denscomp(list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend, ylim = c(0, 0.2))
cdfcomp (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)

png("Figures/Thu Fitdist PosValues GAMETIC.png",width=1000,height=1000)
par(mfrow=c(2,2))
plot.legend = c("Weibull", "Gamma", "Exponential","Gaussian", "Log Normal", "Exponential Power")
denscomp(list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend, ylim = c(0, 0.2))
cdfcomp (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_e, fit_gauss, fit_ln, fit_expower), legendtext = plot.legend)
par(mfrow=c(1,1))
dev.off()

#------------------------
# Full model bin + continuous dist
# Recompute logLik, AIC and BIC
(LLik_gauss = logLik(fit_bin) + logLik(fit_gauss))
(LLik_e = logLik(fit_bin) + logLik(fit_e))
(LLik_expower = logLik(fit_bin) + logLik(fit_expower))
(LLik_w = logLik(fit_bin) + logLik(fit_w))
(LLik_g = logLik(fit_bin) + logLik(fit_g))
(LLik_ln = logLik(fit_bin) + logLik(fit_ln))

(AIC_gauss = AIC(LLik_gauss))
(AIC_e = AIC(LLik_e))
(AIC_expower = AIC(LLik_expower))
(AIC_w = AIC(LLik_w))
(AIC_g = AIC(LLik_g))
(AIC_ln = AIC(LLik_ln))

(BIC_gauss = BIC(LLik_gauss))
(BIC_e = BIC(LLik_e))
(BIC_expower = BIC(LLik_expower))
(BIC_w = BIC(LLik_w))
(BIC_g = BIC(LLik_g))
(BIC_ln = BIC(LLik_ln))


# Bootstrapping CI for parameters of the Weibull
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
quantile(boot_shape, threshold)
quantile(boot_shape, 1-threshold)
mean(boot_scale)
quantile(boot_scale, threshold)
quantile(boot_scale, 1-threshold)

# # Bootstrapping CI for parameters of the Lognormal
# boot_shape = numeric()
# boot_scale = numeric()
# nboot = 1000
# for (i in 1:nboot) {
#   print(i)
#   fit_w_boot  = fitdist(sample(dataThu_pos, replace = TRUE), distr = "weibull")
#   boot_shape[i] = fit_w_boot$estimate[1]
#   boot_scale[i] = fit_w_boot$estimate[2]
# }
# mean(boot_shape)
# quantile(boot_shape, threshold)
# quantile(boot_shape, 1 - threshold)
# mean(boot_scale)
# quantile(boot_scale, threshold)
# quantile(boot_scale, 1 - threshold)

paramThu[2,] = c("Shape", mean(boot_shape), quantile(boot_shape, threshold), quantile(boot_shape, 1 - threshold))
paramThu[3,] = c("Scale", mean(boot_scale), quantile(boot_scale, threshold), quantile(boot_scale, 1 - threshold))
paramThu
write.table(paramThu, file = "Tables/ParamThu GAMETIC.txt", quote = F, col.names = T, row.names = F)

#---------------------------------#
# END
#---------------------------------#
