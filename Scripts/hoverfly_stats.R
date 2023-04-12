##%######################################################%########
#                                                                #
####       THE IMPACT OF FARMLAND VEGETATION MANAGEMENT          ####
####      ON PROMOTING A SPECIES AND FUNCTIONALLY DIVERSE        ####
####                    HOVERFLY COMMUNITY                       ####
#### Property of a University of Edinburgh undergraduate student ####
####                    Created on 01/04/2023                    ####
#                                                                #
##%######################################################%########

# set working directory----
setwd("your_filepath")

# load libraries (install the packages if you haven't yet)----
library(tidyverse) # for efficient data manipulation and visualization, contains ggplot2, dplyr and tidyr packages
library(ggpubr) # to create and customize composite figures in unison with ggplot2
library(patchwork) # to help create composite figures
library(plyr) # to aid in data manipulation
library(permute) # needed to run vegan package
library(lattice) # needed to run vegan package
library(vegan) # for diversity analysis and ordination methods
library(lme4) # for mixed effects models
library(mblm) # for non parametric regression
library(RobustLinearReg) # for non parametric Kendall-Theil Sen Siegel regression
library(FSA) # for post hoc comparisons of Kruskal wallis test using Dunn test
library(MASS) # for negative binomial regression
library(glmmTMB) # generate zero-inflated mixed models
library(DHARMa) # allows for zero-inflation testing
library(lsmeans) # for post hoc comparisons of glmmTMB models
library(rcompanion) # for plots using medians and 95% CI


# import data----
hoverflies <- read.csv("Data/diversity_hoverflies_flowers.csv")
hoverflies_hedgerow <- read.csv("Data/hoverfly_hedgerow.csv")

# data manipulation----
# check data
head(hoverflies)
str(hoverflies)
head(hoverflies_hedgerow)
str(hoverflies_hedgerow)

# have hedgerow_type as a factor
hoverflies_hedgerow$hedgerow_type <- as.factor(hoverflies_hedgerow$hedgerow_type)

# create a new variable classifying pan traps as having high or low plant diversity
# this will be a useful classifying variable in carrying out community analysis
hoverflies$flower_diversity <- as.factor(ifelse(hoverflies$flower_shannon >= 0.75, 
                                                "High", "Low"))

# make other variables factors for statistical analyses
hoverflies_hedgerow$Farm <- as.factor(hoverflies_hedgerow$Farm)
hoverflies$ Farm <- as.factor(hoverflies$Farm)

# RQ1 Analysis----
# How does farmland floral diversity affect hoverfly communities?

# Effect of floral community on taxonomical diversity----
# floral diversity vs hoverfly shannon diversity----
# check response variable distribution
hist(hoverflies$hoverfly_shannon)

# not normal, cannot transform data because it has a lot of zeros

# build linear model
shannon_mod<- lm(hoverfly_shannon ~ flower_shannon, data = hoverflies)
summary(shannon_mod)

# visualise results
(prelim_plot <- ggplot(hoverflies, aes(x = flower_shannon, y = hoverfly_shannon)) +
    geom_point() +
    geom_smooth(method = "lm"))

# check assumptions
plot(shannon_mod, which = 1)  
plot(shannon_mod, which = 2) 

# ASSUMPTIONS VIOLATED

# see if we have independent data points
boxplot(hoverfly_shannon ~ Farm, data = hoverflies)  # could be something going on (grouped per farm)

# split data by farm
(split_plot <- ggplot(aes(flower_shannon, hoverfly_shannon), data = hoverflies) + 
    geom_point() + 
    facet_wrap(~ Farm) + 
    xlab("Floral Shannon") + 
    ylab("Hoverfly Shannon"))

# we can see that different farms have different floral diversities and hoverfly diversities,
# indicating that sites are not independent and we must include transect as a nested random effect

# mixed effects model
# create variable that combines hedgerow transects with farms, to account for experimental design
hoverflies$Farm <- as.factor(hoverflies$Farm)
hoverflies$Transect <- as.factor(hoverflies$Transect)

# create new variable combining farm and transect so we can account for the nested design
hoverflies <- within(hoverflies, hedgerow <- factor(Farm:Transect))

# build model
shannon.lmer <- lmer(hoverfly_shannon ~ flower_shannon + (1|Farm) + (1|hedgerow), data = hoverflies)
summary(shannon.lmer)

# check model assumptions with diagnostic plots
plot(shannon.lmer)
qqnorm(resid(shannon.lmer))
qqline(resid(shannon.lmer)) 

# ASSUMPTIONS VIOLATED

# do non parametric alternative--> Kendal Theill test
hovshannon.mod <- mblm(hoverfly_shannon ~ flower_shannon, dataframe=hoverflies)
summary(hovshannon.mod)

hovshannon.mod <- theil_sen_regression(hoverfly_shannon ~ flower_shannon, data=hoverflies)
summary(hovshannon.mod)

# spearman correlation as a non parametric alternative because Kendal Theill test can't cope with the zeros
hov.shannon_correlation <- cor.test(hoverflies$hoverfly_shannon, hoverflies$flower_shannon, method = "spearman", conf.level = 0.95, exact = FALSE)
hov.shannon_correlation

# can do kruskal test if floral diversity is classed as high or low (Shannon>0.5)
kruskal.test(hoverfly_shannon ~ flower_diversity, data = hoverflies)
boxplot(hoverfly_shannon ~ flower_diversity, data = hoverflies)

# floral abundance vs hoverfly shannon diversity----
# VIOLATES ASSUMPTIONS OF LINEAR MODELS AND LINEAR MIXED MODELS (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs flower shannon)
# do non parametric alternative--> Kendal Theill test
hovshannon.mod2 <- mblm(hoverfly_shannon ~ flower_abundance, dataframe=hoverflies)
summary(hovshannon.mod2)

hovshannon.mod2 <- theil_sen_regression(hoverfly_shannon ~ flower_abundance, data=hoverflies)
summary(hovshannon.mod2)

# spearman correlation as a non parametric alternative because Kendal Theill test can't cope with the zeros
hov.shannon_correlation2 <- cor.test(hoverflies$hoverfly_shannon, hoverflies$flower_abundance, method = "spearman", conf.level = 0.95, exact = FALSE)
hov.shannon_correlation2

# ZERO INFLATED MODELS----
# floral diversity vs hoverfly abundance----
# view distribution
hist(hoverflies$hoverfly_abundance)

# data is non negative count data (whole numbers) with a right skew and A LOT of zeros
# Deduce proportion of zeros in the data 
sum(hoverflies$hoverfly_abundance == 0)/nrow(hoverflies) 
# 40% of our observations are zeros

# see how response variable is distributed according to explanatory variable
(first_plot <- ggplot(hoverflies, aes(x = flower_shannon, y = hoverfly_abundance)) +
    geom_point() +
    geom_smooth(method = "lm"))

# can see that there is no equal variances, and data/residuals are not normally distributed

# start modelling
# fit a basic poisson regression into data
poisson_model <- glm(hoverfly_abundance~flower_shannon, data = hoverflies, family = poisson)
summary(poisson_model)
# seems like flower diversity has an effect on hoverfly abundance, but need to check model assumptions

# check model assumptions
# checking for overdispersion
mean(hoverflies$hoverfly_abundance) %>% round(4) # round to 4 decimal places
var(hoverflies$hoverfly_abundance) %>% round (4) # round to 4 decimal places
var(hoverflies$hoverfly_abundance)/mean(hoverflies$hoverfly_abundance) # calculate ratio between variance and mean 
# THERE IS CLEAR OVERDISPERSION OF THE DATA
# note: this could have also been checked looking at the summary output (residual deviance vs residual df) and checking if ratio is greater than 1
# Poisson model results cannot be trusted as it poorly fits the data
# if overdispersion ratio is greater than 2 then we should try to build another model to fit the data better

# try negative binomial response
# The glm.nb function characterizes the data based on the negative binomial error distribution 
# it deals with the same kind of data (non-negative count data) but doesn't assume equal mean and variance
# actually assumes higher dispersion of the count data
# good to deal with a lot of zeros but we may have too many zeros, we'll see
shannonabundance.nb <- glm.nb(hoverfly_abundance~flower_shannon, data = hoverflies)
summary(shannonabundance.nb)
# no significant effect of flower shannon on hoverfly abundance
# theta is the overdispersion parameter. Larger values means less variance and poisson may be more effective

# compare this model fit to poisson using AIC
AIC(poisson_model, shannonabundance.nb)
# nb model fits data better than poisson

# nb fits data better and account for overdispersion
# however, overdispersion is being caused by the large amount of zeros in the data
# thus, consider zero inflated model and assess model fit later

# build model and plot residuals to assess if we have a case of zero inflation
mod1 <- glmmTMB(hoverfly_abundance~flower_shannon, family = "poisson", data = hoverflies) # build model
simulationOutput <- simulateResiduals(fittedModel = mod1) # generate residuals
plot(simulationOutput) # visualize residuals
# focus on QQ plot
# make QQ plot neater and easier to interpret
plotQQunif(simulationOutput = simulationOutput, 
           testDispersion = FALSE,
           testUniformity = FALSE,
           testOutliers = FALSE)
# we can see that residuals are not correctly distributed following the 1-1 line, verify using zero inflation test

# test for zero inflation using the zero inflation test from the DHARMa package
testZeroInflation(simulationOutput)

# a value of ratioObsSim > 1 indicates zero inflation of the data (OUR CASE, YAY)
# p value indicates there is zero inflation present
# this advocates the use of a zero inflated model
# now we can build zero inflated model and compare it to the previous models we have constructed

# sampling zeros in our dataset are not due to chance, they are due to the fluctuation of hoverfly abundance at different floral diversities
# zero inflated models only work for sampling zeros, not structural

# build model
shannonabundance.zi<- glmmTMB(hoverfly_abundance ~ flower_shannon, ziformula = ~flower_shannon, family = "nbinom2", data = hoverflies)

# use negative binomial model due to overdispersion (argument 'family =' )
# nbiom2 means quadratic increases of the variance with respect to the mean
# look at zero model output
summary(shannonabundance.zi)

# conditional model is response in the absence of zero inflation
# zi model shows the probability of generating a structural zero not accounted for by the conditional model


# assess model fit
AIC(poisson_model, shannonabundance.nb, shannonabundance.zi)

# seems like negative binomial model accounted for the data better than the zero-inflated model

# check if data points are independent
boxplot(hoverfly_abundance ~ Farm, data = hoverflies)  # could be something going on (grouped per farm)

# split data by farm
(split_plot <- ggplot(aes(flower_shannon, hoverfly_abundance), data = hoverflies) + 
    geom_point() + 
    facet_wrap(~ Farm) + 
    xlab("Floral Shannon") + 
    ylab("Hoverfly Abundance"))

# we can see that different farms have different floral diversities and hoverfly abundance,
# indicating that sites are not independent and we must include farm and transect as nested random effects

# need to add sampling design into models (note: hedgerow is Farm + transect which we previously combined earlier wehn trying to fit a mixed effects model)

# in the negative binomial model
shannonabundance_mixed.nb <- glmmTMB(hoverfly_abundance ~ flower_shannon + (1|Farm) + (1|hedgerow), family = "nbinom2", data = hoverflies)

# in the zero inflated model
shannonabundance_mixed.zi <- glmmTMB(hoverfly_abundance ~ flower_shannon + (1|Farm) + (1|hedgerow), ziformula = ~flower_shannon, family = "nbinom2", data = hoverflies)
# this model doesn't converge, so ignore

# also add random effects to see if they affect the generation process of structural zeros
shannonabundance_mixed2.zi <- glmmTMB(hoverfly_abundance ~ flower_shannon + (1|Farm) + (1|hedgerow), ziformula = ~flower_shannon + (1|Farm) + (1|hedgerow), family = "nbinom2", data = hoverflies)

# assess model fit
AIC(shannonabundance_mixed.nb, shannonabundance.nb, shannonabundance_mixed2.zi, shannonabundance.zi)
# Mixed effect negative binomial model seems to be the best. 
# random effects have a strong influence on the generation of both sampling and structural zeros (site characteristics), but model fits worse than negative binomial
# structural zero means site characteristics do not allow for the proliference of hoverflies

# code negative binomial model using glmmTMB rather than glm.nb just to carry out likelihood ratio test
# to confirm that random effects are needed in the model
shannonabundance.nb <- glmmTMB(hoverfly_abundance ~ flower_shannon, family = "nbinom2", data = hoverflies)

# likelihood ratio test (fixed effects are held constant)
anova(shannonabundance_mixed.nb, shannonabundance.nb)
# indeed, random effects do improve model fit

# see findings using best model- mixed effect negative binomial
summary(shannonabundance_mixed.nb)
# found that flower diversity does not have a significant effect on hoverfly abundance
# can look at estimates to determine the ecological implications of our findings (need to do exp(estimate) to undo log transformation in each term (intercept and fixed effect) and compare)
# visualize
(first_plot <- ggplot(hoverflies, aes(x = flower_shannon, y = hoverfly_abundance)) +
    geom_point() +
    geom_smooth(method = "lm"))

# floral abundance vs hoverfly abundance----
# view distribution
hist(hoverflies$hoverfly_abundance)

# data is non negative count data (whole numbers) with a right skew and A LOT of zeros
# Deduce proportion of zeros in the data 
sum(hoverflies$hoverfly_abundance == 0)/nrow(hoverflies) 
# 40% of our observations are zeros

# see how response variable is distributed according to explanatory variable
(first_plot <- ggplot(hoverflies, aes(x = flower_abundance, y = hoverfly_abundance)) +
    geom_point() +
    geom_smooth(method = "lm"))

# can see that there is no equal variances, and data/residuals are not normally distributed

# start modelling
# fit a basic poisson regression into data
poisson_model2 <- glm(hoverfly_abundance~flower_abundance, data = hoverflies, family = poisson)
summary(poisson_model2)
# seems like flower diversity has an effect on hoverfly abundance, but need to check model assumptions

# check model assumptions
# checking for overdispersion
mean(hoverflies$hoverfly_abundance) %>% round(4) # round to 4 decimal places
var(hoverflies$hoverfly_abundance) %>% round (4) # round to 4 decimal places
var(hoverflies$hoverfly_abundance)/mean(hoverflies$hoverfly_abundance) # calculate ratio between variance and mean 
# THERE IS CLEAR OVERDISPERSION OF THE DATA
# note: this could have also been checked looking at the summary output (residual deviance vs residual df) and checking if ratio is greater than 1
# Poisson model results cannot be trusted as it poorly fits the data
# if overdispersion ratio is greater than 2 then we should try to build another model to fit the data better

# try negative binomial response
# The glm.nb function characterizes the data based on the negative binomial error distribution 
# it deals with the same kind of data (non-negative count data) but doesn't assume equal mean and variance
# actually assumes higher dispersion of the count data
# good to deal with a lot of zeros but we may have too many zeros, we'll see
abundance.nb <- glm.nb(hoverfly_abundance~flower_abundance, data = hoverflies)
summary(abundance.nb)
# no significant effect of flower shannon on hoverfly abundance
# theta is the overdispersion parameter. Larger values means less variance and poisson may be more effective

# compare this model fit to poisson using AIC
AIC(poisson_model, abundance.nb)
# nb model fits data better than poisson

# nb fits data better and account for overdispersion
# however, overdispersion is being caused by the large amount of zeros in the data
# thus, consider zero inflated model and assess model fit later

# build model and plot residuals to assess if we have a case of zero inflation
mod2 <- glmmTMB(hoverfly_abundance~flower_shannon, family = "poisson", data = hoverflies) # build model
simulationOutput2 <- simulateResiduals(fittedModel = mod2) # generate residuals
plot(simulationOutput2) # visualize residuals
# focus on QQ plot
# make QQ plot neater and easier to interpret
plotQQunif(simulationOutput = simulationOutput2, 
           testDispersion = FALSE,
           testUniformity = FALSE,
           testOutliers = FALSE)
# we can see that residuals are not correctly distributed following the 1-1 line, verify using zero inflation test

# test for zero inflation using the zero inflation test from the DHARMa package
testZeroInflation(simulationOutput2)

# a value of ratioObsSim > 1 indicates zero inflation of the data (OUR CASE, YAY)
# p value indicates there is zero inflation present
# this advocates the use of a zero inflated model
# now we can build zero inflated model and compare it to the previous models we have constructed

# sampling zeros in our dataset are not due to chance, they are due to the fluctuation of hoverfly abundance at different floral diversities
# zero inflated models only work for sampling zeros, not structural

# build model
abundance.zi<- glmmTMB(hoverfly_abundance ~ flower_abundance, ziformula = ~flower_abundance, family = "nbinom2", data = hoverflies)

# use negative binomial model due to overdispersion (argument 'family =' )
# nbiom2 means quadratic increases of the variance with respect to the mean
# look at zero model output
summary(abundance.zi)

# conditional model is response in the absence of zero inflation
# zi model shows the probability of generating a structural zero not accounted for by the conditional model


# assess model fit
AIC(poisson_model2, abundance.nb, abundance.zi)

# seems like negative binomial model accounted for the data better than the zero-inflated model

# check if data points are independent
boxplot(hoverfly_abundance ~ Farm, data = hoverflies)  # could be something going on (grouped per farm)

# split data by farm
(split_plot <- ggplot(aes(flower_abundance, hoverfly_abundance), data = hoverflies) + 
    geom_point() + 
    facet_wrap(~ Farm) + 
    xlab("Floral Abundance") + 
    ylab("Hoverfly Abundance"))

# we can see that different farms have different floral abundances and hoverfly abundances,
# indicating that sites are not independent and we must include farm and transect as nested random effects

# need to add sampling design into models (note: hedgerow is Farm + transect which we previously combined earlier wehn trying to fit a mixed effects model)

# in the negative binomial model
abundance_mixed.nb <- glmmTMB(hoverfly_abundance ~ flower_abundance + (1|Farm) + (1|hedgerow), family = "nbinom2", data = hoverflies)

# in the zero inflated model
abundance_mixed.zi <- glmmTMB(hoverfly_abundance ~ flower_abundance + (1|Farm) + (1|hedgerow), ziformula = ~flower_abundance, family = "nbinom2", data = hoverflies)
# this model doesn't converge, so ignore

# also add random effects to see if they affect the generation process of structural zeros
abundance_mixed2.zi <- glmmTMB(hoverfly_abundance ~ flower_abundance + (1|Farm) + (1|hedgerow), ziformula = ~flower_abundance + (1|Farm) + (1|hedgerow), family = "nbinom2", data = hoverflies)
# failed to converge

# changing the optimiser to help converge
abundance_mixed2.zib <- update(abundance_mixed2.zi,
                               control=glmmTMBControl(optimizer=optim,
                                                      optArgs=list(method="BFGS")))
# converges now

# assess model fit
AIC(abundance_mixed.nb, abundance.nb, abundance.zi, abundance_mixed.zi, abundance_mixed2.zib)
# Mixed effect negative binomial model seems to be the best. 
# random effects do NOT have a strong influence on the generation of both sampling and structural zeros (site characteristics)

# code negative binomial model using glmmTMB rather than glm.nb just to carry out likelihood ratio test
# to confirm that random effects are needed in the model
abundance.nb <- glmmTMB(hoverfly_abundance ~ flower_abundance, family = "nbinom2", data = hoverflies)

# likelihood ratio test (fixed effects are held constant)
anova(abundance_mixed.nb, abundance.nb)
# indeed, random effects do improve model fit

# see findings using best model- mixed effect negative binomial
summary(abundance_mixed.nb)
# found that flower abundance affects hoverfly abundance
# can look at estimates to determine the ecological implications of our findings (need to do exp(estimate) to undo log transformation in each term (intercept and fixed effect) and compare)

# visualize
(first_plot <- ggplot(hoverflies, aes(x = flower_abundance, y = hoverfly_abundance)) +
    geom_point() +
    geom_smooth(method = "lm"))

# floral diversity vs hoverfly richness----
# view distribution
hist(hoverflies$hoverfly_richness)

# data is non negative count data (whole numbers) with a right skew and A LOT of zeros
# Deduce proportion of zeros in the data 
sum(hoverflies$hoverfly_richness == 0)/nrow(hoverflies) 
# 40% of our observations are zeros

# see how response variable is distributed according to explanatory variable
(first_plot <- ggplot(hoverflies, aes(x = flower_shannon, y = hoverfly_richness)) +
    geom_point() +
    geom_smooth(method = "lm"))

# can see that there is no equal variances, and data/residuals are not normally distributed

# start modelling
# fit a basic poisson regression into data
poisson_model3 <- glm(hoverfly_richness~flower_shannon, data = hoverflies, family = poisson)
summary(poisson_model3)
# seems like flower diversity has an effect on hoverfly abundance, but need to check model assumptions

# check model assumptions
# checking for overdispersion
mean(hoverflies$hoverfly_richness) %>% round(4) # round to 4 decimal places
var(hoverflies$hoverfly_richness) %>% round (4) # round to 4 decimal places
var(hoverflies$hoverfly_richness)/mean(hoverflies$hoverfly_richness) # calculate ratio between variance and mean 
# THERE IS CLEAR OVERDISPERSION OF THE DATA
# note: this could have also been checked looking at the summary output (residual deviance vs residual df) and checking if ratio is greater than 1
# Poisson model results cannot be trusted as it poorly fits the data
# we do not have overdispersion ratio greater than 2, so not awful but we can still find a better model fit

# try negative binomial response
# The glm.nb function characterizes the data based on the negative binomial error distribution 
# it deals with the same kind of data (non-negative count data) but doesn't assume equal mean and variance
# actually assumes higher dispersion of the count data
# good to deal with a lot of zeros but we may have too many zeros, we'll see
shannonrichness.nb <- glm.nb(hoverfly_richness~flower_shannon, data = hoverflies)
summary(shannonrichness.nb)
# no significant effect of flower shannon on hoverfly abundance
# theta is the overdispersion parameter. Larger values means less variance and poisson may be more effective

# compare this model fit to poisson using AIC
AIC(poisson_model3, shannonrichness.nb)
# nb model fits data better than poisson

# nb fits data better and account for overdispersion
# however, overdispersion is being caused by the large amount of zeros in the data
# thus, consider zero inflated model and assess model fit later

# build model and plot residuals to assess if we have a case of zero inflation
mod3 <- glmmTMB(hoverfly_richness~flower_shannon, family = "poisson", data = hoverflies) # build model
simulationOutput3 <- simulateResiduals(fittedModel = mod3) # generate residuals
plot(simulationOutput) # visualize residuals
# focus on QQ plot
# make QQ plot neater and easier to interpret
plotQQunif(simulationOutput = simulationOutput3, 
           testDispersion = FALSE,
           testUniformity = FALSE,
           testOutliers = FALSE)
# we can see that residuals are not correctly distributed following the 1-1 line, verify using zero inflation test

# test for zero inflation using the zero inflation test from the DHARMa package
testZeroInflation(simulationOutput3)

# a value of ratioObsSim > 1 indicates zero inflation of the data (OUR CASE, but not by much this time)
# p value indicates there is zero inflation present
# this advocates the use of a zero inflated model
# now we can build zero inflated model and compare it to the previous models we have constructed

# sampling zeros in our dataset are not due to chance, they are due to the fluctuation of hoverfly abundance at different floral diversities
# zero inflated models only work for sampling zeros, not structural

# build model
shannonrichness.zi<- glmmTMB(hoverfly_richness ~ flower_shannon, ziformula = ~flower_shannon, family = "nbinom2", data = hoverflies)

# changing the optimiser to help converge
shannonrichness.zib <- update(shannonrichness.zi,
                               control=glmmTMBControl(optimizer=optim,
                                                      optArgs=list(method="BFGS")))
# converges now

# use negative binomial model due to overdispersion (argument 'family =' )
# nbiom2 means quadratic increases of the variance with respect to the mean
# look at zero model output
summary(shannonrichness.zib)

# conditional model is response in the absence of zero inflation
# zi model shows the probability of generating a structural zero not accounted for by the conditional model


# assess model fit
AIC(poisson_model3, shannonrichness.nb, shannonrichness.zib)

# seems like negative binomial model accounted for the data better than the zero-inflated model

# check if data points are independent
boxplot(hoverfly_richness ~ Farm, data = hoverflies)  # could be something going on (grouped per farm)

# split data by farm
(split_plot <- ggplot(aes(flower_shannon, hoverfly_richness), data = hoverflies) + 
    geom_point() + 
    facet_wrap(~ Farm) + 
    xlab("Floral Shannon") + 
    ylab("Hoverfly Richness"))

# we can see that different farms have different floral diversities and hoverfly richness,
# indicating that sites are not independent and we must include farm and transect as nested random effects

# need to add sampling design into models (note: hedgerow is Farm + transect which we previously combined earlier wehn trying to fit a mixed effects model)

# in the negative binomial model
shannonrichness_mixed.nb <- glmmTMB(hoverfly_richness ~ flower_shannon + (1|Farm) + (1|hedgerow), family = "nbinom2", data = hoverflies)

# in the zero inflated model
shannonrichness_mixed.zi <- glmmTMB(hoverfly_richness ~ flower_shannon + (1|Farm) + (1|hedgerow), ziformula = ~flower_shannon, family = "nbinom2", data = hoverflies)

# change optimizer to help convergence
shannonrichness_mixed.zib <- update(shannonrichness_mixed.zi,
                              control=glmmTMBControl(optimizer=optim,
                                                     optArgs=list(method="BFGS")))
# still doesn't converge so just ignore

# also add random effects to see if they affect the generation process of structural zeros
shannonrichness_mixed2.zi <- glmmTMB(hoverfly_richness ~ flower_shannon + (1|Farm) + (1|hedgerow), ziformula = ~flower_shannon + (1|Farm) + (1|hedgerow), family = "nbinom2", data = hoverflies)

shannonrichness_mixed2.zib <- update(shannonrichness_mixed2.zi,
                                    control=glmmTMBControl(optimizer=optim,
                                                           optArgs=list(method="BFGS")))
# still doesn't converge so just ignore

# assess model fit
AIC(shannonrichness_mixed.nb, shannonrichness.nb, shannonrichness.zib)
# Mixed effect negative binomial model seems to be the best. 

# code negative binomial model using glmmTMB rather than glm.nb just to carry out likelihood ratio test
# to confirm that random effects are needed in the model
shannonrichness.nb <- glmmTMB(hoverfly_richness ~ flower_shannon, family = "nbinom2", data = hoverflies)

# likelihood ratio test (fixed effects are held constant)
anova(shannonrichness_mixed.nb, shannonrichness.nb)
# indeed, random effects do improve model fit

# see findings using best model- mixed effect negative binomial
summary(shannonrichness_mixed.nb)
# found that flower diversity does not have a significant effect on hoverfly richness
# can look at estimates to determine the ecological implications of our findings (need to do exp(estimate) to undo log transformation in each term (intercept and fixed effect) and compare)
# visualize
(first_plot <- ggplot(hoverflies, aes(x = flower_shannon, y = hoverfly_richness)) +
    geom_point() +
    geom_smooth(method = "lm"))
# floral abundance vs hoverfly richness----
# view distribution
hist(hoverflies$hoverfly_richness)

# data is non negative count data (whole numbers) with a right skew and A LOT of zeros
# Deduce proportion of zeros in the data 
sum(hoverflies$hoverfly_richness == 0)/nrow(hoverflies) 
# 40% of our observations are zeros

# see how response variable is distributed according to explanatory variable
(first_plot <- ggplot(hoverflies, aes(x = flower_abundance, y = hoverfly_richness)) +
    geom_point() +
    geom_smooth(method = "lm"))

# can see that there is no equal variances, and data/residuals are not normally distributed

# start modelling
# fit a basic poisson regression into data
poisson_model4 <- glm(hoverfly_richness~flower_abundance, data = hoverflies, family = poisson)
summary(poisson_model4)
# seems like flower diversity has an effect on hoverfly abundance, but need to check model assumptions

# check model assumptions
# checking for overdispersion
mean(hoverflies$hoverfly_richness) %>% round(4) # round to 4 decimal places
var(hoverflies$hoverfly_richness) %>% round (4) # round to 4 decimal places
var(hoverflies$hoverfly_richness)/mean(hoverflies$hoverfly_richness) # calculate ratio between variance and mean 
# THERE IS CLEAR OVERDISPERSION OF THE DATA
# note: this could have also been checked looking at the summary output (residual deviance vs residual df) and checking if ratio is greater than 1
# Poisson model results cannot be trusted as it poorly fits the data
# we do not have overdispersion ratio greater than 2, so not awful but we can still find a better model fit

# try negative binomial response
# The glm.nb function characterizes the data based on the negative binomial error distribution 
# it deals with the same kind of data (non-negative count data) but doesn't assume equal mean and variance
# actually assumes higher dispersion of the count data
# good to deal with a lot of zeros but we may have too many zeros, we'll see
abundancerichness.nb <- glm.nb(hoverfly_richness~flower_abundance, data = hoverflies)
summary(abundancerichness.nb)
# no significant effect of flower shannon on hoverfly abundance
# theta is the overdispersion parameter. Larger values means less variance and poisson may be more effective

# compare this model fit to poisson using AIC
AIC(poisson_model4, abundancerichness.nb)
# nb model fits data better than poisson

# nb fits data better and account for overdispersion
# however, overdispersion is being caused by the large amount of zeros in the data
# thus, consider zero inflated model and assess model fit later

# build model and plot residuals to assess if we have a case of zero inflation
mod4 <- glmmTMB(hoverfly_richness~flower_abundance, family = "poisson", data = hoverflies) # build model
simulationOutput4 <- simulateResiduals(fittedModel = mod4) # generate residuals
plot(simulationOutput3) # visualize residuals
# focus on QQ plot
# make QQ plot neater and easier to interpret
plotQQunif(simulationOutput = simulationOutput4, 
           testDispersion = FALSE,
           testUniformity = FALSE,
           testOutliers = FALSE)
# we can see that residuals are not correctly distributed following the 1-1 line, verify using zero inflation test

# test for zero inflation using the zero inflation test from the DHARMa package
testZeroInflation(simulationOutput4)

# a value of ratioObsSim > 1 indicates zero inflation of the data (OUR CASE, but not by much this time)
# p value indicates there is zero inflation present
# this advocates the use of a zero inflated model
# now we can build zero inflated model and compare it to the previous models we have constructed

# sampling zeros in our dataset are not due to chance, they are due to the fluctuation of hoverfly abundance at different floral diversities
# zero inflated models only work for sampling zeros, not structural

# build model
abundancerichness.zi<- glmmTMB(hoverfly_richness ~ flower_abundance, ziformula = ~flower_abundance, family = "nbinom2", data = hoverflies)

# use negative binomial model due to overdispersion (argument 'family =' )
# nbiom2 means quadratic increases of the variance with respect to the mean
# look at zero model output
summary(abundancerichness.zi)

# conditional model is response in the absence of zero inflation
# zi model shows the probability of generating a structural zero not accounted for by the conditional model


# assess model fit
AIC(poisson_model4, abundancerichness.nb, abundancerichness.zi)

# seems like negative binomial model accounted for the data better than the zero-inflated model

# check if data points are independent
boxplot(hoverfly_richness ~ Farm, data = hoverflies)  # could be something going on (grouped per farm)

# split data by farm
(split_plot <- ggplot(aes(flower_shannon, hoverfly_richness), data = hoverflies) + 
    geom_point() + 
    facet_wrap(~ Farm) + 
    xlab("Floral Abundance") + 
    ylab("Hoverfly Richness"))

# we can see that different farms have different floral diversities and hoverfly richness,
# indicating that sites are not independent and we must include farm and transect as nested random effects

# need to add sampling design into models (note: hedgerow is Farm + transect which we previously combined earlier wehn trying to fit a mixed effects model)

# in the negative binomial model
abundancerichness_mixed.nb <- glmmTMB(hoverfly_richness ~ flower_abundance + (1|Farm) + (1|hedgerow), family = "nbinom2", data = hoverflies)

# in the zero inflated model
abundancerichness_mixed.zi <- glmmTMB(hoverfly_richness ~ flower_abundance + (1|Farm) + (1|hedgerow), ziformula = ~flower_abundance, family = "nbinom2", data = hoverflies)

# also add random effects to see if they affect the generation process of structural zeros
abundancerichness_mixed2.zi <- glmmTMB(hoverfly_richness ~ flower_abundance + (1|Farm) + (1|hedgerow), ziformula = ~flower_abundance + (1|Farm) + (1|hedgerow), family = "nbinom2", data = hoverflies)

# change optimizer to aid convergence
abundancerichness_mixed2.zib <- update(abundancerichness_mixed2.zi,
                                     control=glmmTMBControl(optimizer=optim,
                                                            optArgs=list(method="BFGS")))
# still doesn't converge so just ignore

# assess model fit
AIC(abundancerichness_mixed.nb, abundancerichness.nb, abundancerichness_mixed.zi, abundancerichness.zi)
# Mixed effect negative binomial model seems to be the best. 

# code negative binomial model using glmmTMB rather than glm.nb just to carry out likelihood ratio test
# to confirm that random effects are needed in the model
abundancerichness.nb <- glmmTMB(hoverfly_richness ~ flower_abundance, family = "nbinom2", data = hoverflies)

# likelihood ratio test (fixed effects are held constant)
anova(abundancerichness_mixed.nb, abundancerichness.nb)
# indeed, random effects do improve model fit

# see findings using best model- mixed effect negative binomial
summary(abundancerichness_mixed.nb)
# found that flower diversity does not have a significant effect on hoverfly richness
# can look at estimates to determine the ecological implications of our findings (need to do exp(estimate) to undo log transformation in each term (intercept and fixed effect) and compare)
# visualize
(first_plot <- ggplot(hoverflies, aes(x = flower_shannon, y = hoverfly_richness)) +
    geom_point() +
    geom_smooth(method = "lm"))
# Effect of floral diversity on functional diversity----

# floral diversity vs bodylength fd----
# VIOLATES ASSUMPTIONS OF LINEAR MODELS AND LINEAR MIXED MODELS (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs flower shannon)
# do non parametric alternative--> Kendal Theill test
bodylengthfd.mod <- mblm(bodylength_fd ~ flower_shannon, dataframe=hoverflies)
summary(bodylengthfd.mod)

bodylengthfd.mod <- theil_sen_regression(bodylength_fd ~ flower_shannon, data=hoverflies)
summary(bodylengthfd.mod)

# spearman correlation as a non parametric alternative because Kendal Theill test can't cope with the zeros
bodylengthfd_correlation <- cor.test(hoverflies$bodylength_fd, hoverflies$flower_shannon, method = "spearman", conf.level = 0.95, exact = FALSE)
bodylengthfd_correlation

# can do kruskal test if floral diversity is classed as high or low (Shannon>0.5)
kruskal.test(bodylength_fd ~ flower_diversity, data = hoverflies)
boxplot(bodylength_fd ~ flower_diversity, data = hoverflies)

# floral diversity vs bodylength cwm----
# VIOLATES ASSUMPTIONS OF LINEAR MODELS AND LINEAR MIXED MODELS (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs flower shannon)
# do non parametric alternative--> Kendal Theill test
bodylengthcwm.mod <- mblm(bodylength_cwm ~ flower_shannon, dataframe=hoverflies, repeated = T)
summary(bodylengthcwm.mod)
# repeated = TRUE is for siegel emthod, FALSE is for Theil Sen method

# significant when considering Siegel method but not significant when using Then method 

bodylengthcwm.mod <- siegel_regression(bodylength_cwm ~ flower_shannon, data=hoverflies)
summary(bodylengthcwm.mod)
# USE THIS TEST
# siegel is better as it is less sensitive to outliers in the data than theil sen

bodylengthcwm.lm <- lm(bodylength_cwm ~ flower_shannon, data=hoverflies)

# see siegel regression vs linear regression
plot(bodylength_cwm ~ flower_shannon, data=hoverflies)
abline(bodylengthcwm.mod,col='blue')
abline(bodylengthcwm.lm,col='red')

# spearman correlation 
bodylengthcwm_correlation <- cor.test(hoverflies$bodylength_cwm, hoverflies$flower_shannon, method = "spearman", conf.level = 0.95, exact = FALSE)
bodylengthcwm_correlation

# can do kruskal test if floral diversity is classed as high or low (Shannon>0.5)
kruskal.test(bodylength_cwm ~ flower_diversity, data = hoverflies)
boxplot(bodylength_cwm ~ flower_diversity, data = hoverflies)

# floral diversity vs wing:body ratio fd----
# VIOLATES ASSUMPTIONS OF LINEAR MODELS AND LINEAR MIXED MODELS (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs flower shannon)
# do non parametric alternative--> Kendal Theill test
wingbodyratiofd.mod <- mblm(wing_body_ratio_fd ~ flower_shannon, dataframe=hoverflies, repeated = TRUE)
summary(wingbodyratiofd.mod)

wingbodyratiofd.mod <- theil_sen_regression(wing_body_ratio_fd ~ flower_shannon, data=hoverflies)
summary(wingbodyratiofd.mod)

# spearman correlation as a non parametric alternative because Kendal Theill test can't cope with the zeros
wingbodyratiofd_correlation <- cor.test(hoverflies$wing_body_ratio_fd, hoverflies$flower_shannon, method = "spearman", conf.level = 0.95, exact = FALSE)
wingbodyratiofd_correlation

# can do kruskal test if floral diversity is classed as high or low (Shannon>0.5)
kruskal.test(wing_body_ratio_fd ~ flower_diversity, data = hoverflies)
boxplot(wing_body_ratio_fd ~ flower_diversity, data = hoverflies)

# floral diversity vs wing:body ratio cwm----
# VIOLATES ASSUMPTIONS OF LINEAR MODELS AND LINEAR MIXED MODELS (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs flower shannon)
# do non parametric alternative--> Kendal Theill test
wingbodyratiocwm.mod <- mblm(wing_body_ratio_cwm ~ flower_shannon, dataframe=hoverflies, repeated = TRUE)
summary(wingbodyratiocwm.mod)

wingbodyratiocwm.mod <- siegel_regression(wing_body_ratio_cwm ~ flower_shannon, data=hoverflies)
summary(wingbodyratiocwm.mod)
# could trust siegel regression?, if not then do spearman correlation

plot(wing_body_ratio_cwm~flower_shannon, data = hoverflies)
abline(wingbodyratiocwm.mod)

# spearman correlation as a non parametric alternative if Kendal Theill test can't cope with the zeros
wingbodyratiocwm_correlation <- cor.test(hoverflies$wing_body_ratio_cwm, hoverflies$flower_shannon, method = "spearman", conf.level = 0.95, exact = FALSE)
wingbodyratiocwm_correlation

# can do kruskal test if floral diversity is classed as high or low (Shannon>0.5)
kruskal.test(wing_body_ratio_cwm ~ flower_diversity, data = hoverflies)
boxplot(wing_body_ratio_cwm ~ flower_diversity, data = hoverflies)
# floral diversity vs sex fd----
# VIOLATES ASSUMPTIONS OF LINEAR MODELS AND LINEAR MIXED MODELS (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs flower shannon)
# do non parametric alternative--> Kendal Theill test
sexfd.mod <- mblm(sex_fd ~ flower_shannon, dataframe=hoverflies, repeated = TRUE)
summary(sexfd.mod)

sexfd.mod <- theil_sen_regression(sex_fd ~ flower_shannon, data=hoverflies)
summary(sexfd.mod)

# spearman correlation as a non parametric alternative because Kendal Theill test can't cope with the zeros
sexfd_correlation <- cor.test(hoverflies$sex_fd, hoverflies$flower_shannon, method = "spearman", conf.level = 0.95, exact = FALSE)
sexfd_correlation

# can do kruskal test if floral diversity is classed as high or low (Shannon>0.5)
kruskal.test(sex_fd ~ flower_diversity, data = hoverflies)
boxplot(sex_fd ~ flower_diversity, data = hoverflies)
# floral diversity vs sex cwm----
# VIOLATES ASSUMPTIONS OF LINEAR MODELS AND LINEAR MIXED MODELS (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs flower shannon)
# do non parametric alternative--> Kendal Theill test
sexcwm.mod <- mblm(sex_cwm ~ flower_shannon, dataframe=hoverflies, repeated = TRUE)
summary(sexcwm.mod)

sexcwm.mod <- theil_sen_regression(sex_cwm ~ flower_shannon, data=hoverflies)
summary(sexcwm.mod)

# spearman correlation as a non parametric alternative because Kendal Theill test can't cope with the zeros
sexcwm_correlation <- cor.test(hoverflies$sex_cwm, hoverflies$flower_shannon, method = "spearman", conf.level = 0.95, exact = FALSE)
sexcwm_correlation

# can do kruskal test if floral diversity is classed as high or low (Shannon>0.5)
kruskal.test(sex_cwm ~ flower_diversity, data = hoverflies)
boxplot(sex_cwm ~ flower_diversity, data = hoverflies)

# Floral diversity NMDS----
# comparing high vs low floral diversity hoverfly communities
# import data
hoverflies_species <- read.csv("data/species_matrix.csv")
View(hoverflies_species)

# data manipulation
# view dataset structure
head(hoverflies_species)
str(hoverflies_species)

# add previously calculated column of high vs low plant diversity
hoverflies_species$flower_diversity <- hoverflies$flower_diversity

# remove rows (pan traps) with zero hoverfly abundance as this will lead to error in NMDS
hoverflies_species <- hoverflies_species %>% 
  filter(., Abundance != 0)

# make community matrix by extracting abundance values from inverts dataframe
hoverfly_community <- hoverflies_species[,7:41]

# turn data frame into matrix so it functions in the vegan package
hoverfly_matrix <- as.matrix(hoverfly_community)

# ordination
# now we will proceed to create the ordination
hov.NMDS <- metaMDS(hoverfly_matrix, distance = "bray", k = 3, autotransform = TRUE, trymax=100) 

# Bray-Curtis distance is chosen because it is not affected by zero values. 
# k represents the number of dimensions we want and is used to reduce stress.
# autotransform arguments ensures the data does a sqrt transformation
# trymax is the number of iterations the algorithm will run

# check stress value
hov.NMDS$stress # ideally stress value is < 0.2

# data visualization
# using vegan package and Base R
plot(hov.NMDS) # circles show different communities/sites, crosses show different species
ordiplot(hov.NMDS, type = "n") # create blank ordination plot
orditorp(hov.NMDS, display = "species", col="red", air = 0.1) # add order names in red
orditorp(hov.NMDS, display = "sites", cex = 1.25, air = 0.1) # add site numbers in black

# by distance
# ellipse plot
ordiplot(hov.NMDS) # plot shows communities (circles) and species (crosses)
ordiellipse(hov.NMDS, hoverflies_species$flower_diversity, label = FALSE, 
            col=c("darkorchid1", "darkslategray1"), 
            draw = "polygon", alpha=120)
legend("topright", title="Floral Diversity",
       c("High","Low"), fill=c("darkorchid1", "darkslategray1"), horiz=FALSE, cex=.9) # adding a legend

# save plot
png("your_filepath/floral_nmds_basic.png", width=6, height=5, units = "in", res = 300)
ordiplot(hov.NMDS)
ordiellipse(hov.NMDS, hoverflies_species$flower_diversity, label = FALSE, 
            col=c("darkorchid1", "darkslategray1"), 
            draw = "polygon", alpha=120)
legend("topright", title="Distance (m)",
       c("High","Low"), fill=c("darkorchid1", "darkslategray1"), horiz=FALSE, cex=.9)
dev.off()

# polygon plot
ordiplot(hov.NMDS) #plot shows communities (circles) and species (crosses)
ordihull(hov.NMDS, groups = hoverflies_species$flower_diversity, draw="polygon", col="grey90", label = TRUE) # adding polygons to the plot, grouping by distance (inverts$Distance)

# spider plot
ordiplot(hov.NMDS) #plot shows communities (circles) and species (crosses)
ordispider(hov.NMDS, groups = hoverflies_species$flower_diversity, label = TRUE) # adding spider plot, grouping by distance (inverts$Distance)


# using ggplot to make an ellipse plot
# make new dataframe with by extracting NMDS scores
nmds.scores <- as.data.frame(scores(hov.NMDS)$sites) # for newest version of vegan package

# if you have a version of the vegan package <2.6-2 then the following code should also work:
# nmds.scores = as.data.frame(scores(nmds))

# add data from your original dataframe into new NMDS dataframe
# useful if you want to group based on different criteria
nmds.scores <- nmds.scores %>% 
  mutate(flower_diversity = as.factor(hoverflies_species$flower_diversity))

# check dataframe to ensure all changes have taken place
head(nmds.scores) 
str(nmds.scores)

# define hidden vegan function that finds coordinates for drawing a covariance ellipse
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# create empty dataframe to combine NMDS data with ellipse data
ellipse_df <- data.frame()

# adding data for ellipse, in this case using distance as a grouping factor
for(g in levels(nmds.scores$flower_diversity)){
  ellipse_df <- rbind(ellipse_df, cbind(as.data.frame(with(nmds.scores[nmds.scores$flower_diversity==g,],
                                                           veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                  wt=rep(1/length(NMDS1),length(NMDS1)))$cov,
                                                                           center=c(mean(NMDS1),mean(NMDS2)))))
                                        ,flower_diversity=g))
}

# create ggplot
(flower_NMDS_plot <- ggplot(data = nmds.scores, aes(NMDS1, NMDS2)) + 
    geom_point(aes(color = flower_diversity, shape = flower_diversity)) + # adding different colours and shapes for points at different distances
    geom_path(data=ellipse_df, aes(x=NMDS1, y=NMDS2, colour=flower_diversity), linewidth=1) + # adding covariance ellipses according to distance # use size argument if ggplot2 < v. 3.4.0
    guides(color = guide_legend(override.aes = list(linetype=c(NA,NA)))) + # removes lines from legend
    theme_bw() + # adding theme
    theme(panel.grid = element_blank()) + # remove background grid
    scale_color_manual(name = "Flower diversity", # legend title
                       labels = c("High", "Low"), # adjusting legend labels
                       values = c("goldenrod1", "black")) + # customising colours
    scale_shape_manual("Flower diversity", # legend title
                       labels = c("High", "Low"), # adjusting legend labels
                       values = c(17, 15))) # customising shapes
# save plot
ggsave(filename = "your_filepath/flower_NMDS.png", flower_NMDS_plot, device = "png")

# data analysis
# using a PERMANOVA to test the differences in community composition
# This is a PERmutational Multivariate ANalysis Of VAriance and tests the differences between groups, like an ANOVA, but with lots of variables.
# it is essentially a multivariate analysis of variance used to compare groups of objects
flower_permanova <- adonis2(as.matrix(hoverflies_species [,7:41]) ~ flower_diversity, 
                            hoverflies_species, permutations = 999, method = "bray")
# permutations is the number of possible arrangements the data could take
# using permutations = 999 is standard practice in science and across the literature
# can increase it if you want to get a more thorough analysis
# use method = bray as it is what we previously used to calculate pairwise distances

# look at model output
flower_permanova

# check model assumptions
# check for multivariate homogeneity of group variances

# generate distance matrix from invertebrate commumunity matrix
flower.dist <- vegdist(hoverfly_matrix, method = "bray")

# use betadisper test to check for multivariate homogeneity of group variances 
flower.dispersion <- betadisper(flower.dist, group=hoverflies_species$flower_diversity)
permutest(flower.dispersion)
# CAN TRUST PERMANOVA RESULTS

# can also visually assess the assumptions
plot(flower.dispersion, hull = FALSE, ellipse=TRUE) # remove hulls and add ellipses


# simper analysis
# this is a SIMilarity PERcentage analysis and compares community differences 
# and reports what species (orders in this scenario)are driving those differences.

flower_simper <- with(hoverflies_species[,c(1,2,3,4,5,6,42)], simper(as.matrix(hoverflies_species[,7:41]), flower_diversity, permutations = 100))
# with argument indicates what grouping variables to use to conduct the analysis 
# the number of permutations required to carry out the analysis

# see most influential species (orders in this case) contributing to community differences
flower_simper

# a more detailed summary showing the contributions of all species (or orders)
summary(flower_simper)

# shows species driving differences in community composition between different distances

# looks like most abundant species (marmalade, migrant, common snout) dominate low diversity
# plant communities. High diversity plant communities allow for other species of hoverflies 
# to come in and use the resources, however, both communities are dominated by these abundant 
# hoverflies and they are not significantly different.

# RQ2 Analysis----
# How does farmland hedgerow quality affect hoverfly communities?

# Effect of hedgerow type on hoverfly taxonomical diversity----
# hedgerow type vs hoverfly shannon----
# check response variable distribution
hist(hoverflies_hedgerow$hoverfly_shannon)

# not normal, cannot transform data because it has a lot of zeros

# build linear model
shannonhedge_mod<- lm(hoverfly_shannon ~ hedgerow_type, data = hoverflies_hedgerow)
anova(shannonhedge_mod)
summary(shannonhedge_mod)

# visualise results
boxplot(hoverfly_shannon ~ hedgerow_type, data = hoverflies_hedgerow)

# check assumptions
plot(shannonhedge_mod, which = 1)  
plot(shannonhedge_mod, which = 2) 
shapiro.test(hoverflies_hedgerow$hoverfly_shannon) 
lm_resids <- resid(shannonhedge_mod)
shapiro.test(lm_resids) # shows that residuals are non normally distributed
bartlett.test(hoverfly_shannon ~ hedgerow_type, data = hoverflies_hedgerow) # shows homoskedasticity (GOOD)

# ASSUMPTIONS VIOLATED

# see if we have independent data points
boxplot(hoverfly_shannon ~ Farm, data = hoverflies_hedgerow)  # could be something going on (grouped per farm)

# split data by farm
(split_plot <- ggplot(aes(hedgerow_type, hoverfly_shannon), data = hoverflies_hedgerow) + 
    geom_boxplot() + 
    facet_wrap(~ Farm) + 
    xlab("Hedgerow Type") + 
    ylab("Hoverfly Shannon"))

# we can see that different farms have different hedgerow types and hoverfly diversities,
# indicating that sites are not independent and we must include farm as a random effect because
# it affects the diversity~hedgerow type relationship

# mixed effects model
# include farm as a random effect because different farms have different hedgerows and different hoverfly
# communities
# build model
shannonhedge.lmer <- lmer(hoverfly_shannon ~ hedgerow_type + (1|Farm), data = hoverflies_hedgerow)
summary(shannonhedge.lmer)

# check model assumptions with diagnostic plots
plot(shannonhedge.lmer)
qqnorm(resid(shannonhedge.lmer))
qqline(resid(shannonhedge.lmer)) 
shapiro.test(resid(shannonhedge.lmer))

# ASSUMPTIONS VIOLATED- not sure, seems quite close to qqline, we could maybe use it

# do non parametric alternative--> Kruskal Wallis Test
hovshannonhedge.mod <- kruskal.test(hoverfly_shannon ~ hedgerow_type, data=hoverflies_hedgerow)
hovshannonhedge.mod

# post hoc comparisons
dunnTest(hoverfly_shannon ~ hedgerow_type, data=hoverflies_hedgerow, method = "holm") 

# visualise results
boxplot(hoverfly_shannon ~ hedgerow_type, data = hoverflies_hedgerow)

# p values were adjusted with holm method, just as good as bonferroni but more powerful (and no additional assumptions)
# holm method has strong control for family-wise error rate

# ZERO INFLATED MODELS----
# hedgerow type vs hoverfly abundance----
# view distribution
hist(hoverflies_hedgerow$hoverfly_abundance)

# data is non negative count data (whole numbers) with a right skew and A LOT of zeros
# Deduce proportion of zeros in the data 
sum(hoverflies_hedgerow$hoverfly_abundance == 0)/nrow(hoverflies_hedgerow) 
# 18% of our observations are zeros

# see how response variable is distributed according to explanatory variable
(first_boxplot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, hoverfly_abundance)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    ylab("Hoverfly Abundance\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factors

# can see that there is no equal variances, and data/residuals are not normally distributed

# start modelling
# fit a basic poisson regression into data
poisson_model5 <- glm(hoverfly_abundance~hedgerow_type, data = hoverflies_hedgerow, family = poisson)
summary(poisson_model5)
# seems like flower diversity has an effect on hoverfly abundance, but need to check model assumptions

# check model assumptions
# checking for overdispersion
mean(hoverflies_hedgerow$hoverfly_abundance) %>% round(4) # round to 4 decimal places
var(hoverflies_hedgerow$hoverfly_abundance) %>% round (4) # round to 4 decimal places
var(hoverflies_hedgerow$hoverfly_abundance)/mean(hoverflies_hedgerow$hoverfly_abundance) # calculate ratio between variance and mean 
# THERE IS CLEAR OVERDISPERSION OF THE DATA, cannot trust resulrs
# note: this could have also been checked looking at the summary output (residual deviance vs residual df) and checking if ratio is greater than 1
# Poisson model results cannot be trusted as it poorly fits the data
# if overdispersion ratio is greater than 2 then we should try to build another model to fit the data better

# try negative binomial response
# The glm.nb function characterizes the data based on the negative binomial error distribution 
# it deals with the same kind of data (non-negative count data) but doesn't assume equal mean and variance
# actually assumes higher dispersion of the count data
# good to deal with a lot of zeros but we may have too many zeros, we'll see
hedgeabundance.nb <- glm.nb(hoverfly_abundance~hedgerow_type, data = hoverflies_hedgerow)
summary(hedgeabundance.nb)
# no significant effect of flower shannon on hoverfly abundance
# theta is the overdispersion parameter. Larger values means less variance and poisson may be more effective

# compare this model fit to poisson using AIC
AIC(poisson_model5, hedgeabundance.nb)
# nb model fits data better than poisson

# nb fits data better and account for overdispersion
# however, overdispersion is being caused by the large amount of zeros in the data
# thus, consider zero inflated model and assess model fit later

# build model and plot residuals to assess if we have a case of zero inflation
mod5 <- glmmTMB(hoverfly_abundance~hedgerow_type, family = "poisson", data = hoverflies_hedgerow) # build model
simulationOutput5 <- simulateResiduals(fittedModel = mod5) # generate residuals
plot(simulationOutput5) # visualize residuals
# focus on QQ plot
# make QQ plot neater and easier to interpret
plotQQunif(simulationOutput = simulationOutput5, 
           testDispersion = FALSE,
           testUniformity = FALSE,
           testOutliers = FALSE)
# we can see that residuals are not correctly distributed following the 1-1 line, verify using zero inflation test

# test for zero inflation using the zero inflation test from the DHARMa package
testZeroInflation(simulationOutput5)

# a value of ratioObsSim > 1 indicates zero inflation of the data (OUR CASE, YAY)
# p value indicates there is zero inflation present
# this advocates the use of a zero inflated model
# now we can build zero inflated model and compare it to the previous models we have constructed

# sampling zeros in our dataset are not due to chance, they are due to the fluctuation of hoverfly abundance at different floral diversities
# zero inflated models only work for sampling zeros, not structural

# build model
hedgeabundance.zi<- glmmTMB(hoverfly_abundance ~ hedgerow_type, ziformula = ~hedgerow_type, family = "nbinom2", data = hoverflies_hedgerow)

# use negative binomial model due to overdispersion (argument 'family =' )
# nbiom2 means quadratic increases of the variance with respect to the mean
# look at zero model output
summary(hedgeabundance.zi)

# conditional model is response in the absence of zero inflation
# zi model shows the probability of generating a structural zero not accounted for by the conditional model


# assess model fit
AIC(poisson_model5, hedgeabundance.nb, hedgeabundance.zi)

# seems like negative binomial model accounted for the data better than the zero-inflated model

# check if data points are independent
boxplot(hoverfly_abundance ~ Farm, data = hoverflies_hedgerow)  # could be something going on (grouped per farm)

# split data by farm
(split_plot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, hoverfly_abundance)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    facet_wrap(~Farm) +
    ylab("Hoverfly Abundance\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factors


# we can see that different farms have different hedgerow types and hoverfly abundances,
# indicating that sites are not independent and we must include farm as a nested random effect

# need to add sampling design into models

# in the negative binomial model
hedgeabundance_mixed.nb <- glmmTMB(hoverfly_abundance ~ hedgerow_type + (1|Farm), family = "nbinom2", data = hoverflies_hedgerow)

# in the zero inflated model
hedgeabundance_mixed.zi <- glmmTMB(hoverfly_abundance ~ hedgerow_type + (1|Farm), ziformula = ~hedgerow_type, family = "nbinom2", data = hoverflies_hedgerow)

# also add random effects to see if they affect the generation process of structural zeros
hedgeabundance_mixed2.zi <- glmmTMB(hoverfly_abundance ~ hedgerow_type + (1|Farm), ziformula = ~hedgerow_type + (1|Farm), family = "nbinom2", data = hoverflies_hedgerow)
# change optimizer to aid convergence
hedgeabundance_mixed2.zib <- update(hedgeabundance_mixed2.zi,
                                       control=glmmTMBControl(optimizer=optim,
                                                              optArgs=list(method="BFGS")))
# converges now

# assess model fit
AIC(hedgeabundance.nb, hedgeabundance_mixed.nb, hedgeabundance_mixed.zi, hedgeabundance_mixed2.zib, hedgeabundance.zi)
# Mixed effect negative binomial model seems to be the best. 
# random effects also have an effect on the data in the zero inflated model, but negative binomial fits data better

# code negative binomial model using glmmTMB rather than glm.nb just to carry out likelihood ratio test
# to confirm that random effects are needed in the model
hedgeabundance.nb <- glmmTMB(hoverfly_abundance ~ hedgerow_type, family = "nbinom2", data = hoverflies_hedgerow)

# likelihood ratio test (fixed effects are held constant)
anova(hedgeabundance_mixed.nb, hedgeabundance.nb)
# indeed, random effects do improve model fit

# see findings using best model- mixed effect negative binomial
summary(hedgeabundance_mixed.nb)
# found that flower diversity does not have a significant effect on hoverfly abundance
# can look at estimates to determine the ecological implications of our findings (need to do exp(estimate) to undo log transformation in each term (intercept and fixed effect) and compare)
# visualize
(plot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, hoverfly_abundance)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    ylab("Hoverfly Abundance\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factors

# hedgerow type  vs hoverfly richness----
# view distribution
hist(hoverflies_hedgerow$hoverfly_richness)

# data is non negative count data (whole numbers) with a right skew and A LOT of zeros
# Deduce proportion of zeros in the data 
sum(hoverflies_hedgerow$hoverfly_richness == 0)/nrow(hoverflies_hedgerow) 
# 18% of our observations are zeros

# see how response variable is distributed according to explanatory variable
(first_boxplot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, hoverfly_richness)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    ylab("Hoverfly Richness\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factors

# can see that there is no equal variances, and data/residuals are not normally distributed

# start modelling
# fit a basic poisson regression into data
poisson_model6 <- glm(hoverfly_richness~hedgerow_type, data = hoverflies_hedgerow, family = poisson)
summary(poisson_model6)
# seems like flower diversity has an effect on hoverfly abundance, but need to check model assumptions

# check model assumptions
# checking for overdispersion
mean(hoverflies_hedgerow$hoverfly_richness) %>% round(4) # round to 4 decimal places
var(hoverflies_hedgerow$hoverfly_richness) %>% round (4) # round to 4 decimal places
var(hoverflies_hedgerow$hoverfly_richness)/mean(hoverflies_hedgerow$hoverfly_richness) # calculate ratio between variance and mean 
# data is slightly overdispersed
# note: this could have also been checked looking at the summary output (residual deviance vs residual df) and checking if ratio is greater than 1
# Poisson model results cannot be trusted as it poorly fits the data
# not incredibly overdispersed but we can try find a better model fit

# try negative binomial response
# The glm.nb function characterizes the data based on the negative binomial error distribution 
# it deals with the same kind of data (non-negative count data) but doesn't assume equal mean and variance
# actually assumes higher dispersion of the count data
# good to deal with a lot of zeros but we may have too many zeros, we'll see
hedgerichness.nb <- glm.nb(hoverfly_richness~hedgerow_type, data = hoverflies_hedgerow)
summary(hedgerichness.nb)
# no significant effect of flower shannon on hoverfly abundance
# theta is the overdispersion parameter. Larger values means less variance and poisson may be more effective

# compare this model fit to poisson using AIC
AIC(poisson_model6, hedgerichness.nb)
# nb model fits data better than poisson

# nb fits data better and account for overdispersion
# however, overdispersion is being caused by the large amount of zeros in the data
# thus, consider zero inflated model and assess model fit later

# build model and plot residuals to assess if we have a case of zero inflation
mod6 <- glmmTMB(hoverfly_richness~hedgerow_type, family = "poisson", data = hoverflies_hedgerow) # build model
simulationOutput6 <- simulateResiduals(fittedModel = mod6) # generate residuals
plot(simulationOutput6) # visualize residuals
# focus on QQ plot
# make QQ plot neater and easier to interpret
plotQQunif(simulationOutput = simulationOutput6, 
           testDispersion = FALSE,
           testUniformity = FALSE,
           testOutliers = FALSE)
# we can see that residuals are not correctly distributed following the 1-1 line, verify using zero inflation test

# test for zero inflation using the zero inflation test from the DHARMa package
testZeroInflation(simulationOutput6)

# a value of ratioObsSim > 1 indicates zero inflation of the data (OUR CASE, but not by much)
# p value indicates there is zero inflation present
# this advocates the use of a zero inflated model
# now we can build zero inflated model and compare it to the previous models we have constructed

# sampling zeros in our dataset are not due to chance, they are due to the fluctuation of hoverfly abundance at different floral diversities
# zero inflated models only work for sampling zeros, not structural

# build model
hedgerichness.zi<- glmmTMB(hoverfly_richness ~ hedgerow_type, ziformula = ~hedgerow_type, family = "nbinom2", data = hoverflies_hedgerow)

# use negative binomial model due to overdispersion (argument 'family =' )
# nbiom2 means quadratic increases of the variance with respect to the mean
# look at zero model output
summary(hedgerichness.zi)

# conditional model is response in the absence of zero inflation
# zi model shows the probability of generating a structural zero not accounted for by the conditional model


# assess model fit
AIC(poisson_model6, hedgerichness.nb, hedgerichness.zi)

# seems like negative binomial model accounted for the data better than the zero-inflated model, but just about, both models equally good

# check if data points are independent
boxplot(hoverfly_richness ~ Farm, data = hoverflies_hedgerow)  # could be something going on (grouped per farm)

# split data by farm
(split_plot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, hoverfly_richness)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    facet_wrap(~Farm) +
    ylab("Hoverfly Richness\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factors


# we can see that different farms have different hedgerow types and hoverfly abundances,
# indicating that sites are not independent and we must include farm as a nested random effect

# need to add sampling design into models

# in the negative binomial model
hedgerichness_mixed.nb <- glmmTMB(hoverfly_richness ~ hedgerow_type + (1|Farm), family = "nbinom2", data = hoverflies_hedgerow)
# change optimizer to aid convergence
hedgerichness_mixed.nb2 <- update(hedgerichness_mixed.nb,
                                    control=glmmTMBControl(optimizer=optim,
                                                           optArgs=list(method="BFGS")))
# converges now

# in the zero inflated model
hedgerichness_mixed.zi <- glmmTMB(hoverfly_richness ~ hedgerow_type + (1|Farm), ziformula = ~hedgerow_type, family = "nbinom2", data = hoverflies_hedgerow)

# change optimizer to aid convergence
hedgerichness_mixed.zib <- update(hedgerichness_mixed.zi,
                                  control=glmmTMBControl(optimizer=optim,
                                                         optArgs=list(method="BFGS")))
# converges now

# also add random effects to see if they affect the generation process of structural zeros
hedgerichness_mixed2.zi <- glmmTMB(hoverfly_richness ~ hedgerow_type + (1|Farm), ziformula = ~hedgerow_type + (1|Farm), family = "nbinom2", data = hoverflies_hedgerow)
# change optimizer to aid convergence
hedgerichness_mixed2.zib <- update(hedgerichness_mixed2.zi,
                                    control=glmmTMBControl(optimizer=optim,
                                                           optArgs=list(method="BFGS")))
#doesn't converge so ignore

# assess model fit
AIC(hedgerichness.nb, hedgerichness_mixed.nb2, hedgerichness.zi, hedgerichness_mixed.zib)
# Mixed effect negative binomial model seems to be the best, but just about
# random effects also have an effect on the data in the zero inflated model, but negative binomial fits data better

# code negative binomial model using glmmTMB rather than glm.nb just to carry out likelihood ratio test
# to confirm that random effects are needed in the model
hedgerichness.nb <- glmmTMB(hoverfly_richness ~ hedgerow_type, family = "nbinom2", data = hoverflies_hedgerow)

# likelihood ratio test (fixed effects are held constant)
anova(hedgerichness_mixed.nb2, hedgerichness.nb)
# indeed, random effects do improve model fit

# see findings using best model- mixed effect negative binomial
summary(hedgerichness_mixed.nb2)

# post hoc comparison using ls means
lsmeans(hedgerichness_mixed.nb2, pairwise ~ hedgerow_type)

# found that hedgerow type has a significant effect on hoverfly richness
# can look at estimates to determine the ecological implications of our findings (need to do exp(estimate) to undo log transformation in each term (intercept and fixed effect) and compare)
# visualize
(plot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, hoverfly_richness)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    ylab("Hoverfly Richness\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factors
# Effect of hedgerow type on hoverfly functional diversity----

# hedgerow type vs bodylength fd----
# ASSUMPTIONS VIOLATED(code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs hedgerow type)
# do non parametric alternative--> Kruskal Wallis Test
bodylengthfdhedge.mod <- kruskal.test(bodylength_fd ~ hedgerow_type, data=hoverflies_hedgerow)
bodylengthfdhedge.mod

# post hoc comparisons
dunnTest(bodylength_fd ~ hedgerow_type, data=hoverflies_hedgerow, method = "holm") 

# visualise results
boxplot(bodylength_fd ~ hedgerow_type, data = hoverflies_hedgerow)

# hedgerow type vs bodylength cwm----
# ASSUMPTIONS VIOLATED (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs hedgerow type)
# do non parametric alternative--> Kruskal Wallis Test
bodylengthcwmhedge.mod <- kruskal.test(bodylength_cwm ~ hedgerow_type, data=hoverflies_hedgerow)
bodylengthcwmhedge.mod

# post hoc comparisons
dunnTest(bodylength_cwm ~ hedgerow_type, data=hoverflies_hedgerow, method = "holm") 

# visualise results
boxplot(bodylength_cwm ~ hedgerow_type, data = hoverflies_hedgerow)

# hedgerow type vs wing:body ratio fd----
# ASSUMPTIONS VIOLATED (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs hedgerow type)
# do non parametric alternative--> Kruskal Wallis Test
wingbodyratiofdhedge.mod <- kruskal.test(wing_body_ratio_fd ~ hedgerow_type, data=hoverflies_hedgerow)
wingbodyratiofdhedge.mod

# post hoc comparisons
dunnTest(wing_body_ratio_fd ~ hedgerow_type, data=hoverflies_hedgerow, method = "holm") 

# visualise results
boxplot(wing_body_ratio_fd ~ hedgerow_type, data = hoverflies_hedgerow)

# hedgerow type vs wing:body ratio cwm----
# ASSUMPTIONS VIOLATED (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs hedgerow type)
# do non parametric alternative--> Kruskal Wallis Test
wingbodyratiocwmhedge.mod <- kruskal.test(wing_body_ratio_cwm ~ hedgerow_type, data=hoverflies_hedgerow)
wingbodyratiocwmhedge.mod

# post hoc comparisons
dunnTest(wing_body_ratio_cwm ~ hedgerow_type, data=hoverflies_hedgerow, method = "holm") 

# visualise results
boxplot(wing_body_ratio_cwm ~ hedgerow_type, data = hoverflies_hedgerow)

# hedgerow type vs sex fd----
# ASSUMPTIONS VIOLATED (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs hedgerow type)
# do non parametric alternative--> Kruskal Wallis Test
sexfdhedge.mod <- kruskal.test(sex_fd ~ hedgerow_type, data=hoverflies_hedgerow)
sexfdhedge.mod

# post hoc comparisons
dunnTest(sex_fd ~ hedgerow_type, data=hoverflies_hedgerow, method = "holm") 

# visualise results
boxplot(sex_fd ~ hedgerow_type, data = hoverflies_hedgerow)

# hedgerow type vs sex cwm----
# ASSUMPTIONS VIOLATED (code not included for practicality, 
# but if interested, it is the same as the one used to check hoverfly shannon vs hedgerow type)
# do non parametric alternative--> Kruskal Wallis Test
sexcwmhedge.mod <- kruskal.test(sex_cwm ~ hedgerow_type, data=hoverflies_hedgerow)
sexcwmhedge.mod

# post hoc comparisons
dunnTest(sex_cwm ~ hedgerow_type, data=hoverflies_hedgerow, method = "holm") 

# visualise results
boxplot(sex_cwm ~ hedgerow_type, data = hoverflies_hedgerow)

# Hedgerow NMDS----
# comparing hoverfly communities across different hedgerow types
# import data
hedgerow_species <- read.csv("Data/species_matrix_transect.csv")
View(hedgerow_species)

# data manipulation
# view dataset structure
head(hedgerow_species)
str(hedgerow_species)

# add hedgerow type to the dataset
hedgerow_species$hedgerow_type <- as.factor(hoverflies_hedgerow$hedgerow_type)

# remove rows (pan traps) with zero hoverfly abundance as this will lead to error in NMDS
hedgerow_species <- hedgerow_species %>% 
  filter(., Abundance != 0)

# make community matrix by extracting abundance values from inverts dataframe
hedgerow_community <- hedgerow_species[,4:39]

# turn data frame into matrix so it functions in the vegan package
hedgerow_matrix <- as.matrix(hedgerow_community)

# ordination
# now we will proceed to create the ordination
hedge.NMDS <- metaMDS(hedgerow_matrix, distance = "bray", k = 3, autotransform = TRUE, trymax=100) 

# Bray-Curtis distance is chosen because it is not affected by zero values. 
# k represents the number of dimensions we want and is used to reduce stress.
# autotransform arguments ensures the data does a sqrt transformation
# trymax is the number of iterations the algorithm will run

# check stress value
hedge.NMDS$stress # ideally stress value is < 0.2

# data visualization
# using vegan package and Base R
plot(hedge.NMDS) # circles show different communities/sites, crosses show different species
ordiplot(hedge.NMDS, type = "n") # create blank ordination plot
orditorp(hedge.NMDS, display = "species", col="red", air = 0.1) # add order names in red
orditorp(hedge.NMDS, display = "sites", cex = 1.25, air = 0.1) # add site numbers in black

# by distance
# ellipse plot
ordiplot(hedge.NMDS)
ordiellipse(hedge.NMDS, hedgerow_species$hedgerow_type, label = FALSE, 
            col=c("darkorchid1", "darkslategray1"), 
            draw = "polygon", alpha=120)
legend("topright", title="Hedgerow Type",
       c("Well Managed","Overgrown", "Overtrimmed"), fill=c("darkorchid1", "darkslategray1", "bisque1"), horiz=FALSE, cex=.9)

# save plot
png("your_filepath/hedgerow_nmds_basic.png", width=6, height=5, units = "in", res = 300)
ordiplot(hedge.NMDS)
ordiellipse(hedge.NMDS, hedgerow_species$hedgerow_type, label = FALSE, 
            col=c("darkorchid1", "darkslategray1"), 
            draw = "polygon", alpha=120)
legend("topright", title="Hedgerow Type",
       c("Well Managed","Overgrown", "Overtrimmed"), fill=c("darkorchid1", "darkslategray1", "bisque1"), horiz=FALSE, cex=.9)
dev.off()

# polygon plot
ordiplot(hedge.NMDS) #plot shows communities (circles) and species (crosses)
ordihull(hedge.NMDS, groups = hedgerow_species$hedgerow_type, draw="polygon", col="grey90", label = TRUE) # adding polygons to the plot, grouping by distance (inverts$Distance)

# spider plot
ordiplot(hedge.NMDS) #plot shows communities (circles) and species (crosses)
ordispider(hedge.NMDS, groups = hedgerow_species$hedgerow_type, label = TRUE) # adding spider plot, grouping by distance (inverts$Distance)


# using ggplot to make an ellipse plot
# make new dataframe with by extracting NMDS scores
nmds.scores1 <- as.data.frame(scores(hedge.NMDS)$sites) # for newest version of vegan package

# if you have a version of the vegan package <2.6-2 then the following code should also work:
# nmds.scores = as.data.frame(scores(nmds))

# add data from your original dataframe into new NMDS dataframe
# useful if you want to group based on different criteria
nmds.scores1 <- nmds.scores1 %>% 
  mutate(hedgerow_type = as.factor(hedgerow_species$hedgerow_type))

# check dataframe to ensure all changes have taken place
head(nmds.scores1) 
str(nmds.scores1)

# define hidden vegan function that finds coordinates for drawing a covariance ellipse
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# create empty dataframe to combine NMDS data with ellipse data
ellipse_df1 <- data.frame()

# adding data for ellipse, in this case using distance as a grouping factor
for(r in levels(nmds.scores1$hedgerow_type)){
  ellipse_df1 <- rbind(ellipse_df1, cbind(as.data.frame(with(nmds.scores1[nmds.scores1$hedgerow_type==r,],
                                                           veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                  wt=rep(1/length(NMDS1),length(NMDS1)))$cov,
                                                                           center=c(mean(NMDS1),mean(NMDS2)))))
                                        ,hedgerow_type=r))
}

# create ggplot
(hedgerow_NMDS_plot <- ggplot(data = nmds.scores1, aes(NMDS1, NMDS2)) + 
    geom_point(aes(color = hedgerow_type, shape = hedgerow_type)) + # adding different colours and shapes for points at different distances
    geom_path(data=ellipse_df1, aes(x=NMDS1, y=NMDS2, colour=hedgerow_type), linewidth=1) + # adding covariance ellipses according to distance # use size argument if ggplot2 < v. 3.4.0
    guides(color = guide_legend(override.aes = list(linetype=c(NA,NA,NA)))) + # removes lines from legend
    theme_bw() + # adding theme
    theme(panel.grid = element_blank()) + # remove background grid
    scale_color_manual(name = "Hedgerow Type", # legend title
                       labels = c("Well Managed","Overgrown", "Overtrimmed"), # adjusting legend labels
                       values = c("#FFC125", "#006400", "#8B4513")) + # customising colours
    scale_shape_manual("Hedgerow Type", # legend title
                       labels = c("Well Managed","Overgrown", "Overtrimmed"), # adjusting legend labels
                       values = c(17, 15, 3))) # customising shapes
# save plot
ggsave(filename = "your_filepath/hedgerow_NMDS.png", hedgerow_NMDS_plot, device = "png")

# data analysis
# using a PERMANOVA to test the differences in community composition
# This is a PERmutational Multivariate ANalysis Of VAriance and tests the differences between groups, like an ANOVA, but with lots of variables.
# it is essentially a multivariate analysis of variance used to compare groups of objects
hedgerow_permanova <- adonis2(as.matrix(hedgerow_species [,5:39]) ~ hedgerow_type, 
                            hedgerow_species, permutations = 999, method = "bray")
# permutations is the number of possible arrangements the data could take
# using permutations = 999 is standard practice in science and across the literature
# can increase it if you want to get a more thorough analysis
# use method = bray as it is what we previously used to calculate pairwise distances

# look at model output
hedgerow_permanova

# check model assumptions
# check for multivariate homogeneity of group variances

# generate distance matrix from invertebrate commumunity matrix
hedgerow.dist <- vegdist(hedgerow_matrix, method = "bray")

# use betadisper test to check for multivariate homogeneity of group variances 
hedgerow.dispersion <- betadisper(hedgerow.dist, group=hedgerow_species$hedgerow_type)
permutest(hedgerow.dispersion)
# CAN TRUST PERMANOVA RESULTS

# can also visually assess the assumptions
plot(hedgerow.dispersion, hull = FALSE, ellipse=TRUE) # remove hulls and add ellipses


# simper analysis
# this is a SIMilarity PERcentage analysis and compares community differences 
# and reports what species (orders in this scenario)are driving those differences.

hedgerow_simper <- with(hedgerow_species[,c(1,2,3,4,40)], simper(as.matrix(hedgerow_species[,5:39]), hedgerow_type, permutations = 100))
# with argument indicates what grouping variables to use to conduct the analysis 
# the number of permutations required to carry out the analysis

# see most influential species (orders in this case) contributing to community differences
hedgerow_simper

# a more detailed summary showing the contributions of all species (or orders)
summary(hedgerow_simper)

# shows species driving differences in community composition between different distances

# Most abundant species (Episyrphus balteatus, Eupeodes corollae, Neoascia podagrica)
# dominate all hedgerows equally and are not affected by hedgerow quality. However, rarer
# species such as Helophilus pendulus appear to have hedgerow type preferences.
# Different hedgerow preferences by rarer species, but do not affect the overall community


# data visualization----
# section for all code to produce dissertation figures
# import hoverfly data by species count
hoverfly_count <- read.csv("Data/species_count.csv")
# Create a bar plot of initial distribution----
# remove row containing total count
hoverfly_count <- hoverfly_count [-36,]

# plot
(hoverfly_barplot <- ggplot(hoverfly_count, aes(x = reorder(Species, -Abundance), y = Abundance)) +
    geom_bar(position = position_dodge(), stat = "identity", colour = "black", fill = "#FFC125") +
    geom_text(aes(label = Abundance), size = 4, angle = 45, vjust = -0.5, hjust = 0, fontface = "bold", color = "black") +
    theme_bw() +
    ylim(0, 830) + 
    ylab("Abundance\n") +                             
    xlab("Hoverfly Species")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1, face = "italic"),  # Angled labels, so text doesn't overlap
          axis.text.y = element_text(size = 12, color = "black"),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),                      
          panel.grid = element_blank(),                                          
          plot.margin = unit(c(1,1,1,1), units = , "cm")))
# save plot
ggsave(filename = "your_filepath/hoverfly_barplot.png", hoverfly_barplot, device = "png")

# plots for hoverfly shannon vs floral shannon, floral abundance and hedgerow type----
# hoverfly shannon vs floral shannon
(shannon_plot <- ggplot(hoverflies, aes (flower_shannon, hoverfly_shannon)) +
  geom_point(col="#000000") +    
  geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
  theme_bw() +
  ylab("Hoverfly Shannon Diversity\n") +                             
  xlab("\nFlower Shannon Diversity")  +
  theme(axis.line = element_line(),
         axis.text.x = element_text(size = 12),     
         axis.text.y = element_text(size = 12),
         panel.border = element_blank(),
         axis.title = element_text(size = 14, face = "bold"),                        
         panel.grid = element_blank(),                                   # Removing the background grid lines               
         plot.margin = unit(c(1,1,1,1), units = , "cm")))

# hoverfly shannon vs floral abundance
(shannon_plot2 <- ggplot(hoverflies, aes (flower_abundance, hoverfly_shannon)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Hoverfly Shannon Diversity\n") +                             
    xlab("\nFlower Abundance")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),     
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                        
          panel.grid = element_blank(),                                   # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

# hoverfly shannon vs hedgerow_type
(shannon_boxplot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, hoverfly_shannon)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    ylab("Hoverfly Shannon Diversity\n") +                             
    xlab("\nHedgerow Type")  +
    ylim(0, 2.5) +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

# combine all plots
shannon_combined <- (shannon_plot|shannon_plot2)/shannon_boxplot + plot_annotation(tag_levels = "A") & 
                    theme(plot.tag = element_text(face = "bold"))
shannon_combined

# save plot
ggsave(filename = "your_filepath/hoverfly_shannon.png", shannon_combined, device = "png")

# plots for hoverfly abundance and richness vs flower diversity, flower abundance, hedgerow type----
# hoverfly abundance vs floral shannon
(abundance_plot <- ggplot(hoverflies, aes (flower_shannon, hoverfly_abundance)) +
   geom_point(col="#000000") +    
   geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
   theme_bw() +
   ylab("Hoverfly Abundance\n") +                             
   xlab("\nFlower Shannon Diversity")  +
   theme(axis.line = element_line(),
         axis.text.x = element_text(size = 12),     
         axis.text.y = element_text(size = 12),
         panel.border = element_blank(),
         axis.title = element_text(size = 14, face = "bold"),                        
         panel.grid = element_blank(),                                   # Removing the background grid lines               
         plot.margin = unit(c(1,1,1,1), units = , "cm")))

# hoverfly abundance vs floral abundance
(abundance_plot2 <- ggplot(hoverflies, aes (flower_abundance, hoverfly_abundance)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Hoverfly Abundance\n") +                             
    xlab("\nFlower Abundance")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),     
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                        
          panel.grid = element_blank(),                                   # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

# hoverfly abundance vs hedgerow_type
(abundance_boxplot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, hoverfly_abundance)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    ylab("Hoverfly Abundance\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

# hoverfly richness vs floral shannon
(richness_plot <- ggplot(hoverflies, aes (flower_shannon, hoverfly_richness)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Hoverfly Richness\n") +                             
    xlab("\nFlower Shannon Diversity")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),     
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                        
          panel.grid = element_blank(),                                   # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

# hoverfly richness vs floral abundance
(richness_plot2 <- ggplot(hoverflies, aes (flower_abundance, hoverfly_richness)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Hoverfly Richness\n") +                             
    xlab("\nFlower Abundance")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),     
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                        
          panel.grid = element_blank(),                                   # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

# hoverfly richness vs hedgerow_type
(richness_boxplot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, hoverfly_richness)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    ylab("Hoverfly Richness\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

richabundance_combined <- abundance_plot + abundance_plot2 + abundance_boxplot + richness_plot + richness_plot2 + richness_boxplot +
                          plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))
richabundance_combined
                          
# save plot
ggsave(filename = "your_filepath/hoverfly_richabundance.png", richabundance_combined, device = "png")



# plots for functional diversity indices for all traits vs flower diversity----
# bodylength fd vs floral shannon
(bodyfd_plot <- ggplot(hoverflies, aes (flower_shannon, bodylength_fd)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Body Length FD\n") +                             
    xlab("\nFlower Shannon Diversity")  +
    theme(axis.line = element_line(),
         axis.text.x = element_text(size = 12),     
         axis.text.y = element_text(size = 12),
         panel.border = element_blank(),
         axis.title = element_text(size = 14, face = "bold"),                        
         panel.grid = element_blank(),                                   # Removing the background grid lines               
         plot.margin = unit(c(1,1,1,1), units = , "cm")))
# bodylength cwm vs floral shannon
(bodycwm_plot <- ggplot(hoverflies, aes (flower_shannon, bodylength_cwm)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Body Length CWM (mm)\n") +                             
    xlab("\nFlower Shannon Diversity")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),     
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                        
          panel.grid = element_blank(),                                   # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

# wing:body ratio fd vs floral shannon
(wingfd_plot <- ggplot(hoverflies, aes (flower_shannon, wing_body_ratio_fd)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Wing:Body Ratio FD\n") +                             
    xlab("\nFlower Shannon Diversity")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),     
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                        
          panel.grid = element_blank(),                                   # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

# wing:body ratio cwm vs floral shannon
(wingcwm_plot <- ggplot(hoverflies, aes (flower_shannon, wing_body_ratio_cwm)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Wing:Body Ratio CWM\n") +                             
    xlab("\nFlower Shannon Diversity")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),     
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                        
          panel.grid = element_blank(),                                   # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

# sex fd vs floral shannon
(sexfd_plot <- ggplot(hoverflies, aes (flower_shannon, sex_fd)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Sex FD\n") +                             
    xlab("\nFlower Shannon Diversity")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),     
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                        
          panel.grid = element_blank(),                                   # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

# sex cwm vs floral shannon
(sexcwm_plot <- ggplot(hoverflies, aes (flower_shannon, sex_cwm)) +
    geom_point(col="#000000") +    
    geom_smooth(method = "lm", col="#FFC125", fill="#FFC125") +
    theme_bw() +
    ylab("Sex CWM\n") +                             
    xlab("\nFlower Shannon Diversity")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 12),     
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                        
          panel.grid = element_blank(),                                   # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

fd_combined <- bodycwm_plot + wingcwm_plot + sexcwm_plot + bodyfd_plot + wingfd_plot + sexfd_plot + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))
fd_combined

# save plot
ggsave(filename = "your_filepath/hoverfly_fd.png", fd_combined, device = "png")


# plots for functional diversity indices for all traits vs hedgerow type----
# body length fd vs hedgerow type
# calculate medians and 95% CI
bodyfd_median <- groupwiseMedian(bodylength_fd ~ hedgerow_type,
                      data = hoverflies_hedgerow,
                      bca=FALSE, percentile=TRUE)

bodyfd_median

# plot
(bodyfd_boxplot <- ggplot(bodyfd_median, aes(hedgerow_type, Median, color = hedgerow_type)) + 
    geom_point(shape = 15, size = 4) +
    geom_errorbar(aes(ymin  =  Percentile.lower, ymax  =  Percentile.upper), width =  0.2, linewidth  =  0.7) +
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    theme_bw() +
    ylim(0,3) +
    ylab("Body Length FD\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

(bodyfd_boxplot <- ggplot(hoverflies_hedgerow, aes(hedgerow_type, bodylength_fd)) + 
    geom_boxplot(aes(fill = hedgerow_type)) +
    theme_bw() +
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    ylab("Body Length FD\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

(bodyfd_barplot <- ggplot(bodyfd_median, aes(hedgerow_type, Median, color = hedgerow_type)) + 
    geom_bar(aes(fill = hedgerow_type), position = position_dodge(), stat = "identity") +
    geom_errorbar(aes(ymin  =  Percentile.lower, ymax  =  Percentile.upper), width =  0.2, linewidth  =  0.5, colour = "black") +
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    theme_bw() +
    ylim(0,2) +
    ylab("Body Length FD\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

# body length cwm vs hedgerow type
# calculate medians and 95% CI
bodycwm_median <- groupwiseMedian(bodylength_cwm ~ hedgerow_type,
                                 data = hoverflies_hedgerow,
                                 bca=FALSE, percentile=TRUE)

bodycwm_median

# plot
(bodycwm_barplot <- ggplot(bodycwm_median, aes(hedgerow_type, Median, color = hedgerow_type)) + 
    geom_bar(aes(fill = hedgerow_type), position = position_dodge(), stat = "identity") +
    geom_errorbar(aes(ymin  =  Percentile.lower, ymax  =  Percentile.upper), width =  0.2, linewidth  =  0.5, colour = "black") +
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    theme_bw() +
    ylim(0,15) +
    ylab("Body Length CWM (mm)\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

# wing:body ratio fd vs hedgerow type
# calculate medians and 95% CI
wingfd_median <- groupwiseMedian(wing_body_ratio_fd ~ hedgerow_type,
                                  data = hoverflies_hedgerow,
                                  bca=FALSE, percentile=TRUE)

wingfd_median

# plot
(wingfd_barplot <- ggplot(wingfd_median, aes(hedgerow_type, Median, color = hedgerow_type)) + 
    geom_bar(aes(fill = hedgerow_type), position = position_dodge(), stat = "identity") +
    geom_errorbar(aes(ymin  =  Percentile.lower, ymax  =  Percentile.upper), width =  0.2, linewidth  =  0.5, colour = "black") +
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    theme_bw() +
    ylab("Wing:Body Ratio FD\n") +                             
    xlab("\nHedgerow Type")  +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

# wing:body ratio fd vs hedgerow type
# calculate medians and 95% CI
wingcwm_median <- groupwiseMedian(wing_body_ratio_cwm ~ hedgerow_type,
                                 data = hoverflies_hedgerow,
                                 bca=FALSE, percentile=TRUE)

wingcwm_median

# plot
(wingcwm_barplot <- ggplot(wingcwm_median, aes(hedgerow_type, Median, color = hedgerow_type)) + 
    geom_bar(aes(fill = hedgerow_type), position = position_dodge(), stat = "identity") +
    geom_errorbar(aes(ymin  =  Percentile.lower, ymax  =  Percentile.upper), width =  0.2, linewidth  =  0.5, colour = "black") +
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    theme_bw() +
    ylab("Wing:Body Ratio CWM\n") +                             
    xlab("\nHedgerow Type")  +
    ylim (0,1) +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

# sex fd vs hedgerow type
# calculate medians and 95% CI
sexfd_median <- groupwiseMedian(sex_fd ~ hedgerow_type,
                                  data = hoverflies_hedgerow,
                                  bca=FALSE, percentile=TRUE)

sexfd_median

# plot
(sexfd_barplot <- ggplot(sexfd_median, aes(hedgerow_type, Median, color = hedgerow_type)) + 
    geom_bar(aes(fill = hedgerow_type), position = position_dodge(), stat = "identity") +
    geom_errorbar(aes(ymin  =  Percentile.lower, ymax  =  Percentile.upper), width =  0.2, linewidth  =  0.5, colour = "black") +
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    theme_bw() +
    ylab("Sex FD\n") +                             
    xlab("\nHedgerow Type")  +
    ylim(0,0.6) +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

# sex cwm vs hedgerow type
# calculate medians and 95% CI
sexcwm_median <- groupwiseMedian(sex_cwm ~ hedgerow_type,
                                  data = hoverflies_hedgerow,
                                  bca=FALSE, percentile=TRUE)

sexcwm_median

# plot
(sexcwm_barplot <- ggplot(sexcwm_median, aes(hedgerow_type, Median, color = hedgerow_type)) + 
    geom_bar(aes(fill = hedgerow_type), position = position_dodge(), stat = "identity") +
    geom_errorbar(aes(ymin  =  Percentile.lower, ymax  =  Percentile.upper), width =  0.2, linewidth  =  0.5, colour = "black") +
    scale_colour_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_fill_manual(values = c("#FFC125", "#006400", "#8B4513")) + # Adding custom colours
    scale_x_discrete(labels= c("Well-managed", "Overgrown", "Overtrimmed")) +
    theme_bw() +
    ylab("Sex CWM\n") +                             
    xlab("\nHedgerow Type")  +
    ylim(0, 1) +
    theme(axis.line = element_line(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          panel.border = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),                     
          panel.grid = element_blank(), # Removing the background grid lines               
          plot.margin = unit(c(1,1,1,1), units = , "cm"), # Adding a margin
          legend.position = "none")) # Removing legend - not needed with only 2 factorsg legend - not needed with only 2 factors

hedgefd_combined <- bodycwm_barplot + wingcwm_barplot + sexcwm_barplot + bodyfd_barplot + wingfd_barplot + sexfd_barplot + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))
hedgefd_combined

# save plot
ggsave(filename = "your_filepath/hoverflyhedege_fd.png", hedgefd_combined, device = "png")

