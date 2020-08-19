library(tidyverse)
library(brms)
library(bayestestR)
studyA <- read.csv(file="studyAresults.csv", header=TRUE, sep=",")

library(ggplot2)

set.seed(112)
data <- matrix(sample(1:30,15) , nrow=3)
data.frame(studyA$id:studyA$frustration)
data <- matrix(studyA$mental:studyA$physical)
colnames(data) <- c("mental","physical")
rownames(data) <- studyA$used_system

# Grouped barplot
barplot(data, 
        col=colors()[c(23,89,12)] , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=rownames(data), 
        xlab="group", 
        font.lab=2)


interactionTime <- brm(formula = interaction_time ~ used_system, 
                  data = studyA, 
                  family = lognormal(),
                  prior = set_prior('normal(0, 10)'))
interactionTime
conditional_effects(interactionTime)

navigationTime <- brm(formula = navigation_time ~ used_system, 
                       data = studyA, 
                       family = lognormal(),
                       prior = set_prior('normal(0, 10)'))
navigationTime
conditional_effects(navigationTime)

bayesfactor_parameters(navigationTime)

aggregate(interaction_time~used_system, studyA, mean)

m1 <- brm(preference ~ 1, data=studyA, family=bernoulli)
m1
# The estimate (.87) is on a different scale, you need to transform it:
inv_logit_scaled(.87)
#...and you get the same estimate (0.7), as you had in the frequentist tests above.

# But the real power of this framework comes into play when you want to test the influence
# of independent variables (predictors). 
# For example, did SBSOD score influence the preference?
studyA <- studyA %>% mutate(SBSOD.sum = rowSums(.[7:21])) #ugly code
m2 <- brm(preference ~ SBSOD.sum, data=studyA, family=bernoulli)
m2
conditional_effects(m2)
# no it didn't. But do you see how easy it is to add predictors?

# How do we test if this predictor is significant? Use bayestestR.
library(bayestestR)
bayesfactor_parameters(m3)

# Oh no!
# In order to formally test the hypothesis that SBSOD influences preference, 
# we need to specify some priors.
# let's go for any prior for now, just to demonstrate.

m3 <- brm(preference ~ SBSOD.sum, data=studyA, family=bernoulli,
          prior = set_prior("normal(0,1)", class="b"),
          sample_prior = TRUE, save_all_pars = TRUE)

# In fact this model has two priors (one was added automatically), you can see them here:
prior_summary(m3)

# And this package gives you easy access to testing various hypotheses, 
# e.g. with the use of ROPE, Bayes Factors, and others
# https://github.com/easystats/bayestestR
library(bayestestR)
bayesfactor_parameters(m3)
rope(m3)

# This is the summary of my workflow:
# 1. Fit a model using brms
# 2. Test with appropriate method using bayestestR.

