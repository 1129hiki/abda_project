############### Dependency ############### 
# install.packages("brms")
# install.packages("ggplot2")
# install.packages("tidyverse")
# install.packages("bayesplot")
# install.packages("parallel")

library(brms)
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(parallel)

d <- read.csv("data/winequality-white.csv", sep = ";")

############### M1 (categorical classification) ############### 
f <- quality ~ citric.acid + residual.sugar +
  total.sulfur.dioxide + free.sulfur.dioxide + 
  chlorides + density + pH + sulphates + alcohol +
  fixed.acidity + volatile.acidity
p1 <- prior(normal(0, 5), class = "b", dpar = "mu4") + 
  prior(normal(0, 5), class = "b", dpar = "mu5") + 
  prior(normal(0, 5), class = "b", dpar = "mu6") + 
  prior(normal(0, 5), class = "b", dpar = "mu7") + 
  prior(normal(0, 5), class = "b", dpar = "mu8") + 
  prior(normal(0, 5), class = "b", dpar = "mu9") + 
  prior(normal(0, 5), class = "Intercept", dpar = "mu4") +
  prior(normal(0, 5), class = "Intercept", dpar = "mu5") +
  prior(normal(0, 5), class = "Intercept", dpar = "mu6") +
  prior(normal(0, 5), class = "Intercept", dpar = "mu7") +
  prior(normal(0, 5), class = "Intercept", dpar = "mu8") +
  prior(normal(0, 5), class = "Intercept", dpar = "mu9")

fit1 <- brm(f, 
            data = d, 
            family = categorical(link = "logit"),
            prior = p1,
            iter = 10000,
            save_pars = save_pars(all = TRUE),
            cores = 4)

############### M2 (Gaussian likelihood) ############### 

p2 <- prior(normal(0, 5), class = "b") +
  prior(normal(0, 5), class = "Intercept") +
  prior(normal(0, 5), class = "sigma")

fit2 <- brm(f, 
            data = d,
            family = gaussian(),
            prior = p2,
            iter = 10000,
            save_pars = save_pars(all = TRUE), 
            cores = 4)

############### M3 (Ordinal Regression) ############### 

p3 <- prior(normal(0, 5), class = "b") +
  prior(normal(0, 1), class = "Intercept")

fit3 <- brm(f, 
            data = d,
            family = cumulative("probit"),
            prior = p3,
            save_pars = save_pars(all = TRUE),
            iter = 10000,
            cores = 4,
            init = 0)

############### Model Comparison ############### 
fit1 <- add_criterion(fit1, "loo", ndraws = 20000)
fit2 <- add_criterion(fit2, "loo", ndraws = 20000)
fit3 <- add_criterion(fit3, "loo", ndraws = 20000)

# loo
loo_compare(fit1, fit2, fit3)

# model probability
post_prob(fit1, fit2, fit3)


saveRDS(fit1, "results/fit1.rds")
saveRDS(fit2, "results/fit2.rds")
saveRDS(fit3, "results/fit3.rds")
