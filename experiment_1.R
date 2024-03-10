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
library(MASS)
library(cmdstanr)

d <- read.csv("data/winequality-white.csv", sep = ";")
attach(d)

############### MODEL USED FOR REPORT ############### 
###############  short version (more interpretable) ###############
###############  linear regression ############### 

f_short <- quality ~ citric.acid + volatile.acidity + residual.sugar + sulphates + chlorides + total.sulfur.dioxide + alcohol

p10 <- prior("normal(0,0.36)", class = "b", coef = "alcohol") +
  prior("normal(0,20.268)", class = "b", coef = "chlorides") +
  prior("normal(0,3.659)", class = "b", coef = "citric.acid") +
  prior("normal(0,0.087)", class = "b", coef = "residual.sugar") +
  prior("normal(0,3.88)", class = "b", coef = "sulphates") +
  prior("normal(0,0.01)", class = "b", coef = "total.sulfur.dioxide") +
  prior("normal(0,4.393)", class = "b", coef = "volatile.acidity") +
  prior(normal(6, 5), class = "Intercept") +
  prior(normal(0, 5), class = "sigma")

fit10 <- brm(f_short,
             data = d,
             prior = p10,
             save_pars = save_pars(all=TRUE),
             iter = 4000,
             chains = 4,
             cores = 4)

fit10 <- add_criterion(fit10, "loo", ndraws = 4000)

saveRDS(fit10, "results/short_liner_reg.rds")

###############  cumulative model ############### 

p7 <- prior("normal(0, 4.132)", class = "b", coef = "citric.acid") +
  prior("normal(0,0.099)", class = "b", coef = "residual.sugar") +
  prior("normal(0,3.88)", class = "b", coef = "sulphates") +
  prior("normal(0,22.885)", class = "b", coef = "chlorides") +
  prior("normal(0,0.012)", class = "b", coef = "total.sulfur.dioxide") +
  prior("normal(0,0.406)", class = "b", coef = "alcohol") +
  prior("normal(0,4.961)", class = "b", coef = "volatile.acidity") +
  prior(normal(-2, 1), class = "Intercept", coef = "1") +
  prior(normal(-1.43, 1), class = "Intercept", coef = "2") +
  prior(normal(-0.86, 1), class = "Intercept", coef = "3") +
  prior(normal(-0.29, 1), class = "Intercept", coef = "4") +
  prior(normal(0.29, 1), class = "Intercept", coef = "5") +
  prior(normal(0.86, 1), class = "Intercept", coef = "6") +
  prior(normal(1.43, 1), class = "Intercept", coef = "7") +
  prior(normal(2, 1), class = "Intercept", coef = "8")


fit7 <- brm(f_short,
            data = d,
            family = cumulative("probit"),
            prior = p7,
            init = 0,
            save_pars = save_pars(all=TRUE),
            chains = 4,
            cores = 4,
            iter = 4000)


plot(conditional_effects(fit8, effects = "residual.sugar", categorical = TRUE, plot = FALSE))[[1]]
plot(conditional_effects(fit8, effects = "residual.sugar", method = "posterior_linpred", plot = FALSE))[[1]]
plot(conditional_effects(fit8, effects = "total.sulfur.dioxide", method = "posterior_linpred", plot = FALSE))[[1]]

linear_reg <- readRDS("results/short_liner_reg.rds")
cumlat <- readRDS("results/short_cumulative.rds")
cumlat_s <- readRDS("results/short_cumulative_with_spline.rds")

bmc_obj <- metabmc::metabmc(linear_reg, cumlat, cumlat_s)

###############  cumulative model with spline ############### 

p8 <- prior("normal(0, 4.132)", class = "b", coef = "citric.acid") +
  prior("normal(0,3.88)", class = "b", coef = "sulphates") +
  prior("normal(0,22.885)", class = "b", coef = "chlorides") +
  prior("normal(0,0.406)", class = "b", coef = "alcohol") +
  prior("normal(0,4.961)", class = "b", coef = "volatile.acidity") +
  prior("normal(0, 3)", class = "b", coef = "sresidual.sugar_1") +
  prior("normal(0, 3)", class = "b", coef = "stotal.sulfur.dioxide_1") +
  prior(normal(-2, 1), class = "Intercept", coef = "1") +
  prior(normal(-1.43, 1), class = "Intercept", coef = "2") +
  prior(normal(-0.86, 1), class = "Intercept", coef = "3") +
  prior(normal(-0.29, 1), class = "Intercept", coef = "4") +
  prior(normal(0.29, 1), class = "Intercept", coef = "5") +
  prior(normal(0.86, 1), class = "Intercept", coef = "6") +
  prior(normal(1.43, 1), class = "Intercept", coef = "7") +
  prior(normal(2, 1), class = "Intercept", coef = "8")


f_short_s <- quality ~ s(residual.sugar) + s(total.sulfur.dioxide) + citric.acid + volatile.acidity + sulphates + chlorides + alcohol

fit8 <- brm(f_short_s,
            data = d,
            family = cumulative("probit"),
            prior = p8,
            init = 0,
            save_pars = save_pars(all=TRUE),
            iter = 4000,
            chains = 4,
            cores = 4)

fit7 <- add_criterion(fit7, "loo", ndraws = 4000)
fit8 <- add_criterion(fit8, "loo", ndraws = 4000)

saveRDS(fit7, "results/short_cumulative.rds")
saveRDS(fit8, "results/short_cumulative_with_spline.rds")



make_prior_string <- function(x, tau = 0.5, sd_y =  sd(quality)){
  return (paste("normal(0,", round(tau*sd_y/sd(x), digits = 3), ")", sep = ""))
}
###############  NOT USED FOR REPORT ###############


f <- quality ~ citric.acid + residual.sugar +
  total.sulfur.dioxide + free.sulfur.dioxide + 
  chlorides + density + pH + sulphates + alcohol + 
  fixed.acidity + volatile.acidity

############### categorical classification ############### 
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
            chains = 4,
            save_pars = save_pars(all = TRUE),
            cores = 4)

############### Gaussian likelihood ############### 

# informativeness to be tau
tau <- 0.5

sd_y <- sd(quality)

make_prior_string <- function(x, tau = 0.5, sd_y =  sd(quality)){
  return (paste("normal(0,", round(tau*sd_y/sd(x), digits = 3), ")", sep = ""))
}

p2 <- prior("normal(0,0.36)", class = "b", coef = "alcohol") +
  prior("normal(0,20.268)", class = "b", coef = "chlorides") +
  prior("normal(0,3.659)", class = "b", coef = "citric.acid") +
  prior("normal(0,0.525)", class = "b", coef = "fixed.acidity") +
  prior("normal(0,0.026)", class = "b", coef = "free.sulfur.dioxide") +
  prior("normal(0,2.933)", class = "b", coef = "pH") +
  prior("normal(0,148.055)", class = "b", coef = "density") +
  prior("normal(0,0.087)", class = "b", coef = "residual.sugar") +
  prior("normal(0,3.88)", class = "b", coef = "sulphates") +
  prior("normal(0,0.01)", class = "b", coef = "total.sulfur.dioxide") +
  prior("normal(0,4.393)", class = "b", coef = "volatile.acidity") +
  prior(normal(6, 5), class = "Intercept") +
  prior(normal(0, 5), class = "sigma")

fit2 <- brm(f, 
            data = d,
            family = gaussian(),
            prior = p2,
            iter = 10000,
            save_pars = save_pars(all = TRUE), 
            cores = 4)

############### cumulative model ############### 

p3 <- prior("normal(0,0.406)", class = "b", coef = "alcohol") +
  prior("normal(0,22.885)", class = "b", coef = "chlorides") +
  prior("normal(0,4.132)", class = "b", coef = "citric.acid") +
  prior("normal(0,0.593)", class = "b", coef = "fixed.acidity") +
  prior("normal(0,0.029)", class = "b", coef = "free.sulfur.dioxide") +
  prior("normal(0,3.311)", class = "b", coef = "pH") +
  prior("normal(0,167.173)", class = "b", coef = "density") +
  prior("normal(0,0.099)", class = "b", coef = "residual.sugar") +
  prior("normal(0,4.381)", class = "b", coef = "sulphates") +
  prior("normal(0,0.012)", class = "b", coef = "total.sulfur.dioxide") +
  prior("normal(0,4.961)", class = "b", coef = "volatile.acidity") +
  prior(normal(-2, 1), class = "Intercept", coef = "1") +
  prior(normal(-1.43, 1), class = "Intercept", coef = "2") +
  prior(normal(-0.86, 1), class = "Intercept", coef = "3") +
  prior(normal(-0.29, 1), class = "Intercept", coef = "4") +
  prior(normal(0.29, 1), class = "Intercept", coef = "5") +
  prior(normal(0.86, 1), class = "Intercept", coef = "6") +
  prior(normal(1.43, 1), class = "Intercept", coef = "7") +
  prior(normal(2, 1), class = "Intercept", coef = "8")

# intercept is equi-distanced between -2 and 2 (95% interval for standard normal)

fit3 <- brm(f, 
            data = d,
            family = cumulative("probit"),
            prior = p3,
            save_pars = save_pars(all = TRUE),
            iter = 10000,
            cores = 4,
            init = 0)

pp_check(fit3, ndraws = 50)

############### Model Comparison ############### 
# fit1 <- add_criterion(fit1, "loo", ndraws = 20000)
fit2 <- add_criterion(fit2, "loo", ndraws = 20000)
fit3 <- add_criterion(fit3, "loo", ndraws = 20000)

# loo
loo_compare(fit2, fit3)

# model probability
post_prob(fit2, fit3)


# saveRDS(fit1, "results/fit1.rds")
saveRDS(fit2, "results/fit2.rds")
saveRDS(fit3, "results/fit3.rds")


############### Sequential model ############### 

p4 <- prior(normal(0, 5), class = "b") +
  prior(normal(0, 1), class = "Intercept")

fit4 <- brm(f,
            data = d,
            sratio("cloglog"),
            prior = p4,
            save_pars = save_pars(all = TRUE),
            iter = 10000,
            cores = 4)

pl4_1 <- plot(conditional_effects(fit4, effects = "alcohol", categorical = TRUE, plot = FALSE))[[1]]
pl4_2 <- plot(conditional_effects(fit4, effects = "alcohol", method = "posterior_linpred", plot = FALSE))[[1]]

############### Sequential model with category specific effect ###############

p5 <- prior(normal(0, 5), class = "b") +
  prior(normal(0, 5), class = "Intercept")

f_cs <- quality ~ 1 + cs(sulphates + alcohol + residual.sugar + citric.acid)
fit5 <- brm(f_cs,
            data = d,
            sratio("probit"),
            prior = p5,
            save_pars = save_pars(all = TRUE),
            iter = 10000,
            chains = 4,
            cores = 4)

############### Cumulative model with spline ############### 

f_spline <- quality ~ s(citric.acid) + s(residual.sugar) +
                      s(sulphates) + s(alcohol) +
                      s(total.sulfur.dioxide) + free.sulfur.dioxide +
                      density + pH + chlorides +
                      fixed.acidity + volatile.acidity

p_spline <- prior(normal(0, 1), class = "Intercept") +
  prior(normal(0, 1), class = "sds") + 
  prior(normal(0, 3), class = "b")

fit6 <- brm(f_spline,
            data = d,
            family = cumulative("probit"),
            prior = p_spline,
            iter = 4000,
            init = 0,
            chains = 4,
            cores = 4
            )

summary(fit6)


############### GP ############### 

f_gp <- quality ~ gp(citric.acid, residual.sugar, sulphates, alcohol,
  total.sulfur.dioxide, free.sulfur.dioxide, pH, chlorides,
  fixed.acidity, volatile.acidity)

fit6 <- brm(f_gp,
            data = d,
            family = gaussian(),
            iter = 4000,
            init = 0,
            chains = 4,
            cores = 4)

summary(fit6)

fit3_tmp <- fit3
fit3 <- readRDS("results/fit3.rds")
saveRDS(fit3_tmp, "results/fit3_tmp.rds")




###############  sequential model with category specific effect ############### 

f_cat <- quality ~ cs(citric.acid + volatile.acidity + residual.sugar + sulphates + chlorides + total.sulfur.dioxide + alcohol)


p9 <- prior("normal(0, 4.132)", class = "b", coef = "citric.acid") +
      prior("normal(0,0.099)", class = "b", coef = "residual.sugar") +
      prior("normal(0,3.88)", class = "b", coef = "sulphates") +
      prior("normal(0,22.885)", class = "b", coef = "chlorides") +
      prior("normal(0,0.012)", class = "b", coef = "total.sulfur.dioxide") +
      prior("normal(0,0.406)", class = "b", coef = "alcohol") +
      prior("normal(0,4.961)", class = "b", coef = "volatile.acidity") +
      prior(normal(0, 1), class = "Intercept")

fit9 <- brm(f_cat,
            data = d,
            family = sratio("cloglog"),
            prior = p9,
            chains = 4,
            cores = 4,
            warmup = 2000,
            iter = 10000)



