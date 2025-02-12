## Read the dataset
oncology_dat <- read.csv2("oncology_ss.txt")

## explore the data
head(oncology_dat)
str(oncology_dat)
hist(as.numeric(oncology_dat$Age))

## fit a linear model
oncology_dat$Age <- as.numeric(oncology_dat$Age)
lm_fit <- lm(Rel_reduction ~ factor(Treatment) + factor(Gender) + factor(Chemo) + Age,
             data = oncology_dat)
summary(lm_fit)
AIC(lm_fit)

## Extend the previous model to check for effect modification of the treatment effect by gender (include an interaction effect between treatment and gender in the model
lm_fit2 <- lm(Rel_reduction ~ factor(Treatment) + factor(Gender) + factor(Chemo) +
                Age + factor(Treatment)*factor(Gender),
              data = oncology_dat)
summary(lm_fit2)
AIC(lm_fit2)

## Bayesian linear regression
library(rstanarm)
oncology_dat$Rel_reduction <- as.numeric(oncology_dat$Rel_reduction)
bayes_lm_fit0 <- stan_glm(Rel_reduction ~ factor(Treatment) + factor(Gender) +
                            factor(Chemo) + Age, data = oncology_dat,
                          family = gaussian(link = "identity"), iter = 10000, 
                          prior = NULL, prior_intercept = NULL, warmup = 5000,
                          chains = 1)
summary(bayes_lm_fit0, digits = 5)
prior_summary(bayes_lm_fit0)

## credible interval
ci <- posterior_interval(bayes_lm_fit0, prob = 0.95, pars = "factor(Treatment)2")
round(ci, 4)

## fit another model with default settings for prior distribution specification
bayes_lm_fit <- stan_glm(Rel_reduction ~ factor(Treatment) + factor(Gender) +
                           factor(Chemo) + Age,
                         data = oncology_dat,
                         family = gaussian(link = "identity"),
                         iter = 10000, warmup = 5000, chains = 1)
bayes_lm_fit
summary(bayes_lm_fit, digits = 5)

## Obtain the MCMC chain and construct a trace plot for the treatment effect
library(bayesplot)
bayesplot::color_scheme_set("viridis")
plot(bayes_lm_fit, plotfun = "trace", pars = "factor(Treatment)2")

plot(bayes_lm_fit, plotfun = "combo", pars = "factor(Treatment)2")

## plot posterior distribution for the age effect
mcmc_dens(bayes_lm_fit, pars = c("Age")) +
  vline_at(0, col="red")
plot(bayes_lm_fit)

## checking the convergence of the MCMC chain(s)
summary(bayes_lm_fit)

## Posterior inference
library(bayestestR)
describe_posterior(bayes_lm_fit)

## change the default settings. for prior distributions
testPriors <- normal(location = c(13, 0, 0, 0),
                     scale = c(2, 100, 100, 100), autoscale = FALSE)
bayes_lm_fit2 <- stan_glm(Rel_reduction ~ factor(Treatment) + factor(Gender) +
                            factor(Chemo) + Age,
                          data = oncology_dat,
                          family = gaussian(link = "identity"),
                          prior = testPriors,
                          iter = 10000, warmup = 5000,
                          chains = 1)
summary(bayes_lm_fit2, digits = 5)
plot(bayes_lm_fit2, fill_color = "skyblue4", est_color = "maroon")

## Fit now a Bayesian linear regression model with 4 chains of 10,000 iterations with a burn-in of 5,000 iterations. Use the previous settings for the specification of the prior distributions (including the informative prior for the treatment effect) and use the default MCMC algorithm.
bayes_lm_fit3 <- stan_glm(Rel_reduction ~ factor(Treatment) + factor(Gender) +
                            factor(Chemo) + Age,
                          data = oncology_dat,
                          family = gaussian(link = "identity"),
                          prior = testPriors,
                          iter = 10000, warmup = 5000,
                          chains = 4)
summary(bayes_lm_fit3)
plot(bayes_lm_fit3)

### Posterior inference for the 4 chains
describe_posterior(bayes_lm_fit3)

posterior_summary <- describe_posterior(bayes_lm_fit3, subset = 1:4)
summary(posterior_summary)
describe_posterior(posterior_summary)

## using the Student t-distribution
bayes_lm_fit4 <- stan_glm(
  Rel_reduction ~ factor(Treatment) + factor(Gender) + factor(Chemo) + Age,
  data = oncology_dat,
  family = gaussian(link = "identity"),
  prior_intercept = student_t(df = 2, location = 13, scale = 2),
  prior = student_t(df = 3, location = 0, scale = 100),
  iter = 10000, warmup = 5000,
  chains = 1
)

summary(bayes_lm_fit4)
plot(bayes_lm_fit4)
describe_posterior(bayes_lm_fit4)

## Bayesian linear regression with interaction effect
bayes_lm_fit_int <- stan_glm(Rel_reduction ~ factor(Treatment) + factor(Gender) +
                               factor(Chemo) + Age + factor(Treatment)*factor(Gender),
                             data = oncology_dat, family = gaussian(link = "identity"),
                             iter = 10000, warmup = 5000, chains = 1)
summary(bayes_lm_fit_int)
describe_posterior(bayes_lm_fit_int)

## how about the clustering of observations within health centers?
### extend to linear mixed model-- a new package 'mmrm' in addition to lme4

fit_re <- stan_glmer(Rel_reduction ~ factor(Treatment) + factor(Gender) +
       factor(Chemo) + Age + factor(Treatment) + 1|Center, data = oncology_dat, family = gaussian(link = "identity"),
       iter = 10000, warmup = 5000, chains = 1)
summary(fit_re)

