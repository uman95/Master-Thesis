library(R2WinBUGS, quietly = TRUE)
library(ggplot2)
library(rjags)
library(runjags)
library(MCMCpack)
library(lattice)
library(mcmcplots)

Dsets <- read.csv("Data_Abdul.csv", header = TRUE)
Dsets$rate <- (Dsets$Deaths/Dsets$Exposure)
Dsets$Year_f <- (2000-Dsets$Year)
Dsets$sex <- factor(Dsets$Gender, labels = c(0,1))
View(Dsets)
str(Dsets)

#--------------------------------------------------------------------------------------------------------------#
## Data Exploration
ggplot(Dsets, aes(x=Age, y=Deaths, colour = Gender)) + geom_point() +facet_wrap(~Year)
ggplot(Dsets, aes(x=Age, y=rate, colour = Gender)) + geom_point() +facet_wrap(~Year)
ggplot(Dsets, aes(x=Age, y=log(rate), colour = Gender)) + geom_point() +facet_wrap(~Year)



## GLM Approach

Model.glm1 <- glm(log(Deaths) ~ offset(log(Exposure)) + Gender + Age + Year_f, data=Dsets,
                  family = gaussian(link = "identity"))
summary(Model.glm1)

par(mfrow=c(1,1))

BIC(Model.glm1)

sink("Model.glm1.txt")
summary(Model.glm1)
sink()
pred_prob1 <- Model.glm1$fitted.values
eta_hat1 <- Model.glm1$linear.predictors
dev_res1 <- residuals(Model.glm1, c="deviance")
qqnorm(dev_res1)
qqline(dev_res1)


Model.glm2 <- glm(cbind(Deaths, Exposure-Deaths) ~ Gender + Age + Year_f, data=Dsets,
                  family = binomial(link = "logit"))

summary(Model.glm2)
BIC(Model.glm2)

sink("Model.glm2.txt")
summary(Model.glm2)
sink()

pred_prob2 <- Model.glm2$fitted.values
eta_hat2 <- Model.glm2$linear.predictors
dev_res2 <- residuals(Model.glm2, c="deviance")
qqnorm(dev_res2)
qqline(dev_res2)

Model.glm3 <- glm(Deaths ~  offset(log(Exposure)) + Gender + Age + Year_f,
                  data = Dsets,  family = poisson(link = "log") )

summary(Model.glm3)
BIC(Model.glm3)

sink("Model.glm3.txt")
summary(Model.glm3)
sink()

pred_prob3 <- Model.glm3$fitted.values
eta_hat3 <- Model.glm3$linear.predictors
dev_res3 <- residuals(Model.glm3, c="deviance")
qqnorm(dev_res3)
qqline(dev_res3)

anova(Model.glm3, test= "Chi")

##MCMC Approach

Dsetslist <- list(n = dim(Dsets)[1], Deaths = Dsets$Deaths, Exposure = Dsets$Exposure,
                  deathprop=log(Dsets$Deaths/Dsets$Exposure),sex = Dsets$Gender,
                  age=Dsets$Age, year_f = Dsets$Year_f)


Model.mcmc1 <- "
#Family = Normal
model {
#Likelihood
for(i in 1:n) {
deathprop[i] ~ dnorm(mu[i], tau)
mu[i] <- Intercept + Age*age[i] + Year*year_f[i] + GenderM*sex[i]
}
#Prior Distributions of parameters
Intercept ~ dnorm(0, 0.01)
Age ~ dnorm(0, 0.01)
Year ~ dnorm(0, 0.01)
GenderM ~ dnorm(0, 0.01)
sigma <- 1/sqrt(tau)
tau ~ dgamma(0.1,0.1)

}"

Model.mcmc1.jags <- file.path(getwd(), "Model.mcmc1.jags")
writeLines(Model.mcmc1, con = "Model.mcmc1.jags")

Model.mcmc1.inits <- list( Intercept = 0, Age = 0, Year = 0, GenderM = 0)

Model.mcmc1.jags <- jags.model("Model.mcmc1.jags",Dsetslist, 
                               inits = Model.mcmc1.inits)
Model.mcmc1.vars <- c("Intercept","Age", "Year", "GenderM")
Model.mcmc1.sim <- coda.samples(Model.mcmc1.jags, Model.mcmc1.vars,
                                n.iter = 10000)
summary(Model.mcmc1.sim)

sink("Model.mcmc1.txt")
summary(Model.mcmc1.sim)
sink()
xyplot(Model.mcmc1.sim)
densityplot(Model.mcmc1.sim)

plot(Model.mcmc1.sim)

Model.mcmc2 <- "
#Family = Poisson with log link function
model{
#Likelihood
for(i in 1:n){
Deaths[i] ~ dpois(lambda[i])
log(lambda[i]) <- log(Exposure[i]) + Intercept + Age*age[i] + Year*year_f[i]
                    + GenderM*sex[i]
}
#Prior distributions of parameters
Intercept ~ dnorm(0, 0.0001)
Age ~ dnorm(0, 0.0001)
Year ~ dnorm(0, 0.0001)
GenderM ~ dnorm(0, 0.0001)
}"

Model.mcmc2.jags <- file.path(getwd(), "Model.mcmc2.jags")
writeLines(Model.mcmc2, con = "Model.mcmc2.jags")

Model.mcmc2.inits <- list( Intercept = 0, Age = 0, Year = 0, GenderM = 0)

jags
Model.mcmc2.jags <- jags.model("Model.mcmc2.jags",Dsetslist, 
                               inits = Model.mcmc1.inits)
Model.mcmc2.vars <- c("Intercept","Age", "Year", "GenderM")
Model.mcmc2.sim <- coda.samples(Model.mcmc2.jags, Model.mcmc2.vars,
                                n.iter = 10000)
coda
sink("Model.mcmc2.txt")
summary(Model.mcmc2.sim)
sink()


summary(Model.mcmc2.sim)

plot(Model.mcmc2.sim)

Model.mcmc3 <- "
#Family = Binomial with logit link function
model{
#Likelihood
for(i in 1:n){
Deaths[i] ~ dbin(pi[i], Exposure[i])
probit(pi[i]) <- Intercept + Age*age[i] + Year*year_f[i] + GenderM*sex[i]
}
#Prior parameter distribution
Intercept ~ dnorm(0, 0.0001)
Age ~ dnorm(0, 0.0001)
Year ~ dnorm(0, 0.0001)
GenderM ~ dnorm(0, 0.0001)
}"

Model.mcmc3.jags <- file.path(getwd(), "Model.mcmc3.jags")
writeLines(Model.mcmc3, con = "Model.mcmc3.jags")

Model.mcmc3.inits <- list( Intercept = 0, Age = 0, Year = 0, GenderM = 0)

Model.mcmc3.jags <- jags.model("Model.mcmc3.jags",Dsetslist, 
                               inits = Model.mcmc1.inits)
Model.mcmc3.vars <- c("Intercept","Age", "Year", "GenderM")
Model.mcmc3.sim <- coda.samples(Model.mcmc2.jags, Model.mcmc3.vars,
                                n.iter = 10000)
sink("Model.mcmc3.txt")
summary(Model.mcmc3.sim)
sink()


summary(Model.mcmc3.sim)
#summary(jags.samples(Model.mcmc2.jags, Model.mcmc3.vars,
#                    n.iter = 10000))

traceplot(Model.mcmc3.sim)
densplot(Model.mcmc2.sim)
xyplot(Model.mcmc3.sim)

geweke.plot(Model.mcmc3.sim)
densityplot(Model.mcmc3.sim)

