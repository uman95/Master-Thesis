###############--Library needed--#################################################
library(R2jags)
library(ggplot2)
library(MCMCpack)
library(lattice)
library(mcmcplots)
library(hydroGOF)
options(digits = 3)

#Importing the data 

Dsets <- read.csv("Data_Abdul.csv", header = TRUE)
Dsets$rate <- (Dsets$Deaths/Dsets$Exposure)
Dsets$Year_f <- (2000-Dsets$Year)
Dsets$sex <- factor(Dsets$Gender, labels = c(0,1))
View(Dsets)
mse
#--------------------------------------------------------------------------------------------------------------#
## Data Exploration
ggplot(Dsets, aes(x=Age, y=Deaths, colour = Gender)) + geom_point() +
  facet_wrap(~Year) + theme(text = element_text(size = 20),
                            axis.text = element_text(size = 12),
                            title = element_text(size = 18))

ggplot(Dsets, aes(x=Age, y=log(rate), colour = Gender)) + geom_point() +
  facet_wrap(~Year) + theme(text = element_text(size = 20),
                            axis.text = element_text(size = 12),
                            title = element_text(size = 18))


## GLM Approach

##########################--NORMAL MODEL--##################################
Model.glm1 <- glm(log(Deaths) ~ offset(log(Exposure)) + Age*Gender + 
                    I(Age^2) + Year_f, data=Dsets,
                  family = gaussian(link = "identity"))
S1 <- summary(Model.glm1)

# Mean Squared Error
mse(exp(Model.glm1$fitted.values)/Dsets$Exposure, Dsets$rate)

# Akaike Information Criteria
S1$aic

# Bayesian Information Criteria
BIC(Model.glm1)

# Output save
sink("Model.glm1.txt")
summary(Model.glm1)$coef
sink()

# Quantile-Quantile plot
dev_res1 <- residuals(Model.glm1, c="deviance")
qqnorm(dev_res1, cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
qqline(dev_res1)

##########################--BINOMIAL MODEL--##################################
Model.glm2 <- glm(cbind(Deaths, Exposure-Deaths) ~ Age*Gender + 
                    I(Age^2) + Year_f, data=Dsets,
                  family = binomial(link = "logit"))

S2 <- summary(Model.glm2)

# Mean Squared Error
mse(Model.glm2$fitted.values, Dsets$rate)

# Akaike Information Criteria
S2$aic

# Bayesian Information Criteria
BIC(Model.glm2)

# Output save
sink("Model.glm2.txt")
summary(Model.glm2)$coef
sink()

# Quantile-Quantile plot
dev_res2 <- residuals(Model.glm2, c="deviance")
qqnorm(dev_res2, cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
qqline(dev_res2)


##########################--POISSON MODEL--##################################
Model.glm3 <- glm(Deaths ~  offset(log(Exposure)) + Age*Gender + 
                    I(Age^2) + Year_f, data = Dsets,  
                  family = poisson(link = "log") )

S3 <- summary(Model.glm3)
Model.glm3$fitted.values
mean((Dsets$Deaths - Model.glm3$fitted.values)^2)
# Mean Squared Error
mse(Model.glm3$fitted.values/Dsets$Exposure, Dsets$rate)

# Akaike Information Criteria
S3$aic

# Bayesian Information Criteria
BIC(Model.glm3)

# Output save
sink("Model.glm3.txt")
summary(Model.glm3)$coef
sink()

# Quantile-Quantile plot
dev_res3 <- residuals(Model.glm3, c="deviance")
qqnorm(dev_res3, cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
qqline(dev_res3)



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##MCMC Approach

####--DATA LIST
Dsetslist <- list(n = dim(Dsets)[1], Deaths = Dsets$Deaths, Exposure = Dsets$Exposure,
                  deathprop=log(Dsets$Deaths/Dsets$Exposure),sex = Dsets$sex,
                  age=Dsets$Age, year_f = Dsets$Year_f)


###################---NORMAL MCMC MODEL---###################################
Model.mod1 <-
#Family = Normal
function() {

for(i in 1:n) {
#Likelihood
deathprop[i] ~ dnorm(mu[i], tau)
mu[i] <- Intercept + Age*age[i] + Age_squared*pow(age[i],2) + Year*year_f[i]
+ (GenderM + Age_GenderM*age[i])*sex[i]

#Residuals
N.residuals[i] <- deathprop[i] - mu[i]
Squared_residuals[i] <- pow(N.residuals[i], 2)
}
# Mean Squared Error
N.MSE <- sum(Squared_residuals[])/n


#Prior Distributions of parameters
Intercept ~ dnorm(0, 0.01)
Age ~ dnorm(0, 0.01)
Age_squared ~ dnorm(0,0.01)
Year ~ dnorm(0, 0.01)
GenderM ~ dnorm(0, 0.01)
Age_GenderM ~ dnorm(0, 0.01)
sigma <- 1/sqrt(tau)
tau ~ dgamma(0.1,0.1)

}

#Parameters initial value
Model.inits1 <- function() {list( Intercept = 0, Age = 0, Age_squared = 0,
                                  Year = 0, GenderM = 0, Age_GenderM = 0)}

#Parameters to save
Model.vars1 <- c("Intercept","Age", "Age_squared", "Year", "GenderM",
                  "Age_GenderM")

#Model fitting
Model.fit1 <- jags(data = Dsetslist, inits= Model.inits1,
                   parameters.to.save = Model.vars1, n.chains = 4,
                   n.iter = 100000,n.burnin = 1000,model.file = Model.mod1)


#Criteria
Model.fit1$BUGSoutput$mean$N.MSE
Model.fit1$BUGSoutput$DIC
Model.fit1$BUGSoutput$pD


sink("Model.fit1.txt")
KK <- print(Model.fit1)
sink()

pdf("Normalplots_summary.pdf")
plot(Model.fit1)
dev.off()


Model.mcmc.fit1 <- as.mcmc(Model.fit1)
LL1 <- summary(Model.mcmc.fit1)

#Posterior mean save
sink("Model.mcmc.fit1.txt")
LL1$statistics
sink()


pdf("Normalplots.pdf")
plot(Model.mcmc.fit1)
dev.off()

#Trace plots of parameters
pdf("Normal_xyplot.pdf")
xyplot(Model.mcmc.fit1, layout = c(2,3), aspect = "fill")
dev.off()

#Density plots of parameters
pdf("Normal_densityplot.pdf")
densityplot(Model.mcmc.fit1, layout = c(2,3), aspect = "fill")
dev.off()



###################---BINOMIAL MCMC MODEL---###################################
Model.mod2 <- 
#Family = Binomial with logit link function
function(){

for(i in 1:n){
#Likelihood
Deaths[i] ~ dbin(pi[i], Exposure[i])
logit(pi[i]) <- Intercept + Age*age[i] + Age_squared*pow(age[i],2) + 
  Year*year_f[i] + (GenderM + Age_GenderM*age[i])*sex[i]

#Pearson Residuals
B.residuals[i] <- (Deaths[i] - (Exposure[i]*pi[i]))/ sqrt(Exposure[i]*pi[i]*(1-pi[i]))
Squared_residuals[i] <- pow(B.residuals[i], 2)
}
#Mean Squared Error
B.MSE <- sum(Squared_residuals[])/n
  
#Prior parameter distribution
Intercept ~ dnorm(0, 0.01)
Age ~ dnorm(0, 0.01)
Age_squared ~ dnorm(0, 0.01)
Year ~ dnorm(0, 0.01)
GenderM ~ dnorm(0, 0.01)
Age_GenderM ~ dnorm(0, 0.01)
}

#Parameters initial value
Model.inits2 <- function() {list(Intercept = -9.2,Age = -0.033,
                                 Age_squared = 0.00190, Year = 0.008,GenderM = 1.37,
                                 Age_GenderM = -0.021)}

#Parameters to save
Model.vars2 <- c("Intercept","Age","Age_squared", "Year", "GenderM",
                 "Age_GenderM")

#Model fitting
Model.fit2 <- jags(data = Dsetslist, inits= Model.inits2,
                   parameters.to.save = Model.vars2, n.chains = 4,
                   n.iter = 100000, n.burnin = 1000,model.file = Model.mod2)

#Criteria
Model.fit2$BUGSoutput$mean$B.MSE
Model.fit2$BUGSoutput$DIC
Model.fit2$BUGSoutput$pD


sink("Model.fit2.txt")
print(Model.fit2)
sink()

pdf("Binomialplots_summary.pdf")
plot(Model.fit2)
dev.off()


Model.mcmc.fit2 <- as.mcmc(Model.fit2)
LL2 <- summary(Model.mcmc.fit2)

#Posterior means save
sink("Model.mcmc.fit2.txt")
LL2$statistics
sink()


pdf("Binomialplots.pdf")
plot(Model.mcmc.fit2)
dev.off()

#Trace plots of parameters
pdf("Binomial_xyplot.pdf")
xyplot(Model.mcmc.fit2, layout = c(2,3), aspect = "fill")
dev.off()

#Density plots of parameters
pdf("Binomial_densityplot.pdf")
densityplot(Model.mcmc.fit2, layout = c(2,3), aspect = "fill")
dev.off()


###################---POISSON MCMC MODEL---###################################
Model.mod3 <-
#Family = Poisson with log link function
function(){
  
for(i in 1:n){
#Likelihood
Deaths[i] ~ dpois(lambda[i])
log(lambda[i]) <- log(Exposure[i]) + Intercept + Age*age[i] +
  Age_squared*pow(age[i],2) + Year*year_f[i] + (GenderM + Age_GenderM*age[i])*sex[i]

#Pearson Residuals
P.residuals[i] <- (Deaths[i] - lambda[i])/ sqrt(lambda[i])
Squared_residuals[i] <- pow(P.residuals[i], 2)
}

#Mean Squared Error
P.MSE <- sum(Squared_residuals[])/n
  
#Prior distributions of parameters
Intercept ~ dnorm(0, 0.01)
Age ~ dnorm(0, 0.01)
Age_squared ~ dnorm(0,0.01)
Year ~ dnorm(0, 0.01)
GenderM ~ dnorm(0, 0.01)
Age_GenderM ~ dnorm(0, 0.01)
}

#Parameters initial value
Model.inits3 <- function() {list( Intercept = 0, Age = 0, Age_squared = 0,
                                  Year = 0, GenderM = 0, Age_GenderM = 0)}

#Parameters to save
Model.vars3 <- c("Intercept","Age", "Age_squared", "Year", "GenderM",
                  "Age_GenderM")

#Model fitting
Model.fit3 <- jags(data = Dsetslist, inits= Model.inits3,
                   parameters.to.save = Model.vars3, n.chains = 4,
                   n.iter = 100000, n.burnin = 1000,model.file = Model.mod3)

#Criteria
Model.fit3$BUGSoutput$mean$P.MSE
Model.fit3$BUGSoutput$DIC
Model.fit3$BUGSoutput$pD


sink("Model.fit3.txt")
print(Model.fit3)
sink()

pdf("Poissonplots_summary.pdf")
plot(Model.fit3)
dev.off()


Model.mcmc.fit3 <- as.mcmc(Model.fit3)
LL3 <- summary(Model.mcmc.fit3)

#Posterior means save
sink("Model.mcmc.fit3.txt")
LL3$statistics
sink()

pdf("Poissonplots.pdf")
plot(Model.mcmc.fit3)
dev.off()

#Trace plots of parameters
pdf("Poisson_xyplot.pdf")
xyplot(Model.mcmc.fit3, layout = c(2,3), aspect = "fill")
dev.off()

#Density plots of parameters
pdf("Poisson_densityplot.pdf")
densityplot(Model.mcmc.fit3, layout = c(2,3), aspect = "fill")
dev.off()

