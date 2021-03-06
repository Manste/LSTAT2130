library(coda)
library(rjags)
library(ggplot2)
library(R2WinBUGS)

set.seed(33)

# Prepare the Dataset
data <- read.table(file = "cannabis.txt", header = TRUE, sep="")
n <- nrow(data)
y <- sum(data$y)
head(data)

#*************************************
# Question 3a
#*************************************

# the logarithm posterior of pi
lpost<-function(pi, n, y) {
  gama = (1 +  35*pi) / 36
  return (dbinom(y, n, gama, log = TRUE))
}

#Check logarithm posterior function
lpost(0.02, n, y)


M <- 11000
sd_prop <- 0.03
burnin <- 1000

metropolis_algorithm <- function(M, init_pi, sd_prop, n, y, burnin) {
  n_accept = 0
  pi = numeric(M)
  pi[1] = init_pi # init the starting value
  
  # iteration loop
  for (i in 2:M){
    pi_prop = pi[i-1] + rnorm(1, 0, sd_prop)
    prob = min(1, exp(lpost(pi_prop, n, y) - lpost(pi[i-1], n, y)))
    accept = (runif(1) <= prob)
    if (accept){
      pi[i] = pi_prop
      n_accept = n_accept + 1
    } else {
      pi[i] = pi[i-1]
    }
  }
  
  accept_rate = round(n_accept/(M-1),digits=2)
  pi = pi[-c(1:burnin)] # exclude the burnin
  
  return (list(pi = pi, accept_rate= accept_rate))
}

# Run the metroplis algorithm
metropolis <- metropolis_algorithm(M, 0.4, sd_prop, n, y, burnin)
metropolis_mcmc = mcmc(metropolis$pi)
cat("Accpetance rate:", metropolis$accept_rate)

#trace plot
traceplot(mcmc(metropolis$pi))

# others runs
metropolis2 <- metropolis_algorithm(M, 0.05, sd_prop, n, y, burnin)
metropolis3 <- metropolis_algorithm(M, 0.25, sd_prop, n, y, burnin)
metropolis4 <- metropolis_algorithm(M, 0.5, sd_prop, n, y, burnin)

group = list(metropolis_mcmc, 
             mcmc(metropolis2$pi), 
             mcmc(metropolis3$pi), 
             mcmc(metropolis4$pi))

traceplot(group)

# Gelman diagnostic
gelman.plot(group)
gelman.diag(group)

# Geweke diagnostic
geweke.plot(metropolis_mcmc, nbins = 100)
geweke.diag(metropolis_mcmc)

# Question 3b
cat("The 95% credible interval of pi:")
HPDinterval(metropolis_mcmc)
cat("\n\n")

cat("The quantile of pi: ")
quantile(metropolis$pi, probs = c(.025,.05,.5,.95,.975))
cat("\n\n")

# Question 3c
prob = mean(metropolis$pi > .1)
cat("Probability of the proportion of recent cannabis users:", prob)
cat("\n\n")


#*************************************
# Question 4
#*************************************


# A. For men
#*************************************
sub_men <- subset(data, male == 1)
n_men <- nrow(sub_men)
y_men <- sum(sub_men$y)
# Run the metroplis algorithm
metropolis_men <- metropolis_algorithm(M, 0.05, sd_prop, n_men, y_men, burnin)
metropolis_men_mcmc = mcmc(metropolis_men$pi)
cat("Accpetance rate:", metropolis_men$accept_rate)
cat("\n")

#trace plot
traceplot(mcmc(metropolis_men$pi))


# Geweke diagnostic
geweke.plot(metropolis_men_mcmc, nbins = 100)
geweke.diag(metropolis_men_mcmc)

# 3b
cat("The 95% credible interval of pi for men: ")
HPDinterval(metropolis_men_mcmc)
cat("The quantile of pi for men: ")
quantile(metropolis_men$pi, probs = c(.025,.05,.5,.95,.975))
cat("\n")

# 3c
prob_men = mean(metropolis_men$pi > .1)
cat("Probability of the proportion of recent cannabis male users :", prob_men)
cat("\n\n")


# A. For women
#*************************************
sub_women <- subset(data, male == 0)
n_women <- nrow(sub_women)
y_women <- sum(sub_women$y)
# Run the metroplis algorithm
metropolis_women <- metropolis_algorithm(M, .2, .025, n_women, y_women, burnin)
metropolis_women_mcmc = mcmc(metropolis_women$pi)
cat("Accpetance rate:", metropolis_women$accept_rate)
cat("\n\n")

#trace plot
traceplot(mcmc(metropolis_women$pi))


# Geweke diagnostic
geweke.plot(metropolis_women_mcmc, nbins = 100)
geweke.diag(metropolis_women_mcmc)

# 3b
cat("The 95% credible interval of pi for women:")
HPDinterval(metropolis_women_mcmc)
cat("\n\n")
cat("The quantile of pi for women:")
quantile(metropolis_women$pi, probs = c(.025,.05,.5,.95,.975))
cat("\n\n")

# 3c
#generate random values from the posterior function data
prob_women = mean(metropolis_women$pi > .1)
cat("Probability of the proportion of recent cannabis female users:",prob_women)
cat("\n\n")

# C. Comparison between the proportion of men and women
#*************************************
mymodel = function(){
  y_men ~ dbin((1 + 35*pi_men)/36, n1)
  y_women ~ dbin((1 + 35*pi_women)/36, n2)  
  
  pi_men ~ dbeta(1.,1.) # Uniform prior
  pi_women ~ dbeta(1.,1.) # Uniform prior
  
  # Quantities to monitor
  delta <- pi_men - pi_women
  odds1 <- pi_men / (1 - pi_men)
  odds2 <- pi_women / (1 - pi_women)
  gam <- odds1 / odds2
}
model.file = "models/proportionsComparison.bug"
write.model(mymodel,model.file)

# Data and initial values for model parameters
mydata = list(n1=n_men, y_men=y_men, n2=n_women, y_women=y_women)
inits = list(list(pi_men=.2, pi_women=.2))

# Model specification
foo = jags.model(file=model.file, data=mydata, inits=inits,
                 n.chains=1)
update(foo, 1000)

# Generation of the chain
out = coda.samples(model=foo, variable.names=c("delta","gam"), n.iter=10000)

# Inspection of the generated chains
plot(out)
summary(out)
HPDinterval(out)

#*************************************
# Question 5
#*************************************

# Question 5a

#We create the model function
logistic_model <- function(){
  #generate sample
  for(i in 1:n) {
    y[i] ~ dbern((1 + 35*pi[i])/36)
    logit(pi[i]) = (alpha0 + alpha1 * (age[i] - 40)) * male[i] + 
      (beta0 + beta1 * (age[i] - 40)) * (1 - male[i])
  }  
  tau = 1.0E-06
  # non informative prior
  alpha0 ~ dnorm(0, tau)
  alpha1 ~ dnorm(0, tau)
  beta0 ~ dnorm(0, tau)
  beta1 ~ dnorm(0, tau) 
  
  # difference between alpha1 and beta1
  delta = alpha1 - beta1
  
  # odds of 25 years old men
  odds_25 = alpha0 + alpha1 * (25 - 40)
}
model2.file = "models/logisticmodel.bug"
write.model(logistic_model, model2.file)

# Data and initial values for model parameters
mydata2 = list(male=data$male, n=nrow(data), age=data$age, y=data$y)
inits2 = list(list(alpha0=.2, alpha1=.2, beta0=.2, beta1=.2))

# Model specification
foo2 = jags.model(file=model2.file, data=mydata2, inits=inits2,
                 n.chains=1)
update(foo2, 1000)

# Generation of the chain
out2 <- coda.samples(model=foo2, 
      variable.names=c("alpha0","alpha1", "beta0", "beta1", "delta", "odds_25"), 
      n.iter=10000)
out2.matrix <- as.matrix(out2)

# Question 5b
# Inspection of the generated chains
summary(out2)
HPDinterval(out2)

# Question 5c
post_dist_25 = sapply(out2.matrix[, "odds_25"], function(i) 1/(1 - exp(- i)))
plot(density(post_dist_25),
     main = "Density of posterior distribution for a 25 year old male.",
     xlab = parse(text = paste0('~ post_dist_25')))
summary(post_dist_25)
