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

group = list(metropolis_mcmc, mcmc(metropolis2$pi), mcmc(metropolis3$pi), mcmc(metropolis4$pi))

traceplot(group)

# Gelman diagnostic
gelman.plot(group)
gelman.diag(group)

# Geweke diagnostic
geweke.plot(metropolis_mcmc, nbins = 100)
geweke.diag(metropolis_mcmc)

# Question 3b
HPD_pi <- HPDinterval(metropolis_mcmc)
CI_pi <- quantile(metropolis$pi, probs = c(.025,.05,.5,.95,.975))
cat("The 95% credible interval of pi: [", HPD_pi, "]" )
cat("The quantile of pi: [", CI_pi, "]" )

# Question 3c
sub_data <- subset(data, age >= 20 & age <= 50)
tmp_n <- nrow(sub_data)

#generate random values from the posterior function data
N <- 10000
gama <- (1 + 35*0.09 ) / 36
pi_post <- rbinom(N, tmp_n, gama)
prob = mean(pi_post > 0.1)
cat("The probability that the proportion of recent cannabis users is at least 10% amoung 20-59 years olds is:", prob)


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

#trace plot
traceplot(mcmc(metropolis_men$pi))


# Geweke diagnostic
geweke.plot(metropolis_men_mcmc, nbins = 100)
geweke.diag(metropolis_men_mcmc)

# 3b
HPD_pi_men <- HPDinterval(metropolis_men_mcmc)
CI_pi_men <- quantile(metropolis_men$pi, probs = c(.025,.05,.5,.95,.975))
cat("The 95% credicle interval of pi for men: [", HPD_pi_men, "]" )
cat("The quantile of pi for men: [", CI_pi_men, "]" )

# 3c
sub_men_data <- subset(sub_men, age >= 20 & age <= 50)
tmp_men_n <- nrow(sub_men_data)

#generate random values from the posterior function data
pi_post_mean = mean(metropolis_men_mcmc) #posterior mean
gama_men <- (1 + 35*pi_post_mean ) / 36
pi_post_men <- rbinom(N, tmp_men_n, gama_men)
prob_men = mean(pi_post_men > 0.1)
cat("The probability that the proportion of recent cannabis male users is at least 10% amoung 20-59 years olds is:", prob_men)


# A. For women
#*************************************
sub_women <- subset(data, male == 0)
n_women <- nrow(sub_women)
y_women <- sum(sub_women$y)
# Run the metroplis algorithm
metropolis_women <- metropolis_algorithm(M, 0.2, 0.025, n_women, y_women, burnin)
metropolis_women_mcmc = mcmc(metropolis_women$pi)
cat("Accpetance rate:", metropolis_women$accept_rate)

#trace plot
traceplot(mcmc(metropolis_women$pi))


# Geweke diagnostic
geweke.plot(metropolis_women_mcmc, nbins = 100)
geweke.diag(metropolis_women_mcmc)

# 3b
HPD_pi_women <- HPDinterval(metropolis_women_mcmc)
CI_pi_women <- quantile(metropolis_women$pi, probs = c(.025,.05,.5,.95,.975))
cat("The 95% credicle interval of pi for women: [", HPD_pi_women, "]" )
cat("The quantile of pi for women: ", CI_pi_women, "]" )

# 3c
sub_women_data <- subset(sub_women, age >= 20 & age <= 50)
tmp_women_n <- nrow(sub_women_data)

#generate random values from the posterior function data
pi_post_mean = mean(metropolis_women_mcmc) #posterior mean
gama_women <- (1 + 35*pi_post_mean) / 36
pi_post_women <- rbinom(N, tmp_women_n, gama_women)
prob_women = mean(pi_post_women > 0.1)
cat("The probability that the proportion of recent cannabis female users is at least 10% amoung 20-59 years olds is:", prob_women)


# C. Comparison between the proportion of men and women
#*************************************
mymodel = function(){
  y1 ~ dbin(pi_men, n1)
  y2 ~ dbin(pi_women, n2)  
  pi_men ~ dbeta(1.,1.) # Uniform prior
  pi_women ~ dbeta(1.,1.) # Uniform prior
  
  # Quantities to monitor
  delta <- pi_men - pi_women
  odds1 <- pi_men / (1 - pi_men)
  odds2 <- pi_women / (1 - pi_women)
  gam <- odds1 / odds2
}
model.file = "proportionsComparison.bug"
write.model(mymodel,model.file)

# Data and initial values for model parameters
mydata = list(n1=n_men, y1=y_men, n2=n_women, y2=y_women)
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
  for(i in 1:N) {
    y[i] 
  }
  tau = 1.0E-06
  # non informative prior
  alpha0 ~ dnorm(0, tau)
  alpha1 ~ dnorm(0, tau)
  beta0 ~ dnorm(0, tau)
  beta1 ~ dnorm(0, tau) 
}

# Question 5b


# Question 5c