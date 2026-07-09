#
# Simulation for studying the effect of the number
# of catch trials on a simple Yes/No model in which 
# the probability of a positive response is:
#
# P(R = 1) = pnorm(theta_1 + s / theta_2)
#
# Here theta_1 correspond to what is commonly referred to
# as the criterion and theta_2 to the magnitude of noise; 
# s is the signal level.
#
# A catch trial is selected with a set probability, other
# stimuli are selected using a bayesian adaptive procedure
# that aims to minimize posterior entropy.
#
# Error in estimates is calculated using (approximate) mean
# absolute error (MAE) between generating values and the
# posterior distribution (essentially the MAE between the 
# generating value and samples from the posterior). 
#

setwd("~/Desktop/GITRepositories/PsychometricFunctionSimulations/")
source("Lib/sis_r6.R")

# Information gain, ig
# INPUT
# s     : stimulus level for which ig is evaluated for
# theta : samples from the current posterior
# w     : weights for the aforementioned samples
# OUTPUT
# Information gain for the stimulus level s.
ig = function(s, theta, w){
  p_corr = pnorm(theta[,1] + s / theta[,2])
  
  p_mat = cbind(p_corr, 1 - p_corr)
  E_H = sum(-(p_mat[,1] * log(p_mat[,1]) + 
                p_mat[,2] * log(p_mat[,2])) * w, na.rm = TRUE)
  
  p_vec = c(sum(p_mat[,1] * w), sum(p_mat[,2] * w))
  M_H = -sum(p_vec * log(p_vec), na.rm = TRUE)
  
  return(M_H - E_H)
}

# INPUT
# A catch trial is selected with the probability 
# p_catch_trial, if a catch trial is not selected
# then a stimulus that minimises expected posterior
# entropy is selected. If optimisation fails, stimulus
# level is placed at the current posterior mean for
# theta_2.
#
# sis : a sequential importance sampling object with 
#       properties theta and w (posterior samples and 
#       corresponding weights)
#
getNextStimulus = function(sis, p_catch_trial = 0.2){
  if(runif(1) < p_catch_trial) return(0)
  
  s = tryCatch({
    optimise(function(x, theta, w) -ig(x, theta, w), 
             theta = sis$theta, 
             w = sis$w,
             lower = 0.1, upper = 100)$minimum
  }, error = function(e){
    sis$get_marginal_mus()[2]
  })
  
  return(s)
}

# Functions for the sequential importance sampling
# class.
draw_from_prior = function(n){
  cbind(rnorm(n, -1, 0.5), rgamma(n, 2, 0.5))
}

prior = function(theta){
  dnorm(theta[,1], -1, 0.5) * dgamma(theta[,2], 2, 0.5)
}

likelihood = function(y, theta){
  p_corr = pnorm(theta[,1] + y[1,1] / theta[,2])
  return(dbinom(y[1,2], 1, p_corr))
}

#### Simulation ####

n_simulations = 100
n_trials = 100
n_particles = 1000
simulations = list()

for(i in 1:n_simulations){
  cat("Running simulation", i, "\n")
  
  p_catch_trial = runif(1)
  gen_theta = runif(2, c(-2, 0.5), c(-1, 5.0))
  sis = SIS$new(
    draw_from_prior,
    prior,
    likelihood,
    n_particles,
    opt = list(
      logging = TRUE
    )
  )
  
  for(j in 1:n_trials){
    s = getNextStimulus(sis, p_catch_trial)
    r = rbinom(1, 1, pnorm(gen_theta[1] + s/gen_theta[2]))
    
    sis$add_observation(cbind(s, r))
  }
  
  simulations[[i]] = list(
    gen_theta = gen_theta,
    sis = sis,
    p_catch_trial = p_catch_trial
  )
}

mse = function(gen_theta, theta){
  e_mat = theta
  
  for(i in 1:length(gen_theta)){
    e_mat[,i] = (e_mat[,i] - gen_theta[i])^2
  }
  
  return(mean(e_mat))
}

theta_mse = matrix(NaN, ncol = n_trials, nrow = n_simulations)
sigma_mae = matrix(NaN, ncol = n_trials, nrow = n_simulations)
crit_mae = matrix(NaN, ncol = n_trials, nrow = n_simulations)
sigma_bias = matrix(NaN, ncol = n_trials, nrow = n_simulations)

for(i in 1:n_simulations){
  for(j in 1:n_trials){
    theta = simulations[[i]]$sis$log$particle_set[[j]]$theta
    w = simulations[[i]]$sis$log$particle_set[[j]]$w
    inds = sample(1:length(w), length(w), TRUE, w)
    theta = theta[inds,]
    sigma_mae[i,j] = mean(abs(simulations[[i]]$gen_theta[2] - theta[,2]))
    crit_mae[i,j] = mean(abs(simulations[[i]]$gen_theta[1] - theta[,1]))
    sigma_bias[i,j] = mean(theta[,2]) - simulations[[i]]$gen_theta[2]
    theta_mse[i,j] = mse(simulations[[i]]$gen_theta, theta)
  }
}

x11()
par(mfrow = c(1, 3))
plot(lapply(simulations, function(x) x$p_catch_trial),
     theta_mse[,n_trials], col = rgb(0, 0, 0, 0.5),
     xlab = "P(catch trial)", 
     ylab = "MAE[Theta]")
plot(lapply(simulations, function(x) x$p_catch_trial),
     sigma_mae[,n_trials], col = rgb(0, 0, 0, 0.5),
     xlab = "P(catch trial)",
     ylab = "MAE[theta_2]")
plot(lapply(simulations, function(x) x$p_catch_trial),
     crit_mae[,n_trials], col = rgb(0, 0, 0, 0.5),
     xlab = "P(catch trial)",
     ylab = "MAE[theta_1]")

