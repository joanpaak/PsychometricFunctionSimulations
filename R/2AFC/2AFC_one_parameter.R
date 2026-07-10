#
# Simulation of a 2AFC task with one parameter.
# Probability of a correct response is...
# 
#   P(R = 1) = pnorm(s / theta_1 / sqrt(2))
#
# ...where theta_1 is taken as the threshold.
# Stimuli are selected randomly between 0.1 and theta_1 * 2,
# implying some prior knowledge of the threshold.
#
# Trial-by-trial mean absolute error between generating
# value and posterior samples is plotted at the end.

setwd("~/Desktop/GITRepositories/PsychometricFunctionSimulations//")
source("Lib/sis_r6.R")

draw_from_prior = function(n){
  cbind(rgamma(n, 2, 0.5))
}

prior = function(theta){
  dgamma(theta[,1], 2, 0.5)
}

likelihood = function(y, theta){
  p_corr = pnorm((y[1,1] / theta[,1]) / sqrt(2))
  return(dbinom(y[1,2], 1, p_corr))
}

n_simulations = 100
n_trials = 100
n_particles = 1000
simulations = list()

for(i in 1:n_simulations){
  cat("Running simulation", i, "\n")
  
  gen_theta = draw_from_prior(1)
  s = runif(n_trials, 0.1, gen_theta[1,1] * 2)
  r = rbinom(n_trials, 1, pnorm((s / gen_theta[1,1]) / sqrt(2)))
  y = cbind(s, r)
  
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
    sis$add_observation(y[j,])
  }
  
  simulations[[i]] = list(
    y = y,
    gen_theta = gen_theta,
    sis = sis
  )
}

sigma_mae = matrix(NaN, ncol = n_trials, nrow = n_simulations)

for(i in 1:n_simulations){
  for(j in 1:n_trials){
    theta = simulations[[i]]$sis$log$particle_set[[j]]$theta
    w = simulations[[i]]$sis$log$particle_set[[j]]$w
    inds = sample(1:length(w), length(w), TRUE, w)
    theta = theta[inds,]
    sigma_mae[i,j] = mean(abs(simulations[[i]]$gen_theta[1] - theta))
  }
}

matplot(t(sigma_mae), type = "l", col = rgb(0, 0, 0, 0.2), lty = 1,
        xlab = "Trial", ylab = "MAE[theta_1]")
