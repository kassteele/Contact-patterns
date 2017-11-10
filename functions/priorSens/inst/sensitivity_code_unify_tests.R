source("~/projects/prior_sensitivity/priorSens/inst/old_code.R")
load_all("~/projects/prior_sensitivity/priorSens/")

# Example 1
#------------

## Checking Gaussian grid
##------------------------

mm0<-70
ll0<-0.5

require(rootSolve)

grid_epsilons <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15) # epsilon used in the grid exploration around prior values.

for (grid_epsilon in grid_epsilons){
  
  gaussian_grid = old.compute_grid_gaussian(number_axis_points = 100, # generate approx. 2*number_axis_points
                                            prior_mean = mm0, 
                                            prior_precision = ll0, 
                                            grid_epsilon = grid_epsilon) # Creating grid for the Gaussian using old code
  
  gaussian_grid_polar = compute_grid_polar(number_axis_points = 200, 
                                            objective_function = "gaussian",
                                            hyperpar = c(mean = mm0, prec = ll0),
                                            grid_epsilon = grid_epsilon)
  
  
  x11()
  par(mfrow = c(1,2))
  plot(gaussian_grid[,2], gaussian_grid[,1])
  plot(gaussian_grid_polar$cartesian, main = paste("epsilon", grid_epsilon, sep=""))
    
}

## Complete example with sensitivity
##-----------------------------------

library(INLA)

mm0<-70
ll0<-0.5
nn<-4
kk<-1
mh<-76

n <-nn
const<-rep(1,n)
kappa<-kk
sd.kappa<-sqrt(1/kappa)
mu.dat<-mh
set.seed(100)
y = rnorm(n,mu.dat,sd.kappa)
y<-y-mean(y)+mu.dat
ex1.data<-data.frame(y,const)

# Running INLA
#--------------

ex1.m0l0 = inla(y ~ const-1, data = ex1.data,
                control.fixed=list(mean=list(const=mm0),prec=list(const=ll0)), 
                control.family = list(hyper = list(prec = list(initial=log(kk),fixed=TRUE))),
                num.threads=1)

inla_marginal_posterior <- ex1.m0l0$"marginals.fixed"$"const" # marginal posterior density at (mm0,ll0)

# Computing sensibility
#-----------------------

grid_epsilon <- 0.05 # epsilon used in the grid exploration around prior values.
number_axis_points <- 200  # How many points for a given parameter in the grid exploration, 
# remember total number of points in the grid will be (approximately)
# 2*(axis_points)

old.gaussian_grid = old.compute_grid_gaussian(number_axis_points = number_axis_points,
                                          prior_mean = mm0, 
                                          prior_precision = ll0, 
                                          grid_epsilon = grid_epsilon) # Creating grid for the Gaussian

gaussian_grid_polar = compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                         objective_function = "gaussian",
                                         hyperpar = c(mean = mm0, prec = ll0),
                                         grid_epsilon = grid_epsilon)


old.inla_sens_gaussian = old.compute_sensitivity_gaussian(grid_values = old.gaussian_grid, 
                                                      prior_mean = mm0, 
                                                      prior_precision = ll0, 
                                                      inla_marginal_posterior = inla_marginal_posterior, 
                                                      grid_epsilon = grid_epsilon, 
                                                      integration_limits = TRUE) # Computing sensibility for each point in the grid

inla_sens_gaussian = compute_sensitivity(grid = gaussian_grid_polar,
                                         inla_marginal_posterior = inla_marginal_posterior, 
                                         integration_limits = TRUE)


old.max_inla_sens_gaussian = max(old.inla_sens_gaussian[,"sensibility"])
max_inla_sens_gaussian = max(inla_sens_gaussian[,"sensitivity"])
print(old.max_inla_sens_gaussian); print(max_inla_sens_gaussian)

# Example 2
#------------

## Checking gamma grid
##------------------------

grid_epsilons <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15) # epsilon used in the grid exploration around prior values.

for (grid_epsilon in grid_epsilons){
  
  number_axis_points <- 50  # How many points for a given parameter in the grid exploration, on my current implementation
  # the number of points on the grid can actually be larger that 5 times number_axis_points, so 50
  # is a high enough number
  prior_shape = 0.5
  prior_rate = 0.01
  
  gamma_grid1 = old.compute_grid_gamma(number_axis_points = number_axis_points, 
                                   prior_shape = prior_shape, 
                                   prior_rate = prior_rate, 
                                   grid_epsilon = grid_epsilon, analytic = TRUE) # This is OK.
  
  gamma_grid_polar = compute_grid_polar(number_axis_points = 200, 
                                        objective_function = "gamma",
                                        hyperpar = c(shape = prior_shape, rate = prior_rate),
                                        grid_epsilon = grid_epsilon)
  
  
  # If I use 'analytic = FALSE' above we get unstable results, needs further improvements. For now,
  # use only 'analytic = TRUE'.
  x11()
  par(mfrow = c(1,2))
  plot(gamma_grid1)
  plot(gamma_grid_polar$cartesian, main = paste("epsilon", grid_epsilon, sep=""))
  
}

## Complete example with sensitivity
##-----------------------------------

# Running INLA
#--------------

require(INLA)

## Parameters for the model in Ex2
aa0<-0.5
bb0<-0.01
n <- 10
m <- 2
sigma2.ML<-0.035
ntsigma2.ML<-n*sigma2.ML
set.seed(100)
y <- rnorm(n, mean = m, sd = sqrt(bb0 / aa0))
y <- m + sqrt(sigma2.ML) * (y - mean(y)) / (sqrt((n - 1) / n) * sd(y))
stopifnot(all.equal(ntsigma2.ML,
                    sum((y - m)^2)))
ex2.data <- data.frame(y, const=1)
ex2.a0b0 <- inla(y ~ -1 + const,
                 data = ex2.data,
                 control.fixed=list(mean=list(const=m),prec=list(const=1e12)),
                 control.family=list(hyper =list(prec =list(prior="loggamma",param=c(aa0, bb0)))),
                 control.compute=list(hyperpar=TRUE),
                 control.predictor = list(initial = 15),
                 num.threads=1)
ex2.a0b0.hyp <- inla.hyperpar(ex2.a0b0, dz = 0.25, diff.logdens = 15)

ex2.post.log.a0b0<-ex2.a0b0.hyp$"internal.marginals.hyperpar"$`Log precision for the Gaussian observations`
exact_posterior <- cbind(ex2.post.log.a0b0[,1], dgamma(exp(ex2.post.log.a0b0[,1]), shape = aa0 + n/2, rate = bb0 + 0.5 * ntsigma2.ML) * exp(ex2.post.log.a0b0[,1]))

curve(inla.dmarginal(x, marginal=ex2.post.log.a0b0), 0, 5)
curve(inla.dmarginal(x, marginal=exact_posterior), add=T, col = "red")

# Computing sensibility
#-----------------------

number_axis_points = 200
grid_epsilon = 0.05

old.gamma_grid = old.compute_grid_gamma(number_axis_points = number_axis_points, 
                                     prior_shape = aa0, 
                                     prior_rate = bb0, 
                                     grid_epsilon = grid_epsilon, 
                                     analytic = TRUE) # This is OK.

gamma_grid_polar = compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                         objective_function = "gamma",
                                         hyperpar = c(shape = aa0, rate = bb0),
                                         grid_epsilon = grid_epsilon)


old.inla_sens_gamma_analytical = old.compute_analytical_sensitivity_gamma_ex2(grid_values = gamma_grid_polar$cartesian, 
                                                                      prior_shape = aa0, 
                                                                      prior_rate = bb0, 
                                                                      sample_size = n,
                                                                      sum_squares = ntsigma2.ML,
                                                                      grid_epsilon = grid_epsilon, 
                                                                      integration_limits = TRUE)

inla_sens_gamma = compute_sensitivity(grid = gamma_grid_polar,
                                      inla_marginal_posterior = ex2.post.log.a0b0, 
                                      integration_limits = TRUE)





old.max_inla_sens_gamma_analytical = max(old.inla_sens_gamma_analytical[,3])
max_inla_sens_gamma = max(inla_sens_gamma[,"sensitivity"])


print(old.max_inla_sens_gamma_analytical)
print(max_inla_sens_gamma)

# Compare analytical and numerical solutions
x = cbind(inla_sens_gamma, old.inla_sens_gamma_analytical)
plot(x[, 5], x[, 8])

# Looking at the plot we see a systematic deviation of some points, lets separate good from bad cases
index = ((x[, 5] - x[, 8]) > 0.04)
problematic_cases = cbind(x[index, ], diff = x[index, 5] - x[index, 8])
good_cases = cbind(x[!index, ], diff = x[!index, 8] - x[!index, 8])

#### plot the prior ratios
cases = good_cases[,3:4]

i=1
compute_ratio_priors_gamma(hyperpar_new = c(shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), 
                           hyperpar_old = c(shape = aa0, rate = bb0), inla_internal_marginal = ex2.post.log.a0b0)
plot(compute_ratio_priors_gamma(hyperpar_new = c(shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), 
                                  hyperpar_old = c(shape = aa0, rate = bb0), inla_internal_marginal = ex2.post.log.a0b0), ylab="", ylim = c(0.8, 1.2))  
for (i in 2:79){
print(i)
  compute_ratio_priors_gamma(hyperpar_new = c(shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), 
                           hyperpar_old = c(shape = aa0, rate = bb0), inla_internal_marginal = ex2.post.log.a0b0)
points(compute_ratio_priors_gamma(hyperpar_new = c(shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), 
                                  hyperpar_old = c(shape = aa0, rate = bb0), inla_internal_marginal = ex2.post.log.a0b0))  
Sys.sleep(1)
}

cases = problematic_cases[,3:4]

for (i in 1:nrow(cases)){
  print(i)
  compute_ratio_priors_gamma(hyperpar_new = c(shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), 
                             hyperpar_old = c(shape = aa0, rate = bb0), inla_internal_marginal = ex2.post.log.a0b0)
  points(compute_ratio_priors_gamma(hyperpar_new = c(shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), 
                                    hyperpar_old = c(shape = aa0, rate = bb0), inla_internal_marginal = ex2.post.log.a0b0), col = "red")  
  Sys.sleep(1)
}

cases = good_cases[,3:4]

for (i in 80:nrow(cases)){
  print(i)
  compute_ratio_priors_gamma(hyperpar_new = c(shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), 
                             hyperpar_old = c(shape = aa0, rate = bb0), inla_internal_marginal = ex2.post.log.a0b0)
  points(compute_ratio_priors_gamma(hyperpar_new = c(shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), 
                                    hyperpar_old = c(shape = aa0, rate = bb0), inla_internal_marginal = ex2.post.log.a0b0), col = "blue")  
  Sys.sleep(1)
}



#### plot the densities
f <- function(x, shape, rate) dgamma(exp(x), shape = shape, rate = rate, log = TRUE) + x
curve(f(x, shape = aa0, rate = bb0), -10, 10)

cases = good_cases[,3:4]

for (i in 1:79){
  print(i)
  curve(f(x, shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), add = T, col = "blue")
  #Sys.sleep(1)
}

cases = problematic_cases[,3:4]

for (i in 1:nrow(cases)){
  print(i)
  curve(f(x, shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), add = T, col = "red")
  #Sys.sleep(1)
}

cases = good_cases[,3:4]

for (i in 80:nrow(cases)){
  print(i)
  curve(f(x, shape = as.numeric(cases[i, 1]), rate = as.numeric(cases[i, 2])), add = T, col = "green")
  #Sys.sleep(1)
}

SHAPE = good_cases[1,3]
RATE = good_cases[1,4]

# SHAPE = problematic_cases[1,3]
# RATE = problematic_cases[1,4]

# posterior and corrected posterior for a problematic case
hyperpar_new = c(shape = as.numeric(SHAPE), rate = as.numeric(RATE))
hyperpar_old = c(shape = aa0, rate = bb0)
inla_marginal_posterior = ex2.post.log.a0b0

prior_shape = aa0
prior_rate = bb0 
new_shape = as.numeric(SHAPE)
new_rate = as.numeric(RATE)

sample_size = n
sum_squares = ntsigma2.ML

corrected_posterior_marginal = compute_corrected_posteriors(prior_type = "gamma", 
                                                            hyperpar_new = hyperpar_new,
                                                            hyperpar_old = hyperpar_old,
                                                            inla_marginal_posterior = inla_marginal_posterior, 
                                                            integration_limits = TRUE)
# posterior
f <- function(x) dgamma(exp(x), shape = prior_shape + sample_size/2, rate = prior_rate + 0.5 * sum_squares) * exp(x)
seqq = inla_marginal_posterior[,1]
inla_style_post = matrix(c(seqq, f(seqq)), length(seqq), 2)
       
curve(f(x), 1.5, 5)
curve(inla.dmarginal(x, marginal=inla_style_post), add = T, col = "red")
curve(inla.dmarginal(x, marginal=inla_marginal_posterior), add = T, col = "red")

compute_hellinger_distance_inla(inla_posterior_1 = inla_marginal_posterior, 
                                inla_posterior_2 = inla_style_post,
                                integration_limits = TRUE)

# corrected posterior
f <- function(x) dgamma(exp(x), shape = new_shape + sample_size/2, rate = new_rate + 0.5 * sum_squares) * exp(x)
seqq = corrected_posterior_marginal[,1]
inla_style_correc_post = matrix(c(seqq, f(seqq)), length(seqq), 2)

curve(f(x), 1.5, 5)
curve(inla.dmarginal(x, marginal=inla_style_correc_post), add = T, col = "red")
curve(inla.dmarginal(x, marginal=corrected_posterior_marginal), add = T, col = "red")

compute_hellinger_distance_inla(inla_posterior_1 = corrected_posterior_marginal, 
                                inla_posterior_2 = inla_style_correc_post,
                                integration_limits = TRUE)

sensitivity_vector[i] <- compute_hellinger_distance_inla(inla_posterior_1 = inla_marginal_posterior, 
                                                         inla_posterior_2 = corrected_posterior_marginal,
                                                         integration_limits = TRUE)/grid_epsilon

# analytical
sensibility_vector[i] <- old.compute_hellinger_gammas(prior_shape = prior_shape + sample_size/2, 
                                                      prior_rate = prior_rate + 0.5 * sum_squares, 
                                                      new_shape = new_shape + sample_size/2, 
                                                      new_rate = new_rate + 0.5 * sum_squares, 
                                                      analytic = TRUE)/grid_epsilon


