root = "~/projects/prior_sensitivity/priorSens/inst/profiling/"
require(devtools)
require(ggplot2)
require(profr)
require(proftools)
load_all("~/projects/prior_sensitivity/priorSens/")

# devtools::install_github("lineprof")
# devtools::install_github("shiny-slickgrid", "wch")
require(lineprof)
require(shinySlickgrid)
#---------------------------------------------------------------------

mm0<-70; ll0<-0.5; nn<-4; kk<-1; mh<-76

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

Rprof(file.path(root, "test1_profile"))
gaussian_grid_polar = compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                         log_prior_density = "gaussian",
                                         hyperpar = c(mean = mm0, prec = ll0),
                                         grid_epsilon = grid_epsilon, parallel = TRUE)
Rprof(NULL)
summaryRprof(file.path(root, "test1_profile"))

count <- 0
x = compute_grid_polar(number_axis_points = 2*number_axis_points, 
                   log_prior_density = "gaussian",
                   hyperpar = c(mean = mm0, prec = ll0),
                   grid_epsilon = grid_epsilon, parallel = FALSE,
                       method = "uniroot.all")
count

Rprof(file.path(root, "test1_profile"))
replicate(n=20, compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                         log_prior_density = "gaussian",
                                         hyperpar = c(mean = mm0, prec = ll0),
                                         grid_epsilon = grid_epsilon, parallel = FALSE))
Rprof(NULL)
summaryRprof(file.path(root, "test1_profile"))

Rprof(file.path(root, "test1_profile"))
replicate(n=20, compute_grid_polar2(number_axis_points = 2*number_axis_points, 
                                   log_prior_density = "gaussian",
                                   hyperpar = c(mean = mm0, prec = ll0),
                                   grid_epsilon = grid_epsilon, parallel = FALSE))
Rprof(NULL)
summaryRprof(file.path(root, "test1_profile"))

x2 = compute_grid_polar2(number_axis_points = 2*number_axis_points, 
                    log_prior_density = "gaussian",
                    hyperpar = c(mean = mm0, prec = ll0),
                    grid_epsilon = grid_epsilon, parallel = FALSE)

Rprof(file.path(root, "test1_profile"))
gaussian_grid_polar = compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                         log_prior_density = "gaussian",
                                         hyperpar = c(mean = mm0, prec = ll0),
                                         grid_epsilon = grid_epsilon, parallel = FALSE)
Rprof(NULL)
summaryRprof(file.path(root, "test1_profile"))

x<-lineprof(compute_grid_polar(number_axis_points = 2*number_axis_points, 
                               log_prior_density = "gaussian",
                                         hyperpar = c(mean = mm0, prec = ll0),
                                         grid_epsilon = grid_epsilon, parallel = FALSE), torture = FALSE)
shine(x)

time1 = Sys.time()
x=compute_grid_polar2(number_axis_points = 2*number_axis_points, 
                     log_prior_density = "gaussian",
                     hyperpar = c(mean = mm0, prec = ll0),
                     grid_epsilon = grid_epsilon, parallel = FALSE)
time2 = Sys.time()
time2 - time1


Rprof(file.path(root, "test1_profile"))
time1 = Sys.time()
x=compute_grid_polar(number_axis_points = 2*number_axis_points, 
                     log_prior_density = "gaussian",
                                  hyperpar = c(mean = mm0, prec = ll0),
                                  grid_epsilon = grid_epsilon, parallel = TRUE)
time2 = Sys.time()
time2 - time1
Rprof(NULL)
summaryRprof(file.path(root, "test1_profile"))

ex1_grid_profile = profr(compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                            objective_function = "gaussian",
                                            hyperpar = c(mean = mm0, prec = ll0),
                                            grid_epsilon = grid_epsilon))
ggplot(ex1_grid_profile)

# Intial test - 29 sec.
## 'coord_correc_func': Change from data.frame to c(), went to 5.84 sec.
## Disable check_hyperpar: went to 4.36 sec.
## Turns out uniroot is unstable. Need to use uniroot.all - check how to do it faster.
# Use uniroot instead of unitroot.all, need to check with Goscha if every 
#  thing stays the same, later delete the old part: went to 2.02
## Later see f used in uniroot and compute_hellinger
## Also, explore that r should be around value 1.

Rprof(file.path(root, "test1_sens_profile"))
inla_sens_gaussian = compute_sensitivity(grid = gaussian_grid_polar,
                                         inla_marginal_posterior = inla_marginal_posterior, 
                                         integration_limits = TRUE)
Rprof(NULL)
summaryRprof(file.path(root, "test1_sens_profile"))

# propose a simpler integration rule than integrate

max_inla_sens_gaussian = max(inla_sens_gaussian$sensitivity[,"sensitivity"])
expect_equal(max_inla_sens_gaussian, 1.970, tolerance = 1e-2)

sensitivity_obj = prior_sensitivity(objective_function = "gaussian",
                                    hyperpar = c(mean = mm0, prec = ll0),
                                    inla_marginal_posterior = inla_marginal_posterior, 
                                    number_axis_points = 2*number_axis_points,
                                    grid_epsilon = grid_epsilon,
                                    search_interval = c(-100, 40),
                                    integration_limits = TRUE)

expect_identical(sensitivity_obj$polar_grid, gaussian_grid_polar)
expect_identical(sensitivity_obj$sensitivity, inla_sens_gaussian$sensitivity)
expect_equal(max(sensitivity_obj$sensitivity[,"sensitivity"]), 1.970, tolerance = 1e-2)

#-------------------------------------------------------------------------
  
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

# Computing sensibility
#-----------------------

number_axis_points = 200
grid_epsilon = 0.05

Rprof(file.path(root, "test2_profile"))
gamma_grid_polar = compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                      objective_function = "gamma",
                                      hyperpar = c(shape = aa0, rate = bb0),
                                      grid_epsilon = grid_epsilon)
Rprof(NULL)
summaryRprof(file.path(root, "test2_profile"))

  inla_sens_gamma = compute_sensitivity(grid = gamma_grid_polar,
                                        inla_marginal_posterior = ex2.post.log.a0b0, 
                                        integration_limits = TRUE)
  
  max_inla_sens_gamma = max(inla_sens_gamma$sensitivity[,"sensitivity"])
  
  
 expect_equal(max_inla_sens_gamma, 0.2135015, tolerance = 1e-3)
