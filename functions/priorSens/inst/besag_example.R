require(devtools)
require(testthat) # this one is useful for running the tests (see below)

load_all("~/projects/prior_sensitivity/priorSens/")
source("~/Dropbox/Hallo/code_unification_tests/package_plots_MR_20130613.R")

# Run INLA
#-----------
data(Scotland)
Scotland$Region2=Scotland$Region
const<-rep(1,length(Scotland$Region))
Scotland<-data.frame(Scotland,const)
g = system.file("demodata/scotland.graph", package="INLA")
mm0<-0.1
ll0<-0.01
ee0<-1
ff0<-.2/.59

op <- options("warn")
options(warn = -1)

ex4m2.m0l0e0f0<-inla(Counts ~ const + 
                       f(Region2, model="besag", graph=g, hyper=list(theta=list(prior="loggamma",param=c(ee0,ff0))),constr=TRUE)-1,
                     data=Scotland, family="poisson", E=E,
                     control.fixed=list(mean=list(const=mm0),prec=list(const=ll0)),
                     control.compute=list(hyperpar=TRUE),
                     num.threads=1)

ex4m2.m0l0e0f0.hyp<-inla.hyperpar(ex4m2.m0l0e0f0)
ef.log.m0l0e0f0<-ex4m2.m0l0e0f0.hyp$internal.marginals.hyperpar$`Log precision for Region2`

options(op)

# Grid
#----------

grid_epsilon <- 0.05
ex4_ef_number_axis_points <- 200  
ex4_ef_gamma_grid_polar = compute_grid_polar(number_axis_points = 2*ex4_ef_number_axis_points, 
                                             objective_function = "gamma",
                                             hyperpar = c(shape = ee0, rate = ff0),
                                             grid_epsilon = grid_epsilon)

# Run sensitivity
#------------------

# with integration_limits = TRUE
ex4_ef_log_m0l0e0f0_inla_sens_gamma_polar_t = compute_sensitivity(grid = ex4_ef_gamma_grid_polar,
                                                                  inla_marginal_posterior = ef.log.m0l0e0f0, 
                                                                  integration_limits = TRUE)
besag_sensitivity_t <- max(ex4_ef_log_m0l0e0f0_inla_sens_gamma_polar_t[,"sensitivity"])
print(besag_sensitivity_t) 

# with integration_limits = FALSE
ex4_ef_log_m0l0e0f0_inla_sens_gamma_polar_f = compute_sensitivity(grid = ex4_ef_gamma_grid_polar,
                                                                  inla_marginal_posterior = ef.log.m0l0e0f0, 
                                                                  integration_limits = FALSE)
besag_sensitivity_f <- max(ex4_ef_log_m0l0e0f0_inla_sens_gamma_polar_f[,"sensitivity"])
print(besag_sensitivity_f)

# Plot
#---------
par(mfrow = c(2, 2),pty = "s")

plot_polar_sens(ex4_ef_log_m0l0e0f0_inla_sens_gamma_polar_t, xname = expression(nu), yname = expression(xi), angle_shift = (pi/4), modulus_scale = 0.9)
title("TRUE")
plot_line_sens(ex4_ef_log_m0l0e0f0_inla_sens_gamma_polar_t)
title("TRUE")

plot_polar_sens(ex4_ef_log_m0l0e0f0_inla_sens_gamma_polar_f, xname = expression(nu), yname = expression(xi), angle_shift = (pi/4), modulus_scale = 0.9)
title("FALSE")
plot_line_sens(ex4_ef_log_m0l0e0f0_inla_sens_gamma_polar_f)
title("FALSE")

dev.copy2pdf(file="~/Downloads/temp.pdf")