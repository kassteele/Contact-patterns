#
# Examine prior sensitivity
#

# Load packages
library(INLA)
devtools::load_all("scripts/functions/priorSens")

# Read model output
load("results/polymod.mod.bin")

# Improved estimation of hyper parameters
hyperpar <- inla.hyperpar(polymod.mod)

# INLA log marginal posterior needed for sensitivity computations
hyperpar.log.prec.node <- hyperpar$internal.marginals.hyperpar$"Log precision for node.id"
hyperpar.log.size.negb <- hyperpar$internal.marginals.hyperpar$"log size for the nbinomial observations (overdispersion)"

# Compute grid
grid.prec.node <- compute_grid_polar(number_axis_points = 400, grid_epsilon = 0.00354,
  log_prior_density = "gamma", hyperpar = c(shape = 1, rate = 0.0001))
grid.size.negb <- compute_grid_polar(number_axis_points = 400, grid_epsilon = 0.00354,
  log_prior_density = "gaussian", hyperpar = c(mean = 1, prec = 0.001))

# Compute sensitivity
sens.prec.node <- compute_sensitivity(grid = grid.prec.node,
  inla_marginal_posterior = hyperpar.log.prec.node, integration_limits = TRUE)
sens.size.negb <- compute_sensitivity(grid = grid.size.negb,
  inla_marginal_posterior = hyperpar.log.size.negb, integration_limits = TRUE)
max(sens.prec.node$sensitivity[, "sensitivity"])
max(sens.size.negb$sensitivity[, "sensitivity"])

plot(sens.prec.node$sensitivity[, "sensitivity"])
plot(sens.size.negb$sensitivity[, "sensitivity"])
