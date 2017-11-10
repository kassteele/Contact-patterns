# profiling compute_grid_polar

require(ggplot2)
require(microbenchmark)

m <- microbenchmark(compute_grid_polar(number_axis_points = number_points, 
                                       log_prior_density = log_loggamma,
                                       grid_epsilon = grid_epsilon,
                                       hyperpar = c(log(shape), log(rate)),
                                       integration_method = "integrate",
                                       domain = c(-Inf, Inf), 
                                       parallel = FALSE),
                    compute_grid_polar(number_axis_points = number_points, 
                                       log_prior_density = log_loggamma,
                                       grid_epsilon = grid_epsilon,
                                       hyperpar = c(log(shape), log(rate)),
                                       integration_method = "integrate",
                                       domain = "compute_domain",
                                       parallel = FALSE),
                    compute_grid_polar(number_axis_points = number_points, 
                                       log_prior_density = log_loggamma,
                                       grid_epsilon = grid_epsilon,
                                       hyperpar = c(log(shape), log(rate)),
                                       integration_method = "simpson13rule",
                                       domain = "compute_domain",
                                       parallel = FALSE), times = 10)
levels(m[,"expr"]) <- c("integrate_Inf", "integrate_quantiles", "simpson_quantiles")
autoplot(m)
