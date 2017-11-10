# Profiling compute_hellinger_numerically_grid_pertubation

root = "~/projects/prior_sensitivity/priorSens/inst/profiling/"
require(ggplot2)
require(microbenchmark)


log_loggamma <- function(x, hyperpar){
  dgamma(exp(x), shape = exp(hyperpar[1]), rate = exp(hyperpar[2]), log = TRUE) + x
}
# f<-function(x, hyperpar) exp(log_loggamma(x, hyperpar))
# curve(f(x, hyperpar = c(log(0.5), log(0.01))), -10,10)
# g <- function(x) integrate(f=f,lower=-Inf, upper = x, hyperpar = c(log(0.5), log(0.01)))
# 
# 
# find_quantile(quantile = 1e-6, log_dist = log_loggamma, hyperpar = c(log(0.5), log(0.01)))
# find_quantile(quantile = 1 - 1e-6, log_dist = log_loggamma, hyperpar = c(log(0.5), log(0.01)))

quantile1 = find_quantile(quantile = 1e-6, log_dist = log_loggamma, hyperpar = c(log(0.5), log(0.01)))
quantile2 = find_quantile(quantile = 1 - 1e-6, log_dist = log_loggamma, hyperpar = c(log(0.5), log(0.01)))

compute_hellinger_numerically_grid_pertubation(log_dist = log_loggamma, 
                                               hyperpar1 = c(log(0.5), log(0.01)), 
                                               hyperpar2 = c(log(0.51), log(0.011)), 
                                               method = "integrate",
                                               domain = c(-Inf, Inf))

compute_hellinger_numerically_grid_pertubation(log_dist = log_loggamma, 
                                               hyperpar1 = c(log(0.5), log(0.01)), 
                                               hyperpar2 = c(log(0.51), log(0.011)), 
                                               method = "integrate",
                                               domain = c(quantile1, quantile2))

Rprof(file.path(root, "test1_profile"))
replicate(n=10000, compute_hellinger_numerically_grid_pertubation(log_dist = log_loggamma, 
                                               hyperpar1 = c(log(0.5), log(0.01)), 
                                               hyperpar2 = c(log(0.51), log(0.011)), 
                                               method = "simpson13rule",
                                               domain = c(quantile1, quantile2)))
Rprof(NULL)
summaryRprof(file.path(root, "test1_profile"))


m <- microbenchmark(compute_hellinger_numerically_grid_pertubation(log_dist = log_loggamma, 
                                                                   hyperpar1 = c(log(0.5), log(0.01)), 
                                                                   hyperpar2 = c(log(0.51), log(0.011)), 
                                                                   method = "integrate",
                                                                   domain = c(-Inf, Inf)),
                    compute_hellinger_numerically_grid_pertubation(log_dist = log_loggamma, 
                                                                   hyperpar1 = c(log(0.5), log(0.01)), 
                                                                   hyperpar2 = c(log(0.51), log(0.011)), 
                                                                   method = "integrate",
                                                                   domain = c(quantile1, quantile2)),
                    compute_hellinger_numerically_grid_pertubation(log_dist = log_loggamma, 
                                                                   hyperpar1 = c(log(0.5), log(0.01)), 
                                                                   hyperpar2 = c(log(0.51), log(0.011)), 
                                                                   method = "simpson13rule",
                                                                   domain = c(quantile1, quantile2)), times = 1000)
levels(m[,"expr"]) <- c("integrate_inf", "integrate", "simpson13rule")
autoplot(m)  


gs = compute_hellinger_numerically_grid_pertubation(log_dist = log_loggamma, 
                                                    hyperpar1 = c(log(0.5), log(0.01)), 
                                                    hyperpar2 = c(log(0.51), log(0.011)), 
                                                    method = "integrate",
                                                    domain = c(quantile1, quantile2))

points = 20:100
sr <- numeric(length(points))
times <- numeric(length(points))
for (i in 1:length(points)){

time1 <- Sys.time()  
sr[i] <- compute_hellinger_numerically_grid_pertubation(log_dist = log_loggamma, 
                                               hyperpar1 = c(log(0.5), log(0.01)), 
                                               hyperpar2 = c(log(0.51), log(0.011)), 
                                               method = "simpson13rule",
                                               domain = c(quantile1, quantile2),
                                               number_points = points[i])
time2 <- Sys.time()
times[i] <- time2 - time1

}

