#
# 1957 pandemic flu mortality
#

#
# Init
#

# Load packages
library(mgcv)
library(RColorBrewer)

#
# Data
#

# Read 1957 mortality and population data
mort.data <- read.delim(file = "data/population data/CBS age sex pop deaths 1957.txt")

# Modify mort.data
mort.data <- within(mort.data, {
  sex <- factor(sex, levels = c("man", "vrouw"), labels = c("male", "female"))
  age <- rep(0:99, 2)
  age.cat <- cut(age, breaks = c(0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, Inf), include.lowest = TRUE, right = FALSE)
  agepop <- age*pop
})

# Aggregate mort.data
mort.data.agg <- aggregate(cbind(pop, deaths, agepop) ~ age.cat + sex, FUN = sum, data = mort.data)
mort.data.agg <- within(mort.data.agg, {
  age <- agepop/pop
  rm(agepop)
})

# Set data influenza mortality
# source: www.ntvg.nl/publicatie/influenzasterfte-de-herfst-van-1957/pdf
inf.data <- with(mort.data, expand.grid(age.cat = levels(age.cat), sex = levels(sex)))
inf.data <- within(inf.data, {
  inf <- c(11, 24, 23, 24, 18, 14, 19, 36, 73, 123, 132, 117, # male
    6, 30, 29, 34, 24, 33, 30, 31, 45, 105, 134, 115) # female
})

# Merge inf.data and mort.data.agg
inf.data <- merge(inf.data, mort.data.agg)

#
# Fit models
#

# Models
tot.mod <- with(mort.data, gam(deaths ~ s(age, by = sex, bs = "ps", k = 30, m = c(2, 2)), offset = log(pop), family = poisson))
inf.mod <- with(inf.data,  gam(inf    ~ s(age, by = sex, bs = "ps", k = 12, m = c(2, 2)), offset = log(pop), family = poisson))

# pred.data = subset mort.data with age <=80
pred.data <- subset(mort.data, subset = age <= 80)
n <- 81

# predict simulations
pred.fun <- function(object, newdata, n.sam = 1000) {
  beta <- coef(object); V.beta <- vcov(object)
  beta.sam <- t(MASS::mvrnorm(n.sam, beta, V.beta))
  X <- predict(object, newdata = newdata, type = "lpmatrix")
  mu <- exp(X%*%beta.sam)
  return(mu)
}
mu.tot.sam <- pred.fun(tot.mod, pred.data)
mu.inf.sam <- pred.fun(inf.mod, pred.data)
p.sam <- mu.inf.sam/mu.tot.sam

# Calulate ratio between females and males (=reference)
ratio.p.sam <- p.sam[(n+1):(2*n), ]/p.sam[1:n, ]

#
# Plot
#


# Settings
ticks <- seq(0, 80, 10) # Tick locations
cols.line <- brewer.pal(n = 8, "Accent")[c(5, 6)] # Color for males and females: lines
cols.bord <- adjustcolor(cols.line, alpha = 0.5)

# Calculate statistics for p.sam and difference M-F
p.stat <- t(apply(p.sam, 1, quantile, prob = c(0.025, 0.50, 0.975)))
ratio.p.stat <- t(apply(ratio.p.sam, 1, quantile, prob = c(0.025, 0.5, 0.975)))

# Figure layout
layout(matrix(c(1, 2), 1, 2), widths = c(1, 1), heights = 1, respect = TRUE)
par(mgp = c(2.5, 1, 0), xaxs = "r", yaxs = "r", cex.axis = 0.9)

# Males and Female separate
par(mar = c(4, 4, 0.3, 0.3))
plot.new(); plot.window(xlim = c(0, 80), ylim = c(0, 0.155))
polygon(c(age, rev(age)), c(p.stat[1:n, 1],         rev(p.stat[1:n, 3])),         col = cols.bord[1], border = NA)
polygon(c(age, rev(age)), c(p.stat[(n+1):(2*n), 1], rev(p.stat[(n+1):(2*n), 3])), col = cols.bord[2], border = NA)
matlines(age, cbind(p.stat[1:n, 2], p.stat[(n+1):(2*n), 2]), col = cols.line[1:2], lwd = 1, lty = 1, lend = "butt")
title(xlab = "Age", ylab = "Fraction influenza mort. / total mort.")
axis(1, at = ticks, labels = ifelse(ticks%%20!=0 , NA, ticks), lwd = 0.5)
axis(2, lwd = 0.5)
box(lwd = 0.5)
legend("topright", legend = c("Males", "Females"), fill = cols.bord, border = cols.line, bty = "n", cex = 0.9)
legend("topleft", legend = "a", bty = "n", adj = 2)

# Difference Males - Females
par(mar = c(4, 4, 0.3, 0.3))
plot.new(); plot.window(xlim = c(0, 80), ylim = range(ratio.p.stat))
title(xlab = "Age", ylab = "Mortality ratio females vs. males")
polygon(c(age, rev(age)), c(ratio.p.stat[, 1], rev(ratio.p.stat[, 3])), col = grey(0.8), border = NA)
lines(age, ratio.p.stat[, 2], col = 1, lwd = 1, lty = 1, lend = "butt")
abline(h = 1, lty = 3)
axis(1, at = ticks, labels = ifelse(ticks%%20!=0 , NA, ticks), lwd = 0.5)
axis(2, lwd = 0.5)
box(lwd = 0.5)
legend("topleft", legend = "b", bty = "n", adj = 2)
