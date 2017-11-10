#
# Eigenvectors
#

#
# Init
#

# Load packages
library(RColorBrewer)

# Source functions
source(file = "functions/make.m.matrix.R")

#
# Data
#

# Read data
load(file = "results/polymod.tab.bin")
load(file = "results/linpred.sam.bin")

# Get ages and number of age classes
age <- unique(polymod.tab$part.age)
n.age <- length(age)

# Compute contact rates and intensities
c.sam <- exp(linpred.sam)/1e6
m.sam <- polymod.tab$w*c.sam

# Calculate normalized dominant right eigenvector
v.sam <- apply(X = m.sam, MARGIN = 2, FUN = function(m) {
  M <- make.m.matrix(m)
  eigen.M <- eigen(M)
  v <- with(eigen.M, Re(vectors[, 1]/sum(vectors[, 1])))
  return(v)
})

# Calulate ratio between females and males (=reference)
ratio.v.sam <- v.sam[(n.age+1):(2*n.age), ]/v.sam[1:n.age, ]

# Calculate statistics for v.sam and difference M-F
v.stat <- t(apply(v.sam, 1, quantile, prob = c(0.025, 0.50, 0.975)))
ratio.v.stat <- t(apply(ratio.v.sam, 1, quantile, prob = c(0.025, 0.5, 0.975)))

#
# Eigenvector plot
#

# Settings
ticks <- seq(0, 80, 10) # Tick locations
cols.line <- brewer.pal(n = 8, "Accent")[c(5, 6)] # Color for males and females: lines
cols.bord <- adjustcolor(cols.line, alpha = 0.5)

# Figure layout
layout(matrix(c(1, 2), 1, 2), widths = c(1, 1), heights = 1, respect = TRUE)
par(mgp = c(2.5, 1, 0), xaxs = "r", yaxs = "r", cex.axis = 0.9)

# Males and Female separate
par(mar = c(4, 4, 0.3, 0.3))
plot.new(); plot.window(xlim = c(0, 80), ylim = c(0, max(v.stat)))
polygon(c(age, rev(age)), c(v.stat[1:n.age, 1],             rev(v.stat[1:n.age, 3])),             col = cols.bord[1], border = NA)
polygon(c(age, rev(age)), c(v.stat[(n.age+1):(2*n.age), 1], rev(v.stat[(n.age+1):(2*n.age), 3])), col = cols.bord[2], border = NA)
matlines(age, cbind(v.stat[1:n.age, 2], v.stat[(n.age+1):(2*n.age), 2]), col = cols.line[1:2], lwd = 1, lty = 1, lend = "butt")
title(xlab = "Age", ylab = "Relative infection incidence")
axis(1, at = ticks, labels = ifelse(ticks%%20!=0 , NA, ticks), lwd = 0.5)
axis(2, lwd = 0.5)
box(lwd = 0.5)
legend("topright", legend = c("Males", "Females"), fill = cols.bord, border = cols.line, bty = "n", cex = 0.9)
legend("topleft", legend = "a", bty = "n", adj = 2)

# Difference Males - Females
par(mar = c(4, 4, 0.3, 0.3))
plot.new(); plot.window(xlim = c(0, 80), ylim = range(ratio.v.stat))
title(xlab = "Age", ylab = "Incidence ratio females vs. males")
polygon(c(age, rev(age)), c(ratio.v.stat[, 1], rev(ratio.v.stat[, 3])), col = grey(0.8), border = NA)
lines(age, ratio.v.stat[, 2], col = 1, lwd = 1, lty = 1, lend = "butt")
abline(h = 1, lty = 3)
axis(1, at = ticks, labels = ifelse(ticks%%20!=0 , NA, ticks), lwd = 0.5)
axis(2, lwd = 0.5)
box(lwd = 0.5)
legend("topleft", legend = "b", bty = "n", adj = 2)

#
# Plot for book chapter
#

# Settings
ticks <- seq(0, 80, 10) # Tick locations
cols.line <- c(brewer.pal(n = 9, "Set1")[2], "gold2") # Color for males and females: lines
cols.bord <- adjustcolor(cols.line, alpha = 0.5)

# Open connection to eps
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
cairo_ps(file = "results/fig_relative_infection_incidence.eps", width = 7, height = 7)

# Males and Female separate
par(mar = c(4, 4, 0.3, 0.3), mgp = c(2.5, 1, 0), xaxs = "r", yaxs = "r")
plot.new(); plot.window(xlim = c(0, 80), ylim = c(0, max(v.stat)))
polygon(c(age, rev(age)), c(v.stat[1:n.age, 1],               rev(v.stat[1:n.age, 3])),               col = cols.bord[1], border = NA)
polygon(c(age, rev(age)), c(v.stat[(n.age + 1):(2*n.age), 1], rev(v.stat[(n.age + 1):(2*n.age), 3])), col = cols.bord[2], border = NA)
matlines(age, cbind(v.stat[1:n.age, 2], v.stat[(n.age + 1):(2*n.age), 2]), col = cols.line[1:2], lwd = 1, lty = 1, lend = "butt")
title(xlab = "Age", ylab = "Relative infection incidence")
axis(1, at = ticks, labels = ifelse(ticks %% 20!=0 , NA, ticks), lwd = 0.5)
axis(2, lwd = 0.5)
box(lwd = 0.5)
legend("topright", legend = c("Men", "Women"), fill = cols.bord, border = cols.line, bty = "n")

# Close connection
dev.off()
