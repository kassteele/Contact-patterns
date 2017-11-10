#
# Cross sections of c for specific ages
#

# Read data
load(file = "results/polymod.tab.bin")
load(file = "results/linpred.sam.bin")

# Here: male participant at age 10 and female participant at age 40
index.M10 <- with(polymod.tab, which(part.age == 10 & part.sex == "male"))
index.F40 <- with(polymod.tab, which(part.age == 40 & part.sex == "female"))

# Linear predictors with 95% CI
linpred.M10 <- t(apply(X = linpred.sam[index.M10, ], MARGIN = 1, FUN = function(x) c(mean(x), quantile(x, prob = c(0.025, 0.975)))))
linpred.F40 <- t(apply(X = linpred.sam[index.F40, ], MARGIN = 1, FUN = function(x) c(mean(x), quantile(x, prob = c(0.025, 0.975)))))

#
# Plot
#

# Settings
euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3))) # "euro" levels
ticks <- seq(from = 0, to = 80, by = 10)             # Tick locations

# Set range
y.range <- range(exp(c(linpred.M10, linpred.F40)))

# Figure layout
layout(matrix(c(3, 4, 1, 2), 2, 2, byrow = T), widths = c(1, 1), heights = c(1, 1), respect = T)
par(mgp = c(0, 1, 0), cex.axis = 1)

# M10 - M
par(mar = c(4, 4, 0.3, 0.3)); plot.new()
plot.window(xlim = range(age), ylim = y.range, log = "y")
grid(lty = 3, lwd = 0.5, col = grey(0.8))
polygon(c(age, rev(age)), exp(c(linpred.M10[1:n, 2], rev(linpred.M10[1:n, 3]))), col = grey(0.8), border = NA)
lines(age, exp(linpred.M10[1:n, 1]), lwd = 1, lend = "butt")
axis(1, at = ticks, labels = ifelse(ticks%%20!=0, NA, ticks), lwd = 0.5)
axis(2, at = euro.levs, labels = euro.levs, lwd = 0.5, las = 1)
mtext("10y male to males", side = 1, line = 2.5, adj = 0.5)
box(lwd = 0.5)
mtext("Age", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
mtext(expression(paste("Contact rate per ", 10^6)), side = 2, adj = 0.5, line = -1.5, outer = TRUE)

# F40 - M
par(mar = c(4, 0.3, 0.3, 4)); plot.new()
plot.window(xlim = range(age), ylim = y.range, log = "y")
grid(lty = 3, lwd = 0.5, col = grey(0.8))
polygon(c(age, rev(age)), exp(c(linpred.F40[1:n, 2], rev(linpred.F40[1:n, 3]))), col = grey(0.8), border = NA)
lines(age, exp(linpred.F40[1:n, 1]), lwd = 1, lend = "butt")
axis(1, at = ticks, labels = ifelse(ticks%%20!=0, NA, ticks), lwd = 0.5)
axis(4, at = euro.levs, labels = euro.levs, lwd = 0.5, las = 1)
mtext("40y female to males", side = 1, line = 2.5, adj = 0.5)
box(lwd = 0.5)

# M10 - F
par(mar = c(0.3, 4, 4, 0.3)); plot.new()
plot.window(xlim = range(age), ylim = y.range, log = "y")
grid(lty = 3, lwd = 0.5, col = grey(0.8))
polygon(c(age, rev(age)), exp(c(linpred.M10[(n+1):(2*n), 2], rev(linpred.M10[(n+1):(2*n), 3]))), col = grey(0.8), border = NA)
lines(age, exp(linpred.M10[(n+1):(2*n), 1]), lwd = 1, lend = "butt")
axis(2, at = euro.levs, labels = euro.levs, lwd = 0.5, las = 1)
axis(3, at = ticks, labels = ifelse(ticks%%20!=0, NA, ticks), lwd = 0.5)
mtext("10y male to females", side = 3, line = 2.5, adj = 0.5)
box(lwd = 0.5)

# F40 - F
par(mar = c(0.3, 0.3, 4, 4)); plot.new()
plot.window(xlim = range(age), ylim = y.range, log = "y")
grid(lty = 3, lwd = 0.5, col = grey(0.8))
polygon(c(age, rev(age)), exp(c(linpred.F40[(n+1):(2*n), 2], rev(linpred.F40[(n+1):(2*n), 3]))), col = grey(0.8), border = NA)
lines(age, exp(linpred.F40[(n+1):(2*n), 1]), lwd = 1, lend = "butt")
axis(3, at = ticks, labels = ifelse(ticks%%20!=0, NA, ticks), lwd = 0.5)
axis(4, at = euro.levs, labels = euro.levs, lwd = 0.5, las = 1)
mtext("40y female to females", side = 3, line = 2.5, adj = 0.5)
box(lwd = 0.5)
