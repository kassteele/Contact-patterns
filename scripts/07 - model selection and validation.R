#
# Model choice and validation
#

# Source PIT function
source(file = "scripts/functions/pit.R")

# Load data
load(file = "results/polymod.tab.bin")

# First models: with different priors, e.g. RW1, RW2 and RW3 prior and save results
# Load model
load(file = "results/polymod.mod.bin")

# Get WAIC and p.eff
round(polymod.mod$waic$waic, digits = 1)
round(polymod.mod$waic$p.eff, digits = 1)

# Get PIT values
pit.object <- pit(polymod.mod, polymod.tab$y)

# Plot
par(mar = c(4, 4, 1.5, 0.2), mgp = c(2.5, 1, 0), lwd = 0.5)
plot(pit.object, col = grey(0.8), ylab = "Relative frequency", main = "RW2 prior", lwd = 0.5)
abline(h = 1, lty = 2)
