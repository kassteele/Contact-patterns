#
# Init
#

# Load packages
library(INLA)

# Source functions
source(file = "functions/construct.recordID.R")
source(file = "functions/construct.nodeID.R")
source(file = "functions/construct.Rmat.R")
source(file = "functions/construct.Dmat.R")

#
# Read data
#

# Read
load(file = "results/polymod.tab.bin")

#
# Model contact patterns with INLA
#

# Get ages and number of age classes
age <- unique(polymod.tab$part.age)
n.age <- length(age)

# Add node and record IDs to polymod.tab.list
polymod.tab <- within(polymod.tab, {
  record.id <- with(construct.recordID(n = n.age, sex = TRUE), c( rec.MM,  rec.FM,  rec.MF,  rec.FF))
  node.id   <- with(construct.nodeID  (n = n.age, sex = TRUE), c(node.MM, node.FM, node.MF, node.FF))
})

# Construct structure matrix R
ord <- 2
R <- construct.Rmat(n = n.age, order = ord, sex = TRUE)$R

# Run model
polymod.mod <- inla(
  y ~ 1 + f(node.id, model = "generic0", Cmatrix = R,
    # log-Gamma(1, 0.0001) prior on log-precision
    hyper = list(prec = list(prior = "loggamma", param = c(1, 0.0001))),
    rankdef = 3*ord^2, constr = TRUE, diagonal = 0.001),
  E = U,
  family = "nbinomial",
  data = polymod.tab,
  # Normal(0, 0.001) prior on intercept
  control.fixed = list(mean.intercept = 0, prec.intercept = 0.001),
  # Normal(0, 0.001) prior on log-dispersion parameter
  control.family = list(
    hyper = list(theta = list(prior = "gaussian", param = c(0, 0.001)))),
  # Integration strategy eb for faster computation
  control.inla = list(int.strategy = "eb"),
  # Compute linear predictor
  control.predictor = list(compute = TRUE, link = 1),
  control.compute = list(
    waic = TRUE,    # Compute WAIC
    cpo = TRUE,     # Compute cross-validated predictive measures
    config = TRUE)) # Enable posterior sampling

# Show summary
summary(polymod.mod)

#
# Post processing
#

# Add expected c and m to polymod.tab
polymod.tab <- within(polymod.tab, {
  # Compute expected linear predictor
  linpred <- polymod.mod$summary.linear.predictor[, "mean"]
  # c = exp(linear predictor) = contact rate
  # Divide c by 1e6 to go back to original scale (deflates contact rate c)
  c <- exp(linpred)/1e6
  # m = w*c = contact intensity
  m <- w*c
  # Remove linpred
  rm(linpred)
})

#
# Save results
#

save(polymod.mod, file = "results/polymod.mod.bin")
save(polymod.tab, file = "results/polymod.tab.bin")

#
# Export to contact_matrix_data.txt
#

tmp <- subset(polymod.tab, select = c(record.id, node.id, part.age, cont.age, part.sex, cont.sex, w, m, c))
tmp <- within(tmp, {
  c <- round(c*1e6, digits = 5)
  m <- round(m, digits = 5)
})
rownames(tmp) <- 1:nrow(tmp)
write.table(tmp, file = "results/contact_matrix_data.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, eol = "\n")

#
# Optional: generate samples from approximated joint posterior distribution
#

# Generate samples from approximated joint posterior distribution
n.samples <- 1000
inla.sam <- inla.posterior.sample(n = n.samples, result = polymod.mod)
# Extract samples for linear predictor and put them in a matrix
# Size: nrow(polymod.tab) x n.samples
linpred.sam <- sapply(X = inla.sam, FUN = function(x) x$latent[grepl(pattern = "Predictor", x = rownames(x$latent))])

#
# Save results
#

save(linpred.sam, file = "results/linpred.sam.bin")
