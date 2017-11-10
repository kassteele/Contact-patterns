#
# Aggregation of contact intensities and rates
#

#
# Init
#

# Load packages
library(RColorBrewer)

#
# Read data
#

load(file = "results/polymod.tab.bin")
load(file = "results/pop.data.bin")

#
# Settings
#

# Get ages and number of age classes
age <- unique(polymod.tab$part.age)
n.age <- length(age)

# Define age cateogries
age.cat <- cut(age, breaks = c(0, 1, 5, 10, 20, 45, 65, 80), include.lowest = TRUE, right = FALSE)

#
# Data wrangling
#

# Reshape pop.data in wide format. Sex is put in separate columns
pop.data.wide <- cbind(
  data.frame(age = age),
  matrix(pop.data$w, nrow = n.age, ncol = 2, dimnames = list(NULL, c("wM", "wF"))))
pop.data.wide <- within(pop.data.wide, w <- wM + wF)

# Create a dataframe contact.data with n.age^2 rows. Sex is put in separate columns
contact.data <- cbind(
  expand.grid(part.age = age, cont.age = age),
  matrix(polymod.tab$c, nrow = n.age^2, ncol = 4, dimnames = list(NULL, c("cMM", "cFM", "cMF", "cFF"))),
  matrix(polymod.tab$m, nrow = n.age^2, ncol = 4, dimnames = list(NULL, c("mMM", "mFM", "mMF", "mFF"))))

#
# Aggregate over sexes (equation 6.1, suppl mat AOAS paper)
#

# Aggregate contact intensities over sexes
record.id.part <- match(x = contact.data$part.age, table = pop.data.wide$age)
contact.data <- within(contact.data, {
  m <- (pop.data.wide[record.id.part, "wM"]*(mMM + mMF) + pop.data.wide[record.id.part, "wF"]*(mFM + mFF))/pop.data.wide[record.id.part, "w"]
})

# Calculate contact rates
record.id.cont <- match(x = contact.data$cont.age, table = pop.data.wide$age)
contact.data <- within(contact.data, {
  wM <- pop.data.wide[record.id.cont, "wM"]
  wF <- pop.data.wide[record.id.cont, "wF"]
  w  <- pop.data.wide[record.id.cont, "w"]
  c <- m/w
})

# Reorder columns
contact.data <- contact.data[, c("part.age", "cont.age", "wM", "wF", "w", "mMM", "mFM", "mMF", "mFF", "m", "cMM", "cFM", "cMF", "cFF", "c")]

#
# Aggregate over ages (equation 6.2, suppl mat AOAS paper)
#

# Add age categories to contact.data and pop.data
contact.data <- cbind(contact.data, expand.grid(part.age.cat = age.cat, cont.age.cat = age.cat))
pop.data.wide <- cbind(pop.data.wide, age.cat = age.cat)

# Aggregate population numbers
pop.data.agg <- aggregate(cbind(wM, wF, w) ~ age.cat, FUN = sum, data = pop.data.wide)

# Aggegrate contact intensities over ages
record.id.part <- match(x = contact.data$part.age, table = pop.data.wide$age)
record.id.part.agg <- match(x = contact.data$part.age.cat, table = pop.data.agg$age.cat)
contact.data.agg <- within(contact.data, {
  m   <- pop.data.wide[record.id.part, "w" ]*m  /pop.data.agg[record.id.part.agg, "w" ]
  mMM <- pop.data.wide[record.id.part, "wM"]*mMM/pop.data.agg[record.id.part.agg, "wM"]
  mFM <- pop.data.wide[record.id.part, "wF"]*mFM/pop.data.agg[record.id.part.agg, "wF"]
  mMF <- pop.data.wide[record.id.part, "wM"]*mMF/pop.data.agg[record.id.part.agg, "wM"]
  mFF <- pop.data.wide[record.id.part, "wF"]*mFF/pop.data.agg[record.id.part.agg, "wF"]
})
contact.data.agg <- aggregate(cbind(mMM, mFM, mMF, mFF, m) ~ part.age.cat + cont.age.cat, FUN = sum, data = contact.data.agg)

# Calculate contact rates
record.id.cont.agg <- match(x = contact.data.agg$cont.age.cat, table = pop.data.agg$age.cat)
contact.data.agg <- within(contact.data.agg, {
  wM <- pop.data.agg[record.id.cont.agg, "wM"]
  wF <- pop.data.agg[record.id.cont.agg, "wF"]
  w  <- pop.data.agg[record.id.cont.agg, "w" ]
  cMM <- mMM/wM
  cFM <- mFM/wM
  cMF <- mMF/wF
  cFF <- mFF/wF
  c   <- m  /w
})

# Reorder columns
contact.data.agg <- contact.data.agg[, c("part.age.cat", "cont.age.cat", "wM", "wF", "w", "mMM", "mFM", "mMF", "mFF", "m", "cMM", "cFM", "cMF", "cFF", "c")]

#
# Figure aggregated contact matrix
#

# Little tric to expand aggregated contact rates
n.age.age <- table(age.cat) # Number of ages per age category
n.age.cat <- length(n.age.age)  # Number of age categories
z <- NULL; for (i in 1:n.age.cat) z <- c(z, rep(rep(log(contact.data.agg$c)[(i-1)*n.age.cat + 1:n.age.cat], times = n.age.age), times = n.age.age[i]))

# Settings
cols <- colorRampPalette(colors = brewer.pal(name = "YlOrRd", n = 9))(n = 100) # Colors
ticks <- seq(from = 0, to = 80, by = 10) # Tick locations
z.range <- range(z)

# Figure
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), xaxs = "i", yaxs = "i", cex.axis = 1)
plot.new()
plot.window(xlim = c(min(age), max(age)), ylim = c(min(age), max(age)))
image(age, age, matrix(z[0*n.age^2 + 1:n.age^2], nrow = n.age, ncol = n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
axis(1, at = ticks, labels = ifelse(ticks%%20 != 0, NA, ticks), lwd = 0.5)
axis(2, at = ticks, labels = ifelse(ticks%%20 != 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
title(xlab = "Participants", ylab = "Contacts")

#
# Calculate normalized dominant eigenvector
#

# Eigenvalues and -vectors
eigen.m <- with(contact.data.agg, eigen(matrix(m, nrow = n.age.cat, ncol = n.age.cat)))

# Print
eigen.m

# Normalized dominant eigenvector
with(eigen.m, vectors[, 1]/sum(vectors[, 1]))
