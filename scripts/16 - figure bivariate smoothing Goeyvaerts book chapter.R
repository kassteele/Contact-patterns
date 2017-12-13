#
# Init
#

# Load packages
library(mgcv)
library(RColorBrewer)

#
# Read data
#

# Read
load(file = "results/polymod.data.bin")
load(file = "results/pop.data.bin")

#
# Cross tabulation
#

# Cross tabulate number of participants for each combination of age
# t = number of participants
polymod.part.tab <- with(subset(polymod.data, subset = !duplicated(local.id)),
  as.data.frame(table(part.age = part.age),
    responseName = "t"))

# Cross tabulate number of contacts for each combination of age
# y = number of participants x number of contacts
polymod.cont.tab <- with(polymod.data,
  as.data.frame(table(part.age = part.age, cont.age = cont.age),
    responseName = "y"))

# Make age an integer again (was factor in previous section)
polymod.part.tab <- within(polymod.part.tab, {
  part.age <- as.numeric(as.character(part.age))
})
polymod.cont.tab <- within(polymod.cont.tab, {
  part.age <- as.numeric(as.character(part.age))
  cont.age <- as.numeric(as.character(cont.age))
})

# Sum population numbers by age
pop.data <- aggregate(w ~ cont.age, FUN = sum, data = pop.data)

# Merge polymod.cont.tab, polymod.part.tab and pop.data
# polymod.tab is the tabulated version of polymod.data
polymod.tab <- merge(polymod.part.tab, merge(polymod.cont.tab, pop.data, all = TRUE), all = TRUE)

# Reorder polymod.tab
polymod.tab <- polymod.tab[c("part.age", "cont.age", "y", "t", "w")]
polymod.tab <- with(polymod.tab, polymod.tab[order(cont.age, part.age), ])

# Calculate denominator U = number of participants t
# If t = 0, then set U = 1 and y = NA (record will not contribute to likelihood)
polymod.tab <- within(polymod.tab, {
  U <- t
  y <- ifelse(U == 0, yes = NA, no = y)
  U <- ifelse(U == 0, yes = 1, no = U)
})

#
# Model contact patterns with mgcv
# Goeyvaerts et al. (2009): bivariate smoothing
#

# Get ages and number of age classes
age <- unique(polymod.tab$part.age)
n.age <- length(age)

# Run model
polymod.mod <- bam(y ~ te(part.age, cont.age, bs = "tp", k = 11, fx = F),
  offset = log(U),
  family = nb,
  data = polymod.tab)

# Show summary
summary(polymod.mod)

#
# Post processing
#

# Calculate contact intensties
m <- exp(predict(polymod.mod,
  newdata = cbind(expand.grid(part.age = age, cont.age = age), data.frame(U = 1)),
  type = "link"))
m <- matrix(m, nrow = n.age, ncol = n.age)
w <- matrix(pop.data$w, nrow = n.age, ncol = n.age, byrow = TRUE)

# Take the reciprocal nature of the contacts into account
m <- ((m*t(w) + t(m)*w)/2)/t(w)
c <- m/w

#
# Figure smooth contact matrix
#

# Settings
cols <- colorRampPalette(brewer.pal(name = "YlGn", n = 9))(n = 100) # Colors
euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))                # "euro" levels for contour lines
ticks <- seq(from = 0, to = 80, by = 10)                            # Tick locations

# Set variable and range
z <- log(c)
z.range <- range(z)

# Open connection to eps
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript(file = "results/fig_contact_rate_bivariate_smooth.eps", width = 7*1.14, height = 7)

# Figure layout
layout(matrix(c(1, 2), 1, 2, byrow = TRUE),
  widths = c(1, 0.14),
  heights = 1,
  respect = TRUE)
par(mgp = c(2.5, 1, 0), xaxs = "i", yaxs = "i", cex.axis = 1)

# Plot
par(mar = c(4, 4, 0.3, 0.3))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image(age, age, matrix(z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
title(xlab = "Participants", ylab = "Contacts")

# Legend
z1.range <- 1e6*exp(z.range)
labs <- euro.levs[euro.levs >= z1.range[1] & euro.levs <= z1.range[2]]
par(mar = c(4, 1, 0.3, 3), xpd = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = log(z1.range))
rasterImage(image = cols, xleft = 0, ybottom = log(z1.range[2]), xright = 1, ytop = log(z1.range[1]), interpolate = FALSE)
axis(side = 4, at = log(labs), labels = prettyNum(labs), lwd = 0, lwd.ticks = 0.5, las = 1)

# Close connection
dev.off()
