#
# Init
#

# Load packages
library(RColorBrewer)

# Read data
load(file = "results/polymod.tab.bin")

# Get ages and number of age classes
age <- unique(polymod.tab$part.age)
n.age <- length(age)

#
# Figure crude contact matrix
#

# Settings
cols <- colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(n = 100) # Colors
euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))                  # "euro" levels for contour lines
ticks <- seq(from = 0, to = 80, by = 10)                              # Tick locations

# Set variable and range
z <- with(polymod.tab, log(1 + y/U))
z.range <- range(z, na.rm = TRUE)

# Figure layout
layout(matrix(c(3, 4, 1, 2), 2, 2, byrow = TRUE),
  widths = c(1, 1),
  heights = c(1, 1),
  respect = TRUE)
par(mgp = c(0, 1, 0), xaxs = "i", yaxs = "i", cex.axis = 1)

# MM
par(mar = c(4, 4, 0.3, 0.3))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image(age, age, matrix(z[0*n.age^2 + 1:n.age^2], n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Men", side = 1, line = 2.5, adj = 0.5)
mtext("Men", side = 2, line = 2.5, adj = 0.5)
mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)

# FM
par(mar = c(4, 0.3, 0.3, 4))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image(age, age, matrix(z[1*n.age^2 + 1:n.age^2], n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
axis(side = 4, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Women", side = 1, line = 2.5, adj = 0.5)
mtext("Men",   side = 4, line = 2.5, adj = 0.5)

# MF
par(mar = c(0.3, 4, 4, 0.3))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image(age, age, matrix(z[2*n.age^2 + 1:n.age^2], n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
axis(side = 3, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Women", side = 2, line = 2.5, adj = 0.5)
mtext("Men",   side = 3, line = 2.5, adj = 0.5)

# FF
par(mar = c(0.3, 0.3, 4, 4))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image(age, age, matrix(z[3*n.age^2 + 1:n.age^2], n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
axis(side = 3, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
axis(side = 4, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Women", side = 3, line = 2.5, adj = 0.5)
mtext("Women", side = 4, line = 2.5, adj = 0.5)

#
# Figure smooth contact matrix
#

# Settings
cols <- colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(n = 100) # Colors
euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))                  # "euro" levels for contour lines
ticks <- seq(from = 0, to = 80, by = 10)                              # Tick locations

# Set variable and range
z <- log(polymod.tab$c)
z.range <- range(z)

# Figure layout
layout(matrix(c(3, 4, 1, 2), 2, 2, byrow = TRUE),
  widths = c(1, 1),
  heights = c(1, 1),
  respect = TRUE)
par(mgp = c(0, 1, 0), xaxs = "i", yaxs = "i", cex.axis = 1)

# MM
par(mar = c(4, 4, 0.3, 0.3))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image  (age, age, matrix(        z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
contour(age, age, matrix(1e6*exp(z[0*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Men", side = 1, line = 2.5, adj = 0.5)
mtext("Men", side = 2, line = 2.5, adj = 0.5)
mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)

# FM
par(mar = c(4, 0.3, 0.3, 4))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image  (age, age, matrix(        z[1*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
contour(age, age, matrix(1e6*exp(z[1*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
axis(side = 4, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Women", side = 1, line = 2.5, adj = 0.5)
mtext("Men",   side = 4, line = 2.5, adj = 0.5)

# MF
par(mar = c(0.3, 4, 4, 0.3))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image  (age, age, matrix(        z[2*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
contour(age, age, matrix(1e6*exp(z[2*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
axis(side = 3, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Women", side = 2, line = 2.5, adj = 0.5)
mtext("Men",   side = 3, line = 2.5, adj = 0.5)

# FF
par(mar = c(0.3, 0.3, 4, 4))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image  (age, age, matrix(        z[3*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
contour(age, age, matrix(1e6*exp(z[3*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(side = 3, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
axis(side = 4, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Women", side = 3, line = 2.5, adj = 0.5)
mtext("Women", side = 4, line = 2.5, adj = 0.5)
