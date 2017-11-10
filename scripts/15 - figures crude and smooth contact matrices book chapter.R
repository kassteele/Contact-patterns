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
cols <- colorRampPalette(brewer.pal(name = "YlGn", n = 9))(n = 100) # Colors
euro.levs <- c(0, 1, 2, 5, 10, 20, 50)
ticks <- seq(from = 0, to = 80, by = 10)                              # Tick locations

# Set variable and range
z <- with(polymod.tab, (y/U)^(1/3))
z.range <- range(z, na.rm = TRUE)

# # Open connection to eps
# setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
# postscript(file = "fig_contact_rate_crude.eps", width = 7.7, height = 7)

# Figure layout
layout(matrix(c(3, 4, 5, 1, 2, 5), 2, 3, byrow = TRUE),
  widths = c(1, 1, 0.2),
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
mtext("Participants", side = 1, adj = 1/sum(c(1, 1, 0.2)), line = -1.5, outer = TRUE)
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

# Legend
z1.range <- z.range^3
labs <- euro.levs[euro.levs >= z1.range[1] & euro.levs <= z1.range[2]]
par(mar = c(4, 1, 4, 3), xpd = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = z1.range^(1/3))
rasterImage(image = cols, xleft = 0, ybottom = z1.range[2]^(1/3), xright = 1, ytop = z1.range[1]^(1/3), interpolate = FALSE)
axis(side = 4, at = labs^(1/3), labels = prettyNum(labs), lwd = 0, lwd.ticks = 0.5, las = 1)

# # Close connection
# dev.off()

#
# Figure smooth contact matrix
#

# Settings
cols <- colorRampPalette(brewer.pal(name = "YlGn", n = 9))(n = 100) # Colors
euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))                  # "euro" levels for contour lines
ticks <- seq(from = 0, to = 80, by = 10)                              # Tick locations

# Set variable and range
z <- log(polymod.tab$c)
z.range <- range(z)

# # Open connection to eps
# setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
# postscript(file = "fig_contact_rate_smooth.eps", width = 7.7, height = 7)

# Figure layout
layout(matrix(c(3, 4, 5, 1, 2, 5), 2, 3, byrow = TRUE),
  widths = c(1, 1, 0.2),
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
#contour(age, age, matrix(1e6*exp(z[0*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Men", side = 1, line = 2.5, adj = 0.5)
mtext("Men", side = 2, line = 2.5, adj = 0.5)
mtext("Participants", side = 1, adj = 1/sum(c(1, 1, 0.2)), line = -1.5, outer = TRUE)
mtext("Contacts",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)

# FM
par(mar = c(4, 0.3, 0.3, 4))
plot.new()
plot.window(
  xlim = c(min(age) - 0.5, max(age) + 0.5),
  ylim = c(min(age) - 0.5, max(age) + 0.5))
image  (age, age, matrix(        z[1*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
#contour(age, age, matrix(1e6*exp(z[1*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
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
#contour(age, age, matrix(1e6*exp(z[2*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
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
#contour(age, age, matrix(1e6*exp(z[3*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(side = 3, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
axis(side = 4, at = ticks, labels = ifelse(ticks %% 20 != 0 | ticks == 0, NA, ticks), lwd = 0.5)
box(lwd = 0.5)
mtext("Women", side = 3, line = 2.5, adj = 0.5)
mtext("Women", side = 4, line = 2.5, adj = 0.5)

# Legend
z1.range <- 1e6*exp(z.range)
labs <- euro.levs[euro.levs >= z1.range[1] & euro.levs <= z1.range[2]]
par(mar = c(4, 1, 4, 3), xpd = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = log(z1.range))
rasterImage(image = cols, xleft = 0, ybottom = log(z1.range[2]), xright = 1, ytop = log(z1.range[1]), interpolate = FALSE)
axis(side = 4, at = log(labs), labels = prettyNum(labs), lwd = 0, lwd.ticks = 0.5, las = 1)

# # Close connection
# dev.off()
