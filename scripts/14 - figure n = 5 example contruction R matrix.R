#
# Figure n = 5 example contruction R matrix
#

#
# Init
#

# Load packages
library(RColorBrewer)

# Source functions
source(file = "functions/construct.recordID.R")
source(file = "functions/construct.nodeID.R")

# Circle function
circle <- function(x = 0, y = 0, rx = 1, ry = rx, col = 1, border = NULL, lwd = 0.5) {
  t <- seq(0, 360, 15)*pi/180
  polygon(x+rx*cos(t), y+ry*sin(t), col = col, border = border, lwd = lwd)
}

#
# Create figure
#

# Node and record IDs
n <- 5; n2 <- 500
record.id  <- construct.recordID(n)
node.id    <- construct.nodeID(n)
record.id2 <- construct.recordID(n2)
node.id2   <- construct.nodeID(n2)

# Set colors
pal <- colorRampPalette(brewer.pal(12, "Set3")[c(5, 1, 7, 12, 6, 4, 10, 3)])
cols1 <- pal(max(node.id$node.FF))
cols2 <- pal(max(node.id2$node.FF))
cols3 <- rgb(255 - 0.7*(255 - t(col2rgb(cols2))), max = 255)

# Figure layout
layout(matrix(c(3, 4, 1, 2), 2, 2, byrow = T), widths = c(1, 1), heights = c(1, 1), respect = TRUE)
par(mgp = c(0, 1, 0), xaxs = "i", yaxs = "i", cex.axis = 1)

# MM
par(mar = c(4, 4, 0.3, 0.3)); plot.new(); plot.window(xlim = c(0.5, n+0.5), ylim = c(0.5, n+0.5))
with(node.id2, image(seq(0.5, n+0.5, l = n2), seq(0.5, n+0.5, l = n2), node.MM, zlim = c(1, max(node.FF)),
  col = cols3, useRaster = TRUE, add = TRUE))
segments(3:n, 1, 3:n, 3:n, lwd = 0.5); segments(1:n, 1:(n-2), n, 1:(n-2), lwd = 0.5) # MM
for (i in 1:n) {for (j in 1:n) {
  circle(i, j, rx = 0.35, ry = 0.35*0.66, col = cols1[node.id$node.MM[i, j]], border = 1, lwd = 0.5)
  text(i, j, node.id$node.MM[i, j], cex = 1)
  text(i+0.32, j+0.35, record.id$rec.MM[i, j], cex = 0.8)
}}
axis(1, at = 1:5, lwd = 0.5)
axis(2, at = 1:5, lwd = 0.5)
mtext("Male", side = 1, line = 2.5, adj = 0.5)
mtext("Male", side = 2, line = 2.5, adj = 0.5)
box(lwd = 0.5)
mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
mtext("Contacts", side = 2, adj = 0.5, line = -1.5, outer = TRUE)

# FM
par(mar = c(4, 0.3, 0.3, 4)); plot.new(); plot.window(xlim = c(0.5, n+0.5), ylim = c(0.5, n+0.5))
with(node.id2, image(seq(0.5, n+0.5, l = n2), seq(0.5, n+0.5, l = n2), node.FM, zlim = c(1, max(node.FF)),
  col = cols3, useRaster = TRUE, add = TRUE))
segments(1:n, 1, 1:n, n, lwd = 0.5); segments(1, 1:n, n, 1:n, lwd = 0.5) # MF
for (i in 1:n) {for (j in 1:n) {
  circle(i, j, rx = 0.35, ry = 0.35*0.66, col = cols1[node.id$node.FM[i, j]], border = 1, lwd = 0.5)
  text(i, j, node.id$node.FM[i, j], cex = 1)
  text(i+0.32, j+0.35, record.id$rec.FM[i, j], cex = 0.8)
}}
axis(1, at = 1:5, lwd = 0.5)
axis(4, at = 1:5, lwd = 0.5)
mtext("Female", side = 1, line = 2.5, adj = 0.5)
mtext("Male", side = 4, line = 2.5, adj = 0.5)
box(lwd = 0.5)

# MF
par(mar = c(0.3, 4, 4, 0.3)); plot.new(); plot.window(xlim = c(0.5, n+0.5), ylim = c(0.5, n+0.5))
with(node.id2, image(seq(0.5, n+0.5, l = n2), seq(0.5, n+0.5, l = n2), node.MF, zlim = c(1, max(node.FF)),
  col = cols3, useRaster = TRUE, add = TRUE))
for (i in 1:n) {for (j in 1:n) {
  circle(i, j, rx = 0.35, ry = 0.35*0.66, col = cols1[node.id$node.MF[i, j]], border = 1, lwd = 0.5)
  text(i, j, node.id$node.MF[i, j], cex = 1)
  text(i+0.32, j+0.35, record.id$rec.MF[i, j], cex = 0.8)
}}
axis(2, at = 1:5, lwd = 0.5)
axis(3, at = 1:5, lwd = 0.5)
mtext("Female", side = 2, line = 2.5, adj = 0.5)
mtext("Male", side = 3, line = 2.5, adj = 0.5)
box(lwd = 0.5)

# FF
par(mar = c(0.3, 0.3, 4, 4)); plot.new(); plot.window(xlim = c(0.5, n+0.5), ylim = c(0.5, n+0.5))
with(node.id2, image(seq(0.5, n+0.5, l = n2), seq(0.5, n+0.5, l = n2), node.FF, zlim = c(1, max(node.FF)),
  col = cols3, useRaster = TRUE, add = TRUE))
segments(3:n, 1, 3:n, 3:n, lwd = 0.5); segments(1:n, 1:(n-2), n, 1:(n-2), lwd = 0.5) # FF
for (i in 1:n) {for (j in 1:n) {
  circle(i, j, rx = 0.35, ry = 0.35*0.66, col = cols1[node.id$node.FF[i, j]], border = 1, lwd = 0.5)
  text(i, j, node.id$node.FF[i, j], cex = 1)
  text(i+0.32, j+0.35, record.id$rec.FF[i, j], cex = 0.8)
}}
axis(3, at = 1:5, lwd = 0.5)
axis(4, at = 1:5, lwd = 0.5)
mtext("Female", side = 3, line = 2.5, adj = 0.5)
mtext("Female", side = 4, line = 2.5, adj = 0.5)
box(lwd = 0.5)
