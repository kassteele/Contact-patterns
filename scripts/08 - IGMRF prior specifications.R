#
# Init
#

# Load packages
library(Matrix)
library(RColorBrewer)

#
# Simulate realisations
#

# Simulate data
n <- 50; n.sim <- 1000; age <- 1:n
D <- construct.Rmat(n = n)$D
R <- construct.Rmat(n = n)$R

# Martin Boer trick
Z <- t(D)
Q <- D%*%t(D)%*%D%*%t(D)
Q <- Q+Diagonal(n = nrow(Q), 0.0001)
L <- chol(Q)
b <- solve(L, matrix(rnorm(nrow(L)*n.sim), nrow(L), n.sim))
y <- Z%*%b

# Compute marginal variance
z <- y[,1]

# Reorder into n^2*2*2 vector
node.id <- with(construct.nodeID(n, sex = TRUE), c(node.MM, node.FM, node.MF, node.FF))
z <- z[node.id]

#
# Smooth contact matrices
#

# Settings
cols <- colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(100) # Colors
euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))              # "euro" levels for contour lines
ticks <- seq(0, 80, 10)                                           # Tick locations

# Open file
#pdf(file = "contact patterns.pdf", width = 17.35/2.54, height = 17.35/2.54, pointsize = 12*(17.35/2.54)/7, onefile = TRUE)

# Set variable and range
z.range <- range(z)

# Figure layout
layout(matrix(c(3, 4, 1, 2), 2, 2, byrow = T), widths = c(1, 1), heights = c(1, 1), respect = T)
par(mgp = c(0, 1, 0), xaxs = "i", yaxs = "i", cex.axis = 1)

# MM
par(mar = c(4, 4, 0.3, 0.3)); plot.new(); plot.window(xlim = range(age), ylim = range(age))
image  (age, age, matrix(    z[0*n^2+1:n^2] , n, n), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
contour(age, age, matrix( exp(z[0*n^2+1:n^2]), n, n), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(1, at = ticks, labels = ifelse(ticks%%20!=0, NA, ticks), lwd = 0.5)
axis(2, at = ticks, labels = ifelse(ticks%%20!=0, NA, ticks), lwd = 0.5)
mtext("Male", side = 1, line = 2.5, adj = 0.5)
mtext("Male", side = 2, line = 2.5, adj = 0.5)
box(lwd = 0.5)
mtext("Participants", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
mtext("Contacts", side = 2, adj = 0.5, line = -1.5, outer = TRUE)
#mtext(names(polymod.tab.split)[i], side = 3, adj = 0.5, line = -1.5, outer = TRUE, font = 2)

# FM
par(mar = c(4, 0.3, 0.3, 4)); plot.new(); plot.window(xlim = range(age), ylim = range(age))
image  (age, age, matrix(    z[1*n^2+1:n^2],  n, n), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
contour(age, age, matrix(1e6*exp(z[1*n^2+1:n^2]), n, n), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(1, at = ticks, labels = ifelse(ticks%%20!=0 | ticks==0, NA, ticks), lwd = 0.5)
axis(4, at = ticks, labels = ifelse(ticks%%20!=0, NA, ticks), lwd = 0.5)
mtext("Female", side = 1, line = 2.5, adj = 0.5)
mtext("Male", side = 4, line = 2.5, adj = 0.5)
box(lwd = 0.5)

# MF
par(mar = c(0.3, 4, 4, 0.3)); plot.new(); plot.window(xlim = range(age), ylim = range(age))
image  (age, age, matrix(    z[2*n^2+1:n^2] , n, n), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
contour(age, age, matrix(1e6*exp(z[2*n^2+1:n^2]), n, n), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(2, at = ticks, labels = ifelse(ticks%%20!=0 | ticks==0, NA, ticks), lwd = 0.5)
axis(3, at = ticks, labels = ifelse(ticks%%20!=0, NA, ticks), lwd = 0.5)
mtext("Female", side = 2, line = 2.5, adj = 0.5)
mtext("Male", side = 3, line = 2.5, adj = 0.5)
box(lwd = 0.5)

# FF
par(mar = c(0.3, 0.3, 4, 4)); plot.new(); plot.window(xlim = range(age), ylim = range(age))
image  (age, age, matrix(    z[3*n^2+1:n^2] , n, n), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
contour(age, age, matrix(1e6*exp(z[3*n^2+1:n^2]), n, n), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
axis(3, at = ticks, labels = ifelse(ticks%%20!=0 | ticks==0, NA, ticks), lwd = 0.5)
axis(4, at = ticks, labels = ifelse(ticks%%20!=0 | ticks==0, NA, ticks), lwd = 0.5)
mtext("Female", side = 3, line = 2.5, adj = 0.5)
mtext("Female", side = 4, line = 2.5, adj = 0.5)
box(lwd = 0.5)



matplot(y, type = "l")
