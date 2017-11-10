#
# Standardized measure of global assortativeness - Farrington et al (2009)
#

# Read data
load(file = "results/polymod.tab.bin")
load(file = "results/pop.data.bin")
load(file = "results/linpred.sam.bin")

# Get ages and number of age classes
age <- unique(polymod.tab$part.age)
n.age <- length(age)

# Get population numbers
wM <- pop.data$w[1:n.age]
wF <- pop.data$w[(n.age + 1):(2*n.age)]
w <- wM + wF

# Samples c.sam
c.sam <- exp(linpred.sam)/1e6

# Get population density as in Farrington 2009 paper
f <- w/sum(w) # ...adds up to 1

# Set s and t as in Farrington 2009 paper
s <- age # Age as integer 0:80
t <- age # Idem

# Variance of age distribution as in Farrington 2009 paper
sigma2 <- sum(f*s^2) - sum(s*f)^2

# Set n.samples
n.samples <- ncol(c.sam)

# Allocate I2s vector
I2s <- rep(0, n.samples)

# Loop over samples
for (sample.i in 1:n.samples) {
  # Print progress
  print(sample.i)

  # Sex specific c.sam
  cMM <- c.sam[0*n.age^2 + 1:n.age^2, sample.i]
  cFM <- c.sam[1*n.age^2 + 1:n.age^2, sample.i]
  cMF <- c.sam[2*n.age^2 + 1:n.age^2, sample.i]
  cFF <- c.sam[3*n.age^2 + 1:n.age^2, sample.i]
  # Recalculate m
  mMM <- wM*cMM
  mFM <- wM*cFM
  mMF <- wF*cMF
  mFF <- wF*cFF
  # Sex weighted average of m
  m <- (wM*(mMM + mMF) + wF*(mFM + mFF))/w
  # And back to c
  c <- m/w

  # Beta is matrix version of c as in Farrington et al (2009) paper
  beta <- matrix(c, nrow = n.age, ncol = n.age)

  # Calculate fc matrix = density of contact pairs
  fc <- matrix(0, nrow = n.age, ncol = n.age)
  for (i in 1:n.age) {
    for (j in 1:n.age) {
      fc[i, j] <- f[i]*beta[i, j]*f[j]
    }
  }
  fc <- fc/sum(fc)

  # Absolute disassortativeness
  I2 <- 0.5*sum(outer(s, t, FUN = "-")^2*fc)

  # Standardized measure of global assortativeness
  I2s[sample.i] <- I2/sigma2
}

# Print results
quantile(I2s, prob = c(0.5, 0.025, 0.975))
