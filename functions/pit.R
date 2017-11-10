# Function to get the non-randomized version of the PIT histogram
# Reference: Czado C, Gneiting T, Held L (2009). Predictive Model Assessment for Count Data. Biometrics, 65(4), 1254–1261.
# Original R code by Andrea Riebler: http://www.r-inla.org/?place=msg%2Fr-inla-discussion-group%2FeSwZJ5Iegcc%2FpJ0jX8BvUi4J
pit <- function(result, x, J = 20) {
  # result = inla object
  # x      = vector with count data
  # J      = number of bins
  # Get new improved estimates of the CPO/PIT values from inla object
  result <- inla.cpo(result = result)
  # Compute Px = P(X <= x) = PIT
  Px <- na.omit(result$cpo$pit)
  # Compute Pxm1 = P(X <= x-1) = PIT - CPO
  Pxm1 <- Px - na.omit(result$cpo$cpo)
  # Be sure to avoid negative PITs
  Pxm1 <- pmax(Pxm1, 0)
  # Omit missing values in data
  x <- na.omit(x)
  # Do the thing
  F_u.bar <- sapply(
    # For each u = real number 0 <= u <= 1
    X = (0:J)/J,
    # Replace randomized PIT value by its conditional CDF given the observed count x
    FUN = function(u, x, Px, Pxm1) {
      F_u <- ifelse(u <= Pxm1 , 0, pmin(1, (u - Pxm1)/(Px - Pxm1)))
      F_u[x == 0] <- pmin(1, u/Px)[x == 0] # Needless?
      if (u == 1) F_u <- 1
      if (u == 0) F_u <- 0
      return(mean(F_u))
    },
    # Arguments passed to FUN
    x = x, Px = Px, Pxm1 = Pxm1)
  f_j <- J*diff(F_u.bar)
  # Create PIT histogram
  hist.pit <- list(breaks = (0:J)/J, counts = f_j,
    #density = f_j,
    mids = (0:(J - 1))/J + diff((0:J)/J)/2, xname = "PIT", equidist = TRUE)
  class(hist.pit) <- "histogram"
  # Return output
  return(hist.pit)
}
