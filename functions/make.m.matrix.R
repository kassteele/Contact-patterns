# Function to construct m matrix from output vector
make.m.matrix <- function(output.vector, sex = TRUE) {
  if (sex) {
    # Sex specific
    n <- sqrt(length(output.vector)/4)
    mMM <- matrix(output.vector[        1:n^2], nrow = n, ncol = n)
    mFM <- matrix(output.vector[  n^2 + 1:n^2], nrow = n, ncol = n)
    mMF <- matrix(output.vector[2*n^2 + 1:n^2], nrow = n, ncol = n)
    mFF <- matrix(output.vector[3*n^2 + 1:n^2], nrow = n, ncol = n)
    return((cbind(rbind(mMM, mFM), rbind(mMF, mFF))))
  } else {
    # Not sex specific
    n <- sqrt(length(output.vector))
    m <- matrix(output.vector[1:n^2], nrow = n, ncol = n)
    return(m)
  }
}
