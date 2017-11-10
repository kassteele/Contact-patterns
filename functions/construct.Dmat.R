# Function to construct D-matrix
construct.Dmat <- function(n, order, tri = FALSE) {
  require(Matrix)
  In <- Diagonal(n)            # n x n matrix
  D0 <- diff(In, diff = order) # Differences for single vector
  D1 <- kronecker(In, D0)      # Differences in horizontal direction
  D2 <- kronecker(D0, In)      # Differences in   vertical direction
  if (tri) {
    # For lower triangular matrix (MM and FF)
    ri1.mat <- matrix(1:((n - order)*n), nrow = n - order, ncol = n)
    ri2.mat <- matrix(1:(n*(n - order)), nrow = n, ncol = n - order)
    ci.mat  <- matrix(1:n^2, nrow = n, ncol = n)
    ri1 <- ri1.mat[row(ri1.mat) >=  col(ri1.mat)         ] # Rows to keep in D1
    ri2 <- ri2.mat[row(ri2.mat) >= (col(ri2.mat) + order)] # Rows to keep in D2
    ci  <-  ci.mat[row(ci.mat)  >=  col(ci.mat)          ] # Columns to keep in D1 & D2
    D <- rBind(D1[ri1, ci], D2[ri2, ci])
    return(list(D0 = D0, D1 = D1, D2 = D2, ri1 = ri1, ri2 = ri2, ci = ci, D = D))
  } else {
    # For full matrix (FM)
    D <- rBind(D1, D2)
    return(list(D0 = D0, D1 = D1, D2 = D2, D = D))
  }
}
