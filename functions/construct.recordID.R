# Function to construct record IDs in matrix form
construct.recordID <- function(n, sex = TRUE) {
  if (sex) {
    # Sex specific
    rec.MM <- matrix(        1:n^2, nrow = n, ncol = n)
    rec.FM <- matrix(  n^2 + 1:n^2, nrow = n, ncol = n)
    rec.MF <- matrix(2*n^2 + 1:n^2, nrow = n, ncol = n)
    rec.FF <- matrix(3*n^2 + 1:n^2, nrow = n, ncol = n)
    return(list(rec.MM = rec.MM, rec.FM = rec.FM, rec.MF = rec.MF, rec.FF = rec.FF))
  } else {
    # Not sex specific
    rec <- matrix(1:n^2, nrow = n, ncol = n)
    return(list(rec = rec))
  }
}
