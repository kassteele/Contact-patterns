# Function to construct node IDs in matrix form
construct.nodeID <- function(n, sex = TRUE) {
  if (sex) {
    # Sex specific
    nn1 <- n*(n + 1)/2 # Number of nodes in MM and FF
    nn2 <- n^2         # Number of nodes in FM and MF
    node.MM <- node.FF <- matrix(0, nrow = n, ncol = n) # Empty matrix
    LT <- lower.tri(node.MM, diag = TRUE) # Lower tri of MM and FF with diag
    node.MM[ LT] <- 1:nn1           # Fill lower tri with node IDs
    node.MM[!LT] <- t(node.MM)[!LT] # Fill upper tri with transposed
    node.FM      <- matrix(nn1 + 1:nn2, nrow = n, ncol = n) # Fill entire matrix
    node.MF      <- t(node.FM)                              # Transposed
    node.FF[ LT] <- nn1 + nn2 + 1:nn1 # As node.MM
    node.FF[!LT] <- t(node.FF)[!LT]
    return(list(node.MM = node.MM, node.FM = node.FM, node.MF = node.MF, node.FF = node.FF))
  } else {
    # Not sex specific
    nn <- n*(n + 1)/2
    node <- matrix(0, nrow = n, ncol = n)
    LT <- lower.tri(node, diag = TRUE)
    node[ LT] <- 1:nn
    node[!LT] <- t(node)[!LT]
    return(list(node = node))
  }
}
