.onLoad <- function(libname = find.package("priorSens"), pkgname = "priorSens") {
  
  if (!require(multicore)){
    
    options(priorSens_grid_parallel = FALSE)
    warning("Try to install multicore package to use parallel computation when computing the grid.")
    
  } else {
    
    options(priorSens_grid_parallel = TRUE)
    
  }
  
}