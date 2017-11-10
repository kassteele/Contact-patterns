options(
  # Global options for the function integrate used in the package.
  priorSens_integrate_rel.tol = .Machine$double.eps^0.4,
  priorSens_integrate_abs.tol = .Machine$double.eps^0.4,
  priorSens_integrate_subdivisions = 100L,

  # Global options for computing the grid
  priorSens_grid_number_axis_points = 200,
  priorSens_grid_epsilon = 0.00354, 
  priorSens_grid_root_method = "uniroot",
  priorSens_grid_parallel = TRUE,
  priorSens_grid_coord_explr_method = "nlminb",
  priorSens_grid_coord_explr_method_options = NULL,
  priorSens_grid_HD_integration_method = "simpson13rule",
  priorSens_grid_HD_integration_domain = "compute_domain",
  priorSens_grid_HD_integration_domain_quantile = 1e-18,
  priorSens_grid_HD_integration_number_points = 100,
  
  
  # integration method used to normalize corrected posterior distributions
  priorSens_correc_post_integration_method = "simpson13rule",
  # integration method used to within the 'compute_hellinger_distance_inla'
  priorSens_inla_HD_integration_method = "simpson13rule", 
  
  # Global options for tests
  priorSens_test_computefinal = TRUE
  )

