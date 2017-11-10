# Q_MR_1: In which scale the posterior inla marginal is usually
# normalized?
# Q_MR_2: How (in which case) do we know that the posterior inla
# marginal is already correctly normalized for the hellinger distance
# computation and there is no need for any renormalization?

## A: The inla marginal is returned as a matrix with x and y coordinates.
## The normalization is made using the x-points as the grid for integration.
## So anything you do with the grid returned by INLA you can consider the marginal
## as normalized. But of course this is a numerical discrete representation
## of the marginal and the normalization is only valid on this scale. If you change the
## scale you should in theory renormalize it. But the change is actually very small and 
## doesn't impact most of the computations, it only matter to us because we want to perform
## infinitesimal computations. So, for example, when we use the 'integrate' function it 
## probabably uses a different grid to perform the integration and that is why we had problems 
## that required renormalization.

# Q_MR_3: Are we sure that the domains in inla_posterior_1[,1] and
# inla_posterior_1[,2] are identical? Should we check for it?

## A: Yes, you are right, I have added a check for this 
##   if (!all(inla_posterior_1[,1] == inla_posterior_1[,2])) 
##    stop("compute_hellinger_distance: Domains of the inla marginals should be the same.")

# Q_MR_4: There are cases, when bc_coef>1. Should bc_coef be set to 1
# in such a case?

## A: This is strange and seems pathological, could you give me examples? We should examine this
## before setting bc = 1.

# Q_MR_5: I have replaced "eps" in the code by grid_epsilon in order
# to make the code more consistent. Do you agree?

## A: Yes, I agree.

# Q_MR_6: I have replaced "eps" in the code by grid_epsilon in order
# to make the code more consistent. Do you agree?

## A: Yes, I agree.

# Q_MR_7: Do you agree that return(c(mean_inf, mean_sup, 
# rep(new_precision,2))) is used? This was we will have mean values 
# first and precision values afterwards. For me it is a more natural
# way of parametrisation of the normal distribution.

## A: code has changed

# Q_MR_8: If you agree with Q_MR_7, we should adjust here
# colnames(grid_roots) <- c("mean", "precision"). Is it OK for you? 

## A: code has changed

# Q_MR_9: I am very sorry to hear that. Is this difficulty possibly
# connected to the extremely wide interval = c(0, sup_interval_shape)
# for uniroot.all? 

## A: Not sure.

# Q_MR_10: Maybe we could use the log-Gamma representation leading to
#    g <- function(x, prior_shape, prior_rate, new_shape, new_rate){
#      exp(x + 0.5 *(dgamma(exp(x), shape = prior_shape, rate = prior_rate, log = TRUE) + 
#                  dgamma(exp(x), shape = new_shape, rate = new_rate, log = TRUE))) 
#    }

## A: No, I might be wrong on this, but we cannot use log-gamma distribution instead
## of the Gamma distribution, because if we do we are going to get the Hellinger
## distance between log-gammas, which is different than the hellinger distance between Gammas.
## Does this makes sense?

# Q_MR_11: If you agree with Q_MR_10 the lower integration range
# should be set to lower = -Inf

## A: Yes, you are right.

# Q_MR_12: Is it necessary for the function to be "hidden" by starting
# the name by a dot? What is the purpose of this point in front of the
# function name?

## A: No, it is not necessary, just a choice I have made.

# Q_MR_13: Deals with all interval settings below: Is interval = c(0,
# sup_interval_shape) too wide for uniroot.all to give precise
# results?

## A: I hope that with polar coordinates this will not be an issue anymore.

# Q_MR_14: Should we mention that only internal marginal posteriors
# (in log-scale) can be used here? Do we need to check if the 
# inla_marginal_posterior is in the correct scale?

## A: Yes, we should mention. Not sure if there is a way to check that, any suggestion?


#----------------------------------

QUESTIONS


1. Could you please explain to me, what is your exact motivation for considering

new_shape = prior_shape + (exp(x) * cos(fixed_theta))*prior_shape
new_rate =  prior_rate + (exp(x) * sin(fixed_theta))*prior_rate?

As I have to understand the idea behind the function just from the code I have decided to write down my thoughts explicitly:
  
  The common transformations prior_shape = r * cos(fixed_theta) and prior_rate = r * sin(fixed_theta) would already transform the polar (fixed_theta, r) coordinates to the Cartesian (prior_shapa, prior_scale) coordinates.
However, such parametrisation is not satisfactory in order to have a stable running gamma grid search function.
We need several amendments.

My motivation for for parametrising x = log(r) would be that this way we can choose a more stable search_interval. 
Do you agree?

As I understand considering 

new_shape = prior_shape + (exp(x) * cos(fixed_theta))*prior_shape 
and
new_rate =  prior_rate + (exp(x) * sin(fixed_theta))*prior_rate 

we are actually looking for a scaling factor "((exp(x) * cos(fixed_theta)), (exp(x) * sin(fixed_theta)))" which transforms (prior_shape, prior_rate) into a grid_epsilon-distant (new_shape, new_rate).
Do you agree with my interpretation?

My motivation for such a transformation choice would be that this way we can apply the function to any arbitrary (prior_shape, prior_rate) vector.
Do you agree?

2. I have just tried to provide a preliminary sketch of the idea behind the compute_grid_gamma_polar function, which we could use as a section in the manuscript.
Could you please let me know if you agree with this description.
Your corrections and changes are welcome!
  
  In order to find an epsilon gamma grid a suitable transformation of the cartesian (shape, scale) coordinates to the polar coordinates (theta, r) was used, where theta and r denotesthe angle and modulus, respectively.
For sake of stability of the algorithm log(r) = x was considered.
We looked for a scaling factor "((exp(x) * cos(fixed_theta)), (exp(x) * sin(fixed_theta)))" which transforms (prior_shape, prior_rate) into a grid_epsilon-distant (new_shape, new_rate) by finding the roots of the analytical equation

LaTeX equation corresponding to compute_hellinger_gammas(prior_shape = prior_shape, 
                                                         prior_rate = prior_rate, 
                                                         new_shape = new_shape, 
                                                         new_rate = new_rate,
                                                         analytic = analytic) - grid_epsilon

numerically by means of the uniroot.all function with the interval ranging from -10 to 10.
The gamma_polar_grid scaling factors found were transformed back to cartesian coordoinates by using

LaTeX equation corresponding to prior_shape + (gamma_grid_polar[,"r"] * cos(gamma_grid_polar[,"theta"]))*prior_shape

and

LaTeX equation corresponding to prior_rate + (gamma_grid_polar[,"r"] * sin(gamma_grid_polar[,"theta"]))*prior_rate

Transformations chosen above guarantee that the function can be applied to any arbitrary (prior_shape, prior_rate) vector.


3. As we use the polar mathodology for finding the gamma grid, should we also adjust the methodology for finding the normal grid accordingly?
What do you think?
