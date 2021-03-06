# calculate the incidence rate per person at the final year
incidence <- function( params, fun, tstop, method="rk45dp7", hmin=0.1, hmax=4, maxsteps=500, ... ){
  library(deSolve)
  if( !is.function(fun) )
    stop("Argument fun is not a function!")
  if( is.na(tstop) )
    stop("Please provide the stop time!")
  if( is.null(init_val) )
    stop("Please provide the initial values!")
  if( is.null(params) )
    stop("Please provide the parameter values!")
  
  times <- seq( 0, tstop, 365 ) # outputs every 365 days 
  res <- rk( y=init_val, times=times, func=fun, parms=params, method=method, hmin=hmin, hmax=hmax, maxsteps=maxsteps, ... )  # Integrate ODEs
  res <- res[ , -1 ]
  S <- res[ , index_s ]
  I <- res[ , index_i ]
  R <- res[ , index_r ]
  V <- res[ , index_v ]
  B <- res[ , index_b ]
  CI <- res[ , index_ci ]
  CV <- res[ , index_cv ]

  N <- S + I + R + V
  
  n <- nrow(res)
  annual_inc <- CI[ 2:n, ] - CI[ 1:(n-1), ]
  pyo <- ( N[ 1:(n-1), ] + N[ 2:n, ] )/2 # mean of the two consecutive years  
  # annual_inc_per_person <- annual_inc / pyo
  
  list( params = params,
        fun = fun,
        init = init_val,
        S = S,
        I = I,
        R = R,
        V = V,
        B = B,
        CI = CI,
        CV = CV,
        N = N,
        pyo = pyo ) # annual_inc_per_person has a low
}

inc_model <- function( params, fun, tstop, pyo, ... ){

  res <- incidence( params=params, fun=fun, tstop=tstop, ... )
  # ir_steady <- tail( res$annual_inc_per_person, 1 )
  n <- nrow(res$CI)
  unit_pyo <- tail( res$pyo, 1 )
  inc_pyo <- rep( NA, 4 ) # inc_obs is a global variable
  inc_pyo[1:3] <- (res$CI[n,1:3] - res$CI[(n-1),1:3]) / unit_pyo[1:3] * pyo # incidence by age group over a perid of pyo
  inc_pyo[4] <- (sum(res$CI[n,4:10]) - sum(res$CI[(n-1),4:10]))/sum(unit_pyo[4:10])*pyo
  
  return( inc_pyo )
}

annual_inc_steady <- function( params, ... ){
  library(rootSolve)
  y <- runsteady( y=init_val[1:(4*nag)], func=cholera_sir_smpl_cpp, parms=params, hmin=1e-3, hmax=4 )
  res <- y$y[,-1] # remove time
   
  S <- res[ index_s ]
  I <- res[ index_i ]
  R <- res[ index_r ]
  V <- res[ index_v ]
  N <- S + I + R + V
  foi <- params[1]*sum(I)/sum(N)
  rel_susc <- rep(1,nag)
  rel_susc[1:3] <- params[2:4]
  new_inf <- frac_report*foi*rel_susc*S
  annual_new_inf_pyo <- new_inf*365/N
  return( annual_new_inf_pyo )
}  


compute_R0_water <- function( params ){
  avg_death_rate <- sum(death_rate*age_dist)
  rel_susc <- rep( 1, nag )
  rel_susc[1:3] <- params[2:4]
  R0 <- frac_symptom*(params[1]/(gamma+avg_death_rate))*(rate_excretion/rate_decay)*sum(rel_susc*age_dist)
  
  return( R0 )
}

compute_R0 <- function( params ){
  avg_death_rate <- sum(death_rate*age_dist)
  rel_susc <- rep( 1, nag )
  rel_susc[1:3] <- params[2:4]
  R0 <- (params[1]/(gamma+avg_death_rate))*sum(rel_susc*age_dist)
  
  return( R0 )
}

# # neg_log_lik
# # function to calculate sum of negative log likelihood
# # input is set to have one compoent, par, to ease the use of optim function
# neg_log_lik <- function( params=NULL, data=NULL, ... ){
#   if( is.null(data) ) 
#     data <- inc_rate_obs
#   
#   inc_model <- annual_incidence( params=params, ... )
#   # assume Poisson distribution 
#   (-1)*sum( dpois( data, lambda=inc_model*100000, log=TRUE ) )
# }

neg_log_lik <- function( params, data, fun, tstop, pyo ){
  if( sum(params<0) > 0 | compute_R0(params) <= 1 | compute_R0(params) >= 20 ){
    return ( Inf )
  } 
  inc <- inc_model( params=params, fun=fun, tstop=tstop, pyo=pyo )
  if( sum(is.na(inc)) > 0 ){
    return ( Inf )
  } 
  else {
    return( (-1)*sum( dpois( data, inc, log=TRUE ) ) )
  }
}

# inc_model <- function( params, fun, tstop, pyo, ... ){
#   # ir <- annual_inc_steady( params=params )
#   ir <- annual_inc_tstop( params=params, fun=fun, tstop=tstop, ... )
#   inc_pyo <- ir$annual_inc_pyo*ir$steady_pop/sum(ir$steady_pop)*pyo
#   inc_model <- rep( NA, 4 ) # inc_obs is a global variable
#   inc_model[1] <- inc_pyo[1]
#   inc_model[2] <- inc_pyo[2]
#   inc_model[3] <- inc_pyo[3]
#   inc_model[4] <- sum( inc_pyo[4:10] )
#   
#   return( inc_model )
# }


log_lik <- function ( params, data, fun, tstop, pyo, ...  ) {
  inc <- inc_model( params=params, fun=fun, tstop=tstop, pyo=pyo, ... )
  if( sum(is.na(inc)) > 0 ){
    return( -Inf )
  } 
  else {
    sum( dpois( data, inc, log=TRUE ) )
  }
}

# Prior distribution
log_prior<- function( params ){
  beta_tr = params[1]
  chi_1_tr = params[2]
  chi_2_tr = params[3]
  chi_3_tr = params[4]

  beta_prior = dunif( beta_tr, 1e-6, 1e2, log=T )
  chi_1_prior = dunif( chi_1_tr, 1e-6, 1e2, log=T )
  chi_2_prior = dunif( chi_2_tr, 1e-6, 1e2, log=T )
  chi_3_prior = dunif( chi_3_tr, 1e-6, 1e2, log=T )

  return( beta_prior + chi_1_prior + chi_2_prior + chi_3_prior  )
  
}

log_posterior <- function( params, data, fun, tstop, pyo, ... ){
  # R0 needs to be larger than 1, but not too large  (e.g., 20), in which case ode integrating routine may complain
  # Also, Ro won't be that high and in fact, it will be safe to assume that R0 will be lower than 10
  # parameters have to be positive
  if( sum(params<0) > 0 | compute_R0_water(params) <= 1 | compute_R0_water(params) >= 30 ){
    return ( -Inf )
  } 
  
  else {
    # cat( "log_posterior, ", params[1], ", ", params[2], ", ", params[3], ", ", params[4], "\n" )
    lik <- log_lik( params=params, data=data, fun=fun, tstop=tstop, pyo=pyo, ... )
    prior <- log_prior( params )
    
    return ( lik + prior )
  }
}

run_MCMC <- function( startvalue=c(0.2,1,1,1), iter=100, burnin = round(iter/2), scale=rep(0.1,4), show_progress_bar=TRUE, data, fun, tstop, pyo ){
  if (show_progress_bar) {
    pb <- txtProgressBar( min=0, max=iter, style=3 )
  }
  update_step <- max( 5, floor(iter/100) )
  
  accept <- numeric( iter )
  npar <- length( startvalue )
  chain = matrix( NA, nrow=iter, ncol=(npar+1) ) # unknown parameters + likelihood store
  chain[ 1, 1:npar ] <-  startvalue
  # chain[ 1, (npar+1) ] <- log_posterior( data, chain[1,1:npar] )
  chain[ 1, (npar+1) ] <- log_posterior( data=data, params=chain[1,1:npar], fun=fun, tstop=tstop )
  
  for (i in 2:iter) {
    # proposal function
    if (show_progress_bar && i%%update_step == 0) {
      setTxtProgressBar(pb, i)
    }
    proposal <- rnorm( 4, mean = chain[(i-1),1:4], sd = sqrt(scale) )
    # repeat{
    #   proposal <- rnorm( 4, mean = chain[(i-1),1:4], sd = scale )
    #   if( sum(proposal<0) == 0 & compute_R0(proposal) > 1 ){
    #     break
    #   }
    # }
    # posterior_proposal <- log_posterior( data, proposal )
    posterior_proposal <- log_posterior( data=data, params=proposal, fun=fun )
    alpha <- exp( posterior_proposal - chain[ (i-1), (npar+1) ] )
    # cat( "\niter = ", i, ", log density = ", posterior_proposal, ", log density prev = ", chain[ (i-1), (npar+1) ], ", alpha = ", alpha )
    if( !is.finite(alpha) ){ 
      alpha <- 0
    }
    if( runif(1) < alpha ) {
      chain[ i, 1:npar ] <- proposal
      chain[ i, (npar+1) ] <- posterior_proposal
      accept[ i ] <- 1
      # cat( "\niter = ", i, ", proposal accepted:", chain[ i,],"\n")
    } else { 
      chain[ i, 1:npar ] <- chain[ (i-1), 1:npar ]
      chain[ i, (npar+1)] <- chain[ (i-1), (npar+1) ]
    } 
  } 
  # cat( "\nAcceptance ratio = ", sum(accept[(burnin + 1):iter]) / (iter - burnin), "\n" )
  
  list( theta = chain[ (burnin + 1):iter, (1:npar)],
        loglik_posterior = chain[ (burnin + 1):iter, (npar+1)], 
        acceptance_ratio = sum(accept[(burnin + 1):iter]) / (iter - burnin))
 
}
