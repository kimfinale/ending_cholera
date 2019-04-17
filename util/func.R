# calculate the incidence at the final year 
annual_incidence <- function( tstop=600, dt=365, integ_method="rk45dp7", params=NULL ){
  # if( !is.function(FUN) ) 
  #   stop("Argument FUN is not a function!")
  library(deSolve)
  if( is.na(tstop) )
    stop("Please provide the stop time!")
  if( is.null(init_val) )
    stop("Please provide the initial values!")
  if( is.null(params) )
    stop("Please provide the parameter values!")
  nday <- 365
  times <- seq( 0, tstop*nday, dt )
  # Integrate ODEs
  out <- rk( y=init_val, times=times, func=cholera_sirw_cpp, parms=params, method=integ_method )
  # thin output by year
  thin_row <- seq( 1 , (tstop*nday/dt+1), by=nday/dt )
  out <- out[ thin_row, ] # make sure that the output appears by year
  S <- out[ , index_S ]
  I <- out[ , index_I ]
  R <- out[ , index_R ]
  V <- out[ , index_V ]
  CI <- out[ , index_CI ]
  CV <- out[ , index_CV ]
  N <- S + I + R + V
  annual_CI <- CI[ 2:nrow(CI), ] - CI[1:(nrow(CI)-1),]
  pyo <- (N[1:tstop,] + N[2:(tstop+1),])/2 # beginning of the two consecutive years  
  annual_CI_pyo <- annual_CI / pyo
  
  tail( annual_CI_pyo, 1 )  
}

annual_inc_steady <- function( params, ... ){
  library(rootSolve)
  y <- runsteady( y=init_val[1:(4*nag)], func=cholera_sirw_smpl_cpp, parms=params, hmin=0.2, hmax=4 )
  out <- y$y
  S <- out[ index_S-1 ]
  I <- out[ index_I-1 ]
  R <- out[ index_R-1 ]
  V <- out[ index_V-1 ]
  N <- S + I + R + V
  foi <- params[1]*sum(I)/sum(N)
  rel_susc <- rep(1,nag)
  rel_susc[1:3] <- params[2:4]
  new_inf <- frac_report*foi*rel_susc*S;
  annual_new_inf_pyo <- new_inf*365/N
  return( annual_new_inf_pyo )
}  

compute_R0 <- function( params ){
  rel_susc <- rep(1,nag)
  rel_susc[1:3] <- params[2:4]
  R0 <- (params[1]/gamma)*sum(rel_susc*age_dist)
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

neg_log_lik <- function( params=NULL, data=NULL ){
  if( is.null(data) )
    data <- inc_rate_obs
  
  inc_model <- annual_inc_steady( params=params )
  (-1)*sum( dpois( data, lambda=inc_model*100000, log=TRUE ) )
}


log_lik <- function ( params ) {
  inc_model <- annual_inc_steady( params=params )
  sum( dpois( inc_rate_obs, lambda=inc_model*100000, log=TRUE ) )
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

log_posterior <- function( params ){
  return ( log_lik(params) + log_prior(params) )
}

run_MCMC <- function( startvalue=c(0.2,1,1,1), iter=100, scale=rep(0.1,4), showProgressBar=TRUE ){
  if (showProgressBar) {
    pb <- txtProgressBar(min = 0, max = iter, style = 3)
  }
  update_step <- max(5, floor(iter/100))
  
  naccepted = 1
  npar <- length(startvalue)
  chain = matrix( nrow=(iter+1), ncol=(npar+1) ) # unknown parameters + likelihood store
  chain[ 1, 1:npar ] <-  startvalue
  check_iter <- 0
  chain[ 1, (npar+1) ] <- log_posterior( chain[1,1:npar] )
  
  for (i in 1:iter) {
    # proposal function
    if (showProgressBar && i%%update_step == 0) {
      setTxtProgressBar(pb, i)
    }
    while (check_iter==0) { 
      proposal <- rnorm( 4, mean = chain[i,1:4], sd = scale )
      R0_proposal <- compute_R0( proposal ) 
      if( sum(proposal<0) == 0 & R0_proposal > 0 ) {
        check_iter <- 1
      }
    }
    check_iter <- 0
    posterior_proposal <- log_posterior( proposal )
    probab <- exp( posterior_proposal - chain[ i, (npar+1) ] ) # add syntax to prevent NaN trigger 
    
    if (runif(1) < probab) {
      chain[ i+1, 1:npar ] <-  proposal
      chain[ i+1, (npar+1) ] <- posterior_proposal
      naccepted <- naccepted + 1
      cat( "iter = ", i, ", accepted parameters:", chain[i+1,],"\n")
    } else { 
      chain[ i+1, 1:npar ] <- chain[ i, 1:npar ]
      chain[ i+1, (npar+1)] <- chain[ i, (npar+1) ]
    } 
  } 
  cat( "Acceptance ratio = ", naccepted / iter, "\n" )
  return(chain)
}
