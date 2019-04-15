# neg_log_lik
# function to calculate sum of negative log likelihood
# input is set to have one compoent, par, to ease the use of optim function
neg_log_lik <- function( params=NULL, data=NULL ){
  stopifnot( length(par) == 5, length(data) == 8, !is.null(init_val_java), !is.null(run_java) )
  beta <- par[ 1 ]
  waning_natural <- par[ 2 ]
  frac_report <- par[ 3 ]
  chi <- rep( 1, num_age_grp )
  for( i in 1:2 ) { # chi[ index_age_grp[[ num_index_age_grp ]] ] <- 1, i.e., reference
    chi[ iag[[ i ]] ] <- par[ i+3 ]
  }
  # model runs up to 52 years such that the comparison 1996-9 (47-50) can be done safely
  t_stop <- par
  rJava::.jcall( init_val_java, "V", "setNumOutPoints", as.integer( t_stop + 1 ) )
  annual_incidence <- function( params=NULL ){
    times <- seq( 0, 60*365, 365 )
    source( "params.R" )
    params <- c( beta=params[1], gamma=params[2] )
    yini <-  sample_Y0
    # Integrate ODEs
    out <- ode( y=yini, times=times, func=cholera_model, parms=params )
    CI <- out[ , 4*num_age_grp+2]
    annual_CI <- CI[2:length(CI)] - CI[1:(length(CI)-1)]
    tail(annual_CI,1)
  }  
  
  m <- sapply( out, rJava::.jevalArray ) # out is converted to m
  
  index_yr_surv_start <- nrow( m ) - 11
  # starting index for different infection states, e.g., 0-1 yo
  pop_age <- matrix( NA, nrow=2, ncol=niag )
  inc_age <- matrix( NA, nrow=2, ncol=niag )
  for( j in 1:2 ){
    for( k in 1:niag ){ # sum across age groups
      index1 <- sapply( 1:9, function(x) { is[ x ] + iag[[ k ]]  } )
      index2 <- sapply( 10:11, function(x) { is[ x ] + iag[[ k ]] } )
      pop_age[ j, k ] <- sum( m[ index_yr_surv_start - 1 + j, as.vector( index1 ) ] ) # the first col is time
      inc_age[ j, k ] <- sum( m[ index_yr_surv_start - 1 + j, as.vector( index2 ) ] )
    }
  }
  inc_rate_model <- rep( NA, niag )
  for( j in 1:niag ){
    inc_rate_model[ j ] <- 10^5*(inc_age[2,j] - inc_age[1,j])/((pop_age[2,j]+pop_age[1,j])/2)
  }

  neg_log_lik <- - sum( dpois( data, lambda=inc_rate_model, log=TRUE ) )
}