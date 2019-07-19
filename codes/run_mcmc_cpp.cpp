#include <Rcpp.h>
using namespace Rcpp;

double compute_R0( NumericVector params ){
  double nag = Environment::global_env()["nag"];
  NumericVector rel_susc ( nag );
  for( int i = 0; i < nag; ++i ){
    rel_susc[ i ] = 1;// reference value = 1
  } 
  for( int i = 0; i < (params.size()-1); ++i ){
    rel_susc[i] = params[i+1];
  }
  NumericVector age_dist = Environment::global_env()["age_dist"];
  double gamma = Environment::global_env()["gamma"];
  double R0 = (params[0]/gamma)*sum(rel_susc*age_dist);
  return( R0 );
}

NumericVector inc_model ( NumericVector params, double pyo ){
  Function annual_inc_steady = Environment::global_env()["annual_inc_steady"]; 
  NumericVector ir = annual_inc_steady( params );
  NumericVector init_val = Environment::global_env()["init_val"];
  NumericVector age_dist = Environment::global_env()["age_dist"];
  NumericVector data = Environment::global_env()["inc_obs_Jakarta"];
  NumericVector N  = age_dist*sum(init_val);
  NumericVector inc = ir*N*pyo;
  
  NumericVector inc_model ( data.size() ); 
  inc_model[0] = (inc[0])/(N[0]);
  inc_model[1] = (inc[1])/(N[1]);
  inc_model[2] = (inc[2])/(N[2]);
  inc_model[3] = sum(inc[Rcpp::Range(3,9)])/sum(N[Rcpp::Range(3,9)]);
  
  return( inc_model );
}

// log likelihood
double log_lik( NumericVector params ) {
  IntegerVector data = Environment::global_env()["inc_obs_Jakarta"];
  NumericVector model = inc_model( params, 100000 );
  NumericVector log_lik (data.size());
  for( int i=0; i < data.size(); ++i ){
    log_lik[ i ] = R::dpois( data[i], model[i], true );  
  }
  
  // double log_lik = R::dpois( data[0], model[0], true );
  double sum_log_lik = sum( log_lik );
  return( sum_log_lik );  
}


// Prior distribution
double log_prior( NumericVector params ){
  NumericVector log_density( params.size() );

  for( int i = 0; i < params.size(); ++i ){
    log_density[ i ] = R::dunif( params[i], 1e-6, 1e2, true ); // 1e-6 - 1e2 is assumed to be large enough 
  }
  
  double log_p_sum = sum( log_density );
  return( log_p_sum );

}


double log_posterior( NumericVector params ){
// # R0 needs to be larger than 1, but not too large  (e.g., 10), in which case ode integrating routine may complain
// # Also, Ro won't be that high and in fact, it will be safe to assume that R0 will be lower than 10
// # parameters have to be positive
  // Function compute_R0 = Environment::global_env()["compute_R0"];
  // double R0 = compute_R0( params );
  double nag = Environment::global_env()["nag"];
  NumericVector rel_susc ( nag );
  for( int i = 0; i < nag; ++i ){
    rel_susc[ i ] = 1;// reference value = 1
  } 
  for( int i = 0; i < (params.size()-1); ++i ){
    rel_susc[i] = params[i+1];
  }
  NumericVector age_dist = Environment::global_env()["age_dist"];
  double gamma = Environment::global_env()["gamma"];
  double R0 = (params[0]/gamma)*sum(rel_susc*age_dist);
  
  if( sum(params<0) > 0 | R0 <= 1 | R0 >= 10 ){
    return ( R_NegInf );
  } else {
// # cat( "log_posterior, ", params[1], ", ", params[2], ", ", params[3], ", ", params[4], "\n" )
    return ( log_lik(params) + log_prior(params) );
  }
}

// [[Rcpp::export]]
List run_MCMC( NumericVector data, NumericVector startvalue, int iter ){

  NumericVector accept( iter );
  double npar = startvalue.size();
  NumericMatrix chain( iter, npar+1 ); // unknown parameters + posterior store
  // chain( 0, Rcpp::Range(0, npar-1) ) =  startvalue;
  chain( 0, 0 ) =  startvalue[ 0 ];
  chain( 0, 1 ) =  startvalue[ 1 ];
  chain( 0, 2 ) =  startvalue[ 2 ];
  chain( 0, npar) = log_posterior( data, startvalue );
  // Rcout << chain( 0, Rcpp::Range(0,2));
  // for (i in 2:iter) {
  //   NumericVector proposal = rnorm( 4, mean = chain( (i-1), 0:3 ), sd = sqrt(scale) );
  //   NumericVector   posterior_proposal = log_posterior( data, proposal );
  //   double alpha = exp( posterior_proposal - chain[ (i-1), (npar+1) ] );
  //     if( !is.finite(alpha) ){
  //       alpha = 0
  //     }
  //     if( runif(1) < alpha ) {
  //       chain( i, 0:(npar-1) ) = proposal;
  //       chain( i, (npar+1) ) = posterior_proposal;
  //       accept[ i ] = 1;
  //     } else {
  //       chain( i, 0:(npar-1) ) = chain( (i-1), 0:(npar-1) );
  //       chain( i, npar ) = chain( (i-1), (npar+1) );
  //     }
  // }
    // List list( theta = chain[ (burnin + 1):iter, (1:npar)],
    //       log_posterior = chain[ (burnin + 1):iter, (npar+1)],
    //                            acceptance_ratio = sum(accept[(burnin + 1):iter]) / (iter - burnin))
      
  return ( List::create( chain ) );

}
