#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List cholera_sirw_smpl_cpp( double t, NumericVector y, NumericVector params ){

  NumericVector ag = Environment::global_env()["ag"];
  NumericVector init_val = Environment::global_env()["init_val"];
  NumericVector mu = Environment::global_env()["death_rate"];
  NumericVector vacc_cov_campaign = Environment::global_env()["vacc_cov_campaign"];
  NumericVector vacc_cov_routine = Environment::global_env()["vacc_cov_routine"];
  NumericVector case_fatality = Environment::global_env()["case_fatality"];
  NumericVector rel_susc = Environment::global_env()["rel_susc"];
  
  double dur_campaign = Environment::global_env()["dur_campaign"];
  double birth_rate = Environment::global_env()["birth_rate"];
  double gamma = Environment::global_env()["gamma"];
  int num_age_grp = Environment::global_env()["num_age_grp"];
  
  NumericVector vacc_rate_campaign( 9 );
  vacc_rate_campaign = - log(1-vacc_cov_campaign) / dur_campaign;
   
  NumericVector S = y[ Rcpp::Range(0, (num_age_grp-1)) ];
  NumericVector I = y[ Rcpp::Range(num_age_grp, (2*num_age_grp-1)) ];
  NumericVector R = y[ Rcpp::Range(2*num_age_grp, (3*num_age_grp-1)) ];
  NumericVector V = y[ Rcpp::Range(3*num_age_grp, (4*num_age_grp-1)) ];
     
  NumericVector dSdt( num_age_grp );
  NumericVector dIdt( num_age_grp );
  NumericVector dRdt( num_age_grp );
  NumericVector dVdt( num_age_grp );
    
  double beta = params[ 0 ];
  rel_susc[ 0 ] = params[ 1 ];
  rel_susc[ 1 ] = params[ 2 ];
  rel_susc[ 2 ] = params[ 3 ];

  double births = birth_rate;
  
  // Rprintf("the value of beta : %f \n", beta );
  double N = sum( S ) + sum( I ) + sum( R ) + sum( V );
  double I_tot = sum( I );
  double foi = beta * I_tot / N;
   
  for(int i = 0; i < num_age_grp; ++i){
     if( i == 0 ){
       dSdt[i] = + births - (foi*rel_susc[i] + mu[i] + vacc_rate_campaign[i] + ag[i])*S[i];
       dIdt[i] = + foi*rel_susc[i]*S[i] - (mu[i] + ag[i] + gamma)*I[i];
       dRdt[i] = + (1 - case_fatality[i])*gamma*I[i] - (mu[i] + ag[i] + vacc_rate_campaign[i])*R[i];
       dVdt[i] = + vacc_rate_campaign[i]*(S[i] + R[i]) - (mu[i] + ag[i])*V[i];		
     } else {
       dSdt[i] = + ag[i-1]*S[i-1] - (foi*rel_susc[i] + mu[i] + vacc_rate_campaign[i] + ag[i])*S[i];
       dIdt[i] = + ag[i-1]*I[i-1] + foi*rel_susc[i]*S[i] - (mu[i] + ag[i] + gamma)*I[i];
       dRdt[i] = + ag[i-1]*R[i-1] + (1 - case_fatality[i])*gamma*I[i] - (mu[i] + ag[i] + vacc_rate_campaign[i])*R[i];
       dVdt[i] = + ag[i-1]*V[i-1] + vacc_rate_campaign[i]*(S[i] + R[i]) - (mu[i] + ag[i])*V[i];
     }
     
   }
   
  //output
  NumericVector out( y.size() );
  // out = c( dSdt, dIdt, dRdt, dVdt, dCIdt, dCVdt  );
  for(int i = 0; i < num_age_grp; ++i) {
    out[ i ] = dSdt[i];
    out[ i + num_age_grp ] = dIdt[i];
    out[ i + 2*num_age_grp ] = dRdt[i];
    out[ i + 3*num_age_grp ] = dVdt[i];
  }
  return ( List::create( out ) );
  // return List::create( dSdt, dIdt, dRdt, dVdt, dCIdt, dCVdt );
    
}
