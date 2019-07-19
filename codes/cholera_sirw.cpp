#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List cholera_sirw( double t, NumericVector y, NumericVector params ){

  NumericVector ag = Environment::global_env()["ag"];
  NumericVector init_val = Environment::global_env()["init_val"];
  NumericVector mu = Environment::global_env()["death_rate"];
  NumericVector case_fatality = Environment::global_env()["case_fatality"];
  NumericVector rel_susc = Environment::global_env()["rel_susc"];
  
  double birth_rate = Environment::global_env()["birth_rate"];
  double gamma = Environment::global_env()["gamma"];
  int num_age_grp = Environment::global_env()["num_age_grp"];
  double frac_report = Environment::global_env()["frac_report"];
  double frac_symptom = Environment::global_env()["frac_symptom"];
  double rate_excretion = Environment::global_env()["rate_excretion"];
  double rate_decay = Environment::global_env()["rate_decay"];
  //vaccination
  NumericVector cov_vacc_campaign = Environment::global_env()["cov_vacc_campaign"];
  NumericVector cov_vacc_routine = Environment::global_env()["cov_vacc_routine"];
  NumericVector rate_vacc_campaign = Environment::global_env()["rate_vacc_campaign"];
  
  bool vacc_campaign = Environment::global_env()["vacc_campaign"];
  // bool vacc_routine = Environment::global_env()["vacc_routine"];
  double start_vacc_campaign = Environment::global_env()["start_vacc_campaign"];
  double stop_vacc_campaign = Environment::global_env()["stop_vacc_campaign"];
  double dur_vacc_campaign = Environment::global_env()["dur_vacc_campaign"];
  NumericVector eff_vacc = Environment::global_env()["eff_vacc"];
  
  rate_vacc_campaign = - log(1-cov_vacc_campaign*eff_vacc) / dur_vacc_campaign;
  
  if( t < start_vacc_campaign || t > stop_vacc_campaign || !vacc_campaign ){
    rate_vacc_campaign = rate_vacc_campaign * 0.0;
  }
  
  // Rprintf("the value of rate_excretion %f at time %f\n", rate_excretion, t );
  
  NumericVector S = y[ Rcpp::Range( 0, (num_age_grp-1)) ];
  NumericVector I = y[ Rcpp::Range( num_age_grp, (2*num_age_grp-1)) ];
  NumericVector R = y[ Rcpp::Range( 2*num_age_grp, (3*num_age_grp-1)) ];
  NumericVector V = y[ Rcpp::Range( 3*num_age_grp, (4*num_age_grp-1)) ];
  NumericVector CV = y[ Rcpp::Range( 4*num_age_grp, (5*num_age_grp-1)) ];
  NumericVector CI = y[ Rcpp::Range( 5*num_age_grp, (6*num_age_grp-1)) ];
  double B = y[ 6*num_age_grp ];
     
  NumericVector dSdt( num_age_grp );
  NumericVector dIdt( num_age_grp );
  NumericVector dRdt( num_age_grp );
  NumericVector dVdt( num_age_grp );
  NumericVector dCVdt( num_age_grp );
  NumericVector dCIdt( num_age_grp );
  double dBdt;
  
  // parameter to be estimated  
  double beta = params[ 0 ];
  rel_susc[ 0 ] = params[ 1 ];
  rel_susc[ 1 ] = params[ 2 ];
  rel_susc[ 2 ] = params[ 3 ];

  double births = birth_rate;
  // Rprintf("the value of beta : %f \n", beta );
  double I_symptom = frac_symptom * sum( I );
  double foi = beta * B; 
  for(int i = 0; i < num_age_grp; ++i){
     if( i == 0 ){
       dSdt[i] = + births - (foi*rel_susc[i] + mu[i] + rate_vacc_campaign[i] + ag[i])*S[i];
       dIdt[i] = + foi*rel_susc[i]*S[i] - (mu[i] + ag[i] + gamma)*I[i];
       dRdt[i] = + (1 - case_fatality[i])*gamma*I[i] - (mu[i] + ag[i] + rate_vacc_campaign[i])*R[i];
       dVdt[i] = + rate_vacc_campaign[i]*(S[i] + R[i]) - (mu[i] + ag[i])*V[i];		
     } else {
       dSdt[i] = + ag[i-1]*S[i-1] - (foi*rel_susc[i] + mu[i] + rate_vacc_campaign[i] + ag[i])*S[i];
       dIdt[i] = + ag[i-1]*I[i-1] + foi*rel_susc[i]*S[i] - (mu[i] + ag[i] + gamma)*I[i];
       dRdt[i] = + ag[i-1]*R[i-1] + (1 - case_fatality[i])*gamma*I[i] - (mu[i] + ag[i] + rate_vacc_campaign[i])*R[i];
       dVdt[i] = + ag[i-1]*V[i-1] + rate_vacc_campaign[i]*(S[i] + R[i]) - (mu[i] + ag[i])*V[i];
     }
     dCIdt[i] = frac_report*frac_symptom*foi*rel_susc[i]*S[i]; 
     dCVdt[i] = rate_vacc_campaign[i]*(S[i] + R[i]);
   }
  dBdt = rate_excretion*I_symptom - rate_decay*B;
  
  //output
  NumericVector out( init_val.size() );
  for(int i = 0; i < num_age_grp; ++i) {
    out[ i ] = dSdt[i];
    out[ i + num_age_grp ] = dIdt[i];
    out[ i + 2*num_age_grp ] = dRdt[i];
    out[ i + 3*num_age_grp ] = dVdt[i];
    out[ i + 4*num_age_grp ] = dCIdt[i];
    out[ i + 5*num_age_grp ] = dCVdt[i];
  }
  out[ 6*num_age_grp ] = dBdt;
  
  return ( List::create( out ) );
    
}
