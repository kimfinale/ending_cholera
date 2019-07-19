cholera_sir_r <- function( t, y, params ){
  source( "data/params.R" )
  S <- y[ 1 : num_age_grp ]
  I <- y[ (1+num_age_grp) : (2*num_age_grp) ]
  R <- y[ (1+2*num_age_grp) : (3*num_age_grp) ]
  V <- y[ (1+3*num_age_grp) : (4*num_age_grp) ]
  # 
  dSdt <- rep( 0, num_age_grp )
  dIdt <- rep( 0, num_age_grp )
  dRdt <- rep( 0, num_age_grp )
  dVdt <- rep( 0, num_age_grp )
  
  with(as.list(c(y, params)),{
    beta <- beta
    gamma <- gamma
    rel_susc[1] <- chi_1
    rel_susc[2] <- chi_2
    rel_susc[3] <- chi_3
    rel_susc[4:9] <- chi_4
    # Derived variables
    N <- sum( S ) + sum( I ) + sum( R ) + sum( V )
    mu <- death_rate
    # births <- birth_rate / sum(sample_init_val)
    births <- birth_rate
    
    foi <- beta * sum(I) / N
    dCIdt <- rep( 0, num_age_grp )
    dCVdt <- rep( 0, num_age_grp )
    
    for( i in  1:num_age_grp ) {
      if( i == 1 ){
        dSdt[i] = + births - (foi*rel_susc[i] + mu[i] + vacc_rate_campaign[i] + ag[i])*S[i]
        dIdt[i] = + foi*rel_susc[i]*S[i] - (mu[i] + ag[i] + gamma)*I[i]
        dRdt[i] = + (1 - case_fatality[i])*gamma*I[i] - (mu[i] + ag[i] + vacc_rate_campaign[i])*R[i]
        dVdt[i] = + vacc_rate_campaign[i]*(S[i] + R[i]) - (mu[i] + ag[i])*V[i]		
      } else {
        dSdt[i] = + ag[i-1]*S[i-1] - (foi*rel_susc[i] + mu[i] + vacc_rate_campaign[i] + ag[i])*S[i]
        dIdt[i] = + ag[i-1]*I[i-1] + foi*rel_susc[i]*S[i] - (mu[i] + ag[i] + gamma)*I[i]
        dRdt[i] = + ag[i-1]*R[i-1] + (1 - case_fatality[i])*gamma*I[i] - (mu[i] + ag[i] + vacc_rate_campaign[i])*R[i]
        dVdt[i] = + ag[i-1]*V[i-1] + vacc_rate_campaign[i]*(S[i] + R[i]) - (mu[i] + ag[i])*V[i]
      }
      
      dCIdt[i] = foi*rel_susc[i]*S[i] # cumulative incidence by age group
      dCVdt[i] = vacc_rate_campaign[i]*(S[i] + R[i]) # cumulative number of vaccine doses by age group
    }
    
    dydt <- c( dSdt, dIdt, dRdt, dVdt, dCIdt, dCVdt )
    list( dydt )
  })
}