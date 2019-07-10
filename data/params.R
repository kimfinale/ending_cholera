num_age_grp <- 10L
beta <- 0.2
gamma <- 1/7
ag <- c( 1/365, 1/(4*365), rep( 1/(10*365), 8) ) # <1, 1-4, 5-14, 15-24, ... 65-74, 75-84
rate_excretion <- 1/7 #Rate of Vibrio cholerae from acutely or chronically infecteds
rate_decay <- 1 #Rate of decay of infectious particles from water supply, average duration of viability, 3 weeks
# Vaccination
vacc_cov_campaign <- rep( 0, num_age_grp )
vacc_cov_routine <- rep( 0, num_age_grp )
dur_campaign <- 14.0
vacc_eff <- rep( 0, num_age_grp )
vacc_rate_campaign <- - log(1-vacc_cov_campaign) / dur_campaign
case_fatality <- rep( 0, num_age_grp )
#

rel_susc <- rep( 1, num_age_grp )
susc_age_grp <- list( 1, 2, 3, 4:10 ) # age group 4:9 are assumed to have the same susceptibility

init_val <- c( 0.019129528368344, 0.0674741649748425, 0.126429307827068, 0.0936790567827395, 
                0.0689075046855928, 0.0497120079378464, 0.0339845165599162, 0.0202690060767594, 
                0.00908927231705442, 0.0024549202470792, 1.17987107841576e-05, 
                4.24569334984788e-05, 7.95619657527049e-05, 5.89522483326505e-05, 
                4.33635058955322e-05, 3.12837761159886e-05, 2.1386462779773e-05, 
                1.27552894059669e-05, 5.71988080987387e-06, 1.54488178164627e-06, 
                0.000594187298716229, 0.0109978614296745, 0.0673712858210493, 
                0.0949157764370814, 0.112867457328733, 0.121100057079751, 0.11573977692341, 
                0.0902496016197851, 0.0488325131992378, 0.0146789226129884, rep(0, 30) )
# # rates are given per year
death_rate <-read.csv("data/death_rate.csv")
death_rate <- death_rate[["Indonesia"]]
birth_rate <-read.csv("data/birth_rate.csv")
birth_rate <- birth_rate[["Indonesia"]]

# incidence rates per 100.000 person-years of observation  
# <2 yr, 2-4 yr, 5-9 yr, 10-19 yr, 20-29 yr, 30-39 yr, 40-49 yr, 50-65 yr,
inc_rate_obs <- as.integer( 
  c( 67.20430108, 225.660864, 499.197718, 301.1815584, 99.59332725, 40.67934506, 42.89390906, 11.83712121, 11.83712121 ) )


