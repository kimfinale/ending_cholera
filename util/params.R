# Demographics
death_rate <-read.csv("data/death_rate.csv")
death_rate <- death_rate[["Indonesia"]]
birth_rate <-read.csv("data/birth_rate.csv")
birth_rate <- birth_rate[["Indonesia"]]

nag <- num_age_grp <- 10L # number of age groups
index_S <- (1:nag) + 1 # first index (column) will be time
index_I <- index_S + nag
index_R <- index_S + 2*nag
index_V <- index_S + 3*nag
index_CI <- index_S + 4*nag
index_CV <- index_S + 5*nag

# Transmission-related parameters
beta <- 1.2  # transmission rate will be calibrated  
gamma <- 1/5 # 5 days for the duration of infectiousness Weil et al. (2009) Clin Infect Dis 
omega <- 1/(3*365) # 3 years for the duration of natural immunity Ali et al. (2011) J Infect Dis
ag <- c( 1/365, 1/(4*365), rep( 1/(10*365), 8) ) # <1, 1-4, 5-14, 15-24, ... 65-74, 75-84
rate_excretion <- 1;
rate_decay <- 1/20;
case_fatality <- rep( 0, nag )
rel_susc <- rep( 1, nag )
frac_report <- 1

# Vaccination
vacc_campaign <- FALSE
vacc_routine <- FALSE
start_vacc_campaign <- 61*365
dur_vacc_campaign <- 14.0
stop_vacc_campaign <- start_vacc_campaign + dur_vacc_campaign

cov_vacc_campaign <- rep( 0, num_age_grp )
cov_vacc_routine <- rep( 0, num_age_grp )

eff_vacc <- rep( 0, nag )
eff_vacc[1:2] <- 0.3 #Bi et al.(2017) Lancet Infect Dis
eff_vacc[3:nag] <- 0.64 #Bi et al.(2017) Lancet Infect Dis

rate_vacc_campaign <- - log(1-cov_vacc_campaign*eff_vacc) / dur_vacc_campaign


init_val <- c( 0.019129528368344, 0.0674741649748425, 0.126429307827068, 0.0936790567827395, 
                0.0689075046855928, 0.0497120079378464, 0.0339845165599162, 0.0202690060767594, 
                0.00908927231705442, 0.0024549202470792, 1.17987107841576e-05, 
                4.24569334984788e-05, 7.95619657527049e-05, 5.89522483326505e-05, 
                4.33635058955322e-05, 3.12837761159886e-05, 2.1386462779773e-05, 
                1.27552894059669e-05, 5.71988080987387e-06, 1.54488178164627e-06, 
                0.000594187298716229, 0.0109978614296745, 0.0673712858210493, 
                0.0949157764370814, 0.112867457328733, 0.121100057079751, 0.11573977692341, 
                0.0902496016197851, 0.0488325131992378, 0.0146789226129884, rep(0, 30), 1e-3 )

# age_dist <- c(0.0197355143778444, 0.0785144833380155, 0.19388015561387, 0.188653785468154, 
#               0.181818325520221, 0.170843348793713, 0.149745679946106, 0.11053136298595, 
#               0.0579275053971021, 0.0171353877418492)
# # rates are given per year


# incidence rates per 100.000 person-years of observation  
# <2 yo, 2-4 yo, 5-9 yo, 10-19 yo, 20-29 yo, 30-39 yo, 40-49 yo, 50-65 yo,
# inc_obs <- as.integer( 
#   c( 67.20430108, 225.660864, 499.197718, 301.1815584, 99.59332725, 40.67934506, 42.89390906, 11.83712121, 11.83712121, 11.83712121 ) )

inc_obs_Jakarta <- as.integer( c( 4.01, 1.55, 0.29, 0.27 )*100 )# <1 yo, 1-4 yo, 5-14 yo, 15+ yo per 10^5 person-years
inc_obs <- inc_obs_Jakarta
nag <- num_age_grp
iag <- index_age_all_states <- lapply( 1:10, function(x) c(x, nag+x, 2*nag+x, 3*nag+x) )
age_dist <- sapply( iag, function(x) sum(init_val[x]))/sum(init_val)

