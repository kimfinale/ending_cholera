# Demographics
library(readr)
dat <- read_csv("data/cholera_data.csv")

names(dat) <- c( "country",	"WHO_stratum", "case_fatality_rate",
                 "inc_total",	"inc_<1",	"inc_1-4","inc_5-14", "inc_15+", "crude_birth_rate_2005-2010", 
                 "crude_birth_rate_2010-2015", "crude_birth_rate_2015-2020", "ref_year",
                 "prop_<1",	"prop_1-4",	"prop_5-14",	"prop_15-24",	"prop_25-34",	"prop_35-44",	"prop_45-54",	"prop_55-64",	"prop_65-74", 
                 "prop_75+",	"prop_tot" )
dat <- dat[-1,] 
# death_rate <- death_rate[["Indonesia"]]
# birth_rate <- read.csv("data/birth_rate.csv")
# birth_rate <- birth_rate[["Indonesia"]]
d <- dat[dat$country==country,]

birth_rate <- (d$`crude_birth_rate_2005-2010` + d$`crude_birth_rate_2010-2015`)/2/1000/365# 21.8/1000/365 # per person per day 
# age_dist <- c( 0.02238036, 0.089521439, 0.2510, 0.1951, 0.1446, 0.1108, 0.0798, 0.0568, 0.0356, 0.0145 ) #Nepal
age_dist <- as.double(d[,13:22])
ag <- c( 1/365, 1/(4*365), rep( 1/(10*365), 7), 0 ) # rate of aging <1, 1-4, 5-14, 15-24, ... 65-74, 75+
death_rate <- rep( 0, length(ag)) # death_rate is adjusted according to  
death_rate[1] <- (birth_rate - age_dist[1]*ag[1])/ age_dist[1]
for( i in 2:10 ){
  death_rate[i] <- (age_dist[i-1]*ag[i-1] - age_dist[i]*ag[i])/ age_dist[i]
}
# Transmission-related parameters
beta <- 1.2  # transmission rate will be calibrated  
gamma <- 1/5 # 5 days for the duration of infectiousness Weil et al. (2009) Clin Infect Dis 
omega <- 1/(3*365) # 3 years for the duration of natural immunity Ali et al. (2011) J Infect Dis
rate_excretion <- 1;
rate_decay <- 1/20;
nag <- num_age_grp <- 10 # number of age groups
case_fatality <- rep( 0.03, nag ) # SEAR-D Ali et al. (2015)
rel_susc <- rep( 1, nag )
frac_symptom <- 0.25 # one fourth of the cholera-infected people develop symptoms
frac_report <- 0.2# one fifth of the people with symptoms are reported to health facilities

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


# init_val <- c( 0.019129528368344, 0.0674741649748425, 0.126429307827068, 0.0936790567827395, 
#                 0.0689075046855928, 0.0497120079378464, 0.0339845165599162, 0.0202690060767594, 
#                 0.00908927231705442, 0.0024549202470792, 1.17987107841576e-05, 
#                 4.24569334984788e-05, 7.95619657527049e-05, 5.89522483326505e-05, 
#                 4.33635058955322e-05, 3.12837761159886e-05, 2.1386462779773e-05, 
#                 1.27552894059669e-05, 5.71988080987387e-06, 1.54488178164627e-06, 
#                 0.000594187298716229, 0.0109978614296745, 0.0673712858210493, 
#                 0.0949157764370814, 0.112867457328733, 0.121100057079751, 0.11573977692341, 
#                 0.0902496016197851, 0.0488325131992378, 0.0146789226129884, rep(0, 30), 1e-3 )

init_val <- c( age_dist[1:10]/2, 
               age_dist[1:10]/4, 
               age_dist[1:10]/4, 
               rep(0, 30), 1e-3 )


# # rates are given per year
inc_obs_Jakarta <- as.integer( c( 4.01, 1.55, 0.29, 0.27 )*100 )# <1 yo, 1-4 yo, 5-14 yo, 15+ yo per 10^5 person-years

inc_obs <- as.integer( c( 7.16, 7.01, 2.19, 0.93 )*100 )

# incidence rates per 100.000 person-years of observation  
# <2 yo, 2-4 yo, 5-9 yo, 10-19 yo, 20-29 yo, 30-39 yo, 40-49 yo, 50-65 yo,
# inc_obs <- as.integer( 
#   c( 67.20430108, 225.660864, 499.197718, 301.1815584, 99.59332725, 40.67934506, 42.89390906, 11.83712121, 11.83712121, 11.83712121 ) )

index_s <- 1:nag
index_i <- index_s + nag
index_r <- index_s + 2*nag
index_v <- index_s + 3*nag
index_ci <- index_s + 4*nag
index_cv <- index_s + 5*nag
index_b <- 6*nag + 1

index_states <- lapply( 1:6, function(x) 10*(x-1) + (1:nag) )
iag <- index_age_groups <- lapply( 1:10, function(x) c(x, nag+x, 2*nag+x, 3*nag+x) )
# age_dist <- sapply( iag, function(x) sum(init_val[x]))/sum(init_val)

