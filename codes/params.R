# Demographics
library(readr)
d <- read_csv("data/cholera_data.csv")
names(d) <- c( "country",	"WHO_stratum", "case_fatality_rate",
                 "inc_total",	"inc_<1",	"inc_1-4","inc_5-14", "inc_15+", "crude_birth_rate_2005-2010", 
                 "crude_birth_rate_2010-2015", "crude_birth_rate_2015-2020", "ref_year",
                 "prop_<1",	"prop_1-4",	"prop_5-14",	"prop_15-24",	"prop_25-34",	"prop_35-44",	"prop_45-54",	"prop_55-64",	"prop_65-74", 
                 "prop_75+",	"prop_tot" )

dat <- d[d$country==country,]

# if( is.na(dat$`inc_<1`) ){
#   dat$`inc_<1` <- dat$`inc_1-4`  ## rewrite the code by using the total incidene 
# }
birth_rate <- (dat$`crude_birth_rate_2005-2010` + dat$`crude_birth_rate_2010-2015`)/2/1000/365# 21.8/1000/365 # per person per day 
age_dist <- as.double(dat[,13:22])
ag <- c( 1/365, 1/(4*365), rep( 1/(10*365), 7), 0 ) # rate of aging <1, 1-4, 5-14, 15-24, ... 65-74, 75+
death_rate <- rep( 0, length(ag)) # death_rate is adjusted according to the age distribution 
death_rate[1] <- (birth_rate - age_dist[1]*ag[1])/ age_dist[1]
for( i in 2:10 ){
  death_rate[i] <- (age_dist[i-1]*ag[i-1] - age_dist[i]*ag[i])/ age_dist[i]
}
inc_obs <- as.integer( c( dat$`inc_<1`, dat$`inc_1-4`, dat$`inc_5-14`, dat$`inc_15+` )*100 ) # incidence per 1e5 pyo

# Transmission-related parameters
beta <- 1.2  # transmission rate will be calibrated  
gamma <- 1/5 # 5 days for the duration of infectiousness Weil et al. (2009) Clin Infect Dis 
omega <- 1/(3*365) # 3 years for the duration of natural immunity Ali et al. (2011) J Infect Dis
rate_excretion <- 1
rate_decay <- 1/20
nag <- num_age_grp <- 10 # number of age groups
case_fatality <- rep( 0.03, nag ) # SEAR-D Ali et al. (2015)
rel_susc <- rep( 1, nag )
frac_symptom <- 0.25 # one fourth of the cholera-infected people develop symptoms
frac_report <- 0.2# one fifth of the people with symptoms are reported to health facilities

# Vaccination
vacc_campaign <- FALSE
vacc_routine <- FALSE
start_vacc_campaign <- 61*365
dur_vacc_campaign <- 30
stop_vacc_campaign <- start_vacc_campaign + dur_vacc_campaign

cov_vacc_campaign <- rep( 0, num_age_grp )
cov_vacc_routine <- rep( 0, num_age_grp )
rate_vacc_campaign_susc <- rep( 0, num_age_grp )
rate_vacc_campaign_removed <- rep( 0, num_age_grp ) 
eff_vacc <- rep( 0, nag )
eff_vacc[1:2] <- 0.3 #Bi et al.(2017) Lancet Infect Dis
eff_vacc[3:nag] <- 0.64 #Bi et al.(2017) Lancet Infect Dis

rate_wane_vacc <- 1/(5*365)
rate_wane_nat <- 1/(5*365)

rate_vacc_campaign <- - log(1-cov_vacc_campaign) / dur_vacc_campaign

init_val <- c( age_dist[1:10]/2, 
               age_dist[1:10]/4, 
               age_dist[1:10]/4, 
               rep(0, 30), 1e-3 )

index_s <- 1:nag
index_i <- index_s + nag
index_r <- index_s + 2*nag
index_v <- index_s + 3*nag
index_ci <- index_s + 4*nag
index_cv <- index_s + 5*nag
index_b <- 6*nag + 1

index_states <- lapply( 1:6, function(x) 10*(x-1) + (1:nag) )
iag <- index_age_groups <- lapply( 1:10, function(x) c(x, nag+x, 2*nag+x, 3*nag+x) )

