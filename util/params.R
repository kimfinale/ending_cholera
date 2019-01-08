num_patch <- 4

contact_matrix <- matrix( 0, nrow=num_patch, ncol=num_patch )
contact_matrix[] <- 0.1
diag( contact_matrix ) <- 0.7

R0 <- rep( NA, num_patch )
R0[ 1 ] <- 3
R0[ 2:num_patch ] <- 0.8

gamma <- 1/14 # 14 days of infectious period
beta <- R0 * gamma
rho <- 0 # vaccination rate
ag <- rep( 0, 4 )
ag[ 1 ] <- 1/365
ag[ 2 ] <- 1/(4*365)
ag[ 3 ] <- 1/(10*365)
ag[ 4 ] <- 1/(55*365)

generate_contact_matrix <- function( num_patch=2, prop_own = 0.5, symmetric=TRUE ){
  contact_matrix <- matrix( 0, nrow=num_patch, ncol=num_patch )
  contact_matrix[] <- (1-prop_own)/(num_patch-1)
  diag( contact_matrix ) <- prop_own
}