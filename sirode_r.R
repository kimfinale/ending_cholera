
sir_r <- function( t, y, params ){
  S <- y[ 1 ]
  I <- y[ 2 ]
  R <- y[ 3 ]
  N <- S + I + R
  beta <- params[1]
  gamma <- params[2]
  
  dS <- - beta * S * I / N
  dI <- beta * S * I / N - gamma * I
  dR <- gamma * I
  
  return( list( c( dS, dI, dR ) ) )
  
}

sird_r <- function( t, y, params ){
  S <- y[ 1 ]
  I <- y[ 2 ]
  R <- y[ 3 ]
  N <- S + I + R
  beta <- params[1]
  gamma <- params[2]
  mu <- 1/(60*365)
  dS <- - beta * S * I / N + mu*(N-S) 
  dI <- beta * S * I / N - gamma * I - mu*I
  dR <- gamma * I - mu*R
  
  return( list( c( dS, dI, dR ) ) )
  
}