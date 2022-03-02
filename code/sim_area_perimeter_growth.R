rm(list = ls())

gain <- function(P, R) { 
  # P is perimeter in cm 
  # R is Resource per cm
  # w is conversion rate from gain to area (area per mass)
  
  P*w*R
}


loss <- function( A, L ) { 
  # A is area cm2 
  # L is a loss rate per area
  
  A*L
}


calc_P <- function( A ) { 
  sqrt(A/pi) 
}

years <- 500 
n <- 1
R <- rexp(years, 1/5)
L <- rbeta(years, 1, 2)
w <- 0.2 
A <- 0.25

for ( i in 2:years ) { 
  
  A[i] <- A[i-1] + gain( calc_P(A[i-1]), R[i-1]) - loss(A[i-1], L[i-1]) 
  
}

plot(A, type = 'l')
B <- rnorm( years, A, 1.5)
points(B, type = 'l', col = 'blue')
B[ B <= 0 ] <- min(A)

dat <- data.frame( A = log(A), B = log(B))
dat$A0 <- lag(A)
dat$B0 <- lag(B)
dat$R <- R 
dat$R0 <- lag(R)
dat$R1 <- lag(dat$R0)
dat$S <- rnorm(R)
dat$S0 <- lag(dat$S)
dat$S1 <- lag(dat$S0)

m1 <- lm( data = dat, B ~ B0 + R0 + R1 )
summary( m1)

m2 <- lm( data = dat, A ~ A0 + R0 + R1 )
summary( m2)

m3 <- lm( data = dat, B ~ B0 + S0 + S1 )
summary( m3)

