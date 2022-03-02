rm(list = ls())

library(tidyverse)
library(lme4)
t <- 30
a <- 0.2
s <- 100 
cutoff <- 0
R <- rnorm(t)
beta <- 0.1
m <- rnorm( s, 0.1, 0.2)
t_start <- ceiling( runif(s, 1, t-2) )
t_end <- ceiling( runif( s, t_start, t))
b <- runif( s, 0.2, 0.8)
X <- matrix( NA, s, t )

for( i in 1:s){ 
  X[ i , t_start[i] ] <- m[i]
}

for ( i in 2:t ) { 
  for( j in 1:s ) {   
    temp <- rnorm(1,  a + b[j]*X[j,i-1] + beta*R[i-1], 0.1)
    temp <- ifelse( i > t_end[j] , NA, temp )
    X[j, i] <- ifelse( !is.na(temp), temp, X[j, i])
  }
}

X[X < cutoff] <- cutoff

par(mfrow = c(2, 1))
plot( X[1, ], type = 'l', col = 1, ylim = c(0, 3))

for ( i in 2:s) { 
  points( X[i, ], type = 'l' , col = i )
}

plot( R, type = 'l', lty = 2)
par(mfrow = c(1,1))

dat <- 
  data.frame(X = t(X) ) %>% 
  gather( plant, size, starts_with('X')) %>% 
  mutate( plant = as.numeric(factor( plant ))) %>% 
  group_by( plant) %>% 
  mutate( size0 = lag(size)) %>% 
  mutate( t = row_number())

dat$R0 <- lag(R)
dat$RL <- lag(R, 2)

m1 <- lm( data = dat, size ~ size0 + R0 + RL)
summary(m1)

mer <- lmer( data = dat, size ~ size0 + R0 + RL + (1|t) + (1|plant) )

summary( mer )


