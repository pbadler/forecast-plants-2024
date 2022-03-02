rm(list = ls())
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

age_fun <- function(x, a, b, c, min_size = 0.25){ 
  a*exp(-(x-b)^2/(2*c^2)) + min_size
}


re <- T
P <- 100
a <- 2
b <- 6 
c <- 1.5
d <- 7
min_size <- 0.25
sigma <- 0.5 
a_sig <- 0.5
b_sig <- 0.5

if(re){ 
  #d <- rnbinom( P, mu = d, size = 1)
  b <- rpois(P, b)  + 1
  a <- rlnorm(P, log(a), a_sig) 
}


plants <- data.frame(p = 1:P, 
                     death = d, 
                     b = b, 
                     a = a, 
                     c = c, 
                     max_sizes = a)

res <- list(NA)
i <- 1

for( i in 1:nrow(plants)){ 
  
  b <- plants$b[i]
  death <- plants$death[i] + 1
  a <- plants$a[i]
  c <- plants$c[i]
  
  mu_size <- age_fun(1:500, a, b, c)[1:death]

  size <- rnorm(death, log( mu_size ) , sigma)
   
  # biggest <- max(size)  
  # plot(1:death, size , main = i, ylim = c(0, max(biggest, a)))
  # curve( age_fun(x, a, b, c), 0, death, add = T)

  res[[i]] <- data.frame( p = rep( i, death ), age = 1:death, size = size, a = a, b = b)
  
} 

dat <- do.call(rbind, res )

plot(dat$age, exp(dat$size) )
curve( age_fun(x, mean(dat$a), mean(dat$b), c), 0, death, add = T)

library(tidyverse)
show <- 20

curve_fits <- 
  expand.grid( age = seq(1, 10, by = 0.5), p = unique(dat$p)) %>% 
  left_join(dat %>% distinct(p, a, b), by = 'p') %>% 
  mutate( size = age_fun(age, a = a, b = b, c = c))

dat %>% 
  filter( p < show) %>% 
  ggplot( aes( x = age, y = exp(size) )) + 
  geom_point() + 
  geom_line(data = curve_fits %>% filter( p < show), aes( y = size )) + 
  facet_wrap(~p)


# fit stan models 
Y <- dat$size 
age <- dat$age
p <- dat$p
P <- length(unique(p))
min_size <- min_size
N <- length(Y)

sdat <- list(N = N, P = P, Y = Y, age = age, p = p, min_size = min_size )

# fit <- stan('code/analysis/growth_curve.stan', data = sdat)
# summary(fit, c('a', 'b', 'c', 'sigma'))$summary

library(brms)

hist( sdat$Y )

fit_arma <- brm( Y ~ 1, autocor = cor_arma( ~  age | p , 1, 1), data = sdat)


prior1 <- prior(normal(0, 5), nlpar = 'a', lb = 0) + 
  prior(normal(0, 5), nlpar = 'b', lb = 0) + 
  prior(gamma(2,1), nlpar = 'c', lb = 0)

fit <- 
  brm(bf(Y ~ a*exp(-(age-b)^2/(2*c^2) + min_size), 
        c ~ 1, a + b ~ (1|p), nl = TRUE), 
    data = sdat, prior = prior1, family = 'student', control = list( adapt_delta = 0.99, max_treedepth = 20))

make_stancode(bf(Y ~ log(a*exp(-(age-b)^2/(2*c^2)) + min_size), 
   c ~ 1,  a + b ~ (1|p), nl = TRUE), data = sdat , prior = prior1)

fit
summary( fit )
a_int <- data.frame( ranef(fit))
a_int

plot( a_int$p.Estimate.a_Intercept,  distinct(dat, a, p)$a )
plot( a_int$p.Est.Error.b_Intercept, distinct(dat, b, p)$b )
pp_check(fit)

yhat <- predict( fit )

plot( dat$size,  yhat[, 1] )
abline(0,1)
