rm(list = ls())
library(brms)
library(future)

start <- Sys.time()

fit0 <- brm(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
            data = epilepsy, family = poisson())

t0 <- Sys.time() - start 

start <- Sys.time()

fit1 <- 
  brm(count ~ zAge + zBase + Trt + (1|patient) + (1|obs),
            data = epilepsy, family = poisson(), future = T)

t1 <- Sys.time() - start 

print(t0)
print(t1)

plan(multiprocess) # run kfold with future 
start <- Sys.time()
k1 <- kfold(fit0)
t2 <- Sys.time() - start 

start <- Sys.time()
k2 <- kfold(fit1)
t3 <- Sys.time() - start 

k1
k2 

print(t2)
print(t3)
