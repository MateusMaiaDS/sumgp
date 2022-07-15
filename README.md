# sumgp
A package to sample the sum of GP

## To install the package
```
devtools::install_github("MateusMaiaDS/gpbart")
```

### Code to run the experiments

```r

# Loading packages
library(tidyverse)
library(sumgp)
rm(list=ls())

# Creating a simple simulated model
n_train <-50
n_test <- 200
n_mcmc <- 1500
n_burn <- 500
n_gp <- 10

# Generating the data
x_train <- seq(-pi,pi,length.out = n_train) %>% as.matrix
# x_test <- seq(-pi,pi,length.out = n_test) %>% as.matrix
x_test <- x_train
y_mean <- sin(x_train)
y <- y_mean + rnorm(n = n_train,sd = 0.1)

# Setting parameters
nu <- tau_mu <- 4*4*n_gp/((max(y)-min(y))^2)
tau <- 10
mu <- 0
kappa <- 0.5
phi <- 1


# Running the sum of GP's
main_result <- sumgp(x_train = x_train,
                     y = y,x_test = x_test,
                     n_gp = n_gp,kappa = kappa,phi = phi,
                     nu = nu, tau = tau, mu = mu ,tau_mu = tau_mu,
                     n_mcmc = n_mcmc, n_burn = n_burn)
# Plotting the data
pi_coverage(y = y,y_hat_post = main_result$y_train_hat,
            sd_post = main_result$tau_post^(-1/2),
            prob = 0.5,n_mcmc_replications = 200)

pi_coverage(y = y,y_hat_post = main_result$y_test_hat,
            sd_post = main_result$tau_post^(-1/2),
            prob = 0.5,n_mcmc_replications = 200)

main_result$y_train_hat %>% apply(2,sd) %>% hist(main = "train-sd")
main_result$y_test_hat %>% apply(2,sd) %>% hist(main = "test-sd")

plot(x_test,y)
points(x_train,main_result$y_train_hat %>% colMeans(), pch= 20, col = "blue")
points(x_test,main_result$y_test_hat %>% colMeans(),  pch= 20, col = "orange")


```
