# Main function to do a sample of GP
sumgp <- function(x_train,
                  y,
                  x_test,
                  n_gp, kappa,
                  nu, tau, mu, n_mcmc, n_burn){


  # Scaling the data
  xscale_train <- scale(x_train)
  x_mean <- attr(xscale_train,"scaled:center")
  x_sd <- attr(xscale_train, "scaled:scale")

  x_test <- scale(x_test, center = x_mean,scale = x_sd)

  y_scale <- normalize_bart(y)

  # Creating the predictions train and test matrix
  y_train_hat <- matrix(0, nrow = (n_mcmc-n_burn), ncol = nrow(x_train))
  y_test_hat <- matrix(0, nrow = (n_mcmc - n_burn), ncol = nrow(x_test))

  # Getting pred_train and pred_test
  pred_train <- matrix(0, nrow = n_gp, ncol = nrow(x_train))
  pred_test <- matrix(0, nrow = n_gp, ncol = nrow(x_test))

  tau_post <- numeric()

  # Curr iter
  curr_iter <- 0
  residuals <- numeric()
  # Iterating over all MCMC samples
  for(i in 1:n_mcmc){

    if(i > n_burn){
      curr_iter = curr_iter +1
      y_train_hat[curr_iter,] <- colSums(pred_train)
      y_test_hat[curr_iter,] <- colSums(pred_test)
      tau_post[curr_iter] <- tau
    }

    for(j in 1:n_gp){

      residuals <- y_scale - colSums(pred_train[-j,,drop = FALSE])

      # Calculating the sampled value for mu
      mu <- update_mu(x_train = x_train,res = res,
                      kappa = kappa, tau = tau,
                      tau_mu = tau_mu, nu = nu , phi = phi)

      psi_update <- update_psi(x_train = x_train,
                               res  = res, x_test = x_test,
                               mu = mu,kappa = kappa,tau = tau,
                               tau_mu = tau_mu,nu = nu,phi = phi)
      # Updating the predictions
      pred_train[j,] <- psi_update$sample_train
      pred_test[j,] <- psi_update$sample_test
    }

    tau_old <- tau
    tau <- update_tau_linero(x_train = x_train,y = y,
                             y_hat = colSums(pred_train),curr_tau = tau_old)
  print(i)
  }


  # Returning all the posteriors
  return(list(y_train_hat = y_train_hat,
              y_test_hat = y_test_hat,
              tau_post = tau_post))
}


# Updating mu
update_mu <- function(x_train,res,kappa,tau,tau_mu,nu,phi){

  ones <- matrix(1,nrow = nrow(x_train))
  omega <- kernel_function(squared_distance_matrix = symm_distance_matrix(m1 = x_train),
                           nu = nu, phi = phi)
  v_obj <- crossprod(ones, solve((diag(tau^-1,nrow = nrow(x_train)) + (1-kappa)*omega),ones)) + tau_mu/kappa

  mu_mean <- (v_obj^(-1))*crossprod(res,solve((diag(tau^-1,nrow = nrow(x_train)) + (1-kappa)*omega),ones))
  mu_sd <- v_obj^(-1/2)

  return(rnorm(n = 1,mean = mu_mean,sd = mu_sd))
}


# Update psi

update_psi <- function(x_train,res,x_test,
                       mu, kappa, tau, tau_mu, nu, phi){

  # Getting the mu vec
  mu_train_vec <- matrix(rep(mu,nrow(x_train)))
  mu_test_vec <- matrix(rep(mu,nrow(x_test)))


  omega_train <- kernel_function(squared_distance_matrix = symm_distance_matrix(x_train),
                                nu = nu,phi = phi)

  omega_train_test <- kernel_function(squared_distance_matrix = distance_matrix(x_train,x_test),
                                      nu = nu,phi = phi)

  omega_test <- kernel_function(squared_distance_matrix = symm_distance_matrix(x_test),
                                nu = nu,phi = phi)


  # Calculating mean and variance vector and matrix
  mu_gp_train <- mu_train_vec + (1-kappa)*crossprod(omega_train,
                                                    solve((diag(tau^(-1),nrow = nrow(x_train))+(1-kappa)*omega_train),(res-mu_train_vec)))

  sigma_gp_train <- (1-kappa)*omega_train-((1-kappa)^(2))*crossprod(omega_train,
                                                                  solve((diag(tau^(-1),nrow = nrow(x_train))+(1-kappa)*omega_train),omega_train))

  # Getting the sample train
  sample_train <- rMVN_var(mean = mu_gp_train,Sigma = sigma_gp_train)

  # Sampling now the test observations
  mu_gp_test <- mu_test_vec + (1-kappa)*crossprod(omega_train_test,
                                                solve((diag(tau^(-1),nrow = nrow(x_train))+(1-kappa)*omega_train),(res-mu_train_vec)))

  sigma_gp_test <- (1-kappa)*omega_test-((1-kappa)^(2))*crossprod(omega_train_test,
                                                                   solve((diag(tau^(-1),nrow = nrow(x_train))+(1-kappa)*omega_train),omega_train_test))


  sample_test <- rMVN_var(mean = mu_gp_test,Sigma = sigma_gp_test)


  # Returning the list of those two values
  return(list(sample_train = sample_train,
              sample_test = sample_test))
}
