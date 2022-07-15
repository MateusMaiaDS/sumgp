# Main function to do a sample of GP
#' @export
sumgp <- function(x_train,
                  y,
                  x_test,
                  n_gp, kappa,
                  nu, tau,
                  tau_mu,
                  phi,
                  mu, n_mcmc, n_burn){


  # Scaling the data ==== (NO MORE SCALING AT THE FIRST MOMENT) ===
  # xscale_train <- scale(x_train)
  # x_mean <- attr(xscale_train,"scaled:center")
  # x_sd <- attr(xscale_train, "scaled:scale")
  #
  # xscale_test <- scale(x_test, center = x_mean,scale = x_sd)
  #
  # y_scale <- normalize_bart(y)


  # Function with no scale;
  xscale_train <- x_train
  xscale_test <- x_test
  y_scale <- y

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

    # Storing the posterior
    if(i > n_burn){
      curr_iter = curr_iter +1
      y_train_hat[curr_iter,] <- colSums(pred_train) #unnormalize_bart(colSums(pred_train),a = min(y),b = max(y))
      y_test_hat[curr_iter,] <- colSums(pred_test) #unnormalize_bart(colSums(pred_test),a = min(y),b = max(y))
      tau_post[curr_iter] <- tau #/((max(y)-min(y))^2)
    }

    for(j in 1:n_gp){

      res <- y_scale - colSums(pred_train[-j,,drop = FALSE])

      # Calculating the sampled value for mu
      mu <- update_mu(x_train = xscale_train,res = res,
                      kappa = kappa, tau = tau,
                      tau_mu = tau_mu, nu = nu , phi = phi)

      psi_update <- update_psi(x_train = xscale_train,
                               res  = res, x_test = xscale_test,
                               mu = mu,kappa = kappa,tau = tau,
                               tau_mu = tau_mu,nu = nu,phi = phi)
      # Updating the predictions
      pred_train[j,] <- psi_update$sample_train
      pred_test[j,] <- psi_update$sample_test
    }

    tau_old <- tau
    tau <- update_tau_linero(x_train = xscale_train,y = y_scale,
                             y_hat = colSums(pred_train),curr_tau = tau_old)
  print(i)
  }

  main_result <- list(y_train_hat = y_train_hat,
      y_test_hat = y_test_hat,
      tau_post = tau_post)


  # Returning all the posteriors
  return(main_result)
}


# Updating mu
update_mu <- function(x_train,res,kappa,tau,tau_mu,nu,phi){

  ones <- matrix(1,nrow = nrow(x_train))
  omega <- kernel_function(squared_distance_matrix = symm_distance_matrix(m1 = x_train),
                           nu = nu, phi = phi)
  v_obj <- crossprod(ones, solve((diag(tau^-1,nrow = nrow(x_train)) + (1-kappa)*omega),ones)) + tau_mu/kappa

  mu_mean <- (v_obj^(-1))*crossprod(res,solve((diag(tau^-1,nrow = nrow(x_train)) + (1-kappa)*omega),ones))
  mu_sd <- v_obj^(-1/2)

  return(stats::rnorm(n = 1,mean = mu_mean,sd = mu_sd))
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
  # sample_train <- rMVN_var(mean = mu_gp_train,Sigma = sigma_gp_train)

  # getting the train sample using a package
  sample_train <- mvnfast::rmvn(n = 1,mu = mu_gp_train,
                                sigma = sigma_gp_train+diag(1e-9,nrow = length(mu_gp_train)))

  # Sampling now the test observations
  mu_gp_test <- mu_test_vec + (1-kappa)*crossprod(omega_train_test,
                                                solve((diag(tau^(-1),nrow = nrow(x_train))+(1-kappa)*omega_train),(res-mu_train_vec)))

  sigma_gp_test <- (1-kappa)*omega_test-((1-kappa)^(2))*crossprod(omega_train_test,
                                                                   solve((diag(tau^(-1),nrow = nrow(x_train))+(1-kappa)*omega_train),omega_train_test))

  sample_test <- mvnfast::rmvn(n = 1,mu = mu_gp_test,
                               sigma = sigma_gp_test+diag(1e-9,nrow = length(mu_gp_test)))


  # Our older way to calculate the sample
  # sample_test <- rMVN_var(mean = mu_gp_test,Sigma = sigma_gp_test)

  if(!identical(mu_gp_test,mu_gp_train)){
   print("THERE'S SOMETHING WRONG (MU)")
  }

  if(!identical(sigma_gp_test,sigma_gp_train)){
    print("THERE'S SOMETHING WRONG (SIGMA)")
  }


  # Returning the list of those two values
  return(list(sample_train = sample_train,
              sample_test = sample_test))
}

# Calculating a PI coverage
#' @export
pi_coverage <- function(y, y_hat_post, sd_post,only_post = FALSE, prob = 0.5,n_mcmc_replications = 1000){

  # Getting the number of posterior samples and columns, respect.
  np <- nrow(y_hat_post)
  nobs <- ncol(y_hat_post)

  full_post_draw <- list()

  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_mcmc_replications,
    style = 3, width = 50 )

  # Only post matrix
  if(only_post){
    post_draw <- y_hat_post
  } else {
    for(i in 1:n_mcmc_replications){
      utils::setTxtProgressBar(progress_bar, i)

      full_post_draw[[i]] <-(y_hat_post + replicate(sd_post,n = nobs)*matrix(stats::rnorm(n = np*nobs),
                                                                             nrow = np))
    }
  }

  if(!only_post){
    post_draw<- do.call(rbind,full_post_draw)
  }

  # CI boundaries
  low_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = prob/2)})
  up_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = 1-prob/2)})

  pi_cov <- sum((y<=up_ci) & (y>=low_ci))/length(y)

  return(pi_cov)
}
