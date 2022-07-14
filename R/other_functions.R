# Normalize BART function (Same way as theOdds code)
normalize_bart <- function(y) {

  # Defining the a and b
  a <- min(y)
  b <- max(y)

  # This will normalize y between -0.5 and 0.5
  y  <- (y - a)/(b - a) - 0.5
  return(y)
}

# Now a function to return everything back to the normal scale

unnormalize_bart <- function(z, a, b) {
  # Just getting back to the regular BART
  y <- (b - a) * (z + 0.5) + a
  return(y)
}
# Calculating RMSE
#' @export
rmse <- function(obs, pred) {
  return(sqrt(mean((obs - pred)^2)))
}

rMVN_var <- function(mean, Sigma) {
  if(length(mean) == 1){
    mean <- rep(mean,nrow(Sigma))
  }
  if(is.matrix(Sigma)) {
    drop(mean + crossprod(PD_chol(Sigma), stats::rnorm(length(mean))))
  } else {
    mean + sqrt(Sigma) * stats::rnorm(length(mean))
  }
}

dh_cauchy <- function(x,location,sigma){

  if(x>=location){
    return ((2/(pi*sigma))*(1/(1+((x-location)^2)/(sigma^2))))
  } else
    return(0)
}

is_diag_matrix <- function(m) all(m[!diag(nrow(m))] == 0)

PD_chol  <- function(x, ...) tryCatch(chol(x, ...), error=function(e) {
  d    <- nrow(x)
  eigs <- eigen(x, symmetric = TRUE)
  eval <- eigs$values
  evec <- eigs$vectors
  return(chol(x + evec %*% tcrossprod(diag(pmax.int(1e-8, 2 * max(abs(eval)) * d * .Machine$double.eps - eval), d), evec), ...))
}
)

# Check the appendix of Linero SoftBART for more details
update_tau_linero <- function(x_train,
                              y,
                              y_hat,
                              curr_tau){
  # Getting number of observations
  n <- length(y)
  # Calculating current sigma
  curr_sigma <- curr_tau^(-1/2)

  sigma_naive <- naive_sigma(x = x_train,y = y)

  proposal_tau <- stats::rgamma(n = 1,shape = 0.5*n+1,rate = 0.5*crossprod( (y-y_hat) ))

  proposal_sigma <- proposal_tau^(-1/2)

  acceptance <- exp(log(dh_cauchy(x = proposal_sigma,location = 0,sigma = sigma_naive)) +
                      3*log(proposal_sigma) -
                      log(dh_cauchy(x = curr_sigma,location = 0,sigma = sigma_naive)) -
                      3*log(curr_sigma))

  if(stats::runif(n = 1)<acceptance){
    return(proposal_sigma^(-2))
  } else {
    return(curr_tau)
  }

}

# Naive sigma_estimation
naive_sigma <- function(x,y){

  # Getting the valus from n and p
  n <- length(y)

  # Getting the value from p
  p <- ifelse(is.null(ncol(x)), 1, ncol(x))

  # Adjusting the df
  df <- data.frame(x,y)
  colnames(df)<- c(colnames(x),"y")

  # Naive lm_mod
  lm_mod <- stats::lm(formula = y ~ ., data =  df)

  # Getting sigma
  sigma <- stats::sigma(lm_mod)

  return(sigma)
}
