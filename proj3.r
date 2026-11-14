

## function definitions
## poission_lasso_objective: penalised cost function of the Poisson changepoint
## model
## inputs: beta, the parameter vector of length n, where beta_1 is the
##         intercept value, with beta_k (k>=2) denoting the change in the
##         log-mean function relative to the last time step
##         y, the data vector
##         lambda, the penalty parameter, set a default value of 1
##         eps, the smoothing parameter for approximating the absolute
##         function, set a default value of 1e-6.
## outputs: the value of the penalised cost
## work: compute the mean function from beta, then construct the negative 
##       log-likelihood function, add lambda times the sum of 
##       sqrt(beta_k^2 + eps) over k = 2,...,n, which shrinks beta_k towards 0

poission_lasso_objective <- function(beta, y, lambda = 1, eps = 1e-6){
  eta <- cumsum(beta)  ## linear predictor, eta_k = beta_1 + ... + beta_k
  mu <- exp(eta)  ## mean function, mu_k = exp(eta_k)
  nll <- sum(mu - y*eta + lgamma(y+1))  ## negative log likelihood
  penalty <- lambda * sum(sqrt(beta[-1]^2 + eps))  ## beta_1 is not penalised
  nll + penalty
}

## poission_lasso_grad: first derivatives of the penalised cost function of the
##                      Poisson changepoint model with respect to beta
## inputs: beta, the parameter vector of length n, where beta_1 is the 
##         intercept value, with beta_k (k>=2) denoting the change in the
##         log-mean function relative to the last time step
##         y, the data vector
##         lambda, the penalty parameter, set a default value of 1
##         eps, the smoothing parameter for approximating the absolute function,
##         set a default value of 1e-6
## outputs: the gradient vector of length n
## work: compute the mean function from beta, construct the gradient of the 
##       negative log likelihood by summing mu_k - y_k where k = j,...,n,
##       for each beta_j (j = 1,...,n), add the derivative of the penalty 
##       for j = 2,..., n

poission_lasso_grad <- function(beta, y, lambda = 1, eps = 1e-6){
  eta <- cumsum(beta)  ## linear predictor, eta_k = beta_1 + ... + beta_k
  mu <- exp(eta)  ## mean function, mu_k = exp(eta_k)
  ## need to calculate the sum of mu - y from each position j to the end, cumsum
  ## adds up from the start, so reverse the vector first, then reverse back to
  ## keep the original order
  grad_nll <- rev(cumsum(rev(mu - y)))
  ## derivative of the penalty, beta_1 is not penalised
  grad_penalty <- c(0, (lambda*beta[-1]) / sqrt(beta[-1]^2 + eps))
  grad_nll + grad_penalty 
}

## poission_lasso_hess: second derivatives of the penalised cost function of the
##                      Poisson changepoint model with respect to beta
## inputs: beta, the parameter vector of length n, where beta_1 is the 
##         intercept value, with beta_k (k>=2) denoting the change in the
##         log-mean function relative to the last time step
##         y, the data vector
##         lambda, the penalty parameter, set a default value of 1
##         eps, the smoothing parameter for approximating the absolute function,
##         set a default value of 1e-6
## outputs: the n×n Hessian matrix
## work: compute the mean function from beta, construct the Hessian matrix of
##       the negative log likelihood by summing mu_k where k = max(i,j),...,n,
##       for each (beta_i, beta_j) with i,j = 1,...,n, add the second derivative
##       of the penalty on the diagonal for j = 2,...,n

poission_lasso_hess <- function(beta, y, lambda = 1, eps = 1e-6){
  n <- length(beta)
  eta <- cumsum(beta)  ## linear predictor, eta_k = beta_1 + ... + beta_k
  mu <- exp(eta)  ## mean function, mu_k = exp(eta_k)
  ## need to calculate the sum of mu from each position k to the end
  ## cumsum adds up from the start, so reverse the vector first, then reverse
  ## back to keep the original order
  mu_remaining <- rev(cumsum(rev(mu)))
  ## only when k >= i, mu_k changes as beta_i changes; only when k >= j,
  ## mu_k changes as beta_j changes
  ## so only when k >= max(i, j), mu_k is affected by both beta_i and beta_j
  hess_nll <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      hess_nll[i, j] <- mu_remaining[max(i, j)]
    }
  }
  ## the first derivative of the penalty with respect to beta_j only contains
  ## beta_j, so differentiating again with respect to beta_i will be 0
  ## so only the diagonal (j, j) will be nonzero
  ## beta_1 is not penalised, entry = 0
  hess_penalty <- diag(c(0, (lambda*eps) / (beta[-1]^2 + eps)^1.5))
  hess_nll + hess_penalty
}

## fd_test: using finite difference approximations to check the gradient and
##          Hessian matrix computed by code
## inputs: beta0, the parameter vector used in testing
##         y, the data vector
##         lambda, the penalty parameter, set a default value of 1
##         eps, the smoothing parameter for approximating the absolute function,
##         set a default value of 1e-6
##         h, the finite difference interval, set a default value of 1e-7
## outputs: a vector with two numbers, the largest absolute difference between
##          the result computed by code and the finite difference
##          approximations for gradient and Hessian matrix
## work: increase each beta_i (i = 1,...,n) by h each time, approximate the
##       gradient from the change in the objective and the Hessian from the
##       change in the gradient, then compare with the derivatives computed by
##       code

fd_test <- function(beta0, y, lambda = 1, eps = 1e-6, h = 1e-7){
  n <- length(beta0)
  obj0 <- poission_lasso_objective(beta0, y, lambda, eps)  ## objective at beta0
  grad0 <- poission_lasso_grad(beta0, y, lambda, eps)  ## gradient at beta0
  fd <- rep(0, n)  ## finite difference approximation for the gradient
  ## finite difference approximation for the Hessian matrix
  Hfd <- matrix(0, n, n)
  for (i in 1:n){
    beta1 <- beta0
    beta1[i] <- beta1[i] + h  ## increase beta_i by h
    obj1 <- poission_lasso_objective(beta1, y, lambda, eps)
    grad1 <- poission_lasso_grad(beta1, y, lambda, eps)
    ## derivative of the objective with respect to beta_i
    fd[i] <- (obj1 - obj0) / h
    ## derivative of the gradient with respect to beta_i
    Hfd[i, ] <- (grad1 - grad0) / h
  }
  hess0 <- poission_lasso_hess(beta0, y, lambda, eps)
  ## if the code is correct, the largest absolute difference between the
  ## derivatives computed by code and the approximations will be close to 0
  grad_diff <- max(abs(grad0 - fd))
  hess_diff <- max(abs(hess0 - Hfd))
  c(grad_diff, hess_diff)
}

## fit_changepoint: fit the Poisson changepoint model to the simulated dataset
##                  with a given penalty parameter
## inputs: y, the data vector
##         lambda, the penalty parameter, set a default value of 1
##         eps, the smoothing parameter for approximating the absolute function,
##         set a default value of 1e-6
## outputs: a list returned by optim
## work: initialise the parameter vector at beta_1 = log(mean(y)) and beta_k = 0
##       for all k >= 2, then minimise the objective using the BFGS method in
##       the optim command, and give a warning if optim did not converge

fit_changepoint <- function(y, lambda = 1, eps = 1e-6) {
  n <- length(y)
  ## start with a constant mean, mean(y), at every time point
  beta_start <- c(log(mean(y)), rep(0, n-1)) 
  ## find the beta that minimises the objective, starting from beta_start
  ## increase the maximum number of iterations from the default 100, as 100 is
  ## not enough to converge
  fit <- optim(beta_start, poission_lasso_objective, gr = poission_lasso_grad, 
               y = y, lambda = lambda, eps = eps, method = "BFGS", 
               control = list(maxit = 1000)) 
  if (fit$convergence != 0) warning("optim did not converge") 
  fit
}

## plt: plot the data with the true mean function and true changepoints
## inputs: time, the time points
##         y, the data vector
##         mu_true, the true mean function
##         tau_true, the ordered changepoint locations tau_0 < tau_1 < ...
##         < tau_m < tau_{m+1}, where m is the number of changepoints,
##         tau_0 = 1 and tau_{m+1} = n+1
##         ylim, the y-axis limits
## outputs: a plot
## work: plot the data with grey points, plot the true mean with black step
##       lines, plot the true changepoints with black dashed vertical lines
plt <- function(time, y, mu_true, tau_true, ylim = range(y)) {
  plot(time, y, pch = 16, col = "grey", xlab = "time", ylab = "y", ylim = ylim)
  lines(mu_true, col="black", lwd=2, type ="s")
  ## drop tau_0 and tau_{m+1} and plot changepoints
  abline(v = tau_true[-c(1,length(tau_true))], lty="dashed")
}

## simulate the data: 
## set.seed(12345) is used so that the results can be reproduced
set.seed(12345)
n <- 100
time <- 1:n
## changepoint positions (with tau_0 and tau_{m+1})
tau_true <- c(1, 30, 70, n+1)
m_true <- length(tau_true)-2  ## number of cpts (m)
theta_true <- c(10, 2, 7)  ## segment parameters (theta_1, theta_2, theta_3)               
mu_true <- numeric(n)  ## mean function
for(i in seq_along(theta_true))
  mu_true[tau_true[i]:(tau_true[i+1]-1)] <- theta_true[i]
y <- rpois(n, lambda = mu_true) ## simulate data y
## plot the simulated data, the true mean function, and the true changepoints
plt(time, y, mu_true, tau_true)

## test the derivative functions by finite differencing:
## beta_1 = log(mean(y)), so that mu_1 equals the average of y, and set the
## other parameters small and nonzero
beta0 <- c(log(mean(y)), rep(0.05, n - 1))
## if the derivatives are correct, the two differences will be close to 0
fd_test(beta0, y) 

## preliminary sanity check:
## fit the model with the default penalty parameter of lambda = 1
fit <- fit_changepoint(y, lambda = 1)
## compute the estimated mean function from the estimated beta
mu_lambda1 <- exp(cumsum(fit$par))
## plot the data, the true mean function, the estimated mean function
plt(time, y, mu_true, tau_true)
lines(mu_lambda1, col = "red", lwd = 2, type = "s")  ## add the estimated mean
legend("topright", legend = c("true mean", "estimated mean"),
       col = c("black", "red"), lwd = 2)

## choosing an appropriate lambda by minimising the mBIC:
## search over the grid of log(lambda) values to find the penalty parameter
## that minimises the mBIC
log_lambda <- seq(-1, 5, by = 0.05)
n_grid <- length(log_lambda)
m_tilde <- rep(0, n_grid)  ## the effective number of changepoints
mbic <- rep(0, n_grid)  ## mBIC at each lambda
for (i in seq_along(log_lambda)){
  lambda <- exp(log_lambda[i])
  ## estimate beta using optim (fit_changepoint warns if it did not converge)
  fit <- fit_changepoint(y, lambda = lambda)
  beta_hat <- fit$par
  ## evaluate the Hessians H_lambda and H_0 for the parameter estimates, and
  ## compute the effective degrees of freedom 
  ## EDF(lambda) = trace(H_lambda^{-1} H_0)
  H_lambda <- poission_lasso_hess(beta_hat, y, lambda = lambda)
  H_0 <- poission_lasso_hess(beta_hat, y, lambda = 0)
  edf_lambda <- sum(diag(solve(H_lambda, H_0)))
  ## evaluate the maximum log likelihood, estimate the effective number of
  ## changepoints, and then compute mBIC(lambda)
  ll <- -poission_lasso_objective(beta_hat, y, lambda = 0)
  m_tilde[i] <- edf_lambda - 1 ## effective number of changepoints
  mbic[i] <- 3*m_tilde[i]*log(n) - 2*ll
}

## find the best penalty value lambda_hat that minimises the mBIC
best <- which.min(mbic) ## index of smallest mBIC
lambda_hat <- exp(log_lambda[best]) ## best penalty value
## display the results
lambda_hat
m_tilde[best]

## fit the model with lambda_hat:
## refit the changepoint model with the best penalty value lambda_hat to get
## the final estimate of beta
fit_best <- fit_changepoint(y,lambda = lambda_hat)
beta_best <- fit_best$par
mu_hat <- exp(cumsum(beta_best))  ## estimated mean function
## plot the data, the true mean function, and the estimated mean function
plt(time, y, mu_true, tau_true)
lines(mu_hat, col = "red", lwd = 2, type = "s")  ## add the estimated mean
legend("topright", legend = c("true mean", "estimated mean"),
       col = c("black", "red"), lwd = 2)
## estimate the number of changepoints by rounding down the effective number
## of changepoints
m_hat <- floor(m_tilde[best])
## drop beta_1, order the estimated parameters |beta_k| from largest to smallest
ord <- order(abs(beta_best[-1]), decreasing = TRUE)
## add 1 because beta_1 was dropped, which shifted every position back by 1
tau_hat <- ord[1:m_hat] + 1
m_hat  ## estimated number of changepoints
tau_hat  ## estimated changepoint positions

## assess the uncertainty of the mean function by parametric bootstrapping:
## simulate B = 200 data sets from a Poisson model with mean function mu_hat,
## refit the model to each one (using the best penalty value lambda_hat for
## all simulations) and evaluate the mean function based on the new estimates

B <- 200
mu_bootstrap <- matrix(0, B, n)  ## row b is the estimate from data set b
for (b in 1:B){
  ## simulate a new data set y_sim 
  y_sim <- rpois(n, lambda = mu_hat)  ## the lambda of rpois is the Poisson mean
  ## refit the model to y_sim, fit_changepoint gives a warning if optim did 
  ## not converge
  fit_sim <- fit_changepoint(y_sim, lambda = lambda_hat)
  mu_bootstrap[b, ] <- exp(cumsum(fit_sim$par))
}
mu_bar <- colMeans(mu_bootstrap) ## estimated mean function from bootstrap
## plot the data, the true mean function, the estimated mean function mu_hat,
## and the bootstrap estimate mu_bar
plt(time, y, mu_true, tau_true)
lines(mu_hat, col = "red", lwd = 2, type = "s") ## estimated mean
lines(mu_bar, col = "blue", lwd = 2, type = "s") ## bootstrap mean
legend("topright", legend = c("true mean", "estimated mean",
                              "bootstrap mean"), 
       col = c("black", "red", "blue"), lwd = 2)
## to get the bias corrected mean function, subtract the estimated bias
## mu_bar - mu_hat from mu_hat
mu_bc <- 2*mu_hat - mu_bar ## bias corrected mean function
## evaluate the bias corrected 95% confidence intervals
q <- apply(mu_bootstrap, 2, quantile, probs = c(0.025, 0.975))
qL <- q[1, ]  ## lower limit
qU <- q[2, ]  ## upper limit
qL_bc <- 2*mu_hat - qU ## bias corrected lower limit
qU_bc <- 2*mu_hat - qL ## bias corrected upper limit

## final plot:
## produce a plot showing the original data y, the true mean function mu, the
## estimated mean function mu_hat, the bias corrected mean function mu_bc, and
## the bias corrected 95% confidence limits qL_bc and qU_bc
## ylim is set so that the plot fits the data and the confidence limits
plt(time, y, mu_true, tau_true, ylim = range(c(y, qL_bc, qU_bc)))
lines(time, mu_hat, col = "red", lwd = 2, type = "s")  ## estimated mean
lines(time, mu_bc, col = "blue", lwd = 2, type = "s")  ## bias corrected mean
lines(time, qL_bc, col = "orange", lty = "dashed", type = "s")  ## lower limit
lines(time, qU_bc, col = "orange", lty = "dashed", type = "s")  ## upper limit
legend("topright", legend = c("data", "true mean", "true changepoints",
                              "estimated mean", "bias corrected mean", 
                              "bias corrected 95% limits"),
       col = c("grey", "black", "black", "red", "blue", "orange"),
       pch = c(16, NA, NA, NA, NA, NA), lty = c(NA, 1, 2, 1, 1, 2),
       lwd = c(NA, 2, 1, 2, 2, 1))







