#############################################################################################################
library(MASS)
# Set the data
y <- c(78, 64, 90, 78, 83, 82, 89)
t <- c(14.1, 13.2, 15.4, 14.9, 15.6, 15.2, 16.6)
t

# Define the priors
alpha <- 1
beta <- 0.1
sigma <- 1

# Define the hyperparameters for the priors
a <- 2  # Shape parameter for alpha gamma distribution
b <- 0.5   # Scale parameter for alpha gamma distribution
lambda <- 1   # Rate parameter for beta exponential distribution

# Define the prior distributions
prior_alpha <- function(alpha) {
  dgamma(alpha, shape = a, scale = b)
}
prior_beta <- function(beta) {
  dexp(beta, rate = lambda)
}
prior_sigma <- function(sigma) {
  dgamma(sigma, shape = 0.1, scale = 1)
}

likelihood <- function(alpha, beta, sigma, y, t) {
  lambda <- alpha * t
  log_likelihood <- sum(dpois(y, lambda, log = TRUE))
  return(exp(log_likelihood))
}

# Set up a range of values for each parameter
alpha_seq <- seq(0, 5, length.out = 100)
beta_seq <- seq(0, 5, length.out = 100)
sigma_seq <- seq(0, 2, length.out = 100)

# Compute the prior densities for each parameter
prior_alpha_dens <- prior_alpha(alpha_seq)
prior_beta_dens <- prior_beta(beta_seq)
prior_sigma_dens <- prior_sigma(sigma_seq)

# Plot the prior distributions
par(mfrow = c(1, 3))
plot(alpha_seq, prior_alpha_dens, type = "l", xlab = expression(alpha), ylab = "Density", main = "Prior for Alpha")
plot(beta_seq, prior_beta_dens, type = "l", xlab = expression(beta), ylab = "Density", main = "Prior for Beta")
plot(sigma_seq, prior_sigma_dens, type = "l", xlab = expression(sigma), ylab = "Density", main = "Prior for Sigma")


# Define the joint posterior
joint_posterior <- function(alpha, beta, sigma, y, t) {
  prior_alpha_dens <- prior_alpha(alpha)
  prior_beta_dens <- prior_beta(beta)
  prior_sigma_dens <- prior_sigma(sigma)
  lik <- likelihood(alpha, beta, sigma, y, t)
  
  posterior <- lik * prior_alpha_dens * prior_beta_dens * prior_sigma_dens
  return(posterior)
}

# Compute the conditional posterior distributions
# Theta_j | alpha, beta, y
theta_cond_posterior <- function(j, alpha, beta, sigma, y, t) {
  prior_alpha_dens <- prior_alpha(alpha)
  prior_beta_dens <- prior_beta(beta)
  lik <- likelihood(alpha, beta, sigma, y, t)
  
  marginal_sigma_dens <- integrate(function(log_sigma) {
    sigma <- exp(log_sigma)
    lik * prior_alpha_dens * prior_beta_dens * prior_sigma(sigma)
  }, lower = -Inf, upper = Inf)$value
  
  posterior_unnorm <- function(theta_j) {
    theta <- c(alpha, beta, sigma)
    theta[j] <- theta_j
    lik_j <- likelihood(theta[1], theta[2], theta[3], y, t)
    posterior <- lik_j * prior_alpha_dens * prior_beta_dens * prior_sigma(theta[3])
    return(posterior)
  }
  theta_cond_posterior_dens <- function(theta_j) {
    posterior_unnorm(theta_j) / marginal_sigma_dens
  }
  
  return(theta_cond_posterior_dens)
}

# Alpha | beta, theta, y
alpha_cond_posterior <- function(alpha, beta, sigma, y, t) {
  prior_beta_dens <- prior_beta(beta)
  lik <- likelihood(alpha, beta, sigma, y, t)
  theta <- c(alpha, beta, sigma)
  
  # Compute the normalizing constant (marginal likelihood) for sigma and alpha
  marginal_sigma_alpha_dens <- integrate(function(log_sigma, alpha) {
    sigma <- exp(log_sigma)
    lik * prior_alpha(alpha) * prior_beta_dens * prior_sigma(sigma)
  }, lower = -Inf, upper = Inf, alpha = alpha)$value
  
  # Define the unnormalized posterior density for alpha
  posterior_unnorm <- function(alpha) {
    lik_alpha <- likelihood(alpha, theta[2], theta[3], y, t)
    posterior <- lik_alpha * prior_alpha(alpha) * prior_beta_dens * prior_sigma(theta[3])
    return(posterior)
  }
  
  # Compute the fully conditional posterior density for alpha
  alpha_cond_posterior_dens <- function(alpha) {
    posterior_unnorm(alpha) / marginal_sigma_alpha_dens
  }
  
  return(alpha_cond_posterior_dens)
}

# Beta | alpha, theta, y
beta_cond_posterior <- function(alpha, beta, sigma, y, t) {
  prior_alpha_dens <- prior_alpha(alpha)
  lik <- likelihood(alpha, beta, sigma, y, t)
  theta <- c(alpha, beta, sigma)
  
  # Compute the normalizing constant (marginal likelihood) for sigma and beta
  marginal_sigma_beta_dens <- integrate(function(log_sigma, beta) {
    sigma <- exp(log_sigma)
    lik * prior_alpha_dens * prior_beta(beta) * prior_sigma(sigma)
  }, lower = -Inf, upper = Inf, beta = beta)$value
  
  # Define the unnormalized posterior density for beta
  posterior_unnorm <- function(beta) {
    lik_beta <- likelihood(theta[1], beta, theta[3], y, t)
    posterior <- lik_beta * prior_alpha_dens * prior_beta(beta) * prior_sigma(theta[3])
    return(posterior)
  }
  
  # Compute the fully conditional posterior density for beta
  posterior_norm <- integrate(posterior_unnorm, lower = -Inf, upper = Inf)$value / marginal_sigma_beta_dens
  return(posterior_norm)
}

# Define the unnormalized posterior of alpha

alpha.posterior <- function(a, b, thetas) {# Data has thetas in column 1, betas in column 2
  beta <- b
  dens <- exp(-a)*prod(((thetas^(a-1))*(b^a))/gamma(a))
  return(dens)
}

#install.packages("truncnorm")
library(truncnorm)

# Proposal distribution
prop.dist.alpha <- function(a, prop.var) {
  rtruncnorm(1, mean=a, sd=sqrt(prop.var), a=0)
} 

# Density of proposal
prop.dist.alpha.dens <- function(a, a.mean, prop.var) {
  dtruncnorm(a, mean = a.mean, sd=sqrt(prop.var), a=0)
}

# Metropolis-Hastings Algorithm
metrop <- function(param, thetas, b, alpha.posterior, prop.dist.alpha, prop.dist.alpha.dens, prop.var, n.iter) {
  # Store sampled alpha values
  alphas <- c()
  # Initialize model
  param.t <- param
  for(t in 1:n.iter) {
    # Draw proposed value of alpha
    param.new <- prop.dist.alpha(param.t, prop.var)
    # Calculate acceptance probability
    u <- runif(1, 0, 1)
    prob.accept <- min(1, (alpha.posterior(param.new, b, thetas)*
                             (prop.dist.alpha.dens(param.t, param.new, prop.var)))/(alpha.posterior(param.t, b, thetas)
                              *(prop.dist.alpha.dens(param.new, param.t, prop.var))))
    
    if(u < prob.accept) {
      value <- param.new
    } else {
      value <- param.t
    }
    alphas <- c(alphas, value)
    param.t <- value
  }
  
  # Modification for MH-within-Gibbs sampling --> if only drawing one sample, return the sampled value.
  # If drawing multiple samples, return the list of all samples
  if (length(alphas) == 1) {
    return(alphas[1])
  } else {
    return(alphas)
  }
}

gibbs <- function(initial, y, t, n.iter) {
  # Initialize variables
  J <- length(y)
  l <- length(initial)
  results <- matrix(NA, n.iter, l)
  results[1,] <- initial
  
  for(i in 2:n.iter) {
    thetas <- results[i-1, 1:7] # Stores all 7 theta_j values
    a <- results[i-1,8]
    b <- results[i-1,9]
    
    # Draw theta_j samples 
    for(j in 1:J) {
      
      # Find alpha, beta parameters for theta_j's gamma posterior distribution
      alpha.theta <- y[j] + a
      beta.theta <- t[j] + b
      
      # Store singular theta_j sample using parameters calculated above
      results[i,j] <- rgamma(1, alpha.theta, beta.theta)
    }
    
    # Find alpha, beta parameters for beta's gamma posterior distribution using theta sample
    alpha.beta <- J*a + 0.1
    beta.beta <- 1 + sum(results[i, 1:7])
    
    # Store singular beta sample from its gamma posterior distribution
    results[i, 9] <- rgamma(1, alpha.beta, beta.beta)
    
    # Use Metropolis-Hastings algorithm from above to draw singular alpha sample
    results[i, 8] <- metrop(a, results[i, 1:7], results[i, 9], alpha.posterior, prop.dist.alpha, prop.dist.alpha.dens, 4, 1)
  }
  return(results)
}

set.seed(5000) # Set seed for reproducability


# Store prior data
yhat <- c(78, 64, 90, 78, 83, 82, 89)
t <- c(14.1, 13.2, 15.4, 14.9, 15.6, 15.2, 16.6)


# Initialize first set of sampling values
initial_1 <- c(rep(.1, 7), 1, 1)
# Draw samples
sample_1 <- gibbs(initial_1, yhat, t, 10000)

par(mfrow=c(3,3))
for (i in 1:9) {
  hist(sample_1[,i], breaks= 30, main = paste0("sample ", i))
}

means <- matrix(NA, 9, 1)
for (i in 1:9) {
  means[i,1] <- mean(sample_1[,i])
}
rownames(means) <- c("Theta1", "Theta2", "Theta3", "Theta4", "Theta5","Theta6","Theta7", "Alpha", "Beta")
colnames(means) <- c("Posterior Mean")
means <- as.table(means)
means

#sample_1
head(sample_1, n=10)

par(mfrow = c(3, 3))
par(mar = c(1, 1, 1, 1))
for (i in 1:9) {
  plot(sample_1[,i], type="l")
}


# Set initial values and run MCMC
initial_1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
sample_1 <- gibbs(initial_1, yhat, t, 10000)

initial_2 <- c(rep(5, 7), 0.5, 0.5)
sample_2 <- gibbs(initial_2, yhat, t, 10000)

initial_3 <- c(rep(1, 7), 10, 10)
sample_3 <- gibbs(initial_3, yhat, t, 10000)

initial_4 <- c(rep(0, 7), 1, 1)
sample_4 <- gibbs(initial_4, yhat, t, 10000)

initial_5 <- c(rep(0, 7), 0.1, 0.1)
sample_5 <- gibbs(initial_5, yhat, t, 10000)

# Store MCMC chain
chains <- list(sample_1, sample_2, sample_3, sample_4, sample_5)

# Run Gelman-Rubin diagnostic

n_chains <- length(chains)
n_samples <- nrow(chains[[1]])
n_params <- ncol(chains[[1]])
B <- n_samples * var(sapply(chains, function(chain) rowMeans(chain)))
W <- mean(sapply(chains, function(chain) var(chain))) * n_samples
var_hat <- ((n_samples - 1) / n_samples) * W + (1 / n_samples) * B
R_hat <- sqrt(var_hat / W)
R_hat

#install.packages("MCMCpack")
library(MCMCpack)

chain_1 <- mcmc(sample_1)
chain_2 <- mcmc(sample_2)
chain_3 <- mcmc(sample_3)
chain_4 <- mcmc(sample_4)
chain_5 <- mcmc(sample_5)


# Store MCMC chain
combined.chains <- mcmc.list(chain_1, chain_2, chain_3, chain_4, chain_5)

# Run Gelman-Rubin diagnostic
gelman.rubin <- gelman.diag(combined.chains)
gelman.rubin

# Combine results from all 5 samples
total.samples <- rbind(sample_1, sample_2, sample_3, sample_4, sample_5)
post.info <- matrix(NA, 9, 3)
for (i in 1:9) {
  post.info[i, 1] <- round(mean(total.samples[,i]), 3)
  post.info[i, 2] <- round(quantile(total.samples[,i], 0.025), 3)
  post.info[i, 3] <- round(quantile(total.samples[,i], 0.975), 3)
}
rownames(post.info) <- c("Theta1", "Theta2", "Theta3", 
                         "Theta4", "Theta5","Theta6","Theta7", "Alpha", "Beta")
colnames(post.info) <- c("Posterior Mean", "Lower Bound", "Upper Bound")
post.info <- as.table(post.info)
names(dimnames(post.info)) <- c("Parameter", "95% Posterior Interval Information")
post.info

###############################################################################################################


