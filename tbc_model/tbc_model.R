# MCMC fit of a mixture of poissons to each farm, full model

# Read the data in
data_tbc <- read.csv(file="../data/tbc_data.csv")

# pull out just the data we need (TBC, Farm), and exclude missing data
d <- cbind(data_tbc$TBC, data_tbc$Farm)
d <- data.frame(na.exclude(d))
names(d) <- c("TBC", "Farm")

# sort by farm
d <- d[order(d$Farm),]

# variables
y     <- d$TBC
farm  <- d$Farm

# initial parameter values
n_farms <- length(unique(d$Farm))

# run through the farms and construct before and after
farm_map <- matrix(0, n_farms, 2)
for (i in 1:n_farms)
{
  wch <- which(d$Farm==i)
  farm_map[i,] <- range(wch)
}

# size of dataset
n     <- length(y)

#' Model is
#' 
#' $$
#' TBC_{ij} \sym Poisson(\alpha_j + \beta_j x_j)
#' $$
#' 
#' where
#' 
#' $$
#' \begin{aligned}
#'    x_j &\sym Bernoulli(w_j)\\
#'    w_j &\sym Beta(w_\alpha, w_\beta)\\
#'    \alpha_j &\sym Normal(\mu_\alpha, \tau_\alpha)\\
#'    \beta_j  &\sym Normal(\mu_\beta, \tau_\beta)\\
#'    \mu_\alpha, \mu_\beta &\sym ~ Normal(0, \tau_mu)\\
#'    \tau_\alpha, \tau_\beta &\sym ~ Gamma(\alpha, \beta)
#' \end{aligned}
#'
#' hyper priors
w_alpha    <- 1      # hyperprior for w_j
w_beta     <- 1      # hyperprior for w_j
tau_alpha  <- 0.001  # prior for precision of \alpha_j
tau_beta   <- 0.001  # prior for precision of \beta_j
tau_mu     <- 0.001  # hyperprior, precision of mu_\alpha, \mu_\beta

burnin     <- 20000;
iters      <- 100000;
thinning   <- 50;

# posteriors
w_post     <- rep(0, iters/thinning)
alpha_post <- matrix(0, iters/thinning, n_farms)
beta_post  <- matrix(0, iters/thinning, n_farms)
x_post     <- matrix(0, iters/thinning, n)
mu_a_post  <- rep(0, iters/thinning)
mu_b_post  <- rep(0, iters/thinning)
tau_a_post <- rep(0, iters/thinning)
tau_b_post <- rep(0, iters/thinning)

# acceptance rates
alpha_accept <- rep(0, n_farms)
beta_accept  <- rep(0, n_farms)
swap_accept  <- rep(0, n_farms)

# parameters in our model
w     <- rep(0.1); 	# Probability of high event
alpha <- rep(0, n_farms);   # random effect for farm level (note: we assume this is additive for high/low)
beta  <- rep(0, n_farms);   # random effect for farm level (note: we assume this is additive for high/low)
x     <- rep(0, n);         # Latent variables for high/low
mu_a  <- 0;                 # Mean level for random effects (alpha)
mu_b  <- 0;                 # Mean level for random effects (beta)
tau_a <- 1;                 # Precision of random effects (alpha)
tau_b <- 1;                 # Precision of random effects (beta)

# Sample the latent variables conditional on the parameters and data
SampleX <- function()
{
  for (i in 1:n)
  {
#    p1 <- p * exp((beta[farm[i]])*y[i] - exp(alpha[farm[i]] + beta[farm[i]]))
#    p0 <- (1-p) * exp(-exp(alpha[farm[i]]))
    log_r1 <- log(w)   + beta[farm[i]]*y[i] - exp(alpha[farm[i]] + beta[farm[i]]);
    log_r0 <- log(1-w) - exp(alpha[farm[i]]);
    max_r  <- max(log_r0, log_r1);
    log_r1 <- log_r1 - max_r;
    log_r0 <- log_r0 - max_r;
    p1 <- exp(log_r1)
    p0 <- exp(log_r0)
    if (!is.finite(p1)) {
      x[i] <<- 1;
    } else {
      x[i] <<- rbinom(1, 1, p1/(p0+p1))
      if (is.na(x[i]))
      {
        cat("p1=", p1, "p0=", p0, "p=", p, "beta[i]=",beta[farm[i]], "alpha[i]=", alpha[farm[i]],"\n");
        stop("na x")
      }
    }
  }
}


RunMCMC <- function(iters, burnin, thinning)
{
  # tuning parameters
  alpha_jump <- 1;
  beta_jump  <- 1;
  
  sample <- 0
  
  for (it in 1:(iters+burnin))
  {
    # 1. Gibbs sample our latent variables
    SampleX()
  
    # 2. Sample alpha (dumb, independent M-H)
    for (j in 1:n_farms)
    {
      y_j  <- y[farm_map[j,1]:farm_map[j,2]]
      x_j  <- x[farm_map[j,1]:farm_map[j,2]]
      alpha_n <- rnorm(1, alpha[j], alpha_jump)
      ll_n <- sum(y_j*alpha_n  - exp(alpha_n  + beta[j]*x_j))
      ll   <- sum(y_j*alpha[j] - exp(alpha[j] + beta[j]*x_j))
      lp_n <- -0.5*(alpha_n  - mu_a)^2*tau_a
      lp   <- -0.5*(alpha[j] - mu_a)^2*tau_a
      la   <- ll_n + lp_n - (ll + lp)
      u    <- runif(1)
      if (is.na(la))
      {
        cat("Error in delta update\n")
      }
      if (la > 0 || u < exp(la))
      {
        alpha[j] <<- alpha_n
        if (it > burnin)
          alpha_accept[j] <<- alpha_accept[j] + 1
      }
    }
    # 3. Sample beta (dumb, independent M-H)
    for (j in 1:n_farms)
    {
      y_j  <- y[farm_map[j,1]:farm_map[j,2]]
      x_j  <- x[farm_map[j,1]:farm_map[j,2]]
      beta_n <- rnorm(1, beta[j], beta_jump)
      ll_n <- sum(y_j*x_j*beta_n  - exp(alpha[j] + beta_n *x_j))
      ll   <- sum(y_j*x_j*beta[j] - exp(alpha[j] + beta[j]*x_j))
      lp_n <- -0.5*(beta_n  - mu_b)^2*tau_b
      lp   <- -0.5*(beta[j] - mu_b)^2*tau_b
      la   <- ll_n + lp_n - (ll + lp)
      u    <- runif(1)
      if (is.na(la))
      {
        cat("Error in beta update\n")
      }
      if (la > 0 || u < exp(la))
      {
        beta[j] <<- beta_n
        if (it > burnin)
          beta_accept[j] <<- beta_accept[j]+1
      }
    }
    # 3b. Swapping alpha and beta (doesn't change the likelihood)
    for (j in 1:n_farms)
    {
      # propose to swap alpha and alpha+beta
      alpha_n <- alpha[j] + beta[j];
      beta_n <- -beta[j];
      x_o  <- x[farm_map[j,1]:farm_map[j,2]];
      x_n  <- 1 - x_o;
  
      # likelihood won't change
      lp_n <- -0.5*(beta_n  - mu_b)^2*tau_b - 0.5*(alpha_n  - mu_a)^2*tau_a
      lp_n <- lp_n + sum(x_n)*log(w) + sum(x_o)*log(1-w)
  
      lp   <- -0.5*(beta[j] - mu_b)^2*tau_b - 0.5*(alpha[j] - mu_a)^2*tau_a
      lp   <- lp + sum(x_o)*log(w) + sum(x_n)*log(1-w)
  
      la   <- lp_n - lp
      u    <- runif(1)
      if (is.na(la))
      {
        cat("Error in beta update\n")
      }
      if (la > 0 || u < exp(la))
      {
        beta[j] <<- beta_n
        alpha[j] <<- alpha_n
        x[farm_map[j,1]:farm_map[j,2]] <<- x_n
        if (it > burnin)
          swap_accept[j] <<- swap_accept[j]+1
      }
    }
    # 4. Gibbs update of w
    {
      w <<- rbeta(1, w_alpha+sum(x), w_beta + length(x) - sum(x))
    }
    # 5. Gibbs update of mu
    posterior_tau <- (tau_mu + n_farms*tau_a);
    posterior_mu  <- tau_a*sum(alpha) / posterior_tau;
    mu_a <<- rnorm(1, posterior_mu, 1/sqrt(posterior_tau))
    posterior_tau <- (tau_mu + n_farms*tau_b);
    posterior_mu  <- tau_b*sum(beta) / posterior_tau;
    mu_b <<- rnorm(1, posterior_mu, 1/sqrt(posterior_tau))
  
    # 6. Gibbs update of sigma
    tau_a <<- rgamma(1, tau_alpha + n_farms/2, tau_beta + sum((alpha-mu_a)^2)/2)
    tau_b <<- rgamma(1, tau_alpha + n_farms/2, tau_beta + sum((beta-mu_b)^2)/2)
  
    if (it > burnin && (it-burnin) %% thinning == 0)
    { # sample
      sample <- sample + 1
      w_post[sample]      <<- w
      alpha_post[sample,] <<- alpha
      beta_post[sample,]  <<- beta
      x_post[sample,]     <<- x
      mu_a_post[sample]   <<- mu_a
      mu_b_post[sample]   <<- mu_b
      tau_a_post[sample]  <<- tau_a
      tau_b_post[sample]  <<- tau_b
    }
    if ((it-burnin) %% thinning == 0)
      cat("iteration", it, "of", burnin+iters,"\n")
  }
}

cat("Time taken for", iters+burnin, "iterations:", system.time(RunMCMC(iters, burnin, thinning))[3],"seconds \n")

# write our results out to the output folder
out_folder <- "output"

system(paste("mkdir", out_folder))
# reorder, forcing beta > 0 if necessary
#for (i in 1:n_farms)
#{
#  found <- beta_post[,i] < 0
#  alpha_post[found,i] <- alpha_post[found,i] + beta_post[found,i]
#  beta_post[found,i] <- -beta_post[found,i]
#  p_post[found,i] <- 1 - p_post[found,i]
#  x_post[found,i] <- 1 - x_post[found,i]
#}

accept <- data.frame(alpha=alpha_accept, beta=beta_accept, swap=swap_accept) / iters
write.csv(accept,     file.path(out_folder, "acceptance.csv"), row.names=F)

write.csv(mu_a_post,  file.path(out_folder, "mu_a.csv"),  row.names=F)
write.csv(mu_b_post,  file.path(out_folder, "mu_b.csv"),  row.names=F)
write.csv(tau_a_post, file.path(out_folder, "tau_a.csv"), row.names=F)
write.csv(tau_b_post, file.path(out_folder, "tau_b.csv"), row.names=F)
write.csv(w_post,     file.path(out_folder, "w.csv"),     row.names=F)
write.csv(alpha_post, file.path(out_folder, "alpha.csv"), row.names=F)
write.csv(beta_post,  file.path(out_folder, "beta.csv"),  row.names=F)
write.csv(x_post,     file.path(out_folder, "x.csv"),     row.names=F)
write.csv(farm_map,   file.path(out_folder, "farm_map.csv"), row.names=F)

