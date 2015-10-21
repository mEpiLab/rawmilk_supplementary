# Survey modelling.

# Idea is that we simulate our surveys a heap and see what prevalence
# we observe given the shedding rates, farm and animal level prevalence
# and TBC distribution.

# Helper function to load in the TBC data from our posterior model fit
# and then fit a mixture of normals to these to see how well it does.
# we'll then sample from these normals rather than directly from the
# posterior as the normal mixture fits really well and it's much easier
# to sample from

norm_tbc <- function(abm) {
  # Those where the means are within 1.4 s.d of each other mayaswell
  # be one distribution, so don't count as 'high'
  abm$high <- (abs(abm$m2-abm$m1) / sqrt(abm$s1^2+abm$s2^2)) > 1.4
  for (i in 1:nrow(abm)) {
    if (abm[i,3] > abm[i,4]) {
      abm[i,1:6] <- abm[i,c(2,1,4,3,6,5)]
    }
  }
  abm
}

load_tbc <- function()
{
  # load in our tbc information per farm and categorise into those that have mixtures and those that don't
  library(mixtools)

  data_tbc <- read.csv(file="../data/tbc_data.csv")

  # pull out just the data we need (TBC, Farm), and exclude missing data
  d <- cbind(data_tbc$TBC, data_tbc$Farm)
  d <- data.frame(na.exclude(d))
  names(d) <- c("TBC", "Farm")

  # sort by farm
  d <- d[order(d$Farm),]

  farm  <- d$Farm

  tbc_folder  <- "../tbc_model/output"

  # now do plots etc.
  mu_a_post  <- read.csv(file.path(tbc_folder, "mu_a.csv"))
  mu_b_post  <- read.csv(file.path(tbc_folder, "mu_b.csv"))
  tau_a_post <- read.csv(file.path(tbc_folder, "tau_a.csv"))
  tau_b_post <- read.csv(file.path(tbc_folder, "tau_b.csv"))
  w_post     <- read.csv(file.path(tbc_folder, "w.csv"))
  alpha_post <- read.csv(file.path(tbc_folder, "alpha.csv"))
  beta_post  <- read.csv(file.path(tbc_folder, "beta.csv"))
  x_post     <- read.csv(file.path(tbc_folder, "x.csv"))
  
  n_farms    <- ncol(alpha_post)
  n_samples  <- nrow(alpha_post)
  
  # convert to matrices/numeric
  mu_a_post  <- mu_a_post[,1]
  mu_b_post  <- mu_b_post[,1]
  tau_a_post <- tau_a_post[,1]
  tau_b_post <- tau_b_post[,1]
  w_post     <- w_post[,1]
  alpha_post <- as.matrix(alpha_post)
  beta_post  <- as.matrix(beta_post)
  x_post     <- as.matrix(x_post)
  
  abm <- matrix(0,n_farms,6)
  for (i in 1:n_farms) {
    b = as.numeric(alpha_post[,i] + x_post[,farm==i]*beta_post[,i])
    nm = normalmixEM(b, k=2, maxit=20000, maxrestarts=100)
    abm[i,] <- c(nm$lambda, nm$mu, nm$sigma)
    if (min(nm$lambda*nm$sigma) < 0.001) { 
      # tiny sd and weight - try fitting a 1 component model
      abm[i,] <- c(1,0,mean(b),0,sd(b),1)
    }
  }
  colnames(abm) <- c("w1","w2","m1","m2","s1","s2")
  abm = data.frame(abm)

  return(norm_tbc(abm))
}

# Function to load in the number of samples and positives for a particular
# disease from the survey data
load_pos <- function(type)
{
  pos <- read.csv("../data/pathogen_positives.csv")

  d <- matrix(0, 80, 2)
  for (j in 1:80)
  {
    d[j,1] <- sum(!is.na(pos[pos$Farm == j,type]))
    d[j,2] <- sum(pos[pos$Farm == j,type], na.rm=T)
  }
  data.frame(total=d[,1], counts=d[,2])
}

# main survey sampling function, assuming normal mixture outcomes
sample_survey <- function(farms, tbc, herd_prev_func, cow_prev_func, path_high_func, path_low_func, sensitivity)
{
  dummy <- 0
  n_sim <- 10000

  n <- rep(1:80,farms)
  len <- length(n)

  values <- rep(0, n_sim)
  for (i in 1:n_sim)
  {
    high <- path_high_func(len, dummy)
    low  <- path_low_func(len, dummy)

    tbc_high <- rbinom(len, 1, tbc[n,2])
    tbc_sample <- exp((1-tbc_high)*rnorm(len, tbc[n,3], tbc[n,5]) + tbc_high * rnorm(len, tbc[n,4], tbc[n,6]))

    farm_pos   <- rbinom(len, 1, herd_prev_func(len, dummy))
    cow_pos    <- rbinom(len, 1, cow_prev_func(len, dummy))

    # we have a high-level event (i.e. sample from the "high" distributions)
    # if all of the following hold:
    # * The farm is positive
    # * The cow is positive
    # * The farm can have high level events
    # * The farm has a high level event
    high_event <- tbc_high & farm_pos & cow_pos

    # we have a low-level event (i.e. sample from the "low" distributions)
    # if all of the following hold:
    # * The farm is positive
    # * We don't have a high level event
    low_event  <- !high_event & farm_pos

    out <- rep(0,len)
    cfu_g <- 10^rnorm(len, mean=4.5, sd=0.3)
    out[high_event] <- tbc_sample[high_event] * high[high_event] / cfu_g[high_event]
    out[low_event]  <- tbc_sample[low_event]  * low[low_event] / cfu_g[low_event]

    s <- rpois(len, out * 0.025)

    # see whether we detected this or not...
    values[i] <- sum(s >= sensitivity)
    if (i %% 1000 == 0)
      cat("simulating",i,"of", n_sim, "\n")
  }
  return(values)
}

# counts forward prediction
sample_counts <- function(tbc, herd_prev_func, cow_prev_func, path_high_func, path_low_func)
{
  dummy <- 0
  n_sim <- 100000

  values <- rep(0, n_sim)
  for (i in 1:(n_sim/80))
  {
    n <- 1:80
    len <- length(n)

    high <- path_high_func(len, dummy)
    low  <- path_low_func(len, dummy)

    tbc_high <- rbinom(len, 1, tbc[n,2])
    tbc_sample <- exp((1-tbc_high)*rnorm(len, tbc[n,3], tbc[n,5]) + tbc_high * rnorm(len, tbc[n,4], tbc[n,6]))

    farm_pos   <- rbinom(len, 1, herd_prev_func(len, dummy))
    cow_pos    <- rbinom(len, 1, cow_prev_func(len, dummy))

    # we have a high-level event (i.e. sample from the "high" distributions)
    # if all of the following hold:
    # * The farm is positive
    # * The cow is positive
    # * The farm can have high level events
    # * The farm has a high level event
    high_event <- tbc_high & farm_pos & cow_pos # & tbc[n,7]

    # we have a low-level event (i.e. sample from the "low" distributions)
    # if all of the following hold:
    # * The farm is positive
    # * We don't have a high level event
    low_event  <- !high_event & farm_pos

    out <- rep(0,len)
    cfu_g <- 10^rnorm(len, mean=4.5, sd=0.3)
    out[high_event] <- tbc_sample[high_event] * high[high_event] / cfu_g[high_event]
    out[low_event]  <- tbc_sample[low_event]  * low[low_event] / cfu_g[low_event]
    values[(i-1)*80 + 1:80] <- out
    if (i %% 100 == 0)
      cat("simulating",i,"of", n_sim/80, "\n")
  }

  as.numeric(table(cut(values, c(-1,0,10^(0:6),Inf))))
}

# Righto, read in the parameters for each data set
source("campy.R")
source("ecoli.R")
source("listeria.R")
source("salmonella.R")

pathogens <- list(campy_params(),
                  ecoli_params(),
                  listeria_params(),
                  salmonella_params())

# And read in the TBC data
TBC_post <- load_tbc()

# set random seed for simulations for repeatability
set.seed(3)

# maximum sensitivity (number bugs in 25mL) to simulate
max_sens <- 5

# counts per litre
counts <- matrix(0, 4, 9);
rownames(counts) <- 1:4

for (j in seq_along(pathogens)) {

  pathogen <- pathogens[[j]]

  # load pathogen data
  pos      <- load_pos(pathogen$data_name)
  
  s <- matrix(0, 10000+1, max_sens)
  for (i in 1:max_sens) {
    cat("Processing", pathogen$data_name, "at sensitivity", i, "of", max_sens, "\n")
    s[1:10000,i] <- sample_survey(pos$total,
                                  TBC_post,
                                  pathogen$farm_prev,
                                  pathogen$cow_prev,
                                  pathogen$high_dist,
                                  pathogen$low_dist,
                                  i) / sum(pos$total)
    
  }
  # final value is actual prevalence in survey
  s[10001,] <- sum(pos$count) / sum(pos$total)
  
  # write data out
  write.csv(s, file.path("output", pathogen$out_file), row.names=F)

  # sample the pathogen counts
  counts[j,] <- sample_counts(TBC_post,
                              pathogen$farm_prev,
                              pathogen$cow_prev,
                              pathogen$high_dist,
                              pathogen$low_dist)
  rownames(counts)[j] <- pathogen$data_name
}

write.csv(counts, file.path("output", "pathogen_counts.csv"), row.names=TRUE)
