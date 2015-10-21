# salmonella params...

salmonella_params <- function() {

  high_dist <- function(n, data)
  {
    mu_high <- 0.89
    sd_high <- 0.83666
    10^rnorm(n, mu_high, sd_high)
  }

  low_dist <- function(n, data)
  {
    mu_high <- 0.89
    sd_high <- 0.83666

    mu <- log(10)*mu_high
    sigma <- log(10)*sd_high
    mu_logn <- exp(mu + sigma^2/2)
    var_logn <- (exp(sigma^2)-1)*exp(2*mu+sigma^2)

    p_cow <- 0.4610
    mean <- p_cow*mu_logn
    var  <- p_cow^2*var_logn + mu_logn^2*p_cow*(1-p_cow) + p_cow*(1-p_cow)*var_logn

    r <- rnorm(n, mean=mean, sd=sqrt(var/400))
    r[r < 0] <- 0
    return(r)
  }

  cow_prev_dist <- function(n, data)
  {
    # fit beta distribution to prevalence
    prev <- data.frame(farm  = c(rep(4,5),rep(5,4),rep(6,4)),
                       pos   = c(112,27,94,45,54,66,105,78,88,59,9,3,7),
                       total = c(179,154,145,121,133,107,149,109,111,97,103,101,106))

    p <- NULL
    for (i in 1:nrow(prev)) {
      p <- c(p, rep(prev$pos[i] / prev$total[i], prev$total[i]))
    }

    m <- mean(p)
    v <- var(p)

    rbeta(n, shape1=(m*(1-m)/v-1)*m, shape2=(m*(1-m)/v-1)*(1-m))
  }

  farm_prev_dist <- function(n, data)
  {
    0.438
  }

  return(list(data_name = "Salmonella",
              cow_prev  = cow_prev_dist,
              farm_prev = farm_prev_dist,
              high_dist = high_dist,
              low_dist  = low_dist,
              out_file  = "salmonella.csv"))
}

