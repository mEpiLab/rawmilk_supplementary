# listeria params...
# ecoli params...

listeria_params <- function() {

  # animal-level prevalence
  cow_prev <- function() {
    prev <- c(rep(44/182,182), rep(0.43,1500)) # Esteban, Mohammed  (Mohammed also found 19% prevalence on teat swabs)
    list(mean = mean(prev), var = var(prev))
  }

  # distribution functions
  high_dist <- function(n, data)
  {
  #  mu_high <- data[1]
    mu_high <- 3.5
    sd_high <- 0.5
    10^rnorm(n, mean=mu_high, sd=sd_high)
  }

  low_dist <- function(n, data)
  {
  #  mu_high <- data[1]
  #  sd_high <- 1 # see above
    mu_high <- 3.5
    sd_high <- 0.5

    mu <- log(10)*mu_high
    sigma <- log(10)*sd_high
    mu_logn <- exp(mu + sigma^2/2)
    var_logn <- (exp(sigma^2)-1)*exp(2*mu+sigma^2)

    # combine cow prevalence
    prev <- cow_prev()
    mean <- prev$mean * mu_logn
    var  <- prev$mean^2*var_logn + mu_logn^2*prev$var + prev$var*var_logn

    r <- rnorm(n, mean=mean, sd=sqrt(var/400))
    r[r < 0] <- 0
    return(r)
  }

  farm_prev_dist <- function(n, data)
  {
    herd_prev <- c(rep(38/82,82), rep(50/50,50))
    m <- mean(herd_prev)
    v <- var(herd_prev)
    rbeta(n, shape1=(m*(1-m)/v-1)*m, shape2=(m*(1-m)/v-1)*(1-m))
  }

  cow_prev_dist <- function(n, data)
  {
    prev <- cow_prev()
    m <- prev$mean
    v <- prev$var
    cow_a <- (m*(1-m)/v-1)*m
    cow_b <- (m*(1-m)/v-1)*(1-m)
    rbeta(n, shape1=cow_a, shape2=cow_b)
  }

  return(list(data_name = "ListMono",
              cow_prev  = cow_prev_dist,
              farm_prev = farm_prev_dist,
              high_dist = high_dist,
              low_dist  = low_dist,
              out_file  = "listeria.csv"))
}
