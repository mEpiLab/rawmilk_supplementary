# ecoli params...

ecoli_params <- function() {

  high_dist <- function(n, data)
  {
    mu_high <- 2.032588
    sd_high <- 0.5624717
  #  mu_high <- data[1]
  #  sd_high <- 1
    10^rnorm(n, mean=mu_high, sd=sd_high)
  }

  low_dist <- function(n, data)
  {
    mu_high <- 2.032588
    sd_high <- 0.5624717
  #  mu_high <- data[1]
  #  sd_high <- 1
    mu <- mu_high*log(10)
    sigma <- sd_high*log(10)
    mu_logn <- exp(mu + sigma^2/2)
    var_logn <- (exp(sigma^2)-1)*exp(2*mu+sigma^2)

    p_cow  <- cow_prev_dist()

    mu_low <- p_cow * mu_logn
    sd_low <- sqrt((p_cow^2*var_logn + mu_logn^2*p_cow*(1-p_cow) + p_cow*(1-p_cow)*var_logn)/400)
    r <- rnorm(n, mean=mu_low, sd=sd_low)
    r[r < 0] <- 0
    return(r)
  }

  cow_prev_dist <- function(n, data)
  {
    16/24
  }

  farm_prev_dist <- function(n, data)
  {
    16/136
  }

  return(list(data_name = "Ecoli",
              cow_prev  = cow_prev_dist,
              farm_prev = farm_prev_dist,
              high_dist = high_dist,
              low_dist  = low_dist,
              out_file  = "ecoli.csv"))
}
