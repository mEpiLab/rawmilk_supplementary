# campy params...

campy_params <- function() {
  cow_prev_dist <- function(n, data) {
    0.8417*0.7476
  }

  farm_prev_dist <- function(n, data) {
    0.95
  }

  high_dist_mu <- function(n, data) {
    mu_high <- 2.865612
    sd_high <- sqrt(2.035025)
    10^rnorm(n, mean=mu_high, sd=sd_high)
    # truncate at 10^5.5
    #  h[h > 5.5] <- 5.5 - (h[h > 5.5]-5.5)  # reflect in h=5.5 
    #  return(10^h)
  }
  low_dist_mu <- function(n, data) {
    mu_high <- 2.865612
    sd_high <- sqrt(2.035025)
    mu <- log(10)*mu_high
    sigma <- log(10)*sd_high
    #  mu_logn <- exp(mu + sigma^2/2)
    #  var_logn <- (exp(sigma^2)-1)*exp(2*mu+sigma^2)
    mu_logn <- 16717.25
    var_logn <- (45598.41) ^2  # computed from high_dist_mu truncated dist

    p_cow  <- cow_prev_dist()

    mean <- p_cow*mu_logn
    var  <- p_cow^2*var_logn + mu_logn^2*p_cow*(1-p_cow) + p_cow*(1-p_cow)*var_logn

    r <- rnorm(n, mean=mean, sd=sqrt(var/400))
    r[r < 0] <- 0
    return(r)
  }

  return(list(data_name = "Campy",
              cow_prev  = cow_prev_dist,
              farm_prev = farm_prev_dist,
              high_dist = high_dist_mu,
              low_dist  = low_dist_mu,
              out_file  = "campy.csv"))
}
