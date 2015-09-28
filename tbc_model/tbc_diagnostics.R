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

out_folder  <- "output"
diag_folder <- "diagnostics"

# now do plots etc.
accept     <- read.csv(file.path(out_folder, "acceptance.csv"))
mu_a_post  <- read.csv(file.path(out_folder, "mu_a.csv"))
mu_b_post  <- read.csv(file.path(out_folder, "mu_b.csv"))
tau_a_post <- read.csv(file.path(out_folder, "tau_a.csv"))
tau_b_post <- read.csv(file.path(out_folder, "tau_b.csv"))
w_post     <- read.csv(file.path(out_folder, "w.csv"))
alpha_post <- read.csv(file.path(out_folder, "alpha.csv"))
beta_post  <- read.csv(file.path(out_folder, "beta.csv"))
x_post     <- read.csv(file.path(out_folder, "x.csv"))
farm_map   <- read.csv(file.path(out_folder, "farm_map.csv"))

n_farms    <- ncol(alpha_post)
n_samples  <- nrow(alpha_post)

# hyperpriors
w_alpha    <- 1      # hyperprior for w_j
w_beta     <- 1      # hyperprior for w_j
tau_alpha  <- 0.001  # prior for precision of \alpha_j
tau_beta   <- 0.001  # prior for precision of \beta_j
tau_mu     <- 0.001  # hyperprior, precision of mu_\alpha, \mu_\beta

# convert to matrices/numeric
mu_a_post  <- mu_a_post[,1]
mu_b_post  <- mu_b_post[,1]
tau_a_post <- tau_a_post[,1]
tau_b_post <- tau_b_post[,1]
w_post     <- w_post[,1]
alpha_post <- as.matrix(alpha_post)
beta_post  <- as.matrix(beta_post)
x_post     <- as.matrix(x_post)

# create output folder if necessary
dir.create(diag_folder)

# MCMC traces
pdf(file.path(diag_folder, "traces.pdf"), paper="a4", width=8, height=11)
par(mfrow=c(4,2))
plot(mu_a_post, type="l", xlab="iterations", ylab="", main=expression(mu[alpha]))
acf(mu_a_post)
plot(mu_b_post, type="l", xlab="iterations", ylab="", main=expression(mu[beta]))
acf(mu_b_post)
plot(tau_a_post, type="l", xlab="iterations", ylab="", main=expression(tau[alpha]))
acf(tau_a_post)
plot(tau_b_post, type="l", xlab="iterations", ylab="", main=expression(tau[beta]))
acf(tau_b_post)
plot(w_post, type="l", xlab="iterations", ylab="", main=expression(w))
acf(w_post)
hist(accept$alpha, xlab="Acceptance", ylab="", main=substitute(paste("Acceptance rate ",alpha)))
hist(accept$beta, xlab="Acceptance", ylab="", main=substitute(paste("Acceptance rate ",beta)))
hist(accept$swap, xlab="Acceptance", ylab="", main=substitute(paste("Acceptance rate ",beta, "<->", alpha)))

# plot pair-wise scatters
par(mfrow=c(5,2))
plot(mu_b_post ~ mu_a_post, xlab=expression(mu[alpha]), ylab=expression(mu[beta]))
plot(mu_b_post ~ tau_a_post, xlab=expression(tau[alpha]), ylab=expression(mu[beta]))
plot(mu_b_post ~ tau_b_post, xlab=expression(tau[beta]), ylab=expression(mu[beta]))
plot(mu_b_post ~ w_post, xlab=expression(w), ylab=expression(mu[beta]))
plot(mu_a_post ~ tau_a_post, xlab=expression(tau[alpha]), ylab=expression(mu[alpha]))
plot(mu_a_post ~ tau_b_post, xlab=expression(tau[beta]), ylab=expression(mu[alpha]))
plot(mu_a_post ~ w_post, xlab=expression(w), ylab=expression(mu[alpha]))
plot(tau_b_post ~ tau_a_post, xlab=expression(tau[alpha]), ylab=expression(tau[beta]))
plot(tau_b_post ~ w_post, xlab=expression(w), ylab=expression(tau[beta]))
plot(tau_a_post ~ w_post, xlab=expression(w), ylab=expression(tau[alpha]))

par(mfrow=c(5,4))
for (i in 1:n_farms)
{
  plot(alpha_post[,i], type="l", xlab="iterations", ylab="", main=substitute(alpha[nn], list(nn=i)))
  plot(beta_post[,i], type="l", xlab="iterations", ylab="", main=substitute(beta[nn], list(nn=i)))
  plot(rowSums(x_post[,farm_map[i,1]:farm_map[i,2]])/(farm_map[i,2] - farm_map[i,1]+1), type="l", xlab="iterations", ylim=c(0,1), ylab="", main=substitute(x[nn], list(nn=i)))
  plot(alpha_post[,i], beta_post[,i], xlab=substitute(alpha[nn], list(nn=i)), ylab=substitute(beta[nn], list(nn=i)))
}

dev.off()

pdf(file.path(diag_folder, "density.pdf"), paper="a4", width=8, height=11)
par(mfrow=c(3,2))

plot(density(mu_a_post), xlab="iterations", ylab="", main=expression(mu[alpha]))
prior_x <- seq(min(mu_a_post), max(mu_a_post), length.out=100)
prior_y <- dnorm(prior_x, 0, 1/sqrt(tau_mu)) # TODO: rescale area?
lines(prior_x, prior_y, col="grey60")

plot(density(mu_b_post), xlab="iterations", ylab="", main=expression(mu[beta]))
prior_x <- seq(min(mu_b_post), max(mu_b_post), length.out=100)
prior_y <- dnorm(prior_x, 0, 1/sqrt(tau_mu)) # TODO: rescale area?
lines(prior_x, prior_y, col="grey60")

plot(density(tau_a_post), xlab="iterations", ylab="", main=expression(tau[alpha]))
prior_x <- seq(min(tau_a_post), max(tau_a_post), length.out=100)
prior_y <- dgamma(prior_x, tau_alpha, tau_beta) # TODO: rescale area?
lines(prior_x, prior_y, col="grey60")

plot(density(tau_b_post), xlab="iterations", ylab="", main=expression(tau[beta]))
prior_x <- seq(min(tau_b_post), max(tau_b_post), length.out=100)
prior_y <- dgamma(prior_x, tau_alpha, tau_beta) # TODO: rescale area?
lines(prior_x, prior_y, col="grey60")

w_range <- range(w_post)
plot(density(w_post), xlim=w_range, xlab="iterations", ylab="", main=expression(w))
prior_x <- seq(min(w_post), max(w_post), length.out=100)
prior_y <- dbeta(prior_x, w_alpha, w_beta)
lines(prior_x, prior_y, col="grey60")

par(mfrow=c(5,3))
alpha_range <- range(alpha_post)
beta_range <- range(beta_post)
for (i in 1:n_farms)
{
  plot(density(alpha_post[,i]), xlim=alpha_range, ylab="", main=substitute(alpha[nn], list(nn=i)))
  plot(density(beta_post[,i]), xlim=beta_range, ylab="", main=substitute(beta[nn], list(nn=i)))
  barplot(table(x_post[,farm_map[i,1]:farm_map[i,2]]), main=substitute(x[nn], list(nn=i)))
}
dev.off()




# random effects




l1 <- apply(alpha_post, 2, quantile, 0.025)
m1 <- apply(alpha_post, 2, mean)
h1 <- apply(alpha_post, 2, quantile, 0.975)

l2 <- apply(alpha_post+beta_post, 2, quantile, 0.025)
m2 <- apply(alpha_post+beta_post, 2, mean)
h2 <- apply(alpha_post+beta_post, 2, quantile, 0.975)

# order by m
o <- order(m2)
l1 <- l1[o]
m1 <- m1[o]
h1 <- h1[o]
l2 <- l2[o]
m2 <- m2[o]
h2 <- h2[o]

pdf(file.path(diag_folder, "tbc_ranef.pdf"), width=8, height=6)
par(mai=c(0.5, 0.4, 0.25, 0.2))
plot(NULL, xlim=c(0,3.2*log(10)), ylim=c(1,n_farms), xlab="", ylab="", xaxt="n", yaxt="n")
for (i in 1:n_farms)
{
  lines(c(l2[i],h2[i]), c(i,i), col="red")
  lines(c(l1[i],h1[i]), c(i,i))
  points(m2[i],i, pch=19, col="red", cex=0.5)
  points(m1[i],i, pch=19, cex=0.5)
}
axis(1, at=(-2:4)*log(10), labels=10^(-2:4))
lines(rep(mean(mu_a_post),2),c(0,n_farms+1), col="black")
lines(rep(mean(mu_b_post+mu_a_post),2),c(0,n_farms+1), col="red")
dev.off()

# posterior for mu and tau
pdf(file.path(diag_folder, "tbc_posterior_mu.pdf"), width=6, height=2)
par(mai=c(0.5, 0.4, 0.25, 0.2))
plot(density(log10(exp(mu_a_post))), xlim=c(0.6,1.6), xlab="", ylab="", main="", cex.axis=0.7)
lines(density(log10(exp(mu_a_post + mu_b_post))), col="red")
text(1.55,22, expression(paste(alpha, ", ")), adj=1)
text(1.55,22, expression(beta), col="red", adj=0)
dev.off()

pdf(file.path(diag_folder, "tbc_posterior_tau.pdf"), width=6, height=2)
par(mfrow=c(1,3), mai=c(0.5, 0.4, 0.25, 0.2), omi=c(0,0,0,0))
plot(density(log10(exp(tau_a_post))), xlab="", ylab="", main="")
text(11.1, 0.3, expression(tau[alpha]), cex=1.2)
plot(density(log10(exp(tau_b_post))), xlab="", ylab="", main="", col="red")
text(0.6, 5.5, expression(tau[beta]), cex=1.2, col="red")
plot(density(w_post), xlab="", ylab="", main="")
text(0.165, 22, expression(w), cex=1.2)
dev.off()




# predictive density.  Idea is to average the density at each sample
#plot(h$mids, log10(h$counts/sum(h$counts)/0.2))
#hist_area <- sum(log10(h$counts/sum(h$counts)/0.2)*0.2)
x <- unique(round(10^seq(0,4,length.out=500)))
total_pred <- rep(0, length(x))
for (i in 1:n_samples)
{
  lambda_i <- alpha_post[i,farm] + x_post[i,]*(beta_post[i,farm])
  pred <- rep(0,length(x))
  for (j in 1:length(x))
    pred[j] <- mean(dpois(x[j], exp(lambda_i)))
  total_pred <- total_pred + pred;
  if (i %% 10 == 0)
    cat("up to ", i, "\n")
}
areas_raw <- diff(x)*(total_pred[-1]+total_pred[-length(x)])
areas_log <- diff(log10(x))*(total_pred[-1]+total_pred[-length(x)])
total_pred[-length(x)] <- total_pred[-length(x)] * areas_raw / areas_log
total_pred[!is.finite(total_pred)] <- 0
scale <- 0.5*sum(diff(log10(x))*(total_pred[-1]+total_pred[-length(x)]))

pdf("tbc_pred_overall.pdf", width=6, height=4)
hist(log10(y), freq=F, main="", xaxt="n", xlab="Total bacterial count per uL")
lines(log10(x), total_pred/scale, col="red")
axis(side=1, at=0:4, labels=10^(0:4))
dev.off()


farm_no <- 13

x <- unique(round(10^seq(0,4,length.out=500)))
areas_raw <- diff(x)*(total_pred[-1]+total_pred[-length(x)])
areas_log <- diff(log10(x))*(total_pred[-1]+total_pred[-length(x)])
total_pred[-length(x)] <- total_pred[-length(x)] * areas_raw / areas_log
total_pred[!is.finite(total_pred)] <- 0
scale <- 0.5*sum(diff(log10(x))*(total_pred[-1]+total_pred[-length(x)]))

# this draws the fits for 4 farms with slightly different distributions as an example.
farms <- c(13, 61, 26, 40)
pdf("tbc_pred_farms.pdf", width=6, height=4)
par(mfrow=c(2,2), mai=c(0.5, 0.5, 0.5, 0.5))
for (farm_no in farms)
{
  total_pred <- rep(0, length(x))
  for (i in 1:n_samples)
  {
    lambda_i <- alpha_post[i,farm_no] + x_post[i,farm==farm_no]*(beta_post[i,farm_no])
    pred <- rep(0,length(x))
    for (j in 1:length(x))
      pred[j] <- mean(dpois(x[j], exp(lambda_i)))
    total_pred <- total_pred + pred;
    if (i %% 10 == 0)
      cat("up to ", i, "\n")
  }
  areas_raw <- diff(x)*(total_pred[-1]+total_pred[-length(x)])
  areas_log <- diff(log10(x))*(total_pred[-1]+total_pred[-length(x)])
  total_pred[-length(x)] <- total_pred[-length(x)] * areas_raw / areas_log
  total_pred[!is.finite(total_pred)] <- 0
  scale <- 0.5*sum(diff(log10(x))*(total_pred[-1]+total_pred[-length(x)]))
  h <- hist(log10(y[farm == farm_no]), freq=F, xlim=c(0,4), plot=F)
  hist(log10(y[farm == farm_no]), ylim=range(h$density, total_pred/scale), freq=F, xlim=c(0,4), xaxt="n", xlab="Total bacterial count per uL", ylab="", main="")
  lines(log10(x), total_pred/scale, col="red")
  axis(side=1, at=0:4, labels=10^(0:4))
}
dev.off()






# see how closely a normal fits alpha...
pdf("alpha_post.pdf", paper="a4", width=8, height=11)
am <- matrix(0,n_farms,3)
par(mfrow=c(5,2))
for (i in 1:n_farms)
{
  d = density(alpha_post[,i], adjust=1.3);
  plot(d);
  lines(d$x, dnorm(d$x, mean(alpha_post[,i]), sd(alpha_post[,i])))
}
dev.off()

pdf("beta_post.pdf", paper="a4", width=8, height=11)
par(mfrow=c(5,2))
bm <- matrix(0,n_farms,6)
for (i in 1:n_farms)
{
  d = density(beta_post[,i], adjust=1.3);
  b = beta_post[,i]
  nm = normalmixEM(beta_post[,i], maxit=20000)
  #  b = b[b > quantile(b,0.025) & b < quantile(b, 0.975)]
  plot(d$x, nm$lambda[1] * dnorm(d$x, nm$mu[1], nm$sigma[1]) + nm$lambda[2] * dnorm(d$x, nm$mu[2], nm$sigma[2]), col="red", type="l", ylab="")
  lines(d);
  bm[i,] <- c(nm$lambda, nm$mu, nm$sigma)
  if (min(nm$lambda*nm$sigma) < 0.001)
  { # tiny sd and weight - try fitting a 1 component model
    lines(d$x, dnorm(d$x, mean(b), sd(b)), col="green", type="l")
    bm[i,] <- c(1,0,mean(b),0,sd(b),1)
  }
  title(paste(i, ":", min(nm$lambda*nm$sigma)))
  cat("farm", i, "variance=", nm$sigma^2, "lambda=", nm$lambda, "\n")
}
dev.off()

library(mixtools)

pdf("alphabeta_post1.pdf", paper="a4", width=8, height=11)
par(mfrow=c(5,2))
abm <- matrix(0,n_farms,6)
for (i in 1:n_farms)
{
  b <- as.numeric(alpha_post[,i] + x_post[,farm==i]*beta_post[,i])
  d = density(b, adjust=1.3);
  #  b = alpha_post[,i] + x_post[,i]*beta_post[,i]
  nm = normalmixEM(b, k=2, maxit=20000, maxrestarts=100)
  #  b = b[b > quantile(b,0.025) & b < quantile(b, 0.975)]
  plot(d$x, nm$lambda[1] * dnorm(d$x, nm$mu[1], nm$sigma[1]) + nm$lambda[2] * dnorm(d$x, nm$mu[2], nm$sigma[2]), col="red", type="l", ylab="")
  lines(d);
  abm[i,] <- c(nm$lambda, nm$mu, nm$sigma)
  if (min(nm$lambda*nm$sigma) < 0.001)
  { # tiny sd and weight - try fitting a 1 component model
    lines(d$x, dnorm(d$x, mean(b), sd(b)), col="green", type="l")
    abm[i,] <- c(1,0,mean(b),0,sd(b),1)
  }
  title(paste(i, ":", min(nm$lambda*nm$sigma)))
  cat("farm", i, "variance=", nm$sigma^2, "lambda=", nm$lambda, "\n")
}
dev.off()
colnames(abm) <- c("w1","w2","m1","m2","s1","s2")
write.csv(abm, "TBC_post.csv", row.names=F)

pdf("alphabeta_post2.pdf", paper="a4", width=8, height=11)
par(mfrow=c(5,2))
for (i in 1:n_farms)
{
  b <- as.numeric(alpha_post[,i] + x_post[,farm==i]*beta_post[,i])
  d = density(b, adjust=1.3);
  plot(d);
  lines(d$x, abm[i,1] * dnorm(d$x, abm[i,3], abm[i,5]) + abm[i,2] * dnorm(d$x, abm[i,4], abm[i,6]), col="red", type="l", ylab="")
}
dev.off()
