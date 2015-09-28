# Figure 1: comparison of farm sizes to NZ data
f <- c(10, 22 ,19 , 6, 5, 3, 6 , 5 , 0,  0,  0 , 1)
freq <- c(2.4, 7, 12.2, 14.3, 11, 10.7, 7, 6.4, 4.8, 4.5, 3.4, 3.2, 2.1, 1.8, 1.6, 1.5, 1.0, 0.9, 0.7, 1.2, 0.7, 1.2, 0.7)
size <- c(seq(100,1000,by=50)-25, 1050, 1150, 1350, 1700)
nz <- data.frame(size=size, freq=freq)

pdf("figure1.pdf", width=7, height=5)
barplot(f, space=0, col="white", xlab="Number of cows", ylab="Number of Farms")
d <- density(nz$size, weights=nz$freq/sum(nz$freq), adjust=0.5)
lines((d$x-50)/100, 9000*d$y, col="black", lwd=2)
axis(1, at=c(0:12), labels=(0:12)*100, cex.axis=0.85)
mtext("Fig 1", side=3, at=11.50, adj=0, line=1)
dev.off()

# Figure 2: Posterior distributions of TBC model parameters

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

out_folder  <- "../tbc_model/output"

# now do plots etc.
mu_a_post  <- read.csv(file.path(out_folder, "mu_a.csv"))
mu_b_post  <- read.csv(file.path(out_folder, "mu_b.csv"))
tau_a_post <- read.csv(file.path(out_folder, "tau_a.csv"))
tau_b_post <- read.csv(file.path(out_folder, "tau_b.csv"))
w_post     <- read.csv(file.path(out_folder, "w.csv"))
alpha_post <- read.csv(file.path(out_folder, "alpha.csv"))
beta_post  <- read.csv(file.path(out_folder, "beta.csv"))
x_post     <- read.csv(file.path(out_folder, "x.csv"))

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

pdf("figure2.pdf", width=6.3, height=4.4)

split.screen(c(2, 1))       # split display into two screens
split.screen(c(1, 3), screen = 2) # now split the bottom half into 3
screen(1) # prepare screen 1 for output
par(mai=c(0.8, 0.8, 0.25, 0.2), omi=c(0,0,0,0))
plot(density(log10(exp(mu_a_post))), xlim=c(0.6,1.6), xlab=expression(log[10] * "(TBC)/" * mu *"L"), ylab="Density", main="", xaxt="n")
axis(1, at=seq(0,5,by=0.2))
lines(density(log10(exp(mu_a_post + mu_b_post))), col="black", lty="dashed", lwd=1.5)
text(log10(6.5),22, expression(alpha))
text(log10(21),5, expression(beta), col="black")
mtext("Fig 2", side=3, at=1.6, adj=0, line=0.2)

screen(3) # bottom left
par(mai=c(0.8, 0.8, 0.25, 0.2), omi=c(0,0,0,0))
plot(density(log10(exp(tau_a_post))), xlab=expression(tau[alpha]), ylab="Density", main="")
screen(4) # bottom middle
par(mai=c(0.8, 0.8, 0.25, 0.2), omi=c(0,0,0,0))
plot(density(log10(exp(tau_b_post))), xlab=expression(tau[beta]), ylab="", main="", lty="dashed", lwd=1.5)
screen(5) # bottom right
par(mai=c(0.8, 0.8, 0.25, 0.2), omi=c(0,0,0,0))
plot(density(w_post), xlab=expression(w), ylab="", main="")

close.screen(all.screens=TRUE)
dev.off()


# Figure 3: Random effects plots for alpha and beta (log mean low and high counts)
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

pdf("figure3.pdf", width=8, height=6)
par(mai=c(0.8, 0.8, 0.5, 0.2))
plot(NULL, xlim=c(0,3.2*log(10)), ylim=c(1,n_farms), xlab=expression("TBC per " * mu * "L"), ylab="Farms", xaxt="n", yaxt="n")
for (i in 1:n_farms)
{
  lines(c(l2[i],h2[i]), c(i,i), col="grey50")
  lines(c(l1[i],h1[i]), c(i,i))
  points(m2[i],i, pch=22, col="grey50", cex=0.5, bg="grey50")
  points(m1[i],i, pch=19, cex=0.5)
}
axis(1, at=(-2:4)*log(10), labels=10^(-2:4))
lines(rep(mean(mu_a_post),2),c(0,n_farms+1), col="black")
lines(rep(mean(mu_b_post+mu_a_post),2),c(0,n_farms+1), col="grey50", lwd=1.5)
mtext("Fig 3", side=3, at=3*log(10), adj=0, line=1, cex=1.3)
dev.off()


# Figure 4: Predictive densities overall and for four of the farms

farms <- c(13, 61, 26, 40)

x <- unique(round(10^seq(0,4,length.out=500)))
total_pred <- matrix(0, length(x), 1 + length(farms))
for (i in 1:n_samples)
{
  lambda_i <- alpha_post[i,farm] + x_post[i,]*(beta_post[i,farm])
  pred <- rep(0,length(x))
  for (j in 1:length(x))
    pred[j] <- mean(dpois(x[j], exp(lambda_i)))
  total_pred[,1] <- total_pred[,1] + pred;
  if (i %% 10 == 0)
    cat("up to ", i, "\n")
}
areas_raw <- diff(x)*(total_pred[-1,1]+total_pred[-length(x),1])
areas_log <- diff(log10(x))*(total_pred[-1,1]+total_pred[-length(x),1])
total_pred[-length(x),1] <- total_pred[-length(x),1] * areas_raw / areas_log
total_pred[!is.finite(total_pred[,1]),1] <- 0
scale <- 0.5*sum(diff(log10(x))*(total_pred[-1,1]+total_pred[-length(x),1]))

pdf("figure4.pdf", width=6, height=8)

split.screen(c(2, 1))             # split display into two screens
split.screen(c(2, 2), screen = 2) # now split the bottom half into 4
screen(1) # prepare screen 1 for output
par(mai=c(0.8, 0.8, 0.5, 0.2), omi=c(0,0,0,0), cex = 0.7)

hist(log10(y), freq=F, main="", xaxt="n", xlab=expression("Total bacterial count per " * mu * "L"))
lines(log10(x), total_pred[,1]/scale, col="black", lwd=2)
axis(side=1, at=0:4, labels=10^(0:4))
mtext("Fig 4", side=3, at=3.1, adj=0, line=1, cex=2)

# this draws the fits for some farms with slightly different distributions as an example.
for (farm_it in 1:length(farms))
{
  farm_no <- farms[farm_it]
  for (i in 1:n_samples)
  {
    lambda_i <- alpha_post[i,farm_no] + x_post[i,farm==farm_no]*(beta_post[i,farm_no])
    pred <- rep(0,length(x))
    for (j in 1:length(x))
      pred[j] <- mean(dpois(x[j], exp(lambda_i)))
    total_pred[,1+farm_it] <- total_pred[,1+farm_it] + pred;
    if (i %% 10 == 0)
      cat("up to ", i, "\n")
  }
  areas_raw <- diff(x)*(total_pred[-1,1+farm_it]+total_pred[-length(x),1+farm_it])
  areas_log <- diff(log10(x))*(total_pred[-1,1+farm_it]+total_pred[-length(x),1+farm_it])
  total_pred[-length(x),1+farm_it] <- total_pred[-length(x),1+farm_it] * areas_raw / areas_log
  total_pred[!is.finite(total_pred[,1+farm_it]),1+farm_it] <- 0
  scale <- 0.5*sum(diff(log10(x))*(total_pred[-1,1+farm_it]+total_pred[-length(x),1+farm_it]))
  
  screen(2 + farm_it)
  par(mai=c(0.8, 0.8, 0.25, 0.2), omi=c(0,0,0,0), cex = 0.7)
  h <- hist(log10(y[farm == farm_no]), freq=F, xlim=c(0,4), plot=F)
  hist(log10(y[farm == farm_no]), ylim=range(h$density, total_pred[,1+farm_it]/scale), freq=F, xlim=c(0,4), xaxt="n", xlab=expression("Total bacterial count per " * mu * "L"), ylab="", main="")
  lines(log10(x), total_pred[,1+farm_it]/scale, col="black", lwd=1.5)
  axis(side=1, at=0:4, labels=10^(0:4))
}
close.screen(all.screens=TRUE)
dev.off()


# Figure 5: Predictive pathogen counts

sens_folder <- "../sensitivity_model/output"
d <- read.table(file.path(sens_folder, "pathogen_counts.txt"))

pdf("figure5.pdf", width=10, height=8)
par(mfrow=c(2,2))
par(omi=rep(0.25,4), mai=c(0.9,0.9,0.5,0.5))
for (i in 1:4)
{
  plot(NULL, xlim=c(0,6), ylim=c(0,90), xlab="Bacteria per L", ylab="Frequency", yaxs="i", xaxt="n")
  barplot(as.numeric(d[i,2:7])/1000, space=0,add=T,yaxt="n")
  axis(side=1,at=1:6,labels=c(0,10^(0:4)))
  if (i == 2)
    mtext("Fig 5", side=3, at=6.5, adj=0, line=2)
}
dev.off()

s <- list()
s[[1]] <- read.csv(file.path(sens_folder, "campy_sens.csv"))
s[[2]] <- read.csv(file.path(sens_folder, "ecoli_sens.csv"))
s[[3]] <- read.csv(file.path(sens_folder, "salmonella_sens.csv"))
s[[4]] <- read.csv(file.path(sens_folder, "listeria_sens.csv"))
max_sens <- 0
max_prev <- 0
for (i in 1:4)
{
  max_sens <- max(max_sens, ncol(s[[i]]))
  max_prev <- max(max_prev, max(s[[i]][-nrow(s[[i]]),]))
}

edge_density <- function(x, h)
{
	x0 <- 0
	d <- density(x, h)
	d$int <- d$x > x0
	d$x <- d$x[d$int]
	d$y <- d$y[d$int]
	d$corr <- 1 - pnorm(x0, d$x, d$bw)
	d$y <- d$y / d$corr
	return(d)
}

my_vioplot <- function(dat, bw, border, col, at)
{
  e <- edge_density(dat, bw)
  m <- 0.3/max(e$y)
  incl <- e$y > max(e$y)/100
  polygon(c(at-e$y[incl]*m,rev(at+e$y[incl]*m)), c(e$x[incl], rev(e$x[incl])), col=col, border=border)
}

# Figure 6: Expected prevalence given sensitivity of detection in 25mL
pdf("figure6.pdf", width=10, height=8)
par(mfrow=c(2,2))
par(omi=c(0.25,0.25,0.25,0.25), mai=c(0.9,0.9,0.5,0.5))
for (j in 1:4)
{
  plot(NULL, xlim=c(0.5,max_sens+0.5), ylim=c(0,0.05*1.05), ylab="Prevalence", xlab="Detection limit (per 25mL)", xaxt="n", yaxs="i")
  for (i in 1:max_sens)
    my_vioplot(s[[j]][-nrow(s[[j]]),i], bw=0.0015, border="grey60", col="grey60", at=i)
  axis(1, at=0:(max_sens+1), labels=0:(max_sens+1))
  abline(h=s[[j]][nrow(s[[j]]),1], col="black", lwd=2)
  if (j == 2)
    mtext("Fig 6", side=3, at=max_sens + 1, adj=0, line=2)
}
dev.off()

