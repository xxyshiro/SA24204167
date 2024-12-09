## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
f <- function(n) { cat((1:n)*n, sep = " ", "\n") }
a <- sapply(1:9, f)

## -----------------------------------------------------------------------------
set.seed(101)

m = 500; n = 1000

X <- runif(n, 0, 1)
x1 <- X[1:m]; x2 <- X[(m+1):n]

n1 <- rnorm(n, 0, 1)
n2 <- rnorm(n, 3, 1)

y1 <- (x1<=.5)*n1[1:m] + (x1>.5)*n2[1:m]
y2 <- (x2<=.7)*n1[(m+1):n] + (x2>.7)*n2[(m+1):n]

head(as.data.frame(cbind(y1, y2)))

hist(y1, main = "Sample of Y: Group 1", xlab = "Y")
hist(y2, main = "Sample of Y: Group 2", xlab = "Y")

## -----------------------------------------------------------------------------
library(MASS)

set.seed(123)

n = 200

sigma <- diag(3)
sigma[2, 1] <- sigma[1, 2] <- sigma[3, 2] <- sigma[2, 3] <- -.5
sigma[1, 3] <- sigma[3, 1] <- .5 ## 生成协方差矩阵

x <- mvrnorm(n, mu = rep(0, 3), Sigma = sigma)

head(as.data.frame(x))

## ----warning=FALSE------------------------------------------------------------
pairs(x)

## -----------------------------------------------------------------------------
gen.Rayleigh <- function(n, sigma = 1){
  u = runif(n)
  x = sigma*sqrt(2*log(1/(1-u)))
  return(x)
}

## -----------------------------------------------------------------------------
## Generate Rayleigh sample with different sigma
set.seed(124)
n = 1000
sigma.all = c(1, 5, 10)

x <- array(NA, dim = c(3, n))

for (i in 1:3){
  sigma = sigma.all[i]
  x[i, ] <- gen.Rayleigh(n, sigma)
}

## ----message=FALSE, warning=FALSE---------------------------------------------
## library
library(tidyverse)
library(ggplot2)

## -----------------------------------------------------------------------------
## Build a dataset containing all samples, with their groups labeled
data = data.frame(Group = factor(x = c(rep(1, n), rep(2, n), rep(3, n)),
                                 levels = c(1, 2, 3),
                                 labels = c("Rayleigh(1)", "Rayleigh(5)", "Rayleigh(10)")),
                  value = c(x[1, ], x[2, ], x[3, ]))

## Represent it
p <- data %>%
  ggplot(aes(x=value, fill=Group)) +
    geom_histogram(color="#e9ecef", alpha = 0.7, position = "identity",
                   binwidth = 1) +
    scale_x_continuous(breaks = c(0, 10, 20, 30),
                       limits = c(0, 30)) +
    scale_fill_manual(values=c("#f7ea83", "#69b3a2", "#404080")) +
    geom_vline(xintercept = c(1, 5, 10), color = "#e05c5c", size = .1) + 
    labs(caption = "choose binwidth = 1")

## ----echo=FALSE---------------------------------------------------------------
print(p)

## -----------------------------------------------------------------------------
## function to generate samples

gen_mix2 <- function(p1, n = 1000, mu = c(0, 3)){
  x1 = rnorm(n, mu[1])
  x2 = rnorm(n, mu[2])

  u = sample(c(1, 0), size = n, prob = c(p1, 1-p1), replace = TRUE)

  return(u*x1 + (1-u)*x2)
}

## generate samples with p1 = 0.75
set.seed(123)
x <- gen_mix2(.75)

## -----------------------------------------------------------------------------
## Histogram it
data = data.frame(x = x)

p <- data %>% 
  ggplot(aes(x = x)) + 
  geom_histogram(aes(y = ..density..),color = "#e9ecef", alpha = 0.7, 
                 binwidth = 0.15) + 
  geom_density(color = "#69b3a2", size = 1) +
  labs(title = "p1 = 0.75")

## ----echo = FALSE-------------------------------------------------------------
print(p)

## -----------------------------------------------------------------------------
## Generate samples with different p1s
n = 1000
p1 = 5:9 / 10
x = array(NA, dim = c(5, n))
for (i in 1:5){
  x[i, ] = gen_mix2(p1[i])
}
## Build dataset
dat = data.frame(value = c(x[1, ], x[2, ], x[3, ], x[4, ], x[5, ]),
                 p1 = factor(c(rep(p1[1], n), rep(p1[2], n),
                                       rep(p1[3], n), rep(p1[4], n),
                                       rep(p1[5], n))))

## -----------------------------------------------------------------------------
## Plot them
p <- ggplot(dat) +
  geom_density(aes(x = value, fill = p1), color = NA, alpha = 0.5) + 
  facet_grid(rows = vars(p1))

## ----echo = FALSE-------------------------------------------------------------
print(p)

## -----------------------------------------------------------------------------
## generate time nodes
tau = 0.025; m = 500
t = cumsum(c(0, rep(tau, m)))

## generate increments of Nt between time nodes
lambda = 1
h = rpois(m, lambda*tau)

## generate Nt
Nt = cumsum(c(0, h))

plot(x = t, y = Nt, type = "l",
     main = "Simulation of Poi(1)", ylab = "N(t)")

## -----------------------------------------------------------------------------
alpha = 2; beta = 2
n = max(Nt)
y = rgamma(n, alpha, beta)

## -----------------------------------------------------------------------------
y_cumsum = c(0, cumsum(y))
Xt = NULL
for (i in 0:m){
  Xt <- c(Xt, y_cumsum[Nt[i+1] + 1])
}

plot(x = t, y = Xt, type = "l",
     main = "Simulation of Compound Poi(1)-Gamma(2, 2) process",
     ylab = "X(t)")

## -----------------------------------------------------------------------------
gen_Xt <- function(lambda, alpha, beta){
  ## generate time nodes
  tau = 0.025; m = 500
  t = cumsum(c(0, rep(tau, m)))

  ## generate increments of Nt between time nodes
  h = rpois(m, lambda*tau)

  ## generate Nt
  Nt = cumsum(c(0, h))
  
  ## generate Gamma r.v.s
  n = max(Nt)
  y = rgamma(n, alpha, beta)
  
  ## generate Xt
  y_cumsum = c(0, cumsum(y))
  Xt = NULL
  for (i in 0:m){
    Xt <- c(Xt, y_cumsum[Nt[i+1] + 1])
  }
  return(data.frame(t = t, Xt = Xt))
}

## -----------------------------------------------------------------------------
lambda = c(1, 2, 2)
alpha = c(1, 2, 3)
beta = c(1, 2, 1)

mean_theo <- function(l, a, b, t = 10){
  l*t*a/b
}

var_theo <- function(l, a, b, t = 10){
  l*t*a*(1 + a)/b^2
}

## -----------------------------------------------------------------------------
R = 200 ## number of reputations
mean_all <- NULL
var_all <- NULL
mean_th <- NULL
var_th <- NULL

for (j in 1:3){ ## the i-th coupling of parameters
  X10 <- NULL
  for (i in 1:R){
    X <- gen_Xt(lambda[j], alpha[j], beta[j])
    X10 <- c(X10, X[t == 10, "Xt"])
  }
  
  mean_all <- c(mean_all, mean(X10)) ## sample mean of X(10)
  var_all <- c(var_all, var(X10)) ## sample variance of X(10)
  
  ## theoretical values
  mean_th <- c(mean_th, mean_theo(lambda[j], alpha[j], beta[j]))
  var_th <- c(var_th, var_theo(lambda[j], alpha[j], beta[j]))
}

result = data.frame(lambda = lambda, alpha = alpha, beta = beta, 
                    emp_mean = mean_all, theo_mean = mean_th,
                    emp_variance = var_all, theo_variance = var_th)

knitr::kable(result, align = c("l", "l", "l", "l", "l", "l", "l"))

## -----------------------------------------------------------------------------
## Quick sort algorithm
qs <- function(x){
  n = length(x)
  if (n == 0 || n == 1){
    return(x)
  }
  
  a = x[1]
  y = x[-1]
  lower = y[y < a]
  upper = y[y >= a]
  return(c(qs(lower), a, qs(upper)))
}

## ----eval=FALSE---------------------------------------------------------------
#  n.all = c(1, 2, 4, 6, 8)*1e4
#  an = NULL
#  m = 100 ## repetitions
#  
#  for (i in 1:5){
#    n = n.all[i]
#    am = rep(NA, m)
#    for (j in 1:m){ ## Monte Carlo simulations
#      x = sample(n) ## a permutation of 1:n
#      am[j] = system.time(qs(x))[1]
#    }
#    an = c(an, mean(am))
#  }

## ----include=FALSE------------------------------------------------------------
load("qs.RData")

## -----------------------------------------------------------------------------
tn = n.all*log(n.all)
data = data.frame(an = an, tn = tn)

plot(tn, an)

fit = lm(an~tn, data = data)
abline(fit, col = "#3db30e")

## -----------------------------------------------------------------------------
## order statistic generator
kth.gen <- function(n, k = 3, m = 5){
  u = matrix(runif(n*m), n, m)
  for (i in 1:n){
    u[i, ] = qs(u[i, ])
  }
  return(u[, k])
}

## ecdf generator
beta.ecdf <- function(a, b){
  k = a; m = k + b - 1
  n = 1000 # sample size
  x <- kth.gen(n, k, m)
  return(ecdf(x))
}

## -----------------------------------------------------------------------------
R = 500 ## times of repetitions
x = 1:9/10
y = array(NA, dim = c(R, length(x)))
est = rep(0, length(x))

for (j in 1:R){
  cdf <- beta.ecdf(3, 3)
  y[j, ] = cdf(x)
}

for (k in 1:length(x)){
  est[k] = mean(y[, k])
}

data = data.frame(x = x, Fn = est, F.real = pbeta(x, 3, 3))

## -----------------------------------------------------------------------------
knitr::kable(data, align = c("l", "l", "l"))

## -----------------------------------------------------------------------------
gen.Rayleigh <- function(n, sigma = 1){
  u = runif(n)
  x = sigma*sqrt(2*log(1/(1-u)))
  return(x)
}

## -----------------------------------------------------------------------------
gen.Rayleigh.anti <- function(n, sigma = 1){
  stopifnot(n %% 2 == 0)
  m = n/2
  u = runif(m)
  uu = c(u, 1 - u) ## antithetic variables
  x = sigma*sqrt(2*log(1/(1-uu)))
  return(x)
}

## -----------------------------------------------------------------------------
x1 = gen.Rayleigh.anti(500)
x2 = gen.Rayleigh(500)

f <- function(x){
  m = length(x)
  stopifnot(m %% 2 == 0)
  y = (x[1:(m/2)] + x[(m/2 + 1):m]) / 2
  return(y)
}

y1 = f(x1)
y2 = f(x2)

var.anti = var(y1) ## antithetic method
var.indp = var(y2) ## independent sampling

(var.indp - var.anti)/var.indp ## percentage of variance reduction

## -----------------------------------------------------------------------------
theta.f1 <- function(n){
  k = 0
  u <- NULL
  repeat{
    uu = runif(1)
    if (uu > 1 - exp(-1/2)){ ## F(1) = 1 - exp(-1/2)
      u <- c(u, uu)
      k <- k + 1
    }
    if (k == n) {break}
  }
  x = sqrt(2*log(1/(1-u)))
  return(x)
}

## -----------------------------------------------------------------------------
theta.f2 <- function(n){
  k = 0
  x = NULL
  repeat{
    xx = rnorm(1)
    if (xx > 1) {x = c(x, xx); k = k + 1}
    if (k == n){break}
  }
  return(x^2)
}

## -----------------------------------------------------------------------------
n = 500 ## sample size
c1 = exp(-1/2)/sqrt(2*pi)
c2 = 1 - pnorm(1)

## Monte Carlo estimation
theta1 <-  theta.f1(n)
intgral.hat1 <- c1 * mean(theta1)
var1 <- c1^2*var(theta1)

theta2 <- theta.f2(n)
intgral.hat2 <- c2* mean(theta2)
var2 <- c2^2*var(theta2)

## -----------------------------------------------------------------------------
intgral.hat1 ## f_1
intgral.hat2 ## f_2

var1 ## f_1
var2 ## f_2

## -----------------------------------------------------------------------------
sk <- function(x) {## computes the sample skewness coefficient
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

## -----------------------------------------------------------------------------
n = 500 ## sample size
m = 1000 ## times of simulations
alpha <- c(0.025, 0.05, 0.95, 0.975)
b1 = numeric(m) ## storing sample skewness

for (i in 1:m){
  x = rnorm(n)
  b1[i] = sk(x)
}

q_emp = quantile(b1, alpha) ## empirical quantiles

q_asy = qnorm(alpha, mean = 0, sd = sqrt(6/n)) ## asymptotic quantiles

## -----------------------------------------------------------------------------
rbind(q_emp, q_asy)

## -----------------------------------------------------------------------------
sd = sqrt(alpha*(1 - alpha)/(n * dnorm(q_asy)^2))
rbind(alpha, sd)

## -----------------------------------------------------------------------------
library(MASS)
sim.copula <- function(n, rho, qm1, qm2){
  mu <- rep(0, 2)
  M <- matrix(rho, 2, 2)
  diag(M) <- 1
  dat <- mvrnorm(n, mu, M)
  dat[,1] <- qm1(pnorm(dat[,1]))
  dat[,2] <- qm2(pnorm(dat[,2]))
  return(dat)
}

## -----------------------------------------------------------------------------
simu <- function(n, m, rho.all, alpha = 0.05, qm1, qm2){
  res <- array(NA, dim = c(3, length(rho.all), m))
  ## test results of Pearson, Spearman, Kendall respectively

  for (i in 1:length(rho.all)){
    rho = rho.all[i]
    
    for (j in 1:m){
      dat = sim.copula(n, rho, qm1, qm2)
      
      res[1, i, j] = (cor.test(dat[, 1], dat[, 2])$p.value < alpha)
      res[2, i, j] = 
        (cor.test(dat[, 1], dat[, 2], method = "spearm")$p.value < alpha)
      res[3, i, j] = 
        (cor.test(dat[, 1], dat[, 2], method = "kendall")$p.value < alpha)
    }
  }
  
  return(res)
}

## -----------------------------------------------------------------------------
n = 500 ## sample size
m = 500 ## simulations
alpha = 0.05 ## significance level
rho.all = c(0.05, 0.10, 0.15)

## ----eval=FALSE---------------------------------------------------------------
#  res <- simu(n, m, rho.all, alpha, qnorm, qnorm)

## ----include=FALSE------------------------------------------------------------
load("res1.RData")

## -----------------------------------------------------------------------------
power1 = numeric(length(rho.all))
power2 = numeric(length(rho.all))
power3 = numeric(length(rho.all))

for (i in 1:length(rho.all)){
  power1[i] = mean(res[1, i, ]) ## pearson
  power2[i] = mean(res[2, i, ]) ## spearman
  power3[i] = mean(res[3, i, ]) ## kendall
}

rbind(rho.all, power1, power2, power3)

## ----eval=FALSE---------------------------------------------------------------
#  res <- simu(n, m, rho.all, alpha, qnorm, qexp)

## ----include=FALSE------------------------------------------------------------
load("res2.RData")

## -----------------------------------------------------------------------------
power1 = numeric(length(rho.all))
power2 = numeric(length(rho.all))
power3 = numeric(length(rho.all))

for (i in 1:length(rho.all)){
  power1[i] = mean(res[1, i, ]) ## pearson
  power2[i] = mean(res[2, i, ]) ## spearman
  power3[i] = mean(res[3, i, ]) ## kendall
}

rbind(rho.all, power1, power2, power3)

## -----------------------------------------------------------------------------
N = 1000
M = 950
H = c(rep(0, M), rep(1, N-M)) ## type of hypothesis

p.original <- function(N, M){
  p.null = runif(M)
  p.alter = rbeta(N-M, 0.1, 1)
  return(c(p.null, p.alter))
}

## -----------------------------------------------------------------------------
adj <- function(p.value){
  p.Bon <- p.adjust(p.value, method = "bonferroni")
  p.BH <- p.adjust(p.value, method = "BH")
  return(rbind(p.Bon, p.BH))
}

## -----------------------------------------------------------------------------
alpha = 0.1
m = 10000
TypeIerror <- true_positive <- rejected <-  array(0, dim = c(3, m)) ## frequencies
FWER <- FDR <- TPR <- numeric(2)

for (j in 1:m){
  p.value <- p.original(N, M)
  p.adj <- adj(p.value)
  ## original
  TypeIerror[1, j] = sum(p.value[1:M] < alpha)
  true_positive[1, j] = sum(p.value[-(1:M)] < alpha)
  rejected[1, j] = sum(p.value < alpha)
  ## adjusted
  for (k in 2:3){
    TypeIerror[k, j] = sum(p.adj[k-1, 1:M] < alpha)
    true_positive[k, j] = sum(p.adj[k-1, -(1:M)] < alpha)
    rejected[k, j] <- sum(p.adj[k-1, ] < alpha)
  }
}

for (k in 1:3){
  FWER[k] <- sum(TypeIerror[k, ] >= 1)/m
  FDR[k] <- mean(TypeIerror[k, ]/rejected[k, ])
  TPR[k] <- mean(true_positive[k, ])/(N-M)
}

res <- data.frame(FWER, FDR, TPR, row.names = c("Original", "Bon", "B-H"))
print(res)

## -----------------------------------------------------------------------------
library(boot)
x <- aircondit$hours
n = length(x)

## -----------------------------------------------------------------------------
f <- function(x){
  1/mean(x)
}

theta.hat = f(x)

## -----------------------------------------------------------------------------
B = 50
theta.boot = numeric(B)

for (b in 1:B){
  x.star <- sample(x, size = n, replace = T)
  theta.boot[b] <- f(x.star)
}

res = data.frame(lambda.hat = theta.hat,
                 bias.estimated = mean(theta.boot) - theta.hat,
                 se.estimated = sd(theta.boot), B = B)

print(res)

## -----------------------------------------------------------------------------
theta.hat = mean(x)

## -----------------------------------------------------------------------------
stat.mean = function(x, i){
  mean(x[i])
}

## -----------------------------------------------------------------------------
m = 1000 ## replications
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)

for (i in 1:m){
  de <- boot(x, statistic = stat.mean, R = 2000)
  ci = boot.ci(de, type=c("norm","basic","perc","bca"))
  ci.norm[i,]<-ci$norm[2:3];    ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]; ci.bca[i,]<-ci$bca[4:5]
}

## -----------------------------------------------------------------------------
cat("norm:(", mean(ci.norm[, 1]), ",", mean(ci.norm[, 2]), ")\n")
cat("basic:(", mean(ci.basic[, 1]), ",", mean(ci.basic[, 2]), ")\n")
cat("percentile:(", mean(ci.perc[, 1]), ",", mean(ci.perc[, 2]), ")\n")
cat("bca:(", mean(ci.bca[, 1]), ",", mean(ci.bca[, 2]), ")")

## -----------------------------------------------------------------------------
data(scor, package = "bootstrap")
n = nrow(scor)
pairs(scor)
mcor <- cor(scor)
round(mcor, digits = 2) ## sample correlation

## -----------------------------------------------------------------------------
theta <- function(ind){
  vals <- eigen(var(scor[ind, ]), symmetric = T,
                only.values = T)$values
  vals[1]/sum(vals)
}
theta.hat = theta(1:n)
print(theta.hat)

## -----------------------------------------------------------------------------
#compute the jackknife replicates, leave-one-out estimates
theta.jack <- numeric(n)
for (i in 1:n)
  theta.jack[i] <- theta(-i)
#jackknife estimates of bias and standard error
bias <- (n - 1)*(mean(theta.jack) - theta.hat)
print(bias)#jackknife estimate of bias
sd <- sqrt((n - 1)*
             mean((theta.jack - mean(theta.jack))^2))
print(sd)#jackknife estimate of standard error

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)

# the four models
L1 <- lm(magnetic ~ chemical)
L2 <- lm(magnetic ~ chemical + I(chemical^2))
L3 <- lm(log(magnetic) ~ chemical)
L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))

# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n){
  y <- magnetic[-k]
  x <- chemical[-k]
  ## linear
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2]*chemical[k]
  e1[k] <- magnetic[k] - yhat1
  ## quadratic
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2]*chemical[k] + 
    J2$coef[3]*chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  ## log linear
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2]*chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  ## cubic
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2]*chemical[k] + 
    J4$coef[3]*chemical[k]^2 + J4$coef[4]*chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}
detach(ironslag)

## -----------------------------------------------------------------------------
# MSE
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
# adjusted R squared
summary(L1)$adj.r.squared
summary(L2)$adj.r.squared
summary(L3)$adj.r.squared
summary(L4)$adj.r.squared

## -----------------------------------------------------------------------------
cramer_von_Mises.stat <- function(x, y){
  n = length(x); m = length(y)
  Fn = ecdf(x); Gm = ecdf(y)
  z = c(x, y)
  n*m*mean((Fn(z) - Gm(z))^2)/(n + m)
}

perm.test <- function(x, y, R = 999, test_stat){
  z <- c(x, y) # pooled sample
  N = length(z)
  n = length(x); m = N - n
  K <- 1:N # indices for the pooled sample
  W <- numeric(R) # storage for replicates
  
  for (i in 1:R){
    #generate indices k for the first sample
    k <- sample(K, size = n)
    x1 = z[k]; y1 = z[-k]
    W[i] <- test_stat(x1, y1)
  }
  W0 <- test_stat(x, y)
  p <- mean(c(W0, W) >= W0)
  return(list(W0 = W0, W.perm = W, p.value = p))
}

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

res <- perm.test(x, y, 999, cramer_von_Mises.stat)

res$p.value
hist(res$W.perm, freq = FALSE, xlab = "test_stat_permuted",
     main = NULL)
abline(v = res$W0, col = "red")

## -----------------------------------------------------------------------------
data(scor, package = "bootstrap")
x = scor[, 3]; y = scor[, 4]
z <- c(x, y)

cor_bv <- function(z, ix){
  stopifnot(length(z) %% 2 == 0)
  n = length(z)/2
  z <- z[ix]
  cor(z[1:n], z[-(1:n)], method = "spearman")
}

library(boot)
boot.obj <- boot(z, statistic = cor_bv,
                 sim = "permutation", R = 999)

boot.obj$t0 ## original correlation

tb <- c(boot.obj$t, boot.obj$t0)
p.value <- mean(tb >= boot.obj$t0)
p.value

hist(boot.obj$t, freq = FALSE, main = "", xlim = c(-1, 1),
     xlab = "replicates of Spearman correlations")
abline(v = boot.obj$t0, col = "red")

## -----------------------------------------------------------------------------
options(warn = -1)
cor.test(x, y, method = "spearman")

## -----------------------------------------------------------------------------
## random walk Metropolis
mcmc.rw <- function(N = 5000, x0, f=dcauchy, sigma = 1){
  x <- numeric(N)
  x[1] = x0
  u <- runif(N)
  for (i in 2:N) {
    xt <- x[i-1]
    y <- rnorm(1, sd = sigma) + xt
    num = f(y)
    den = f(xt)
    if (u[i] <= num/den) x[i] <- y
    else x[i] <- xt
  }
  return(x)
}

## -----------------------------------------------------------------------------
N <- 5000
b <- 1000

x = mcmc.rw(N, x0 = rcauchy(1), f = dcauchy)
x = x[(b + 1):N]

plot((b + 1):N, y = x, type = "l", main = "", xlab = "t", ylab = "X")
abline(h = 0, col = "green")

## -----------------------------------------------------------------------------
quantile(x, c(0.25, 0.75))
qcauchy(p = c(0.25, 0.75))

## -----------------------------------------------------------------------------
Gibbs <- function(N = 5000, x0 = c(0, 0),
                  a = 1, b = 1, n = 100) {
  X <- matrix(0, N, 2) # the chain, a bivariate sample
  X[1, ] = x0          # initialize
  
  for (i in 2:N){
    y <- X[i-1, 2]
    X[i, 1] <- rbinom(1, n, y)
    x <- X[i-1, 1]
    X[i, 2] <- rbeta(1, x + a, n - x + b)
  }
  return(X)
}

## -----------------------------------------------------------------------------
N <- 5000            # the length of chain
burn <- 1000         # the burn-in length

####set parameters####
a <- b <- 1
n = 100

X <- Gibbs(N, c(0,0), a, b, n)
X <- X[(burn + 1):N, ]

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j]) for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  m <- nrow(psi)
  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est. 
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
m = 4
N = 25000
burn = 1000

X <- array(dim = c(m, N))

#choose overdispersed initial values
x0 = c(-10, -5, 5, 10)

####generating chains####
for (i in 1:m){
  X[i, ] <- mcmc.rw(N, x0[i])
}

## -----------------------------------------------------------------------------
#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i, ] <- psi[i, ] / (1:ncol(psi))
print(Gelman.Rubin(psi))

## -----------------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
for (i in 1:m){
  plot(psi[i, (burn + 1):N], type = "l",
       xlab = i, ylab = bquote(psi))
}

## -----------------------------------------------------------------------------
par(mfrow = c(1, 1))
R.hat <- rep(0, N)
for (j in (burn + 1):N){
  R.hat[j] <- Gelman.Rubin(psi[, 1:j])
}
plot(R.hat[(burn + 1):N], type = "l", ylim = c(0, 3.5),
     xlab = "", ylab = "R")
abline(h = 1.2, lty = 2)

## -----------------------------------------------------------------------------
m = 4
N = 25000
burn = 1000

X <- array(0, dim = c(m, N, 2))

#choose overdispersed initial values
x0 = matrix(c(0,  25, 50,  75,
              0, .25, .5, .75), 4, 2)

a <- b <- 1
n = 100

####generating chains####
for (i in 1:m){
  X[i, , ] <- Gibbs(N, x0[i, ])
}

## -----------------------------------------------------------------------------
#compute diagnostic statistics for Xt and Yt respectively
psi.x <- t(apply(X[, , 1], 1, cumsum))
for (i in 1:nrow(psi.x))
  psi.x[i, ] <- psi.x[i, ] / (1:ncol(psi.x))
print(Gelman.Rubin(psi.x))

psi.y <- t(apply(X[, , 2], 1, cumsum))
for (i in 1:nrow(psi.y))
  psi.y[i, ] <- psi.y[i, ] / (1:ncol(psi.y))
print(Gelman.Rubin(psi.y))

## -----------------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
for (i in 1:m){
  plot(psi.x[i, (burn + 1):N], type = "l",
       xlab = i, ylab = "psi.x")
}

## -----------------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
for (i in 1:m){
  plot(psi.y[i, (burn + 1):N], type = "l",
       xlab = i, ylab = "psi.y")
}

## -----------------------------------------------------------------------------
par(mfrow = c(1, 1))
R.hat.x <- rep(0, N)
for (j in (burn + 1):N){
  R.hat.x[j] <- Gelman.Rubin(psi.x[, 1:j])
}
plot(R.hat.x[(burn + 1):N], type = "l", ylim = c(0, 3.5),
     xlab = "", ylab = "R of Xt")
abline(h = 1.2, lty = 2)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 1))
R.hat.y <- rep(0, N)
for (j in (burn + 1):N){
  R.hat.y[j] <- Gelman.Rubin(psi.y[, 1:j])
}
plot(R.hat.y[(burn + 1):N], type = "l", ylim = c(0, 3.5),
     xlab = "", ylab = "R of Yt")
abline(h = 1.2, lty = 2)

## -----------------------------------------------------------------------------
fk <- function(k, a, d){
  stopifnot(d >= 1 && length(a) == d)
  A = log(sum(a^2))
  
  signk = (k %% 2 == 0)
  
  signk*exp((k + 1)*A + lgamma((d + 1)/2) + lgamma(k + 3/2) - lgamma(k + d/2 + 1) - lgamma(k + 1) - k*log(2))/(2*k + 1)/(2*k + 2)
}

fk(20, rep(1, 10), 10)

## -----------------------------------------------------------------------------
Sk <- function(k, a, d){
  stopifnot(d >= 1 && length(a) == d)
  A = log(sum(a^2))
  
  K = 0:k
  signk = rep(c(1, -1), (k + 2) %/% 2)
  signk = signk[1:(k + 1)]
  sum(signk*exp((K + 1)*A + lgamma(d/2 + 1/2) + lgamma(K + 3/2) - 
                  lgamma(K + d/2 + 1) - lgamma(K + 1) - K*log(2))/(2*K + 1)/(2*K + 2))
}

Sk(20, rep(1, 10), 10)
Sk(20, rep(1, 10), 10) - Sk(19, rep(1, 10), 10)

## -----------------------------------------------------------------------------
a = c(1, 2); d = 2
K = 0:15
S <- rep(0, length(K))

for (k in K){
  S[k + 1] <- Sk(k, a, d)
}

data = data.frame(k = K, sum = S)

knitr::kable(data, align = c("c", "c"))
plot(data, type = "l")

## -----------------------------------------------------------------------------
f <- function(k, a){
  stopifnot(k + 1 - a^2 > 0)
  ck = sqrt(a^2*k/(k+1-a^2))
  
  b = 2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))
  
  integrand <- function(u){
    (1 + u^2/k)^{-(k+1)/2}
  }
  
  I = integrate(integrand, 0, ck)
  
  b*I$value
}

## -----------------------------------------------------------------------------
K = c(4:25, 100, 500, 1000)
ans = rep(0, length(K))
eps <- .Machine$double.eps^{0.25}

for (i in 1:length(K)){
  k = K[i]
  ans[i] <- uniroot(function(a){
    f(k - 1, a) - f(k, a)
  }, interval = c(1, 1.9))$root
  
}

data = data.frame(k = K, a = ans)
knitr::kable(data, align = c("c", "c"))

