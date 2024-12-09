## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(MASS)
library(slam)
n = 500
data = mvrnorm(n, mu = c(0, 0, 3), Sigma = diag(3))
data = data / row_norms(data)

## -----------------------------------------------------------------------------
library(SphRankSign)
temp = Fn(data, 25, 20, 0)
temp$FnTheta  ## the axis (North Pole)
Fnn = temp$Fn ## empirical distribution
Fnn[1:5, ]

res = RankSign(data, 25, 20, 0)
Rn = res$rank ## ranks
Rn[1:10]
Sn = res$sign ## signs
Sn[1:5, ]
res$FnTheta   ## the axis

## -----------------------------------------------------------------------------
p = length(data[1, ])

UU = array(rnorm(n*p) , dim = c(n,p))
U = diag(sqrt(diag(1/(UU%*%t(UU)))))%*%UU

d = matrix(NA, n, n)  ## transportation cost matrix
for(i in 1:n){
  for(j in 1:n){
    d[i,j] = (acos(t(data[i,])%*%U[j,]))^2/2
  }
}
library(transport)
library(Rcpp)
sourceCpp('../src/KM.cpp')

library(microbenchmark)
tm2 <- microbenchmark(
  ot1 = transport(rep(1, n), rep(1, n), costm = d, method = "networkflow")$to,
  ot2 = opt_coupling(d, rep(1, n))
)

knitr::kable(summary(tm2)[,c(1,3,5,6)])

