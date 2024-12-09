#' @import RiemBase
#' @importFrom Rcpp evalCpp
#' @importFrom transport transport
#' @importFrom stats rnorm
NULL

Fr = function(data){
  da = list()
  size = length(data[, 1])
  for(i in 1:size){
    da[[i]] = data[i,]
  }
  da = RiemBase::riemfactory(da, name="sphere")
  Fr_rob = RiemBase::rbase.robust(da)
  return(Fr_rob$x)
}

geodist = function(x, y){
  return(acos(t(x)%*%y))
}

unifgrid=function(n, p){
  UU = array(rnorm(n*p) , dim = c(n,p))
  U = diag(sqrt(diag(1/(UU%*%t(UU)))))%*%UU
  return(U)
}

Approx_Fr = function(x){
  hat_Fr = Fr(x)
  dist = apply(x, 1, geodist, hat_Fr)
  ind = which.min(dist)
  return(list(index = ind, sample_Approx_Fr = x[ind, ]))
}

opt_assign = function(data, Yn){
  size = length(data[, 1])
  d = matrix(NA, size, size)
  for(i in 1:size){
    for(j in 1:size){
      d[i,j] = (acos(t(data[i,])%*%Yn[j,]))^2/2
    }
  }
  return(transport::transport(rep(1, size), rep(1, size),
                              costm = d, method = "networkflow")$to)
}

Fn_theta = function(x, p){
  size = length(x[,1])
  Grid = unifgrid(size, p)
  Fnn = Grid[opt_assign(x, Grid), ]
  approx_theta = Approx_Fr(x)
  FnTheta = Fnn[approx_theta$index,]
  return(FnTheta)
}

Ebasis = function(data, p, FnTheta){
  e1 = FnTheta
  x = c(1, 0, 0)
  y = c(0, 1, 0)
  u2 = x - sum(e1*x)*e1
  e2 = u2/sqrt(sum(u2^2))
  u3 = y - sum(e1*y)*e1 - sum(e2*y)*e2
  e3 = u3/sqrt(sum(u3^2))
  return(rbind(as.vector(e1), as.vector(e2), as.vector(e3)))
}

EuclidCor = function(zeta, data, p, FnTheta){
  E = Ebasis(data, p, FnTheta)
  eta = rep(NA, p)
  eta[1] = cos(zeta[1])
  eta[2] = sin(zeta[1])*cos(zeta[2])
  eta[3] = sin(zeta[1])*sin(zeta[2])
  eta = matrix(eta, ncol = 1)
  return(solve(E)%*%eta)
}

SphCor = function(X){
  theta = acos(X[1])
  phi = atan2(X[3], X[2])
  return(180*c(theta, phi))
}

Sgrid = function(data, nS, nR, n0, p){
  phi = seq(0, nR/(nR+1), length = (nR+1))
  phi = phi[2:(nR+1)]
  phi = acos(2*phi - 1)
  theta = seq(0, 2*pi, length = (nS+1))
  theta = theta[2:(nS+1)]

  angle = expand.grid(phi, theta, KEEP.OUT.ATTRS = FALSE)

  n = nR*nS+n0
  Y = matrix(NA, p, n)
  FnTheta = Fn_theta(data, p)
  if(n0>0){
    for(i in 1:n0){
      Y[ ,i] = FnTheta
    }
  }
  Y[ ,(n0+1):n] = apply(angle, MARGIN = 1, EuclidCor, data, p, FnTheta)
  Y = t(Y)
  Y = Y / slam::row_norms(Y)
  return(list(Grid = Y, FnTheta = FnTheta))
}

Fn0 = function(data, p){
  size = length(data[,1])
  Grid = unifgrid(size, p)
  Fnn = Grid[opt_assign(data, Grid), ]
  return(Fnn)
}

#' Empirical distribution function for spherical data
#'
#' This function calculates the empirical distribution function for spherical data defined by optimal transport. The axis around which the uniform grids with rotational-symmetry is generated is also given in the results.
#' @param data spherical data
#' @param nS number of signs
#' @param nR number of ranks
#' @param n0 number of sample points transported to the axis
#' @return empirical distribution, axis
#' @keywords distribution spherical OT
#' @export
Fn = function(data, nS, nR, n0){
  p = length(data[1, ])
  stopifnot(p == 3)
  size = length(data[, 1])
  temp = Sgrid(data, nS, nR, n0, p)
  Yn = temp$Grid
  Fnn = matrix(NA, size, p)
  Fnn = Yn[opt_assign(data, Yn), ]
  return(list(Fn = Fnn, FnTheta = temp$FnTheta))
}

#' Ranks and signs for spherical data
#'
#' Calculates the ranks and signs based on empirical distribution for spherical data.
#' @param data spherical data
#' @param nS number of signs
#' @param nR number of ranks
#' @param n0 number of sample points transported to the axis
#' @return ranks, signs, indices for data points that are transported to the axis, the axis
#' @keywords rank sign spherical
#' @export
RankSign = function(data, nS, nR, n0){
  p = length(data[1, ])
  stopifnot(p == 3)
  size = length(data[, 1])
  temp = Fn(data, nS, nR, n0)
  Fnn = temp$Fn
  FnTheta = temp$FnTheta
  prod_FnFr = rep(NA, size)
  TanDecom = matrix(NA, size, p)
  TanDecom_Norm = rep(NA, size)
  for(i in 1:size){
    prod_FnFr[i] = (t(Fnn[i,])%*%FnTheta)
    TanDecom[i, ] = Fnn[i, ] - prod_FnFr[i]*FnTheta
    TanDecom_Norm[i] = sqrt(sum(TanDecom[i,]^2))
  }
  Mn = rank(-prod_FnFr)
  Rn = ceiling((Mn-n0)/nS)*(Mn>n0)

  Sn = matrix(0, size, p)
  j = which(Rn == 0)
  fullsize = 1:size

  if(length(j) == 0){index = fullsize}
  else {index = fullsize[-j]}

  for(i in 1:p){
    Sn[index, i] = TanDecom[index, i]/TanDecom_Norm[index]
  }

  return(list(rank = Rn, sign = Sn, pole_index = j, FnTheta = FnTheta))
}

