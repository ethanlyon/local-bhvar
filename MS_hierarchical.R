#################################
# Fit SVAR by Gibbs Sampler     #
# multi subject hierarchical    #
#                               #
# using auxilliary variables    #
# to break variable dependence  #
#################################




library(MASS)
library(statmod)
library(mvnfast)




MSGibbs <- function(eta,                          # array of dim TT * B * J * N
                      H,
                      HH,
                      eta0,
                      w_MLE,
                      TT,       # a vector of T1, ..., TN for N subjects
                      B,
                      p=2, 
                      N,
                      seed=100, 
                      numIter=2000, 
                      thin=1, 
                      warmup=0,
                      K_inv = diag(1/(B-1), B),      # inverse of scale matrix of Lambda prior
                      nu = 1,                        # df of Lambda prior
                      mu1 = 1,                       # mean of lambda1sq prior
                      nu1 = 0.001,                   # df of lambda1sq prior
                      mu2 = 1,                       # mean of lambda2 prior
                      nu2 = 0.01,                    # df of lambda2 prior
                      k1 = 0.0005,                      # shape prior of omega_v
                      theta1 = 200,                # scale prior of omega_v
                      a = 100,
                      initial=NULL
){               
  
  
  
  # store values
  sim.v_n <- sim.w <- vector(mode="list", length=(numIter-warmup)/thin+1) 
  sim.Lambda <- sim.lambda1sq <- sim.lambda2 <- sim.tausq <- vector(mode="list", length=(numIter-warmup)/thin+1)
  sim.U_v <- sim.alpha <- vector(mode="list", length=(numIter-warmup)/thin+1)
  sim.v_n.star <- vector(mode="list", length=(numIter-warmup)/thin+1)
  sim.omega_v <- vector(mode="list", length=(numIter-warmup)/thin+1)
  
  
  # initialize
  set.seed(seed)
  if(is.null(initial)){
    sim.w[[1]] <- w <- rnorm(B^2*p, 0, 1)
    sim.v_n[[1]] <- v_n <- lapply(vector("list", N), function(x){rnorm(B^2*p, 0, 1)})
    sim.Lambda[[1]] <- Lambda <- Posdef(n=B, ev=runif(B, 0, 0.001))
    sim.lambda1sq[[1]] <- lambda1sq <- runif(B^2*p, 0.1, 1000)
    sim.lambda2[[1]] <- lambda2 <- runif(B^2*p, 0.1, 10)
    sim.tausq[[1]] <- tausq <- runif(B^2*p, 0.1, 10)
    sim.U_v[[1]] <- U_v <- runif(B^2*p, 0.1, 1)
    sim.alpha[[1]] <- alpha <- runif(B^2*p, 0.1, 1)          
    sim.v_n.star[[1]] <- lapply(1:N, function(n){alpha*v_n[[n]]})
    sim.omega_v[[1]] <- U_v/alpha^2
  } else {
    sim.w[[1]] <- w <- initial$w
    sim.v_n[[1]] <- v_n <- initial$v_n
    sim.Lambda[[1]] <- Lambda <- initial$Lambda
    sim.lambda1sq[[1]] <- lambda1sq <- initial$lambda1sq
    sim.lambda2[[1]] <- lambda2 <- initial$lambda2
    sim.tausq[[1]] <- tausq <- initial$tausq
    sim.U_v[[1]] <- U_v <- initial$U_v
    sim.alpha[[1]] <- alpha <- initial$alpha          
    sim.v_n.star[[1]] <- v_n.star <- initial$v_n.star
    sim.omega_v[[1]] <- omega_v <- initial$omega_v
  }
  
  
  
  for(i in 1:numIter){
    
    # recover G, xisq and prec
    temp <- Gamma_recover(tausq=tausq, lambda1sq=lambda1sq, Lambda=Lambda, B=B, p=p)      
    G <- temp$G
    xisq <- temp$xisq
    prec <- prec_get(HH=HH, Lambda=Lambda)
    
    # update parameters
    w <- w_update(prec=prec, w_MLE=w_MLE, v_n=v_n, alpha=alpha, tausq=tausq, lambda2, N=N)
    v_n <- v_n_update(prec=prec, w_MLE=w_MLE, w=w, alpha=alpha, U_v=U_v, N=N)
    Lambda <- Lambda_update(eta0=eta0, w=w, v_n=v_n, alpha=alpha, H=H, G=G, K_inv=K_inv, nu=nu, TT=TT, B=B, p=p, N=N)
    lambda1sq <- lambda1sq_update(tausq=tausq, xisq=xisq, B=B, p=p, mu1=mu1, nu1=nu1)          
    lambda2 <- lambda2_update(w=w, B=B, p=p, mu2=mu2, nu2=nu2)                         
    tausq <- tausq_update(lambda1sq=lambda1sq, w=w, xisq=xisq, B=B, p=p) 
    U_v <- U_v_update(v_n=v_n, k1=k1, theta1=theta1, N=N, B=B, p=p)
    alpha <- alpha_update(prec=prec, w_MLE=w_MLE, w=w, v_n=v_n, a=a, N=N, B=B, p=p)
    v_n.star <- lapply(1:N, function(n){alpha*v_n[[n]]})
    omega_v <- U_v/alpha^2

    
    
    # store thinned sims after warmup
    if((i>warmup) & (i %% thin == 0)){  
      sim.w[[(i-warmup)/thin+1]] <- w
      sim.v_n[[(i-warmup)/thin+1]] <- v_n
      sim.Lambda[[(i-warmup)/thin+1]] <- Lambda
      sim.lambda1sq[[(i-warmup)/thin+1]] <- lambda1sq
      sim.lambda2[[(i-warmup)/thin+1]] <- lambda2
      sim.tausq[[(i-warmup)/thin+1]] <- tausq
      sim.U_v[[(i-warmup)/thin+1]] <- U_v 
      sim.alpha[[(i-warmup)/thin+1]] <- alpha 
      sim.v_n.star[[(i-warmup)/thin+1]] <- v_n.star
      sim.omega_v[[(i-warmup)/thin+1]] <- omega_v
    }
    
    # report iter
    cat("i=", i, "\n", sep="")
  }
  
  
  return(list(sim.w=sim.w, 
              sim.v_n=sim.v_n,
              sim.Lambda=sim.Lambda, 
              sim.lambda1sq=sim.lambda1sq, 
              sim.lambda2=sim.lambda2, 
              sim.tausq=sim.tausq, 
              sim.U_v=sim.U_v,
              sim.alpha=sim.alpha,
              sim.v_n.star=sim.v_n.star,
              sim.omega_v=sim.omega_v))   
}








# generate n-by-n positive definite matrix
Posdef <- function (n, ev = runif(n, 0, 10)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}







# (0.1) recover Gamma, xisq
Gamma_recover <- function(tausq, lambda1sq, Lambda, B, p){
  M <- kronecker(diag(1, B*p), chol2inv(chol(Lambda)))
  l <- B^2*p
  xisq <- c(rep(NA, B-1), M[B,B])                        # for xisq
  for(j in (B-1):1) xisq[j] <- M[j,j]- M[j,(j+1):B] %*% solve(M[(j+1):B,(j+1):B]) %*% M[(j+1):B,j]
  xisq <- rep(xisq, B*p)
  m.all.possible <- vector("list", B-1)                  # for m_j
  for(j in (B-1):1) m.all.possible[[j]] <- M[j,(j+1):B] %*% solve(M[(j+1):B,(j+1):B])
  G <- c(rep(NA, l-1), sqrt(tausq[l]*lambda1sq[l]))      # for gamma
  for(j in (l-1):1){
    k1 <- (l-j)%%B
    k2 <- floor((l-j)/B)
    if(k1>0){
      m <- m.all.possible[[B-k1]]
      G[j] <- sqrt(tausq[j]*lambda1sq[j]) - sum(m*G[(j+1):(j+k1)])
    } else {
      G[j] <- sqrt(tausq[j]*lambda1sq[j])
    }
  }
  return(list(G=G, xisq=xisq))
}



# (0.2) get prec: list of N*J
prec_get <- function(HH, Lambda){lapply(HH, function(x){kronecker(x, Lambda)})}   # a list of length N




# (1) update w
w_update <- function(prec, w_MLE, v_n, alpha, tausq, lambda2, N){
  D_vec <- 2*tausq/(2*lambda2*tausq+1)
  mu_n <- lapply(1:N, function(n){w_MLE[[n]] - alpha*v_n[[n]] })  
  Ssq_w <- chol2inv(chol(Reduce("+", prec) + diag(1/D_vec)))
  M_w <- Ssq_w %*% apply(mapply("%*%", prec, mu_n), 1, sum)
  w <- as.vector(rmvn(n=1, mu=M_w, sigma=chol(Ssq_w), isChol=TRUE,  ncores = 4))
  return(w)
}




# (2) update v_n
v_n_update <- function(prec, w_MLE, w, alpha, U_v, N){        
  v_n <- lapply(1:N, function(n){
    Ssq <- chol2inv(chol(t(t(prec[[n]]*alpha)*alpha) + diag(U_v)))   
    mu <- w_MLE[[n]]-w
    M <- Ssq %*% (alpha*(prec[[n]]%*%mu))
    as.vector(rmvn(n=1, mu=M, sigma=chol(Ssq), isChol=TRUE, ncores=4))
  })
  return(v_n)
}




# (4) update Lambda
Lambda_update <- function(eta0, w, v_n, alpha, H, G, K_inv, nu, TT, B, p, N){
  S <- lapply(1:N, function(n){
    W <- matrix(w + alpha*v_n[[n]], B)
    tcrossprod(eta0[[n]] - W%*%H[[n]])
  })
  S_Lambda <- chol2inv(chol(Reduce('+', S) + K_inv + 2*tcrossprod(matrix(G, nrow=B))))
  Sm <- min(eigen(S_Lambda)$values)
  if(Sm<=0) S_Lambda <- S_Lambda+diag(abs(Sm)+0.1, B)
    Lambda <- rWishart(1, df=sum(TT-p)+2*B*p+nu, Sigma=S_Lambda)[,,1]    
  while(sum(eigen(Lambda)$value<=0) | (!isSymmetric(Lambda))){                                    
    Lambda <- rWishart(1, df=sum(TT-p)+2*B*p+nu, Sigma=S_Lambda)[,,1]
  }
  return(Lambda)
}





# (5) update lambda1sq for j=1,...,B^2*p
lambda1sq_update <- function(tausq, xisq, B, p, mu1, nu1){
  if(length(mu1)==1) mu1 <- rep(mu1, B^2*p)
  if(length(nu1)==1) nu1 <- rep(nu1, B^2*p)
  nu_lambda1sq <- nu1+2
  mu_lambda1sq <- nu_lambda1sq*xisq*mu1/(2*tausq*mu1+nu1*xisq)
  lambda1sq <- abs(rgamma(n=rep(1, B^2*p), shape=nu_lambda1sq/2, scale=2*mu_lambda1sq/nu_lambda1sq))
  return(lambda1sq)
}





# (6) update lambda2 for j=1,...,B^2*p
lambda2_update <- function(w, B, p, mu2, nu2){
  if(length(mu2)==1) mu2 <- rep(mu2, B^2*p)
  if(length(nu2)==1) nu2 <- rep(nu2, B^2*p)
  nu_lambda2 <- nu2+2
  mu_lambda2 <- nu_lambda2*mu2/(w^2*mu2+nu2)
  lambda2 <- abs(rgamma(n=rep(1, B^2*p), shape=nu_lambda2/2, scale=2*mu_lambda2/nu_lambda2))
  return(lambda2)
}





# (7) update tausq for j=1,...,B^2*p
tausq_update <- function(lambda1sq, w, xisq, B, p){
  M <- sqrt(lambda1sq/(w^2)/xisq)
  S <- lambda1sq/xisq
  temp <- abs(rinvgauss(n=rep(1, B^2*p), mean=M, shape=S))
  tausq <- 1/2/temp
  return(tausq)
}




# (8) update U_v
U_v_update <- function(v_n, k1, theta1, N, B, p){
  d <- Reduce("+", lapply(v_n, function(x){x^2}))
  kv <- rep(N/2+k1, B^2*p)
  thetav <- 1/(d/2+1/theta1)     
  U_v <- abs(rgamma(n=rep(1, B^2*p), shape=kv, scale=thetav))  
  return(U_v)
}





# (10) update alpha
alpha_update <- function(prec, w_MLE, w, v_n, a, N, B, p){
  Ssq_alpha <- chol2inv(chol(Reduce("+",lapply(1:N, function(n){
    t(t(prec[[n]]*v_n[[n]])*v_n[[n]])
    })) + diag(a, B^2*p)))
  M_alpha <- Ssq_alpha %*%  Reduce("+", lapply(1:N, function(n){(prec[[n]]%*%(w_MLE[[n]]-w))*v_n[[n]]}))
  alpha <- as.vector(rmvn(n=1, mu=M_alpha, sigma=chol(Ssq_alpha), isChol=TRUE, ncores=4))
  return(alpha)
}





