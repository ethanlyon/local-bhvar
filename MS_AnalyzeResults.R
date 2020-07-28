######################
# anaylize results   #
######################
library(rstan)


w_inference <- function(chains, B, p, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=3){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim+1, by=thin)
  l <- length(sel)
  sw <- array(NA, dim=c(l, n.chain, B^2*p))
  for(i in 1:n.chain) sw[,i,] <- do.call("rbind", chains[[i]]$sim.w[sel])
  mon.w <- rstan::monitor(sw, print=FALSE,probs=c(0.025, 0.975, 0.05, 0.95))
  w.mode <-  apply(sw, 3, function(x){find.mode(x, d=d)}) 
  w.median <- mon.w[,6]
  
  if(plot){
    pdf(filename, width=14)   # trace & hist
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*p)) {
      matplot(sw[,,j], type="l", lty=1, main=paste("w", j, sep=""), ylab="")
      hist(sw[,,j], main=paste("w", j, sep=""), breaks=30, ylab="", xlab="")
      #abline(v= w.median[j], col=2)
      #axis(side=1, at= w.median[j] , labels=round( w.median[j], d))
      abline(v=w.mode[j], col=4)
      axis(side=3, at=w.mode[j], labels=w.mode[j])
    }
    dev.off()
  }
  return(list(sw=sw, mon.w=mon.w, w.mode=w.mode, w.median=w.median))
}



v_inference <- function(chains, B, p,N, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=3){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim+1, by=thin)
  l <- length(sel)
  sv <- array(NA, dim=c(l, n.chain, B^2*p*N))
  for(i in 1:n.chain){sv[,i,] <-  do.call("rbind", lapply(chains[[i]]$sim.v_n.star[sel], unlist))}
  mon.v <- rstan::monitor(sv, print=F, probs=c(0.025, 0.975, 0.05, 0.95))
  v.mode <-  apply(sv, 3, function(x){find.mode(x, d=d)}) 
  v.median <- mon.v[,6]
  
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*p*N)) {
      matplot(sv[,,j], type="l", lty=1, main=paste("v", j, sep=""), ylab="")
      hist(sv[,,j], main=paste("v", j, sep=""), breaks=40, ylab="", xlab="")
      #abline(v=v.median[j], col=2)
      #axis(side=1, at=v.median[j] , labels=round(v.median[j], d))
      abline(v=v.mode[j], col=4)
      axis(side=3, at=v.mode[j], labels=v.mode[j])
    }
    dev.off()
  }
  return(list(sv=sv, mon.v=mon.v, v.mode=v.mode, v.median=v.median))
}








u_inference <- function(chains, B, p, J, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=3){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim, by=thin)
  l <- length(sel)
  su <- array(NA, dim=c(l, n.chain, B^2*p*N*J))
  for(i in 1:n.chain){su[,i,] <-  do.call("rbind", lapply(chains[[i]]$sim.u_nj.star[sel], unlist))}   # changed u_j to u_nj
  mon.u <- rstan::monitor(su, print=F)
  u.mode <-  apply(su, 3, function(x){find.mode(x, d=d)}) 
  u.median <- mon.u[,6]
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*p*J)){
      matplot(su[,,j], type="l", lty=1, main=paste("u", j, sep=""),ylab="")
      hist(su[,,j], main=paste("u", j, sep=""), breaks=40, ylab="", xlab="")
      #abline(v=u.median[j], col=2)
      #axis(side=1, at=u.median[j] , labels=round(u.median[j], d))
      abline(v=u.mode[j], col=4)
      axis(side=3, at=u.mode[j], labels=u.mode[j])
    }
    dev.off()
  }
  return(list(su=su, mon.u=mon.u, u.mode=u.mode, u.median=u.median))
}









L_inference <- function(chains, B, p, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=6){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim+1, by=thin)
  l <- length(sel)
  sL <- array(NA, dim=c(l, n.chain, B^2))
  for(i in 1:n.chain){sL[,i,] <-  do.call("rbind", lapply(chains[[i]]$sim.Lambda[sel], c))}
  mon.L <- rstan::monitor(sL, print=F)
  L.mode <-  apply(sL, 3, function(x){find.mode(x, d=d)}) 
  L.median <- mon.L[,6]
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2)){
      matplot(sL[,,j], type="l", lty=1, main=paste("lambda", j, sep=""),ylab="")
      hist(sL[,,j], main=paste("lambda", j, sep=""), breaks=40, ylab="", xlab="")
      #abline(v=L.median[j], col=2)
      #axis(side=1, at=L.median[j] , labels=round(L.median[j], d))
      abline(v=L.mode[j], col=4)
      axis(side=3, at=L.mode[j], labels=L.mode[j])
    }
    dev.off()
  }
  return(list(sL=sL, mon.L=mon.L, L.mode=L.mode, L.median=L.median))
}







omega_v_inference <- function(chains, B, p, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=3, log=T){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim+1, by=thin)
  l <- length(sel)
  somega_v <- array(NA, dim=c(l, n.chain, B^2*p))
  for(i in 1:n.chain){somega_v[,i,] <-  do.call("rbind", chains[[i]]$sim.omega_v[sel])}
  mon.omega_v <- rstan::monitor(somega_v, print=F)
  # omega_v.mode <-  apply(somega_v, 3, function(x){find.mode(x, d=-1)}) 
  omega_v.mode <-  apply(somega_v, 3, function(x){find.mode(x, d=d, log=log)}) 
  omega_v.median <- mon.omega_v[,6]
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*p)){
      matplot(somega_v[,,j], type="l", lty=1, log="y", main=paste("omega_v", j, sep=""), ylab="")
      hist(log(somega_v[,,j]), main=paste("log(omega_v", j, ")", sep=""), breaks=40, ylab="", xlab="")
      #abline(v=log(omega_v.median[j]), col=2)
      #axis(side=1, at=log(omega_v.median[j]), labels=round(log(omega_v.median[j]), 2))
      abline(v=log(omega_v.mode[j]), col=4)
      axis(side=3, at=log(omega_v.mode[j]), labels=round(log(omega_v.mode[j]), 2))
    }
    dev.off()
  }
  return(list(somega_v=somega_v, 
              mon.omega_v=mon.omega_v, 
              omega_v.mode=omega_v.mode, 
              omega_v.median=omega_v.median))
}





lambda1sq_inference <- function(chains, B, p,  n.sim, warmup=0, thin=1, filename, plot=FALSE){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim+1, by=thin)
  l <- length(sel)
  slambda1sq <- array(NA, dim=c(l, n.chain, B^2*p))
  for(i in 1:n.chain){slambda1sq[,i,] <-  do.call("rbind", chains[[i]]$sim.lambda1sq[sel])}
  mon.lambda1sq <- rstan::monitor(slambda1sq, print=F)
  lambda1sq.median <- mon.lambda1sq[,6]
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*p)){
      matplot(slambda1sq[,,j], type="l", lty=1, log="y", main=paste("lambda1sq", j, sep=""), ylab="")
      hist(log(slambda1sq[,,j]), main=paste("log(lambda1sq", j, ")", sep=""), breaks=40, ylab="", xlab="")
      abline(v=log(lambda1sq.median[j]), col=2)
      axis(side=1, at=log(lambda1sq.median[j]), labels=round(log(lambda1sq.median[j]), 2))
    }
    dev.off()
  }
  return(list(slambda1sq=slambda1sq, mon.lambda1sq=mon.lambda1sq, lambda1sq.median=lambda1sq.median))
}









lambda2_inference <- function(chains, B, p,  n.sim, warmup=0, thin=1, filename, plot=FALSE){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim+1, by=thin)
  l <- length(sel)
  slambda2<- array(NA, dim=c(l, n.chain, B^2*p))
  for(i in 1:n.chain){slambda2[,i,] <-  do.call("rbind", chains[[i]]$sim.lambda2[sel])}
  mon.lambda2<- rstan::monitor(slambda2, print=F)
  lambda2.median <- mon.lambda2[,6]
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*p)){
      matplot(slambda2[,,j], type="l", lty=1, log="y", main=paste("lambda2_", j, sep=""), ylab="")
      hist(log(slambda2[,,j]), main=paste("log(lambda2_", j, ")", sep=""), breaks=40, ylab="", xlab="")
      abline(v=log(lambda2.median[j]), col=2)
      axis(side=1, at=log(lambda2.median[j]), labels=round(log(lambda2.median[j]), 2))
    }
    dev.off()
  }
  return(list(slambda2=slambda2, mon.lambda2=mon.lambda2, lambda2.median=lambda2.median))
}








tausq_inference <- function(chains, B, p,  n.sim, warmup=0, thin=1, filename, plot=FALSE){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim+1, by=thin)
  l <- length(sel)
  stausq <- array(NA, dim=c(l, n.chain, B^2*p))
  for(i in 1:n.chain){stausq[,i,] <-  do.call("rbind", chains[[i]]$sim.tausq[sel])}
  mon.tausq <- rstan::monitor(stausq, print=F)
  tausq.median <- mon.tausq[,6]
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*p)){
      matplot(stausq[,,j], type="l", lty=1, log="y", main=paste("tausq", j, sep=""), ylab="")
      hist(log(stausq[,,j]), main=paste("log(tausq", j, ")", sep=""), breaks=40, ylab="", xlab="")
      abline(v=log(tausq.median[j]), col=2)
      axis(side=1, at=log(tausq.median[j]), labels=round(log(tausq.median[j]), 2))
    }
    dev.off()
  }
  return(list(stausq=stausq, mon.tausq=mon.tausq, tausq.median=tausq.median))
}








find.mode <- function(x, d=3, log=F){
  
  # y <- as.numeric(names(table(round(x, d)))[which.max(table(round(x, d)))])   # using rounding
  
  # hh <- hist(x, breaks=30, plot=F)
  # h <- hh$breaks[-length(hh$breaks)] + diff(hh$breaks)/2   # midpoints
  # y <- h[which.max(hh$counts)]
   
  if(log){
    hh <- density(log(x))
    y <- exp(hh$x[which.max(hh$y)])
    
  } else {
    hh <- density(x)
    y <- hh$x[which.max(hh$y)]
  }

   
   
  return(y)
  
}





# AIC

getAIC <- function(sw=sw, sv=sv, N=N, B=B, p=p, L_hat_2S_Gibbs=L_hat_2S_Gibbs){
  
  w_Bayes <- vector("list", N)
  for(i in 1:N){
    w_Bayes[[i]] <- apply(sw$sw + sv$sv[,,((i-1)*B^2*p+1):(i*B^2*p)], 3, function(x){find.mode(x, d=3)})
  }
  
  W <- vector("list", N)
  for(n in  1:N){W[[n]] <-  matrix(w_Bayes[[n]], B)}
  Lambda <- matrix(L_hat_2S_Gibbs, B, B) 
  
  k <- (N+6)*B^2*p+B^2
  
  2*k + sum(TT-p)*log(2*pi) - sum(TT-p)* log(det(Lambda)) + Reduce("+", lapply(1:N,  function(x){
    e <- eta0[[x]] - W[[x]]%*%H[[x]]    # residual for nj: B*(TT-p)
    sum(apply(e, 2, function(et){sum (et * (Lambda %*% et))}))    # Sum of t(et) %*% Lambda %*% et 
  }))
  
}



