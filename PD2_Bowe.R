##################
# PD 2           #
##################

library(date)
library(stats)
library(imputeTS)
library(vars)




rm(list=ls())

setwd("C:/Columbia/projects/David/PD/PLOS")
source("MS_hierarchical.R")
source("MS_AnalyzeResults.R")




############################## --------- import RData here -----------#######################

# 1. load data
load("PD2_Bowe.Rdata")


p <- 1                              # VAR order
n.test <- 10      # number of test data points: default=5 days
var.sel <- 1:6     # all 6 vars



-

############################## --------- posterior inference -----------#######################


# (1) w
sw <- w_inference(chains=list(sims1, sims2, sims3, sims4), B=B, p=p, n.sim=n.sim, warmup=n.sim/2+1, thin=1, filename="PD2_Bowe_w_trace.pdf", plot=F, d=3)
sw$mon.w                    # all converged
w_hat_2S_Gibbs <- sw$w.mode
sw_mat <-rbind(sw$sw[,1,], sw$sw[,2,], sw$sw[,3,], sw$sw[,4,])



tiff("Fig7s.tiff", height=10, width=10, res=300, units="in")
par(mfrow=c(3, 3))
for(j in 1:9) {
  matplot(sw$sw[,,j], type="l", lty=1, main=paste("w", j, sep=""), ylab="")
}
dev.off()


tiff("Fig8s.tiff", width=10, height=10, res=300, units="in")
par(mfrow=c(6, 6), mar=c(3, 2, 4, 2))
for(i in 1:(9)){
  for(j in 1:4){
    acf(sw$sw[,j,i], type="correlation", plot=T, main=paste("w_", i, " chain ", j, sep=""), xlab="", ylab="")
  }
}
dev.off()



# (2) v
sv <- v_inference(chains=list(sims1, sims2, sims3, sims4), B=B, p=p, N=N, n.sim=n.sim, warmup=n.sim/2+1, thin=1, filename="PD2_Bowe_v_trace.pdf", plot=F, d=3)
sv$mon.v[order(sv$mon.v[,9], decreasing=TRUE),][1:20,]    # all converged
v_hat_2S_Gibbs <- sv$v.mode
sv_mat <-rbind(sv$sv[,1,], sv$sv[,2,], sv$sv[,3,], sv$sv[,4,])



# (3) Lambda
sL <- L_inference(chains=list(sims1, sims2, sims3, sims4), B=B, p=p, n.sim=n.sim, warmup=n.sim/2+1, thin=1, filename="L_trace.pdf", plot=F, d=3)
sL$mon.L[order(sL$mon.L[,10], decreasing=TRUE),]
L_hat_2S_Gibbs <- sL$L.mode
sL_mat <-rbind(sL$sL[,1,], sL$sL[,2,], sL$sL[,3,], sL$sL[,4,])



# (4) omega_v
somega_v <- omega_v_inference(chains=list(sims1, sims2, sims3, sims4), B=B, p=p, n.sim=n.sim, warmup=n.sim/2+1, thin=1, filename="omega_v.pdf", plot=F, d=3)
somega_v$mon.omega_v
omega_v_hat_2S_Gibbs <- somega_v$omega_v.mode
sdv_hat_2S_Gibbs <- apply(sqrt(1/somega_v$somega_v), 3, function(x){find.mode(x, d=3)})




# AIC
getAIC(sw=sw, sv=sv, N=N, B=B, p=p, L_hat_2S_Gibbs=L_hat_2S_Gibbs)

R_hat <- c(sw$mon.w[,9], sv$mon.v[,9], sL$mon.L[,10], somega_v$mon.omega_v[,10])
hist(R_hat, breaks=30)
range(R_hat)  # 0.992814 1.018621


# save image for later use
save.image("PD2_Bowe_fit.Rdata")




############################## --------- prediction -----------#######################

load("pred6_Bowe_p1.RData")

a <- 0.05


##################### ---------  (1) use Bayesian to predict   -----------#######################


# recursive model matrix for Bayes and MLE
eta.pred.Bayes <- vector("list", N)

for(n in 1:N){
  #eta.pred.MLE[[n]] <- rbind(eta[[n]], matrix(NA, ncol=B, nrow=n.test, dimnames=list(NULL, vars[var.sel])))
  eta.pred.Bayes[[n]] <- vector("list", 1000)
  for(j in 1:1000){
    eta.pred.Bayes[[n]][[j]] <- rbind(eta[[n]], matrix(NA, ncol=B, nrow=n.test, dimnames=list(NULL, vars[var.sel])))
  }
}



pred.Bayes <- array(NA, dim=c(n.test, B, 1000, N))       # store all predictive samples
set.seed(1)
for(n in 1:N){
  for(j in 1:1000){     # 1000 repetitions
    temp <- sw_mat[j,] + sv_mat[j,((n-1)*(B^2*p)+1):(n*B^2*p)]     # A^hat^(j)
    Sigma <- solve(matrix(sL_mat[j,], B))                          # Psi^(j)

    for(h in 1:n.test){
      H.regressor.Bayes <- c(t(eta.pred.Bayes[[n]][[j]][(TT[[n]]+h-1):(TT[[n]]+h-p),]))
      linear <- c(matrix(temp, B)%*% H.regressor.Bayes)
      error <- c(rmvn(n=1, mu=rep(0, B), sigma=Sigma, ncores = 4))
      pred <- linear + error
      eta.pred.Bayes[[n]][[j]][TT[[n]]+h,] <- pred
      pred.Bayes[h,,j,n]  <- pred
    }
  }
  cat("n=", n, "\n")
}


pred.Bayes.hat <- vector("list", N)     # point estimate: posterior mean
pred.Bayes.lower <- vector("list", N)   # lower quantile
pred.Bayes.upper <- vector("list", N)   # upper quantile
for(n in 1:N){
  for(h in 1:n.test){
    #pred.B.hat[[n]] <- rbind(pred.B.hat[[n]], apply(pred.Bayes[h,,,n], 1, median))
    pred.Bayes.hat[[n]] <- rbind(pred.Bayes.hat[[n]], apply(pred.Bayes[h,,,n], 1, mean))
    pred.Bayes.lower[[n]] <- rbind(pred.Bayes.lower[[n]], apply(pred.Bayes[h,,,n], 1, function(x){quantile(x, a/2)}))
    pred.Bayes.upper[[n]] <- rbind(pred.Bayes.upper[[n]], apply(pred.Bayes[h,,,n], 1, function(x){quantile(x, 1-a/2)}))
  }
}



pred.Bayes.length <- array(NA, dim=c(n.test, B, N))   # interval length
for(n in 1:N){
  pred.Bayes.length[,,n] <- pred.Bayes.upper[[n]] - pred.Bayes.lower[[n]]
}


# Bayesian: MSE and coverage
MSE1 <- mean(unlist((lapply(1:N, function(n){(pred.Bayes.hat[[n]]-eta.test[[n]])^2}))))
C1 <- mean(unlist(lapply(1:N, function(n){(pred.Bayes.lower[[n]]<=eta.test[[n]])&(pred.Bayes.upper[[n]]>=eta.test[[n]])})))





###################### ---------  (2). use subject specific MLE -----------#######################
# predictives

eta.pred.MLE <- vector("list", N)
for(n in 1:N){
  eta.pred.MLE[[n]] <- rbind(eta[[n]], matrix(NA, ncol=B, nrow=n.test, dimnames=list(NULL, vars[var.sel])))
}
# the same as VAR() function
pred.MLE <- array(NA, dim=c(n.test, B, N))
for(n in 1:N){
  for(h in 1:n.test){
    H.regressor.MLE <- c(t(eta.pred.MLE[[n]][(TT[[n]]+h-1):(TT[[n]]+h-p),]))
    eta.pred.MLE[[n]][TT[[n]]+h,] <- matrix(w_MLE[[n]], B) %*% H.regressor.MLE
    pred.MLE[h,,n] <-  matrix(w_MLE[[n]], B) %*% H.regressor.MLE

  }
}





# use VAR
pred.MLE.hat <- array(NA, dim=c(n.test, B, N))
pred.MLE.lower <- array(NA, dim=c(n.test, B, N))   # lower quantile
pred.MLE.upper <- array(NA, dim=c(n.test, B, N))   # upper quantile
pred.MLE.length <- array(NA, dim=c(n.test, B, N))   # interval length

for(n in 1:N){

  fit <- VAR(eta[[n]], p=p, type="none")
  pred <- predict(fit, n.ahead=10, ci=1-a)

  for(j in 1:length(var.sel)){

    pred.MLE.hat[,j,n] <- pred$fcst[[j]][,1]
    pred.MLE.lower[,j,n] <- pred$fcst[[j]][,2]
    pred.MLE.upper[,j,n] <- pred$fcst[[j]][,3]
    pred.MLE.length[,j,n] <-  pred.MLE.upper[,j,n]- pred.MLE.lower[,j,n]
  }



}



# MSE: MLE & coverage
MSE2 <- mean(unlist(lapply(1:N, function(n){(pred.MLE.hat[,,n] - eta.test[[n]])^2})))
C2 <- mean(unlist(lapply(1:N, function(n){(pred.MLE.lower[,,n]<=eta.test[[n]])&(pred.MLE.upper[,,n]>=eta.test[[n]])})))





####################### ---------   (3) use glmnet  ------------- ##########

pred.glmnet.hat <- array(NA, dim=c(n.test, B, N))
pred.glmnet.lower <- array(NA, dim=c(n.test, B, N))   # lower quantile
pred.glmnet.upper <- array(NA, dim=c(n.test, B, N))   # upper quantile
pred.glmnet.length <- array(NA, dim=c(n.test, B, N))   # interval length


set.seed(12345)
for(i in 1:N){

  temp <- data.frame(eta[[i]])

  train1 <- cbind(temp[-1, 1], temp[-TT[i], ])     # first column is Y=V1,
  train2 <- cbind(temp[-1, 2], temp[-TT[i], ])     # first column is Y=V2,
  train3 <- cbind(temp[-1, 3], temp[-TT[i], ])     # first column is Y=V3,
  train4 <- cbind(temp[-1, 4], temp[-TT[i], ])     # first column is Y=V4,
  train5 <- cbind(temp[-1, 5], temp[-TT[i], ])     # first column is Y=V5,
  train6 <- cbind(temp[-1, 6], temp[-TT[i], ])     # first column is Y=V6,


  glmnet.fit1 <- cv.glmnet(y=train1[, 1], x=as.matrix(train1[, -1]), family="gaussian", nfolds=3, alpha=0.5)
  glmnet.fit2 <- cv.glmnet(y=train2[, 1], x=as.matrix(train2[, -1]), family="gaussian", nfolds=3, alpha=0.5)
  glmnet.fit3 <- cv.glmnet(y=train3[, 1], x=as.matrix(train3[, -1]), family="gaussian", nfolds=3, alpha=0.5)
  glmnet.fit4 <- cv.glmnet(y=train4[, 1], x=as.matrix(train4[, -1]), family="gaussian", nfolds=3, alpha=0.5)
  glmnet.fit5 <- cv.glmnet(y=train5[, 1], x=as.matrix(train5[, -1]), family="gaussian", nfolds=3, alpha=0.5)
  glmnet.fit6 <- cv.glmnet(y=train6[, 1], x=as.matrix(train6[, -1]), family="gaussian", nfolds=3, alpha=0.5)


  # 1-step forecast
  pred.glmnet.hat[1,1,i] <- predict(glmnet.fit1, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
  pred.glmnet.hat[1,2,i] <- predict(glmnet.fit2, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
  pred.glmnet.hat[1,3,i] <- predict(glmnet.fit3, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
  pred.glmnet.hat[1,4,i] <- predict(glmnet.fit4, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
  pred.glmnet.hat[1,5,i] <- predict(glmnet.fit5, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
  pred.glmnet.hat[1,6,i] <- predict(glmnet.fit6, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")



  # h step forecast, h>=2
  for(h in 2:n.test){

    pred.glmnet.hat[h,1,i] <- predict(glmnet.fit1, newx=t(as.matrix(pred.glmnet.hat[h-1,,i])), s="lambda.1se")
    pred.glmnet.hat[h,2,i] <- predict(glmnet.fit2, newx=t(as.matrix(pred.glmnet.hat[h-1,,i])), s="lambda.1se")
    pred.glmnet.hat[h,3,i] <- predict(glmnet.fit3, newx=t(as.matrix(pred.glmnet.hat[h-1,,i])), s="lambda.1se")
    pred.glmnet.hat[h,4,i] <- predict(glmnet.fit4, newx=t(as.matrix(pred.glmnet.hat[h-1,,i])), s="lambda.1se")
    pred.glmnet.hat[h,5,i] <- predict(glmnet.fit5, newx=t(as.matrix(pred.glmnet.hat[h-1,,i])), s="lambda.1se")
    pred.glmnet.hat[h,6,i] <- predict(glmnet.fit6, newx=t(as.matrix(pred.glmnet.hat[h-1,,i])), s="lambda.1se")


  }

}



# use bootstrap to get CI

n.boots <- 1000
pred.glmnet.hat.boots <- array(NA, dim=c(n.test, B, N, n.boots))

for(b in 1:n.boots){


  for(i in 1:N){

    temp <- data.frame(eta[[i]][sample(1:TT[i], replace=T),])     # a bootstrap sample

    train1 <- cbind(temp[-1, 1], temp[-TT[i], ])     # first column is Y=V1, others are X
    train2 <- cbind(temp[-1, 2], temp[-TT[i], ])     # first column is Y=V2, others are X
    train3 <- cbind(temp[-1, 3], temp[-TT[i], ])     # first column is Y=V3, others are X
    train4 <- cbind(temp[-1, 4], temp[-TT[i], ])     # first column is Y=V4, others are X
    train5 <- cbind(temp[-1, 5], temp[-TT[i], ])     # first column is Y=V5, others are X
    train6 <- cbind(temp[-1, 6], temp[-TT[i], ])     # first column is Y=V6, others are X


    glmnet.fit1 <- cv.glmnet(y=train1[, 1], x=as.matrix(train1[, -1]), family="gaussian", nfolds=3, alpha=0.5)
    glmnet.fit2 <- cv.glmnet(y=train2[, 1], x=as.matrix(train2[, -1]), family="gaussian", nfolds=3, alpha=0.5)
    glmnet.fit3 <- cv.glmnet(y=train3[, 1], x=as.matrix(train3[, -1]), family="gaussian", nfolds=3, alpha=0.5)
    glmnet.fit4 <- cv.glmnet(y=train4[, 1], x=as.matrix(train4[, -1]), family="gaussian", nfolds=3, alpha=0.5)
    glmnet.fit5 <- cv.glmnet(y=train5[, 1], x=as.matrix(train5[, -1]), family="gaussian", nfolds=3, alpha=0.5)
    glmnet.fit6 <- cv.glmnet(y=train6[, 1], x=as.matrix(train6[, -1]), family="gaussian", nfolds=3, alpha=0.5)


    # 1-step forecast
    pred.glmnet.hat.boots[1,1,i, b] <- predict(glmnet.fit1, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
    pred.glmnet.hat.boots[1,2,i, b] <- predict(glmnet.fit2, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
    pred.glmnet.hat.boots[1,3,i, b] <- predict(glmnet.fit3, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
    pred.glmnet.hat.boots[1,4,i, b] <- predict(glmnet.fit4, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
    pred.glmnet.hat.boots[1,5,i, b] <- predict(glmnet.fit5, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")
    pred.glmnet.hat.boots[1,6,i, b] <- predict(glmnet.fit6, newx=as.matrix(temp[TT[i], ]), s="lambda.1se")



    # h step forecast, h>=2
    for(h in 2:n.test){

      pred.glmnet.hat.boots[h,1,i,b] <- predict(glmnet.fit1, newx=t(as.matrix(pred.glmnet.hat.boots[h-1,,i, b])), s="lambda.1se")
      pred.glmnet.hat.boots[h,2,i,b] <- predict(glmnet.fit2, newx=t(as.matrix(pred.glmnet.hat.boots[h-1,,i, b])), s="lambda.1se")
      pred.glmnet.hat.boots[h,3,i,b] <- predict(glmnet.fit3, newx=t(as.matrix(pred.glmnet.hat.boots[h-1,,i, b])), s="lambda.1se")
      pred.glmnet.hat.boots[h,4,i,b] <- predict(glmnet.fit4, newx=t(as.matrix(pred.glmnet.hat.boots[h-1,,i, b])), s="lambda.1se")
      pred.glmnet.hat.boots[h,5,i,b] <- predict(glmnet.fit5, newx=t(as.matrix(pred.glmnet.hat.boots[h-1,,i, b])), s="lambda.1se")
      pred.glmnet.hat.boots[h,6,i,b] <- predict(glmnet.fit6, newx=t(as.matrix(pred.glmnet.hat.boots[h-1,,i, b])), s="lambda.1se")


    }

  }

  cat("b=", b, "\n")

}





for(i in 1:N){
  for(h in 1:n.test){
    pred.glmnet.lower[h,,i] <- apply(pred.glmnet.hat.boots[h,,i,], 1, function(x){quantile(x, probs=0.025)})
    pred.glmnet.upper[h,,i] <- apply(pred.glmnet.hat.boots[h,,i,], 1, function(x){quantile(x, probs=0.975)})
    pred.glmnet.length[h,,i] <- pred.glmnet.upper[h,,i] - pred.glmnet.lower[h,,i]
  }
}







# MSE & coverage
MSE3 <- mean(unlist(lapply(1:N, function(n){(pred.glmnet.hat[,,n] - eta.test[[n]])^2})))
C3 <- mean(unlist(lapply(1:N, function(n){(pred.glmnet.lower[,,n]<=eta.test[[n]])&(pred.glmnet.upper[,,n]>=eta.test[[n]])})))





# % reduction in MSE
MSE1
MSE2
MSE3
C1
C2
C3
1-MSE1/MSE2                                      # % reduction in MSE
1-MSE1/MSE3                                      # % reduction in MSE
mean(pred.Bayes.length/pred.MLE.length < 1)      # % shorter interval
mean(pred.Bayes.length/pred.glmnet.length < 1)   # % shorter interval


save.image("pred6_Bowe_p1.RData")

