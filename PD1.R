
library(vars)
library(mvnfast)
library(png)
library(geepack)
library(doBy)
library(vioplot)
library(fields)
library(gplots)
library(lme4)
library(RColorBrewer)
library(colorspace)
library(glmnet)


setwd("C:/Users/ethan/Documents/MEE/Research/Rcode")


source("MS_hierarchical.R")
source("MS_AnalyzeResults.R")

# Specify which columns are ID and time series. The rest should be variable data.
ID_COL <- 1
TIME_COL <- 2


## ============== input data ================== ###
dat <- read.csv("snap_data4.csv")
#dat <- read.csv("diary1.csv")
ID <- read.csv("ID1.csv")
length(unique(dat$ID))
dat_L <- dim(dat)[1]



## ============== preprocessing  ================== ###
# log transform
dat2 <- dat
print(dat2)
col_names <- colnames(dat)
names(dat2)[ID_COL] <- "ID"
names(dat2)[TIME_COL] <- "day"
var_cols <- col_names[-c(ID_COL, TIME_COL)]
dat2
#for(col_n in var_cols){dat2[col_n] <- log(dat[col_n]+1) }

for(col_n in var_cols){dat2[col_n] <- log(dat[col_n] + 1 + runif(dat_L,0, 10^-4)) } # + random data for testing purposes 
                                                                                   # Change max value of unif. dist. if crossprod. matrix isn't posdef.

# demean, detrend, scale

# REFACTOR: Cols 3 to Ncols - Nfeatures + 1
B <- length(var_cols)

Y1 <- split(dat2, dat2$ID) # Split dataframe into a list of data based on user ID
#sapply(Y1, function(x){range(x$day)})
cnt <- 1
eta.full <- lapply(Y1, function(x){apply(x[,3:(3 + B - 1)], 2, function(y){      # demean, detrend, scale 
  #print(x)
  yy <- scale(y) # Demean and scale variance to 1 of each column for each user's data 
  # Fit a linear model from the normalized data time series data
  # Replace columns with their residuals/error with linear trend. 
  # This is to make the data "more" i.i.d. Standard with Autoregression.
  fit <- lm(yy~x$day) 
  cnt <- cnt+1
  return(fit$residuals)
})})


# set training and set data``
N <- length(eta.full) 
eta <- lapply(eta.full, function(x){x[-nrow(x),]})      # training data - Last day of data removed
eta.test <- sapply(eta.full, function(x){x[nrow(x),]})  # test data - Last day of data for each ID

# prepare for modeling fitting
TT <- unlist(lapply(eta, nrow))    # a vector of training days (Number of days for each ID)
summary(TT)
sd(TT)
# REFACTOR: CHANGE B to fit Nfeat
#B <- 3 # Number of features 
p <- 2 # p-day lag VAR Model                     # when p=2, No.14 ID=100518 is nearly singular

# Create the H-matrix organizing the data into time lag-structure.
# B*p rows. TT - p cols. First B rows correspond to data with lag p, next B rows to data with p-1 and so on.
H <- lapply(1:N, function(n){H <- NULL; for(i in 1:(TT[[n]]-p)){H <- cbind(H, c(t(eta[[n]][(i+p-1):i,])))}; H})

H.test <- lapply(1:N, function(n){c(t(eta[[n]][(TT[[n]]):(TT[[n]]-p+1),]))})
eta0 <- lapply(1:N, function(n){t(eta[[n]][(p+1):TT[[n]],])}) # Removes the first p rows from each ID's dataframe
HH <- lapply(H, tcrossprod) # Crossproduct matrix of H
# Calculate the least-squared, maximum likelihood weights for the autoregressive model. i.e. solve: H*w_MLE = eta0 for each ID.
w_MLE <- lapply(1:N, function(n){print(HH[[n]]); c(eta0[[n]] %*% t(H[[n]]) %*% chol2inv(chol(HH[[n]])))}) # 



## ============== fit Bayesian model  ================== ###
# initial values
set.seed(1)
w_hat_2S <- Reduce("+", w_MLE)/N # Averages the vectors of w_MLE 
# Creates v_hat_2s: Subtracts the average of w_MLE from every vector in w_MLE
v_hat_2S <- vector("list", N) 
for(n in 1:N) v_hat_2S[[n]] <- w_MLE[[n]] - w_hat_2S 
# Create the list describing the initial conditions for the Monte Carlo simulation
# Initializes many of the variables used in the prior distributions
initial = list(w=w_hat_2S, v_n=v_hat_2S,
               Lambda=Posdef(n=B, ev=runif(B, 0, 0.001)), #Positive definitematrix
               lambda1sq = runif(B^2*p, 0.1, 1000),
               lambda2=runif(B^2*p, 0.1, 10),
               tausq =runif(B^2*p, 0.1, 10),
               U_v = runif(B^2*p, 0.1, 1),
               alpha = runif(B^2*p, 0.1, 1))
initial$v_n.star <- lapply(1:N, function(n){initial$alpha*initial$v_n[[n]]})
initial$omega_v <- initial$U_v/(initial$alpha^2)


# Create sample distributions using Markov Chain Monte Carlo (MCMC)
# run 4 chains
n.iter <- 10000
n.thin <- 20 # Sample every 20 simulated data points
n.sim <- n.iter/n.thin 

sims1 <- MSGibbs(eta=eta, H=H, HH=HH, eta0=eta0, w_MLE=w_MLE, TT=TT, B=B, p=p, N=N, seed=1800, numIter=n.iter, thin=n.thin, warmup=0, initial=initial)

sims2 <- MSGibbs(eta=eta, H=H, HH=HH, eta0=eta0, w_MLE=w_MLE, TT=TT, B=B, p=p, N=N, seed=2000, numIter=n.iter, thin=n.thin, warmup=0, initial=initial)

sims3 <- MSGibbs(eta=eta, H=H, HH=HH, eta0=eta0, w_MLE=w_MLE, TT=TT, B=B, p=p, N=N, seed=4000, numIter=n.iter, thin=n.thin, warmup=0, initial=initial)

sims4 <- MSGibbs(eta=eta, H=H, HH=HH, eta0=eta0, w_MLE=w_MLE, TT=TT, B=B, p=p, N=N, seed=6000, numIter=n.iter, thin=n.thin, warmup=0, initial=initial)


save.image("PD1.RData")

# P is good up to here?

### =================== posterior inference ====================== ###

load("PD1.RData")

# (1) w

# Posterier inference for population level coeffs
sw <- w_inference(chains=list(sims1, sims2, sims3, sims4), B=B, p=p, n.sim=n.sim, warmup=n.sim/2+1, thin=1, filename="PD1_w_trace.pdf", plot=F, d=3)
# Obtain the weights for the population level coefficients via the posterier modes (max value)
w_hat_Bayes <- sw$w.mode
round(matrix(w_hat_Bayes, B),3)

# diagonostic plots

# REFACTOR: Need to change code to generate correct plots

# Plot showing the simulated sample for each variable for each chain 
tiff("Fig1s.tiff", height=10, width=10, res=300, units="in")
par(mfrow=c(3, 3))
for(j in 1:(B^2*p)) {
  matplot(sw$sw[,,j], type="l", lty=1, main=paste("w", j, sep=""), ylab="")
}
dev.off()

# Plot showing auto corellation for each variable for each sim chain
tiff("Fig2s.tiff", width=10, height=10, res=300, units="in")
par(mfrow=c(6, 6), mar=c(3, 2, 4, 2))
for(i in 1:(9)){
  for(j in 1:4){
    acf(sw$sw[,j,i], type="correlation", plot=T, main=paste("w_", i, " chain ", j, sep=""), xlab="", ylab="")
    }
  }
dev.off()



# (2) 
# Posterier inference for patient/ID level coeffs 
sv <- v_inference(chains=list(sims1, sims2, sims3, sims4), B=B, p=p, N=N, n.sim=n.sim, warmup=n.sim/2+1, thin=1, filename="PD1_v_trace.pdf", plot=F, d=3)
v_hat_Bayes <- sv$v.mode # List of population level coeffs obtained via maximum likelihood of distribution

# EXCLUDED DUE TO LACK OF GENDER/AGE INFO
if(FALSE){
  # distribution of vn by gender and gender
  
  # REFACTOR: No info on gender in SNAPSHOT, maybe remove or view distribution by other variable? 
  
  vn_mat <- as.data.frame(cbind(ID, matrix(round(v_hat_Bayes, 3), nrow=N, byrow=T)))
  vn_mat$age.cat <- as.numeric(vn_mat$age<=21)    # 3 groups: young female, young male and old male
  vn_mat <- vn_mat[order(vn_mat$age.cat),]
  
  # plot 1: by association
  tiff("Fig4.tiff", height=7.5, width=7.5, res=300, units="in")
  main <- c("T=>T", "T=>N", "T=>C", "N=>T", "N=>N", "N=>C", "C=>T", "C=>N", "C=>C")
  par(mfrow=c(3, 3))
  for(i in 4:12){
    if(i==4){ylim=range(vn_mat[,i])} else {ylim=c(-0.022, 0.022)}
    boxplot(list(vn_mat[vn_mat$age.cat==1 & vn_mat$female==1,i],
                 #vn_mat[vn_mat$age.cat==0 & vn_mat$female==1,i],
                 vn_mat[vn_mat$age.cat==1 & vn_mat$female==0,i],
                 vn_mat[vn_mat$age.cat==0 & vn_mat$female==0,i]), ylim=ylim, col=2:5, main=main[i-3])
  
  
  }
  dev.off()
  
  
  # plot 2: by group
  
  
  tiff("Fig4.tiff", height=4, width=9, res=300, units="in")
  
  par(mfrow=c(1, 3), las=2, oma=c(2, 2, 2, 2))
  
  col <- rainbow_hcl(n=9)
  boxplot(vn_mat[vn_mat$age.cat==1 & vn_mat$female==1, 4:12], ylim=c(-0.123, 0.05), col=col, xaxt="n")  
  title( main="Young female", line=3)
  polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
  abline(v=1:9, col="white", lty=3)
  abline(h=0, col="red", lty=1)
  box()
  boxplot(vn_mat[vn_mat$age.cat==1 & vn_mat$female==1, 4:12], col=col, xaxt="n",add=T)
  axis(side=1, at=1:9, labels=c("T=>T", "T=>N", "T=>C", "N=>T", "N=>N", "N=>C", "C=>T", "C=>N", "C=>C"))
  
  
  boxplot(vn_mat[vn_mat$age.cat==1 & vn_mat$female==0, 4:12], ylim=c(-0.13, 0.05), col=col, xaxt="n")
  title( main="Young male", line=3)
  polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
  abline(v=1:9, col="white", lty=3)
  abline(h=0, col="red", lty=1)
  box()
  boxplot(vn_mat[vn_mat$age.cat==1 & vn_mat$female==0, 4:12], col=col, xaxt="n",add=T)
  axis(side=1, at=1:9, labels=c("T=>T", "T=>N", "T=>C", "N=>T", "N=>N", "N=>C", "C=>T", "C=>N", "C=>C"))
  
  
  boxplot(vn_mat[vn_mat$age.cat==0 & vn_mat$female==0, 4:12], ylim=c(-0.13, 0.05), col=col, xaxt="n")
  title( main="Old male", line=3)
  polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
  abline(v=1:9, col="white", lty=3)
  abline(h=0, col="red", lty=1)
  box()
  boxplot(vn_mat[vn_mat$age.cat==0 & vn_mat$female==0, 4:12], col=col, xaxt="n",add=T)
  axis(side=1, at=1:9, labels=c("T=>T", "T=>N", "T=>C", "N=>T", "N=>N", "N=>C", "C=>T", "C=>N", "C=>C"))
  
  dev.off()
  
  
}#EXCLUDED
  
# (3) Lambda
sL <- L_inference(chains=list(sims1, sims2, sims3, sims4), B=B, p=p, n.sim=n.sim, warmup=n.sim/2+1, thin=1, filename="PD1_L_trace.pdf", plot=T, d=3) # d is for decimal rounding (not implemented in code, so value doesn't matter)
L_hat_Bayes <- sL$L.mode
  
# Used for evaluation of simulation variables. Not referenced again. Need explanation of what these vars are exactly.
if(FALSE){  
# (4) omega_v
somega_v <- omega_v_inference(chains=list(sims1, sims2, sims3, sims4), B=B, p=p, n.sim=n.sim, warmup=n.sim/2+1, thin=1, filename="PD1_theta_v_trace.pdf", plot=T, d=3)
sdv <- sqrt(1/somega_v$somega_v)
round(matrix(apply(sdv, 3, function(x){find.mode(x, d=3)}), B), B)   # posterior mode of SD of v_n # p-refactor ?
  
  
  
R_hat <- c(sw$mon.w[,9], sv$mon.v[,9], sL$mon.L[,10], somega_v$mon.omega_v[,10])
hist(R_hat, breaks=30)
range(R_hat)
}


### ===================  prediction ====================== #####
# (1). Bayesian predictive distribution

## REFACTOR: Change the dimensions to accomodate new features

y_hat_post_array <- array(NA, dim=c(1000, B, N))
y_hat_post <- matrix(NA, B, N)                         # posterior mean of y_hat
y_hat_post_quantile <- array(NA, dim=c(2, B, N))       # posterior quantile
coverage_post <- matrix(NA, nrow=N, ncol=B)            # coverage probability
y_hat_post_length <- NULL                              # interval length

for(i in 1:N){
  temp <- sw$sw + sv$sv[,,((i-1)*(p*B^2)+1):(i*p*B^2)]
  #temp <- sw$sw + sv$sv[,,((i-1)*9+1):(i*9)] # Combine population and patient level weights #Refactor to accomodate more than three vars and more than 1 time lag
  # temp's dim's are B x n_samples (from sims) x n_sim chains 
  linear <- apply(temp, c(1,2), function(x){matrix(x, B)%*% H.test[[i]]}) # Create H_n-hat mult by y_n (test samples)
  set.seed(123)
  error <- apply(sL$sL, c(1,2), function(x){Sigma <- chol2inv(chol(matrix(x, B))); rmvn(n=1, mu=rep(0, B), sigma=chol(Sigma), isChol=TRUE,  ncores = 4)}) # Create error terms for ID
  pred <- linear + error # Create predictions 
  
  # Predictions for each variable
  for(m in 1:B){y_hat_post_array[,m,i] <- pred[m,,]}

  # 95% Intervals for each variable
  for(m in 1:B){y_hat_post_quantile[,m,i] <- quantile(pred[m,,], probs=c(0.025, 0.975))}

  # Test to see if test data falls within 95% prediction intervals for each variable
  for(m in 1:B){coverage_post[i,m] <- eta.test[m,i] >=  y_hat_post_quantile[1,m,i] & eta.test[m,i] <=  y_hat_post_quantile[2,m,i]}

  y_hat_post[,i] <- apply(pred, 1, mean)    # posterior mean of y_hat
  #y_hat_post[,i] <- apply(pred, 1, function(x){find.mode(x, d=3)})    # posterior mode of y_hat for each variable and ID
  y_hat_post_length <- rbind(y_hat_post_length, apply(y_hat_post_quantile[,,i], 2, diff)) # Length of 95% confidence intervals

}


MSE_Bayes <- mean((eta.test-y_hat_post)^2)                        # MSE
MSE_Bayes_by_var <- apply((eta.test-y_hat_post)^2, 1, mean)       # MSE by variable
mean(coverage_post)                                               # coverage




# (2). Subject-specific MLE
# For comparison against MSE_Bayes.
y_hat_MLE_quantile <- array(NA,  dim=c(2, B, N))
coverage_MLE <-  matrix(NA, nrow=N, ncol=B)
y_hat_MLE_length <- NULL
for(n in 1:N){
  fit <- VAR(y=eta[[n]], p=p, type="none")
  pred <- predict(fit, n.ahead=p, ci=.95)[[1]]
  for(m in 1:B){y_hat_MLE_quantile[,m,n]<-pred[[m]][2:3]}
  
  for(m in 1:B){coverage_MLE[n,m] <- eta.test[m,n] >=  y_hat_MLE_quantile[1,m,n] & eta.test[m,n] <=  y_hat_MLE_quantile[2,m,n]}
  
  y_hat_MLE_length <- rbind(y_hat_MLE_length, apply(y_hat_MLE_quantile[,,n], 2, diff))
}
y_hat_MLE <- sapply(1:N, function(n){matrix(w_MLE[[n]], B) %*% H.test[[n]]})
MSE_MLE <- mean((eta.test-y_hat_MLE)^2)
MSE_MLE_by_var <- apply((eta.test-y_hat_MLE)^2, 1, mean)
mean(coverage_MLE)

round(do.call("rbind", v_hat_2S), 3)
matrix(round(v_hat_Bayes, 3), nrow=22, byrow=T) #REFACTOR



# Excluded regularized linear model from predictions due to irrelevance and non-generalized code
if(FALSE){

# (3.3)regularized linear model
# Other model for comparisons sake

glmnet.pred1 <- glmnet.pred2 <- glmnet.pred3 <- rep(NA, N)

train1 <- train2 <- train3 <- NULL
test1 <- test2 <- test3 <- NULL
for(i in 1:N){

  temp <- data.frame(eta[[i]])
  temp1 <- as.data.frame(cbind(temp[-1, 1], temp[-TT[i], ], ID=i, day=1:(TT[i]-1)))
  temp2 <- as.data.frame(cbind(temp[-1, 2], temp[-TT[i], ], ID=i, day=1:(TT[i]-1)))
  temp3 <- as.data.frame(cbind(temp[-1, 3], temp[-TT[i], ], ID=i, day=1:(TT[i]-1)))

  names(temp1)[1] <- "smoke"
  names(temp2)[1] <- "negativeaffect"
  names(temp3)[1] <- "cravings"

  train1 <- rbind(train1, temp1)
  train2 <- rbind(train2, temp2)
  train3 <- rbind(train3, temp3)

  test1 <- rbind(test1, unlist(c(smoke=eta.test[1,i], temp[TT[i],], ID=i)))
  test2 <- rbind(test2, unlist(c(negativeaffect=eta.test[2,i], temp[TT[i],], ID=i)))
  test3 <- rbind(test3, unlist(c(cravings=eta.test[3,i], temp[TT[i],], ID=i)))

}


test1 <- as.data.frame(test1)
test2 <- as.data.frame(test2)
test3 <- as.data.frame(test3)

set.seed(12345)

coef_glmnet <- matrix(NA, nrow=N, ncol=9)

for(id in 1:N){

  glmnet.fit1 <- cv.glmnet(y=train1[train1$ID==id, 1], x=as.matrix(train1[train1$ID==id, 2:4]), family="gaussian", nfolds=3, alpha=0.5)
  glmnet.fit2 <- cv.glmnet(y=train2[train1$ID==id, 1], x=as.matrix(train1[train1$ID==id, 2:4]), family="gaussian", nfolds=3, alpha=0.5)
  glmnet.fit3 <- cv.glmnet(y=train3[train1$ID==id, 1], x=as.matrix(train1[train1$ID==id, 2:4]), family="gaussian", nfolds=3, alpha=0.5)

  glmnet.pred1[id] <- predict(glmnet.fit1, newx=as.matrix(test1[test1$ID==id, 2:4]), s="lambda.1se")
  glmnet.pred2[id] <- predict(glmnet.fit2, newx=as.matrix(test2[test2$ID==id, 2:4]), s="lambda.1se")
  glmnet.pred3[id] <- predict(glmnet.fit3, newx=as.matrix(test3[test3$ID==id, 2:4]), s="lambda.1se")
  
  coef_glmnet[id,] <- c(coef(glmnet.fit1, s="lambda.1se")[2:4], coef(glmnet.fit2, s="lambda.1se")[2:4], coef(glmnet.fit3, s="lambda.1se")[2:4])
  
}

round(coef_glmnet, 3)


# use bootstrap for prediction intervals

n.boots <- 1000
set.seed(12345)
glmnet.pred1.boots <- glmnet.pred2.boots <- glmnet.pred3.boots <- matrix(NA, nrow=N, ncol=n.boots)


for(b in 1:n.boots){
  for(id in 1:N){
    train1.boots <- train1[train1$ID==id,][sample(1:(TT[id]-1), replace = T),]
    train2.boots <- train2[train2$ID==id,][sample(1:(TT[id]-1), replace = T),]
    train3.boots <- train3[train3$ID==id,][sample(1:(TT[id]-1), replace = T),]

    try(glmnet.fit1.boots <- cv.glmnet(y=train1.boots[train1.boots$ID==id, 1], x=as.matrix(train1.boots[train1.boots$ID==id, 2:4]), family="gaussian", nfolds=3, alpha=0.5), silent=T)
    try(glmnet.fit2.boots <- cv.glmnet(y=train2.boots[train1.boots$ID==id, 1], x=as.matrix(train1.boots[train1.boots$ID==id, 2:4]), family="gaussian", nfolds=3, alpha=0.5), silent=T)
    try(glmnet.fit3.boots <- cv.glmnet(y=train3.boots[train1.boots$ID==id, 1], x=as.matrix(train1.boots[train1.boots$ID==id, 2:4]), family="gaussian", nfolds=3, alpha=0.5), silent=T)

    try(glmnet.pred1.boots[id, b] <- predict(glmnet.fit1.boots, newx=as.matrix(test1[test1$ID==id, 2:4]), s="lambda.1se"), silent=T)
    try(glmnet.pred2.boots[id, b] <- predict(glmnet.fit2.boots, newx=as.matrix(test2[test2$ID==id, 2:4]), s="lambda.1se"), silent=T)
    try(glmnet.pred3.boots[id, b] <- predict(glmnet.fit3.boots, newx=as.matrix(test3[test3$ID==id, 2:4]), s="lambda.1se"), silent=T)
  }

   cat("b=", b, "\n")
}


y_hat_glmnet_lower1 <- apply(glmnet.pred1.boots, 1, function(x){quantile(x, na.rm=T, probs=0.025)})
y_hat_glmnet_lower2 <- apply(glmnet.pred2.boots, 1, function(x){quantile(x, na.rm=T, probs=0.025)})
y_hat_glmnet_lower3 <- apply(glmnet.pred3.boots, 1, function(x){quantile(x, na.rm=T, probs=0.025)})


y_hat_glmnet_upper1 <- apply(glmnet.pred1.boots, 1, function(x){quantile(x, na.rm=T, probs=0.975)})
y_hat_glmnet_upper2 <- apply(glmnet.pred2.boots, 1, function(x){quantile(x, na.rm=T, probs=0.975)})
y_hat_glmnet_upper3 <- apply(glmnet.pred3.boots, 1, function(x){quantile(x, na.rm=T, probs=0.975)})



MSE_glmnet <- mean(c((glmnet.pred1-test1[,1])^2, (glmnet.pred2-test2[,1])^2, (glmnet.pred3-test3[,1])^2))
MSE_glmnet_by_var <- c(mean((glmnet.pred1-test1[,1])^2),
                      mean((glmnet.pred2-test2[,1])^2),
                      mean((glmnet.pred3-test3[,1])^2))




coverage_glmnet <- mean(c(y_hat_glmnet_lower1 < test1[,1] & y_hat_glmnet_upper1 > test1[,1],
                          y_hat_glmnet_lower2 < test2[,1] & y_hat_glmnet_upper2 > test2[,1],
                          y_hat_glmnet_lower3 < test3[,1] & y_hat_glmnet_upper3 > test3[,1]))


y_hat_glmnet_length <- cbind(y_hat_glmnet_upper1-y_hat_glmnet_lower1, y_hat_glmnet_upper2-y_hat_glmnet_lower2, y_hat_glmnet_upper3-y_hat_glmnet_lower3)

}
save.image("PD1_pred.RData")



# (4) compare Bayesian prediction with the other 2 methods
# with subject-specific MLE
MSE_Bayes
MSE_MLE
#MSE_glmnet

MSE_Bayes_by_var
MSE_MLE_by_var
#MSE_glmnet_by_var


MSE_Bayes/MSE_MLE-1                       # % reduction in MSE
#MSE_Bayes/MSE_glmnet-1 


MSE_Bayes_by_var/MSE_MLE_by_var -1        # % reduction in MSE by variables
#MSE_Bayes_by_var/MSE_glmnet_by_var -1 
 
mean(coverage_post)
mean(coverage_MLE)
#coverage_glmnet


mean(y_hat_post_length<y_hat_MLE_length)  # % shorter interval
#mean(y_hat_post_length<y_hat_glmnet_length)




# Fig 5
tiff("Fig5.tiff", res=300, width=10, height=12, units="in")
par(mfrow=c(5, 5))
for(n in 1:N){
  plot(x=1:3+0.1, y=eta.test[,n], pch="*", col=2, cex=2, xlim=c(0.7, 3.3), ylim=range(c(c(y_hat_MLE_quantile[,,n]), c(y_hat_post_quantile[,,n]), eta.test[,n])), xlab="", ylab="", xaxt="n", type="n", main=paste("subject", n))
  title(main=paste("Tn","=", nrow(eta[[n]]), sep=""), line=0.5)
  abline(a=0, b=0, col="grey", lty=2)
  axis(side=1, at=1:3, labels=c("T", "N", "C"))
  points(x=1:3, y=y_hat_MLE[,n], col=1, cex=1.3)
  points(x=1:3+0.2, y=y_hat_post[,n], pch=2, col=1, cex=1.3)
  arrows(x0=1:3, x1=1:3,  y0=y_hat_MLE_quantile[1,,n], y1=y_hat_MLE_quantile[2,,n], angle=90, length=0, code=3)
  arrows(x0=1:3+0.2, x1=1:3+0.2,  y0=y_hat_post_quantile[1,,n], y1=y_hat_post_quantile[2,,n], angle=90, length=0, code=3)
  points(x=1:3+0.1, y=eta.test[,n], pch="*", col=2, cex=2)

  points(x=1-0.2, y=glmnet.pred1[n], pch=0, col=1, cex=1.3)
  points(x=2-0.2, y=glmnet.pred2[n], pch=0, col=1, cex=1.3)
  points(x=3-0.2, y=glmnet.pred3[n], pch=0, col=1, cex=1.3)
  
  arrows(x0=1-0.2, x1=1-0.2, y0=y_hat_glmnet_lower1[n], y1=y_hat_glmnet_upper1[n], angle=90, length=0, code=3)
  arrows(x0=2-0.2, x1=2-0.2, y0=y_hat_glmnet_lower2[n], y1=y_hat_glmnet_upper2[n], angle=90, length=0, code=3)
  arrows(x0=3-0.2, x1=3-0.2, y0=y_hat_glmnet_lower3[n], y1=y_hat_glmnet_upper3[n], angle=90, length=0, code=3)
  


}
dev.off()
