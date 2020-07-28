##################
# PD 2           #
##################

library(date)
library(stats)
library(imputeTS)
library(vars)
library(fields)
library(gplots)
library(rstan)
library(lme4)
library(reshape2)
library(circlize)
library(tseries)
library(vioplot)
library(RColorBrewer)
library(colorspace)

rm(list=ls())

setwd("C:/Columbia/projects/David/PD/PLOS")

ID <- read.csv("ID2.csv")
ID$age.cat <- as.numeric(ID$age <= 44)

w_mat1 <- NULL
#w_mat2 <- NULL
sig95_mat1 <- NULL
#sig95_mat2 <- NULL
sig90_mat1 <- NULL
#sig90_mat2 <- NULL
sdv_mat1 <- NULL
#sdv_mat2 <- NULL


load("PD2_Head_fit.Rdata")
vn_mat1 <- cbind(ID[ID$name %in% names(eta),], round(matrix(v_hat_2S_Gibbs, nrow=N, byrow=T), 3))
load("PD2_Join_fit.Rdata")
vn_mat2 <- cbind(ID[ID$name %in% names(eta),], round(matrix(v_hat_2S_Gibbs, nrow=N, byrow=T), 3))
load("PD2_Bowe_fit.Rdata")
vn_mat3 <- cbind(ID[ID$name %in% names(eta),], round(matrix(v_hat_2S_Gibbs, nrow=N, byrow=T), 3))
load("PD2_Musc_fit.Rdata")
vn_mat4 <- cbind(ID[ID$name %in% names(eta),], round(matrix(v_hat_2S_Gibbs, nrow=N, byrow=T), 3))


table(vn_mat1$age.cat)
table(vn_mat2$age.cat)
table(vn_mat3$age.cat)
table(vn_mat4$age.cat)



col <- rainbow_hcl(n=36)

tiff("fig8.tiff", height=10, width=10, res=300, units="in")
par(mfrow=c(4, 2), mar=c(3, 3, 4, 3))


boxplot(vn_mat1[vn_mat1$age.cat==1,6:41], ylim=range(vn_mat1[,6:41]), col=col, xaxt="n", main="")
title( main="Head (young)", line=3)
polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
abline(v=1:36, col="white", lty=3)
axis(side=1, at=2:6, labels=c("F", "S", "D", "A", "C"),  col="blue", cex.axis=0.6)
axis(side=1, at=seq(7, 31, 6), labels=c("F", "S", "D", "A", "C"),  col="purple", cex.axis=0.6)
axis(side=3, at=seq(1, 36, 7), labels=c("H", "F", "S", "D", "A", "C"),  col=1, cex.axis=0.6)

box()
boxplot(vn_mat1[vn_mat1$age.cat==1,6:41], ylim=range(vn_mat1[,6:41]), col=col, add=T, xaxt="n", main="")



boxplot(vn_mat1[vn_mat1$age.cat==0,6:41], ylim=range(vn_mat1[,6:41]), main="", col=col, xaxt="n")
title( main="Head (old)", line=3)
polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
abline(v=1:36, col="white", lty=3)
axis(side=1, at=2:6, labels=c("F", "S", "D", "A", "C"),  col="blue", cex.axis=0.6)
axis(side=1, at=seq(7, 31, 6), labels=c("F", "S", "D", "A", "C"),  col="purple", cex.axis=0.6)
axis(side=3, at=seq(1, 36, 7), labels=c("H", "F", "S", "D", "A", "C"),  col=1, cex.axis=0.6)

box()
boxplot(vn_mat1[vn_mat1$age.cat==0,6:41], ylim=range(vn_mat1[,6:41]), main="", col=col, add=T, xaxt="n")


boxplot(vn_mat2[vn_mat2$age.cat==1,6:41], ylim=range(vn_mat2[,6:41]), main="", col=col, xaxt="n")
title( main="Join (young)", line=3)
polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
abline(v=1:36, col="white", lty=3)
axis(side=1, at=2:6, labels=c("F", "S", "D", "A", "C"),  col="blue", cex.axis=0.6)
axis(side=1, at=seq(7, 31, 6), labels=c("F", "S", "D", "A", "C"),  col="purple", cex.axis=0.6)
axis(side=3, at=seq(1, 36, 7), labels=c("J", "F", "S", "D", "A", "C"),  col=1, cex.axis=0.6)

box()
boxplot(vn_mat2[vn_mat2$age.cat==1,6:41], ylim=range(vn_mat2[,6:41]), main="", col=col, add=T, xaxt="n")


boxplot(vn_mat2[vn_mat2$age.cat==0,6:41], ylim=range(vn_mat2[,6:41]), main="", col=col, xaxt="n")
title( main="Join (old)", line=3)
polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
abline(v=1:36, col="white", lty=3)
axis(side=1, at=2:6, labels=c("F", "S", "D", "A", "C"),  col="blue", cex.axis=0.6)
axis(side=1, at=seq(7, 31, 6), labels=c("F", "S", "D", "A", "C"),  col="purple", cex.axis=0.6)
axis(side=3, at=seq(1, 36, 7), labels=c("J", "F", "S", "D", "A", "C"),  col=1, cex.axis=0.6)

box()
boxplot(vn_mat2[vn_mat2$age.cat==0,6:41], ylim=range(vn_mat2[,6:41]), main="", col=col, add=T, xaxt="n")


boxplot(vn_mat3[vn_mat3$age.cat==1,6:41], ylim=range(vn_mat3[,6:41]), main="", col=col, xaxt="n")
title( main="Bowe (young)", line=3)
polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
abline(v=1:36, col="white", lty=3)
axis(side=1, at=2:6, labels=c("F", "S", "D", "A", "C"),  col="blue", cex.axis=0.6)
axis(side=1, at=seq(7, 31, 6), labels=c("F", "S", "D", "A", "C"),  col="purple", cex.axis=0.6)
axis(side=3, at=seq(1, 36, 7), labels=c("B", "F", "S", "D", "A", "C"),  col=1, cex.axis=0.6)

box()
boxplot(vn_mat3[vn_mat3$age.cat==1,6:41], ylim=range(vn_mat3[,6:41]), main="", col=col, add=T, xaxt="n")


boxplot(vn_mat3[vn_mat3$age.cat==0,6:41], ylim=range(vn_mat3[,6:41]), main="", col=col, xaxt="n")
title( main="Bowe (old)", line=3)
polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
abline(v=1:36, col="white", lty=3)
axis(side=1, at=2:6, labels=c("F", "S", "D", "A", "C"),  col="blue", cex.axis=0.6)
axis(side=1, at=seq(7, 31, 6), labels=c("F", "S", "D", "A", "C"),  col="purple", cex.axis=0.6)
axis(side=3, at=seq(1, 36, 7), labels=c("B", "F", "S", "D", "A", "C"),  col=1, cex.axis=0.6)

box()
boxplot(vn_mat3[vn_mat3$age.cat==0,6:41], ylim=range(vn_mat3[,6:41]), main="", col=col, add=T, xaxt="n")


boxplot(vn_mat4[vn_mat4$age.cat==1,6:41], ylim=range(vn_mat4[,6:41]), main="", col=col, xaxt="n")
title( main="Musc (young)", line=3)
polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
abline(v=1:36, col="white", lty=3)
axis(side=1, at=2:6, labels=c("F", "S", "D", "A", "C"),  col="blue", cex.axis=0.6)
axis(side=1, at=seq(7, 31, 6), labels=c("F", "S", "D", "A", "C"),  col="purple", cex.axis=0.6)
axis(side=3, at=seq(1, 36, 7), labels=c("M", "F", "S", "D", "A", "C"),  col=1, cex.axis=0.6)

box()
boxplot(vn_mat4[vn_mat4$age.cat==1,6:41], ylim=range(vn_mat4[,6:41]), main="", col=col, add=T, xaxt="n")

boxplot(vn_mat4[vn_mat4$age.cat==0,6:41], ylim=range(vn_mat4[,6:41]), main="", col=col, xaxt="n")
title( main="Musc (old)", line=3)
polygon(x=c(-100, 100, 100, -100), y=c(-100, -100, 100, 100), col="grey85")
abline(v=1:36, col="white", lty=3)
axis(side=1, at=2:6, labels=c("F", "S", "D", "A", "C"),  col="blue", cex.axis=0.6)
axis(side=1, at=seq(7, 31, 6), labels=c("F", "S", "D", "A", "C"),  col="purple", cex.axis=0.6)
axis(side=3, at=seq(1, 36, 7), labels=c("M", "F", "S", "D", "A", "C"),  col=1, cex.axis=0.6)

box()
boxplot(vn_mat4[vn_mat4$age.cat==0,6:41], ylim=range(vn_mat4[,6:41]), main="", col=col, add=T, xaxt="n")

dev.off()





w <- vector("list", length=4)
sig95 <- vector("list", length=4)
sig90 <- vector("list", length=4)
nam <- vector("list", length=4)
v <- vector("list", length=4)
sdv <- vector("list", length=4)
varnames <- vector("list", length=4)



load("pred6_Head_p1.RData")
w[[1]] <- matrix(w_hat_2S_Gibbs, B)
sig95[[1]] <- matrix(as.numeric(!((sw$mon.w[,4]<0) & ((sw$mon.w[,5]>0)))), B)
sig90[[1]] <- matrix(as.numeric(!((sw$mon.w[,6]<0) & ((sw$mon.w[,7]>0)))), B)
sdv[[1]] <- matrix(sdv_hat_2S_Gibbs, B)
#nam[[1]] <- names(dat.list)
#v[[1]] <- matrix(v_hat_2S_Gibbs, ncol=B, byrow=T)
varnames[[1]] <- vars



load("pred6_Join_p1.RData")
w[[2]] <- matrix(w_hat_2S_Gibbs, B)
sig95[[2]] <- matrix(as.numeric(!((sw$mon.w[,4]<0) & ((sw$mon.w[,5]>0)))), B)
sig90[[2]] <- matrix(as.numeric(!((sw$mon.w[,6]<0) & ((sw$mon.w[,7]>0)))), B)
sdv[[2]] <- matrix(sdv_hat_2S_Gibbs, B)
#nam[[2]] <- names(dat.list)
#v[[2]] <- matrix(v_hat_2S_Gibbs, ncol=B, byrow=T)
varnames[[2]] <- vars


load("pred6_Bowe_p1.RData")
w[[3]] <- matrix(w_hat_2S_Gibbs, B)
sig95[[3]] <- matrix(as.numeric(!((sw$mon.w[,4]<0) & ((sw$mon.w[,5]>0)))), B)
sig90[[3]] <- matrix(as.numeric(!((sw$mon.w[,6]<0) & ((sw$mon.w[,7]>0)))), B)
sdv[[3]] <- matrix(sdv_hat_2S_Gibbs, B)
#nam[[3]] <- names(dat.list)
#v[[3]] <- matrix(v_hat_2S_Gibbs, ncol=B, byrow=T)
varnames[[3]] <- vars


load("pred6_Musc_p1.RData")
w[[4]] <- matrix(w_hat_2S_Gibbs, B)
sig95[[4]] <- matrix(as.numeric(!((sw$mon.w[,4]<0) & ((sw$mon.w[,5]>0)))), B)
sig90[[4]] <- matrix(as.numeric(!((sw$mon.w[,6]<0) & ((sw$mon.w[,7]>0)))), B)
sdv[[4]] <- matrix(sdv_hat_2S_Gibbs, B)
#nam[[4]] <- names(dat.list)
#v[[4]] <- matrix(v_hat_2S_Gibbs, ncol=B, byrow=T)
varnames[[4]] <- vars



# w plot
tiff("Fig6.tiff", height=9, width=10, units="in", res=300)
par(mar=c(4,5,5,5), las=1, mfrow=c(2,2), oma=c(2, 2, 2, 3))

for(k in 1:4){
  image.plot(x=1:B, y=1:B, z=t(w[[k]][B:1,]), col=colorpanel(200, "blue", "white", "red"), breaks=sort(c(seq(-0.5, 0.5, length.out=200), 0)), main="", axes=F, xlab="", ylab="")
  mtext(side=3, font=2, text=paste("(", letters[k], ")", sep=""), line=3, cex=1.5)
  for(i in B:1){
    for(j in 1:B){
      text(x=j, y=i, labels=round(w[[k]][B+1-i, j], 3))
      if(sig90[[k]][B+1-i, j]==1){
        text(x=j, y=i, labels=round(w[[k]][B+1-i, j], 3), font=2)
      }
    }
  }
  axis(side=2, at=B:1, labels=varnames[[k]],  tick=FALSE)
  axis(side=3, at=1:B, labels=varnames[[k]],  tick=FALSE)
  
  polygon(x=c(1:B+0.5, B:1+0.5), y=c(rep(5.5, 6), rep(6.5, 6)), border="purple", lwd=3)
  polygon(x=c(rep(0.5, 6), rep(1.5, 6)), y=c(1:B-0.5, B:1-0.5), border="blue", lwd=3)
  
}
dev.off()






# sdv plot
tiff("Fig7.tiff", height=9, width=10, units="in", res=300)
par(mar=c(4,5,5,5), las=1, mfrow=c(2,2), oma=c(2, 2, 2, 3))

for(k in 1:4){
  image.plot(x=1:B, y=1:B, z=t(sdv[[k]][B:1,]), col=colorpanel(200, "white", "red"), breaks=sort(c(seq(0, 0.3, length.out=201))), main="", axes=F, xlab="", ylab="")
  mtext(side=3, font=2, text=paste("(", letters[k], ")", sep=""), line=3, cex=1.5)
  for(i in B:1){
    for(j in 1:B){
      text(x=j, y=i, labels=round(sdv[[k]][B+1-i, j], 3))
    }
  }
  axis(side=2, at=B:1, labels=varnames[[k]],  tick=FALSE)
  axis(side=3, at=1:B, labels=varnames[[k]],  tick=FALSE)
  
  polygon(x=c(1:B+0.5, B:1+0.5), y=c(rep(5.5, 6), rep(6.5, 6)), border="purple", lwd=3)
  polygon(x=c(rep(0.5, 6), rep(1.5, 6)), y=c(1:B-0.5, B:1-0.5), border="blue", lwd=3)
  
}
dev.off()






