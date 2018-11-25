### =================================================================== ###
###                             APPLICATION                             ###  
### ------------------------------------------------------------------- ###
### Likelihood Based Inference for Censored Linear Regression Models    ###
###         with Scale Mixtures of Skew-Normal Distributions            ###
###  Authors: Thalita do Bem Mattos, Aldo M. Garay and Victor H. Lachos ###
### =================================================================== ###

rm(list=ls(all=TRUE))

source("auxiliary_functions_SMSNCR.R")
source("SAEM_estimation.R")

palette(gray(seq(0,.9,len =30)))

## CASE STUDY I 

library(SMNCensReg)
data(wage.rates)
y <- wage.rates$wage
X <- cbind(1,wage.rates$age,wage.rates$educ,wage.rates$kidslt6,wage.rates$kidsge6)
cc <- c(rep(0,428),rep(1,325))
n <- length(y)

## Skew-Normal

est_SN <- SAEM_EST(y,X,cc,cens="left",LS=NULL,precisao=0.0001,MaxIter=400, M=20, pc=0.40, dist="SN",nu.fixed=FALSE,nu=NULL,show.convergence="TRUE")

EnvelopeRMT(est_SN$teta,y,X,cc,cens="left",dist="SN")

## Skew-t 

est_ST <- SAEM_EST(y,X,cc,cens="left",LS=NULL,precisao=0.0001,MaxIter=400, M=20, pc=0.40, dist="ST",nu.fixed=FALSE,nu=3,show.convergence="TRUE")

EnvelopeRMT(est_ST$teta,y,X,cc,cens="left",dist="ST")

## Skew-Normal Contaminada

est_SCN <- SAEM_EST(y,X,cc,cens="left",LS=NULL,precisao=0.0001,MaxIter=400, M=20, pc=0.40, dist="SCN",nu.fixed=FALSE,nu=c(0.5,0.5),show.convergence="TRUE")

EnvelopeRMT(est_SCN$teta,y,X,cc,cens="left",dist="SCN")

## Skew-Slash

est_SSL <- SAEM_EST(y,X,cc,cens="left",LS=NULL,precisao=0.0001,MaxIter=400, M=20, pc=0.40, dist="SSL",nu.fixed=FALSE,nu=3,show.convergence="TRUE")

EnvelopeRMT(est_SSL$teta,y,X,cc,cens="left",dist="SSL")


## CASE STUDY II 

library(astrodatR)

data(censor_Be)
dados <- censor_Be

y <- dados[,5]
cc <- dados[,4]
x <- dados[,3]/1000
n <- nrow(dados)

for(i in 1:n)
{
  if(cc[i]==1)
  {
    cc[i] <- 0
  }else
  {
    cc[i] <- 1
  }
}

## ------------------

X <- cbind(1,x)

p <- ncol(X)
n <- nrow(X)

## Skew-Normal


SN_model2 <- SAEM_EST(y,X,cc,cens="left",LS=NULL,precisao=0.0001,MaxIter=300, M=20, pc=0.40, dist="SN",nu.fixed=FALSE,nu=NULL,show.convergence="TRUE")

EnvelopeRMT(SN_model2$teta,y,X,cc,cens="left",dist="SN")

## Skew-t 

ST_model2 <- SAEM_EST(y,X,cc,cens="left",LS=NULL,precisao=0.0001,MaxIter=300, M=20, pc=0.30, dist="ST",nu.fixed=TRUE,nu=3,show.convergence="TRUE")

EnvelopeRMT(ST_model2$teta,y,X,cc,cens="left",dist="ST")

## Skew-Normal Contaminada

SCN_model2 <- SAEM_EST(y,X,cc,cens="left",LS=NULL,precisao=0.0001,MaxIter=300, M=20, pc=0.30, dist="SCN",nu.fixed=TRUE,nu=c(0.5,0.1),show.convergence="TRUE")

EnvelopeRMT(SCN_model2$teta,y,X,cc,cens="left",dist="SCN")


## Skew-Slash

SSL_model2 <- SAEM_EST(y,X,cc,cens="left",LS=NULL,precisao=0.0001,MaxIter=300, M=20, pc=0.30, dist="SSL", nu.fixed=TRUE,nu=1.2, show.convergence="TRUE")

EnvelopeRMT(SSL_model2$teta,y,X,cc,cens="left",dist="SSL")

##### 




