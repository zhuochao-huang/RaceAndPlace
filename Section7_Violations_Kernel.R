## This code produces the simulation results when
## there are violations from unmeasured confounders and we use
## kernel smoothing for estimation of densities with
## respect to R

library(mgcv)
library(nnet)
library(betareg)
library(ggplot2)

rm(list=ls())

sigmoid<-function(x){return(1/(1+exp(-x)))}

# number of simulations
nSim = 100

# sample size
n=100000

# number of precincts
n_precinct=35

# R values to consider
r_seq = seq(0.2, 0.8, length=50)

# store estimated and true curves
trueCurve = estCurve = matrix(NA, nSim, length(r_seq))
trueCrr = estCrr = rep(0,nSim)

# Loop through all simulated data sets
for (ni in 1 : nSim) {
  cat(ni, "\r")
  
  ## SIMULATE DATA
  D<-rep(0,n)
  M<-rep(0,n)
  M0<-rep(0,n)
  M1<-rep(0,n)
  Y<-rep(0,n)
  p<-rep(0,n)
  R_value<-rep(0,n)
  Y_full<-rep(0,n)
  Y_D0M1<-rep(0,n)
  Y_D1M1<-rep(0,n)
  Y_1<-rep(0,n)
  Y_0<-rep(0,n)
  
  # Generate covariates and precinct probabilities
  X <- sample(1:4, size = n, replace = TRUE, prob = rep(1,4))
  p <- sample(1:n_precinct, size = n, replace = TRUE, prob = rep(0.1,n_precinct))
  u<-runif(n, min = 0, max = 1)
  
  # Generate D variable
  P_D <- sigmoid(u+2*p/n_precinct-1+rnorm(n))
  for (i in 1:n) {
    D[i] <- rbinom(1,1,P_D[i])
  }
  
  # Generate M variable
  P_M0 <- sigmoid((X-2.5+rnorm(n)))
  P_M1 <- sigmoid((X-2.5+rnorm(n)))
  
  for (i in 1:n) {
    M0[i] <- rbinom(1,1,P_M0[i])
    M1[i] <- rbinom(1,1,P_M1[i])
    if (D[i]==0){M[i]=M0[i]}
    else{M[i]=M1[i]}
  }
  
  # Generate R variable
  R<-rep(0,n_precinct)
  for (i in 1:n_precinct) {
    R[i] <- sum(D[p==i])/sum(p==i)
  }
  for (i in 1:n) {
    R_value[i] <- R[p[i]]
  }
  X_factor<-factor(X)
  
  # Generate potential and observed outcomes
  Y_D0M0<-rep(0,n)
  Y_D1M0<-rep(0,n)
  
  P_Y0<- u*sigmoid(X-0.5 - 2*R_value)
  P_Y1<- u*sigmoid(X-0.5 - 2*R_value)
  for (i in 1:n) {
    Y_D0M1[i]<-rbinom(1,1,P_Y0[i])
    Y_D1M1[i]<-rbinom(1,1,P_Y1[i])
    if (M[i]==0){Y[i]=0}
    else{
      if (D[i]==0){Y[i]=Y_D0M1[i]}
      if (D[i]==1){Y[i]=Y_D1M1[i]}
    }
    if (M0[i]==0){Y_0[i]=0}else{Y_0[i] <- Y_D0M1[i]}
    if (M1[i]==0){Y_1[i]=0}else{Y_1[i] <- Y_D1M1[i]}
  }
  
  ## CALCULATE TRUE CURVE

  # model for potential outcomes that is only available for simulated data
  Y_RX<- glm(Y_1~ R_value + X_factor,family = "binomial")
  
  # make predictions at relevant locations
  newdata <- data.frame(R_value=r_seq,X_factor="3")
  newdataR0 <- data.frame(R_value=0.6,X_factor="3")
  predictGrid = predict(Y_RX, newdata = newdata, type = "response")
  predictR0 = predict(Y_RX, newdata = newdataR0, type = "response")
  trueCurve[ni,] = as.numeric(predictGrid) / as.numeric(predictR0)
  
  ## True CRR 
  trueCrr[ni] <- mean(Y_1[X == 3]) / mean(Y_0[X == 3])
  
  ## Estimated CRR
  term1 <- mean(Y[D == 1 & M==1 & X == 3]) / mean(Y[D == 0 & M==1 & X == 3])
  term2 <- mean(D[X == 3 & M==1]) / (1 - mean(D[X == 3 & M==1]))
  
  #The authors who developed this estimand used secondary data sources 
  # to try and estimate this quantity.
  term3 <- mean(D[X == 3]) / (1 - mean(D[X == 3]))
  estCrr[ni] <- term1*term2/term3
  
  ## CALCULATE ESTIMATED CURVE
  r = r_seq
  R<-rep(0,n_precinct)
  for (i in 1:n_precinct) {
    R[i] <- sum(D[p==i])/sum(p==i)
  }
  P_R<- as.numeric(table(p))/n
  
  # kernel smoothing for marginal distribution of R
  denR_kde = density(R, weights=P_R, bw = 0.25)
  R_kde = function(u) approx(denR_kde$x, denR_kde$y, xout = u)$y
  
  # estimating conditional probability of R given M=1
  P_M1<-rep(0,n_precinct)
  for (i in 1:n_precinct) {
    P_M1[i] <- sum(p==i&M==1)/sum(M==1)
  }
  
  denR_M1_kde = density(R, weights=P_M1, bw = 0.25)
  R_M1_kde=function(u) approx(denR_M1_kde$x, denR_M1_kde$y, xout = u)$y
  
  X_factor <- factor(X)
  Y_M1D1<-Y[M==1&D==1]
  X_factor_M1D1<-X_factor[M==1&D==1]
  for (i in 1:n) {
    R_value[i] <- R[p[i]]
  }
  R_value_M1D1<-R_value[M==1&D==1]
  
  ## Fitting the outcome model 
  Y_RM1D1X<- glm(Y_M1D1~ R_value_M1D1+X_factor_M1D1,family = binomial)

  ## Making predictions at relevant locations
  newdata <-   data.frame(R_value_M1D1=r,X_factor_M1D1="3")
  newdataR0 <-   data.frame(R_value_M1D1=0.6,X_factor_M1D1="3")
  out1 <- as.numeric(predict(Y_RM1D1X, newdata = newdata, type = "response"))
  out1R0 <- as.numeric(predict(Y_RM1D1X, newdata = newdataR0, type = "response"))
  
  ## Fitting model for D given M=1
  D_M1<-D[M==1]
  X_factor_M1<-X_factor[M==1]
  R_value_M1<-R_value[M==1]
  D_RM1X<- glm(D_M1~ R_value_M1 + X_factor_M1,family = "binomial")
  newdata <- data.frame(R_value_M1=r,X_factor_M1="3")
  newdataR0 <- data.frame(R_value_M1=0.6,X_factor_M1="3")
  out2 <- as.numeric(predict(D_RM1X, newdata = newdata, type = "response"))
  out2R0 <- as.numeric(predict(D_RM1X, newdata = newdataR0, type = "response"))
  
  ## Estimating term 3 in identification formula
  X_M1<-X[M==1]
  X_MR <- multinom(X_M1 ~ R_value_M1)
  newdata <- data.frame(R_value_M1=r)
  newdataR0 <- data.frame(R_value_M1=0.6)
  P_X_MR<-as.numeric(predict(X_MR, newdata = newdata, "probs")[3])
  R_M1X<-lm(R_value_M1~X_factor_M1)
  P_X_M1<-sum(X==3&M==1)/sum(M==1)
  out3 <- as.numeric(P_X_MR*  R_M1_kde(r)/ P_X_M1)
  out3R0 <- as.numeric(P_X_MR*  R_M1_kde(0.6)/ P_X_M1)
  
  ## subsetting data to where M=1 and D=1
  X_M1D1<-X[M==1&D==1]
  R_value_M1D1<-R_value[M==1&D==1]
  
  P_MD<-rep(0,n_precinct)
  for (i in 1:n_precinct) {
    P_MD[i] <- sum(p==i&M==1&D==1)/sum(M==1&D==1)
  }
  
  ## kernel smoothed estimate of R given M=1 and D=1
  denR_MD_kde = density(R, weights=P_MD, bw = 0.25)
  R_MD_kde=function(u) approx(denR_MD_kde$x, denR_MD_kde$y, xout = u)$y
  
  
  ## Doing the same for this specific X value of interest in the simulation
  P_MDX<-rep(0,n_precinct)
  for (i in 1:n_precinct) {
    P_MDX[i] <- sum(D==1&p==i&M==1&X==3)/sum(D==1&M==1&X==3)
  }
  
  denR_MDX_kde = density(R, weights=P_MDX, bw = 0.25)
  R_MDX_kde=function(u) approx(denR_MDX_kde$x, denR_MDX_kde$y, xout = u)$y
  
  
  P_R_MDX=R_MDX_kde(r)
  P_R_MDX_R0=R_MDX_kde(0.6)
  
  P_D<-rep(0,n_precinct)
  for (i in 1:n_precinct) {
    P_D[i] <- sum(p==i&D==1)/sum(D==1)
  }

  ## Estimating term 4 in the identification formula
  P_D_R<- gam(D~ s(R_value),family = binomial)
  newdata <- data.frame(R_value=r)
  newdataR0 <- data.frame(R_value=0.6)
  P_R_D=R_kde(r)*  predict(P_D_R, newdata = newdata, type = "response")/(sum(D)/n)
  P_R_D_R0=R_kde(0.6)*  predict(P_D_R, newdata = newdataR0, type = "response")/(sum(D)/n)
  out4 <- as.numeric(P_R_MDX*P_R_D/R_MD_kde(r))
  out4R0 <- as.numeric(P_R_MDX_R0*P_R_D_R0/R_MD_kde(0.6))
  
  ## Overall estimate of curve
  estCurve[ni,] = as.numeric(out1*out2*out3/out4) / as.numeric(out1R0*out2R0*out3R0/out4R0)
  
}

avgTrueCurve = apply(trueCurve, 2, mean, na.rm=TRUE)
avgEstCurve = apply(estCurve, 2, mean, na.rm=TRUE)

# Create a dataframe for the average true and estimated curves
df_avg <- data.frame(r_seq = r_seq,
                     avgTrueCurve = avgTrueCurve,
                     avgEstCurve = avgEstCurve)

df_est <- data.frame(r_seq = rep(r_seq, each = nrow(estCurve)),
                     estCurve = as.vector(estCurve),
                     group = factor(rep(1:nrow(estCurve), length(r_seq))))

colors <- c("Estimated Curve" = "darkgrey","Average Estimated Curve" = "blue", "True Curve" = "red")

p2<-ggplot()+
  geom_line(data = df_est,aes(x = r_seq, y = estCurve, group = group,color = "Est Curve"), 
            linewidth = 0.4, show.legend = FALSE) +
  geom_line(data = df_avg, aes(x = r_seq, y = avgEstCurve,color = "Average Estimated Curve"), 
            linewidth = 1.6, show.legend = FALSE,alpha = 0.55) +
  geom_line(data = df_avg, aes(x = r_seq, y = avgTrueCurve,color = "True Curve"), 
            size = 1, linetype = "dashed", show.legend = FALSE) +
  ylim(range(estCurve, na.rm = TRUE)) +
  xlim(0.2,0.8) +
  labs(x = "r value",color = "",linetype="",title = "Existence of Unmeasured Confounders")+
  scale_color_manual(values = colors)+
  guides(colour = guide_legend(override.aes = list(linewidth=0.6)))+
  theme_set(theme_gray() + theme(legend.key=element_rect(fill = NA)))+
  theme(legend.key.size = unit(1, "cm"))+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 13, face = "bold"))+ 
  theme(plot.title = element_text(size = 20, face = "bold"))+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_text(size = 16, face = "bold"))+ 
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y  = element_text(size = 12))

p2

