## This code produces the simulation results when
## there are no unmeasured confounders and we use
## beta regression for estimation of densities with
## respect to R

library(mgcv)
library(nnet)
library(betareg)
library(ggplot2)

rm(list=ls())

# number of simulations
nSim = 100

# sample size
n=100000

# number of precincts
n_precinct=35

# R values to consider
r_seq = seq(0.2, 0.8, length=50)
r0 = 0.6

# x value of interest
x = c("3")

## expit function for generating probabilities
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

## Calculate beta density for beta regression parameters
betaDensity = function(r, mu, phi) {
  part1 = gamma(phi) / (gamma(mu*phi)*gamma((1 - mu)*phi))
  part2 = r^(mu*phi - 1)
  part3 = (1 - r)^((1 - mu)*phi - 1)
  
  return(part1*part2*part3)
}

## Function to simulate data set
simData = function(n, n_precinct) {
  D<-rep(0,n)
  M<-rep(0,n)
  Y00<-rep(0,n)
  Y10<-rep(0,n)
  Y01<-rep(0,n)
  Y11<-rep(0,n)
  p<-rep(0,n)
  R_value<-rep(0,n)
  Y<-rep(0,n)
  
  ## Generate covariates and precinct probabilities
  X <- sample(1:4, size = n, replace = TRUE, prob = rep(1,4))
  X_factor<-factor(X)
  p <- sample(1:n_precinct, size = n, replace = TRUE, prob = rep(0.1,n_precinct))
  
  # Generate D variable
  P_D <- expit(2*p/n_precinct-1+rnorm(n))
  D <- rbinom(n,1,P_D)
  
  # Generate M variable
  P_M1 <- expit((X-2.5+rnorm(n)))
  P_M0 <- expit((X-2.5+rnorm(n)))
  M1 <- rbinom(n,1,P_M1)
  M0 <- rbinom(n,1,P_M0)
  
  M = M1*(D == 1) + M0*(D == 0)
  
  # Generate R variable
  R<-rep(0,n_precinct)
  for (i in 1:n_precinct) {
    R[i] <- sum(D[p==i])/sum(p==i)
  }
  for (i in 1:n) {
    R_value[i] <- R[p[i]]
  }
  
  ## Simulate potential outcomes
  P_Y01<- expit(X-0.5 - 2*R_value)
  P_Y11<- expit(X-0.5 - 2*R_value)
  Y01 = rbinom(n, 1, P_Y01)
  Y11 = rbinom(n, 1, P_Y11)
  
  Y0 = Y00*(M0 == 0) + Y01*(M0 == 1)
  Y1 = Y10*(M1 == 0) + Y11*(M1 == 1)
  
  # observed outcomes
  Y = Y00*(D == 0 & M == 0) +
    Y01*(D == 0 & M == 1) +
    Y10*(D == 1 & M == 0) +
    Y11*(D == 1 & M == 1)
  
  data = list(Y=Y,
              Y1=Y1,
              Y0=Y0,
              R_value=R_value,
              M=M,
              X_factor=X_factor,
              D=D)
  
  return(data)
}

## Estimate psi(r,x) for a grid of r values and some pre-specified x values
estPsiFunction = function(r_seq, x, r0, Y, D, R_value, M, X_factor) {
  
  term1 = term2 = term3 = term4 = matrix(NA, length(x), length(r_seq))
  
  ## Get subset of data with M=1 and D=1 to fit outcome model to 
  dataYmod = data.frame(Y = Y[D == 1 & M == 1],
                        R_value = R_value[D == 1 & M == 1],
                        X_factor = X_factor[D == 1 & M == 1])
  Yx <- glm(Y ~ R_value + X_factor, family = binomial, data=dataYmod)
  
  ## data frame for predictions at r0
  newdataR0 <- data.frame(R_value=r0, X_factor=x)
  
  ## predictions at x and r0
  out1R0 <- as.numeric(predict(Yx, newdata = newdataR0, type = "response"))
  
  ## Now make predictions at all r values and x values
  for (ii in 1 : length(r_seq)) {
    newdata = data.frame(R_value=r_seq[ii],X_factor=x)
    out1 <- as.numeric(predict(Yx, newdata = newdata, type = "response"))
    
    term1[,ii] = out1/out1R0
  }
  
  ## Now get second term
  dataDmod = data.frame(D = D[M == 1],
                        R_value = R_value[M == 1],
                        X_factor = X_factor[M == 1])
  
  ## D model
  Dx<- glm(D~ R_value + X_factor,family = "binomial", data=dataDmod)
  
  ## Now make predictions at relevant x and r0
  out2R0 <- as.numeric(predict(Dx, newdata = newdataR0, type = "response"))
  
  ## Now make predictions at all r values and x values
  for (ii in 1 : length(r_seq)) {
    newdata = data.frame(R_value=r_seq[ii],X_factor=x)
    out2 <- as.numeric(predict(Dx, newdata = newdata, type = "response"))
    
    term2[,ii] = out2/out2R0
  }
  
  ## Now estimate the third term using beta regression
  dataBeta3mod = data.frame(R_value = R_value[M == 1],
                            X_factor = X_factor[M == 1])
  
  betaModel3 = betareg(R_value ~ X_factor, data=dataBeta3mod, link="logit")
  
  mu3 = predict(betaModel3, newdata = data.frame(X_factor = x))
  phi3 = betaModel3$coefficients$precision
  
  for (ii in 1 : length(r_seq)) {
    term3[,ii] = betaDensity(r = r_seq[ii], mu = mu3, phi = phi3) /
      betaDensity(r = r0, mu = mu3, phi = phi3)    
  }
  
  
  ## Now estimate the fourth term using beta regression
  dataBeta4mod = data.frame(R_value = R_value[D == 1],
                            X_factor = X_factor[D == 1])
  
  betaModel4 = betareg(R_value ~ X_factor, data=dataBeta4mod, link="logit")
  
  mu4 = predict(betaModel4, newdata = data.frame(X_factor = x))
  phi4 = betaModel4$coefficients$precision
  
  for (ii in 1 : length(r_seq)) {
    term4[,ii] = betaDensity(r = r0, mu = mu4, phi = phi4) /
      betaDensity(r = r_seq[ii], mu = mu4, phi = phi4)    
  }
  
  psiMat = term1*term2*term3*term4
  
  return(psiMat)
}

## Estimate true value of psi(r,x) for a grid of r values 
## and some pre-specified x values
truePsiFunction = function(r_seq, x, r0, Y1, R_value, X_factor) {
  trueValue = matrix(NA, length(x), length(r_seq))
  
  ## Fit model of Y1 against R,X to calculate psi(r,x)
  
  #Y_RX<- gam(Y1~ s(R_value) + X_factor,family = "binomial")
  Y_RX<- glm(Y1~ R_value + X_factor,family = "binomial")
  
  newdataR0 <- data.frame(R_value=r0,X_factor=x)
  predictR0 = predict(Y_RX, newdata = newdataR0, type = "response")
  
  for (ii in 1 : length(r_seq)) {
    newdata <- data.frame(R_value=r_seq[ii],X_factor=x)
    predictGrid = predict(Y_RX, newdata = newdata, type = "response")
    
    trueValue[,ii] = predictGrid / predictR0
  }
  
  return(trueValue)
}

## estimate CRR(x)
estCRRfunction = function(x, Y, D, M, X_factor) {
  
  ## Get subset of data with M=1 to fit outcome model to 
  dataYmod = data.frame(Y = Y[M == 1],
                        D = D[M == 1],
                        X_factor = X_factor[M == 1])
  Yx <- glm(Y ~ D + X_factor, family = binomial, data=dataYmod)
  
  ## predictions at x and D=1
  newdataD1 <- data.frame(D=1, X_factor=x)
  outD1 <- as.numeric(predict(Yx, newdata = newdataD1, type = "response"))
  
  ## predictions at x and D=0
  newdataD0 <- data.frame(D=0, X_factor=x)
  outD0 <- as.numeric(predict(Yx, newdata = newdataD0, type = "response"))
  
  term1 = outD1 / outD0
  
  ## Now get second term
  dataDmod = data.frame(D = D[M == 1],
                        X_factor = X_factor[M == 1])
  
  ## D model
  Dx<- glm(D~ X_factor,family = "binomial", data=dataDmod)
  
  ## Now make predictions at relevant x and r0
  out2D1 <- as.numeric(predict(Dx, newdata = newdataD1, type = "response"))
  out2D0 = 1 - out2D1
  
  term2 = out2D1 / out2D0
  
  ## Now estimate the third term assuming we have access to full data
  dataDmodFull = data.frame(D = D,
                            X_factor = X_factor)
  
  DxFull<- glm(D~ X_factor,family = "binomial", data=dataDmodFull)
  
  ## Now make predictions at relevant x and r0
  out2D1Full <- as.numeric(predict(DxFull, newdata = newdataD1, type = "response"))
  out2D0Full = 1 - out2D1Full
  
  term3 = out2D0Full / out2D1Full
  return(term1*term2*term3)
  
}

## Calculate true value of CRR
trueCRRfunction = function(x, Y1, Y0, X_factor) {
  Y1mod <- glm(Y1~ X_factor,family = "binomial")
  Y0mod <- glm(Y0~ X_factor,family = "binomial")
  
  newdata <- data.frame(X_factor=x)
  predictY0 = predict(Y0mod, newdata = newdata, type = "response")
  predictY1 = predict(Y1mod, newdata = newdata, type = "response")
  
  return(predictY1 / predictY0)
}

## Store true and estimated psi values
trueCurve = estCurve = matrix(NA, nSim, length(r_seq))

## Store true and estimated CRR values
trueCRR = estCRR = rep(NA, nSim)

for (ni in 1 : nSim) {
  cat(ni, "\r")
  
  ## Generate data
  data = simData(n = n, n_precinct = n_precinct)
  
  Y = data$Y
  Y1 = data$Y1
  Y0 = data$Y0
  M = data$M
  D = data$D
  X_factor = data$X_factor
  R_value = data$R_value
  
  # ## Estimate psi(r,x)
  estimates = estPsiFunction(r_seq = r_seq,
                             x = x,
                             r0 = r0,
                             Y = Y,
                             D = D,
                             R_value = R_value,
                             M = M,
                             X_factor = X_factor)
  
  ## True value of psi(r,x)
  truth = truePsiFunction(r_seq = r_seq,
                          x = x,
                          r0 = r0,
                          Y1 = Y1,
                          R_value = R_value,
                          X_factor = X_factor)
  
  ## Just store results for first x value for now
  estCurve[ni,] = estimates[1,]
  trueCurve[ni,] = truth[1,]
  
  
  ## Estimate CRR(x)
  estimatesCRR = estCRRfunction(x=x, 
                                Y=Y, 
                                D=D, 
                                M=M, 
                                X_factor=X_factor)
  
  ## True value of CRR
  truthCRR = trueCRRfunction(x=x, 
                             Y1=Y1, 
                             Y0=Y0, 
                             X_factor=X_factor)
  
  estCRR[ni] = estimatesCRR[1]
  trueCRR[ni] = truthCRR[1]
}

## Obtain average estimated curves and truth for psi(r,x=1)
avgTrueCurve = apply(trueCurve, 2, mean, na.rm=TRUE)
avgEstCurve = apply(estCurve, 2, mean, na.rm=TRUE)

## Look at true CRR and distribution of estimates
avgTrueCRR = mean(trueCRR, na.rm=TRUE)
hist(estCRR, col="darkgrey", main = "Estimated CRR values", breaks = 5)
abline(v = avgTrueCRR, lwd=3, col=2)


# Create a dataframe for the average true and estimated curves
df_avg <- data.frame(r_seq = r_seq,
                     avgTrueCurve = avgTrueCurve,
                     avgEstCurve = avgEstCurve)

df_est <- data.frame(r_seq = rep(r_seq, each = nrow(estCurve)),
                     estCurve = as.vector(estCurve),
                     group = factor(rep(1:nrow(estCurve), length(r_seq))))

colors <- c("Estimated Curve" = "darkgrey","Average Estimated Curve" = "blue", "True Curve" = "red")

p1<-ggplot()+
  geom_line(data = df_est,aes(x = r_seq, y = estCurve, group = group,color = "Estimated Curve"), 
            linewidth = 0.4, show.legend = FALSE) +
  geom_line(data = df_avg, aes(x = r_seq, y = avgEstCurve,color = "Average Estimated Curve"), 
            linewidth = 1.6,alpha = 0.55, show.legend = FALSE) +
  geom_line(data = df_avg, aes(x = r_seq, y = avgTrueCurve,color = "True Curve"), 
            linewidth = 1, linetype = "dashed", show.legend = FALSE) +
  ylim(range(estCurve, na.rm = TRUE)) +
  xlim(0.2,0.8) +
  labs(x = "R Value",color = "",linetype="",title = "No violation of Assumption")+
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

p1
