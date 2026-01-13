## This script reproduces the bottom half of Figure 1
## Which looks at the interpretability simulations
## For violations of the race and place ratio assumption
## Note that the code takes some time to run, but can
## be sped up by reducing the sample size, though results
## will be more variable in that case

library(splines)
library(ggplot2)

set.seed(2025)

# sample size and reference R values
n       <- 50000
r0      <- -0.75 
r_pick  <- 0.75

# vector of values for coefficient sizes
a_vals  <- seq(0, 0.5, length.out = 10)

# storing final output
mean_alpha <- numeric(length(a_vals))
var_alpha  <- numeric(length(a_vals))
bias_alpha  <- numeric(length(a_vals))

# Looping through all coefficient sizes
for (i in seq_along(a_vals)) {
  print(i)
  a <- a_vals[i]
  
  # generate X, R, and D
  X     <- rnorm(n)
  R <- rnorm(n, mean=-0.5*X, sd=1)
  
  pD <- plogis(0.3*X - 0.3*R)
  D  <- rbinom(n, 1, pD)
  
  # examine all possible combinations of the signs of the coefficients
  combos <- expand.grid(a1 = c(1, -1), a2 = c(1, -1), 
                        a3 = c(1, -1))
  
  # Store the result generated from each combination of 1 and -1
  mean_combo_alpha <- numeric(nrow(combos))
  var_combo_alpha  <- numeric(nrow(combos))
  bias_combo_alpha  <- numeric(nrow(combos))
  
  # Loop through combinations
  for (j in 1:nrow(combos)){
    ###########################################################
    ###################### Generate data ######################
    ###########################################################
    
    # Generate M
    pM <- plogis(0.3*X - 0.3*R + 0.3*D + 
                   combos$a1[j]*a*D*X + combos$a2[j]*a*R*X)
    M  <- rbinom(n, 1, pM)
    
    # Generate Y
    Yprob <- plogis(-0.6 + 0.5*R - 0.5*X + 0.5*M + 0.5*D + 
                      combos$a3[j]*a*R*X)
    Y <- rbinom(n, 1, Yprob)
    
    ## outcome model probabilities
    Y_estimates = matrix(NA, n, 2)
    Y_estimates[,1] = plogis(-0.6 + 0.5*r0 - 0.5*X + 
                               0.5*1 + 0.5*1 + combos$a3[j]*a*r0*X)
    Y_estimates[,2] = plogis(-0.6 + 0.5*r_pick - 0.5*X + 
                               0.5*1 + 0.5*1 + combos$a3[j]*a*r_pick*X)
    
    ## treatment model
    wM1 = which(M == 1)
    dataDmod = data.frame(D=D, X=X, R=R)
    modD = glm(D ~ X + ns(R, 3), family=binomial, data=dataDmod)
    
    D_estimates = matrix(NA, n, 2)
    D_estimates[,1] = predict(modD, newdata = data.frame(X=X, R = r0), 
                              type="response")
    D_estimates[,2] = predict(modD, newdata = data.frame(X=X, R = r_pick),
                              type = "response")
    
    ## Density of R | M=1, X=x
    grid = seq(min(X)-0.001, max(X)+0.001, length=100)
    gridDens = matrix(NA, length(grid), 2)
    
    for (ni in 1 : length(grid)) {
      nn = order(abs(grid[ni] - X))[1:10000]
      wM1X = which(M == 1 & 1:n %in% nn)
      densM1X = density(R[wM1X])
      gridDens[ni,] = approx(densM1X$x, densM1X$y, xout = c(r0, r_pick))$y
    }
    
    P_M1X = matrix(NA, n, 2)
    
    for (ni in 1 : n) {
      P_M1X[ni,1] = approx(grid, gridDens[,1], xout = X[ni])$y
      P_M1X[ni,2] = approx(grid, gridDens[,2], xout = X[ni])$y
    }
    
    OtherComponents = (Y_estimates[,2] / Y_estimates[,1]) *
      (D_estimates[,2] / D_estimates[,1]) *
      (P_M1X[,2] / P_M1X[,1])
    
    
    ## Density of R | D=1
    wD1 = which(D == 1)
    densD1 = density(R[wD1])
    
    P_D1= approx(densD1$x, densD1$y, xout = c(r0, r_pick))$y
    
    ## Density of R | D=1, M=1
    wD1M1 = which(D == 1 & M==1)
    densD1M1 = density(R[wD1M1])
    
    P_D1M1= approx(densD1M1$x, densD1M1$y, xout = c(r0, r_pick))$y
    
    ## Density of R | D=1, X=x
    grid = seq(min(X)-0.001, max(X)+0.001, length=100)
    gridDens = matrix(NA, length(grid), 2)
    
    for (ni in 1 : length(grid)) {
      nn = order(abs(grid[ni] - X))[1:10000]
      wD1X = which(D == 1 & 1:n %in% nn)
      densD1X = density(R[wD1X])
      gridDens[ni,] = approx(densD1X$x, densD1X$y, xout = c(r0, r_pick))$y
    }
    
    P_D1X = matrix(NA, n, 2)
    
    for (ni in 1 : n) {
      P_D1X[ni,1] = approx(grid, gridDens[,1], xout = X[ni])$y
      P_D1X[ni,2] = approx(grid, gridDens[,2], xout = X[ni])$y
    }
    
    ## Density of R | D=1, M=1, X=x
    grid = seq(min(X)-0.001, max(X)+0.001, length=100)
    gridDens = matrix(NA, length(grid), 2)
    
    for (ni in 1 : length(grid)) {
      nn = order(abs(grid[ni] - X))[1:10000]
      wD1M1X = which(D == 1 & M == 1 & 1:n %in% nn)
      densD1M1X = density(R[wD1M1X])
      gridDens[ni,] = approx(densD1M1X$x, densD1M1X$y, xout = c(r0, r_pick))$y
    }
    
    P_D1M1X = matrix(NA, n, 2)
    
    for (ni in 1 : n) {
      P_D1M1X[ni,1] = approx(grid, gridDens[,1], xout = X[ni])$y
      P_D1M1X[ni,2] = approx(grid, gridDens[,2], xout = X[ni])$y
    }
    
    # remove extreme or missing values
    wKeep = which(!is.na(P_D1M1X[,1]) &
                    !is.na(P_D1M1X[,2]) &
                    !is.na(P_D1X[,1]) & 
                    !is.na(P_D1X[,2]) &
                    P_D1X[,1] > 0.05 &
                    P_D1X[,2] > 0.05 &
                    P_D1M1X[,1] > 0.05 &
                    P_D1M1X[,2] > 0.05)
    
    # take ratio
    wStar = P_D1X[wKeep,1] / P_D1X[wKeep,2]
    
    w = (P_D1M1X[wKeep,1]*P_D1[1]/P_D1M1[1]) / 
      (P_D1M1X[wKeep,2]*P_D1[2]/P_D1M1[2])
    
    
    EstW = w*OtherComponents[wKeep]
    EstWstar = wStar*OtherComponents[wKeep]
    
    alpha_vals = wStar / w
    
    mean_combo_alpha[j]  <- mean(alpha_vals)
    var_combo_alpha[j]   <- var(alpha_vals)
    bias_combo_alpha[j] <- abs(mean(EstW - EstWstar) / mean(EstWstar))
  }
  
  mean_alpha[i]  <- abs(mean_combo_alpha[which.max(abs(mean_combo_alpha - 1))] - 1) + 1
  var_alpha[i]   <- max(var_combo_alpha)
  bias_alpha[i] <- max(abs(bias_combo_alpha))
  
}

# Make plot in GGplot
interpretationData = data.frame(type = rep(c("Mean", "SD", "Percent bias"), 
                                           each=length(a_vals)),
                                value = c(mean_alpha-1, 
                                          sqrt(var_alpha),
                                          bias_alpha),
                                a = rep(a_vals, 3))
interpretationData$type = factor(interpretationData$type,
                                 levels = c("Mean", "SD", "Percent bias"))

gg1 = ggplot(data=interpretationData, aes(x=a, y=value)) +
  geom_line(linewidth=1.4)+
  geom_point(size=2.2) +
  xlab("Coefficient") +
  ylab("") +
  theme(strip.text = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 10),  # x-axis tick labels
        axis.text.y = element_text(size = 10)) +
  facet_wrap(~ type, nrow=1)

gg1



