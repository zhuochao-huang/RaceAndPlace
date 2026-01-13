## This script reproduces the top half of Figure 1
## Which looks at the interpretability simulations
## For violations of the unmeasured confounding 
## Assumption 

library(splines)
library(ggplot2)

rm(list=ls())

set.seed(2025)

# sample size and reference R values
n       <- 50000  
r0      <- 0.42 
r_pick  <- 0.55

# vector of values for coefficient sizes
a_vals  <- seq(0, 1, length.out = 10)

# function for the outcome model
outModel = function(combo_a, R, M, D, X, U) {
  return(plogis(-0.6 + combo_a*a*U + 
                  1*R + 0.5*X + 0.5*M + 0.5*D + 1*R*X))
}

# function for the mediator model
medModel = function(combo_a, R, D, X, U) {
  return(plogis(combo_a * a * U + X + 0.5*D + 1*R + 1*R*X))
}

# storing final output
mean_alpha <- numeric(length(a_vals))
var_alpha  <- numeric(length(a_vals))
bias_alpha  <- numeric(length(a_vals))

mean_beta <- numeric(length(a_vals))
var_beta  <- numeric(length(a_vals))
bias_beta  <- numeric(length(a_vals))

# Looping through all coefficient sizes
for (i in seq_along(a_vals)) {
  print(i)
  a <- a_vals[i]
  
  # generate covariates and random noise
  U <- rbinom(n, 1, p=0.5)
  X     <- rnorm(n)
  eps   <- rnorm(n)
  
  # examine all possible combinations of the signs of the coefficients
  combos <- expand.grid(a1 = c(1, -1), a2 = c(1, -1), 
                        a3 = c(1, -1), a4 = c(1, -1))
  
  # Store the results generated from each combination of 1 and -1
  mean_combo_alpha <- numeric(nrow(combos))
  var_combo_alpha  <- numeric(nrow(combos))
  bias_combo_alpha  <- numeric(nrow(combos))
  
  mean_combo_beta <- numeric(nrow(combos))
  var_combo_beta  <- numeric(nrow(combos))
  bias_combo_beta  <- numeric(nrow(combos))
  
  # Loop through combinations
  for (j in 1:nrow(combos)){
    ###########################################################
    ###################### Generate data ######################
    ###########################################################
    
    # Generate R and normalize
    R <- combos$a4[j]*a*U + X + eps 
    R <- ((R - min(R)) / (1.2*diff(range(R)))) + 0.1
    
    # Generate D
    pD <- plogis(combos$a1[j] * a * U + X)
    D  <- rbinom(n, 1, pD)
    
    # Generate M
    pM1 <- medModel(combo_a = combos$a2[j],
                    R=R, X=X, U=U, D=D)
    pM0 <- pM1
    
    M1  <- rbinom(n, 1, pM1)
    M0  <- rbinom(n, 1, pM0)
    
    M <- D*M1 + (1 - D)*M0
    
    # Generate Y
    Y_prob <- outModel(combo_a = combos$a3[j],
                       M=M, R=R, X=X, U=U, D=D)
    Y      <- rbinom(n, 1, Y_prob)
    
    ###########################################################
    ######################## Calculate W ######################
    ###########################################################
    
    ## First get a model for U | R, X, U, D=1, M=1
    wD1M1 = which(D == 1 & M == 1)
    dataD1M1 = data.frame(R=R[wD1M1], X=X[wD1M1], U=U[wD1M1])
    modD1M1 = glm(U ~ ns(R, 3) + X, data=dataD1M1, family=binomial)
    
    ## We will also need a model for U | R, X, D=1
    wD1 = which(D == 1)
    dataD1 = data.frame(R=R[wD1], X=X[wD1], U=U[wD1])
    modD1 = glm(U ~ ns(R, 3) + X, data=dataD1, family=binomial)
    
    ## Data to make predictions at
    new2_rp <- data.frame(R = r_pick, X = X, U = U)
    new2_r0 <- data.frame(R = r0,     X = X, U = U)
    
    ## now calculate Y portion of W
    numerator_w = outModel(combo_a = combos$a3[j],
                           U = 1, R=new2_rp$R, X=new2_rp$X,
                           M=1, D=1)* 
      predict(modD1M1, new2_rp, type="response") +
      outModel(combo_a = combos$a3[j],
               U = 0, R=new2_rp$R, X=new2_rp$X,
               M=1, D=1)*
      (1 - predict(modD1M1, new2_rp, type="response"))
    
    denominator_w = outModel(combo_a = combos$a3[j],
                             U = 1, R=new2_r0$R, X=new2_r0$X,
                             M=1, D=1)* 
      predict(modD1M1, new2_r0, type="response") +
      outModel(combo_a = combos$a3[j],
               U = 0, R=new2_r0$R, X=new2_r0$X,
               M=1, D=1)*
      (1 - predict(modD1M1, new2_r0, type="response"))
    
    w_y = numerator_w / denominator_w
    
    ## Now calculate M portion of W
    numerator_w = medModel(combo_a = combos$a2[j],
                           U = 1, R=new2_rp$R, X=new2_rp$X,
                           D=1)* 
      predict(modD1, new2_rp, type="response") +
      medModel(combo_a = combos$a2[j],
               U = 0, R=new2_rp$R, X=new2_rp$X,
               D=1)*
      (1 - predict(modD1, new2_rp, type="response"))
    
    denominator_w = medModel(combo_a = combos$a2[j],
                             U = 1, R=new2_r0$R, X=new2_r0$X,
                             D=1)* 
      predict(modD1, new2_r0, type="response") +
      medModel(combo_a = combos$a2[j],
               U = 0, R=new2_r0$R, X=new2_r0$X,
               D=1)*
      (1 - predict(modD1, new2_r0, type="response"))
    
    w_m = numerator_w / denominator_w
    
    w = w_m*w_y
    
    ###########################################################
    ############ Calculate Wstar for estimand 1 ###############
    ###########################################################
    
    w_star_y = outModel(combo_a = combos$a3[j],
                        U = new2_rp$U, R=new2_rp$R, X=new2_rp$X,
                        M=1, D=1) /
      outModel(combo_a = combos$a3[j],
               U = new2_r0$U, R=new2_r0$R, X=new2_r0$X,
               M=1, D=1)
    
    w_star_m = medModel(combo_a = combos$a3[j],
                        U = new2_rp$U, R=new2_rp$R, X=new2_rp$X,
                        D=1) /
      medModel(combo_a = combos$a3[j],
               U = new2_r0$U, R=new2_r0$R, X=new2_r0$X,
               D=1)
    
    wstarEstimand1 = w_star_m*w_star_y  
    
    ###########################################################
    ############ Calculate Wstar for estimand 2 ###############
    ###########################################################
    
    ## For estimand 2, we need a model for U | R=r, X=x, M(1)=1
    wM1 = which(M1 == 1)
    dataM1 = data.frame(R=R[wM1], X=X[wM1], U=U[wM1])
    modM1 = glm(U ~ ns(R, 3) + X, data=dataM1, family=binomial)
    
    ## Also need a model for U | R=r, X=x
    dataAll = data.frame(R=R, X=X, U=U)
    modAll = glm(U ~ ns(R, 3) + X, data=dataAll, family=binomial)
    
    
    # now average over the distribution of U in the numerator and denominator
    numerator = outModel(combo_a = combos$a3[j],
                         U = 1, R=new2_rp$R, X=new2_rp$X,
                         M=1, D=1)* 
      predict(modM1, new2_rp, type="response") +
      outModel(combo_a = combos$a3[j],
               U = 0, R=new2_rp$R, X=new2_rp$X,
               M=1, D=1)*
      (1 - predict(modM1, new2_rp, type="response"))
    
    denominator = outModel(combo_a = combos$a3[j],
                           U = 1, R=new2_r0$R, X=new2_r0$X,
                           M=1, D=1)* 
      predict(modM1, new2_r0, type="response") +
      outModel(combo_a = combos$a3[j],
               U = 0, R=new2_r0$R, X=new2_r0$X,
               M=1, D=1)*
      (1 - predict(modM1, new2_r0, type="response"))
    
    wstarEstimand2_y = numerator / denominator
    
    ## Now do M portion
    numerator = medModel(combo_a = combos$a2[j],
                         U = 1, R=new2_rp$R, X=new2_rp$X,
                         D=1)* 
      predict(modAll, new2_rp, type="response") +
      medModel(combo_a = combos$a2[j],
               U = 0, R=new2_rp$R, X=new2_rp$X,
               D=1)*
      (1 - predict(modAll, new2_rp, type="response"))
    
    denominator = medModel(combo_a = combos$a2[j],
                           U = 1, R=new2_r0$R, X=new2_r0$X,
                           D=1)* 
      predict(modAll, new2_r0, type="response") +
      medModel(combo_a = combos$a2[j],
               U = 0, R=new2_r0$R, X=new2_r0$X,
               D=1)*
      (1 - predict(modAll, new2_r0, type="response"))
    
    wstarEstimand2_m = numerator / denominator
    
    wstarEstimand2 = wstarEstimand2_m*wstarEstimand2_y
    
    ## Store values
    alpha_vals     <- wstarEstimand1 / w
    mean_combo_alpha[j]  <- mean(alpha_vals)
    var_combo_alpha[j]   <- var(alpha_vals)
    bias_combo_alpha[j] <- mean(w - wstarEstimand1) / mean(wstarEstimand1)
    
    beta_vals     <- wstarEstimand2 / w
    mean_combo_beta[j]  <- mean(beta_vals)
    var_combo_beta[j]   <- var(beta_vals)
    bias_combo_beta[j] <- mean(w - wstarEstimand2) / mean(wstarEstimand2)
  }
  
  ## Grab the largest ones across all combinations
  mean_alpha[i]  <- 
    abs(mean_combo_alpha[which.max(abs(mean_combo_alpha - 1))] - 1) + 1
  var_alpha[i]   <- max(var_combo_alpha)
  bias_alpha[i] <- max(abs(bias_combo_alpha))
  
  mean_beta[i]  <- 
    abs(mean_combo_beta[which.max(abs(mean_combo_beta - 1))] - 1) + 1
  var_beta[i]   <- max(var_combo_beta)
  bias_beta[i] <- max(abs(bias_combo_beta))
  
}

## Now make figures in GGplot
interpretationData = data.frame(type = rep(rep(c("Mean", "SD", "Percent bias"), 
                                               each=length(a_vals)), 2),
                                estimand = rep(c("Benchmark estimand", "True estimand"),
                                               each = length(a_vals)*3),
                                value = c(mean_alpha - 1, 
                                          sqrt(var_alpha),
                                          bias_alpha,
                                          mean_beta - 1, 
                                          sqrt(var_beta),
                                          bias_beta),
                                a = rep(a_vals, 6))
interpretationData$type = factor(interpretationData$type,
                                 levels = c("Mean", "SD", "Percent bias"))
interpretationData$estimand = factor(interpretationData$estimand,
                                     levels = c("Benchmark estimand", "True estimand"))

gg1 = ggplot(data=interpretationData, aes(x=a, y=value, color = estimand)) +
  geom_line(linewidth=1.4)+
  geom_point(size = 2.2) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),        # Size of legend labels
        legend.title = element_blank(),       # Size of legend title
        legend.key.size = unit(1.5, "lines"),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 10),  # x-axis tick labels
        axis.text.y = element_text(size = 10)) +
  xlab("Coefficient") +
  ylab("") +
  scale_color_discrete(
    labels = c(
      "Benchmark estimand" =
        expression("Benchmark estimand based on " ~ widetilde(xi)(r, x)),
      "True estimand" = expression("True estimand based on " ~ xi(r, x))
    )
  ) +
  facet_wrap(~ type, nrow=1)
gg1
