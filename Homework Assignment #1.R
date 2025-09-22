#First Homework Assignment
library(tidyverse)
library(R.matlab)
library(splines)
library(glmnet)
library(GGally)

#set seed
set.seed(123456)

#Problem 3
#we are asked to write a computer program to sample a function f with a 
#given analytical form at some given time points in [0,1].  
#part a

#define a function
f1 <- function(t){
  2*sin(2*pi*t) + cos(4*pi*t)
}

n <- 10
t_points <- seq(0, 1, length.out = n)
y_values <- sapply(t_points, f1)

#function to create a design matrix for a given Fourier Basis, assuming the 3
#Fourier basis functions provided, asuming K = odd >= 3
design_matrix <- function(t, K){
  t <- as.numeric(t)
  n <- length(t)
  
  #number of sine/cosine pairs (subtract the intercept and divide by 2)
  #if K = 3, then there is 1, K=5, then there is two, etc
  J <- (K - 1) %/% 2
  
  X <- matrix(0, nrow = n, ncol = K)
  
  #first Fourier basis function is 1
  X[,1] <- 1
  
  col <- 2
  
  #fills sin/cosine columns of Fourier design matrix one at a time
  #J = # of pairs to include
  for (k in 1:J) {
    X[, col]   <- sqrt(2)*sin(2*pi*k*t)
    col <- col + 1
    X[, col]   <- sqrt(2)*cos(2*pi*k*t) 
    col <- col + 1
  }
  
  return(X)
}

#penalty matrix assuming quadtratic penalty and K ?>= 3
fourier_penalty <- function(K, order = 2) {
  #number of sine/cosine pairs (subtract the intercept and divide by 2)
  #if K = 3, then there is 1, K=5, then there is two, etc
  J <- (K - 1) %/% 2
  
  P <- diag(0,K)
  
  #intercept is unpenalized
  #iterate of sin/cos pairs
  for (j in 1:J){
    #term for each "block"
    term <- (2*pi*j)^4
    
    #sin term
    #need the 2j-1 term or else the sin term goes in the cosine slot
    P[1+(2*j - 1), 1+(2*j - 1)] <- term
    
    #cosine term
    #don't need the -1 addeded here
    P[1+(2*j), 1+(2*j)] <- term
  }
  
  return(P)
}

#part b and c
#definding basis functions, creating design matrix and penalty matrix
K <- 5
lambda.grid <- seq(0, 0.0675, length.out = 10)
Phi <- design_matrix(t_points, K)
P <- fourier_penalty(K=K)

#the dense predictions grid
t_dense <- seq(0, 1, length.out = 5*n)
Phi_dense <- design_matrix(t_dense, K)
y_true_dense <- f1(t_dense)

#precompute
XtX <- crossprod(Phi)
Xty <- crossprod(Phi, y_values)

#estimate the function for the various values of lambda
fits <- lapply(lambda.grid, function(lam) {
  A <- XtX + (lam*P)
  c_hat <- solve(A, Xty)
  y_hat <- as.vector(Phi_dense %*% c_hat)
  data.frame(t = t_dense,
             predicted = y_hat,
             lambda = lam,
             lambda_lab = format(lam, scientific = TRUE, digits = 2),
             stringsAsFactors = FALSE)
})

#bind the rows and turns the lambda into labels via scientific notation out to 2 digits
df_pred <- bind_rows(fits) %>%
  mutate(lambda_lab = factor(format(lambda, scientific = TRUE, digits = 2),
                             levels = unique(format(lambda.grid, scientific = TRUE, digits = 2))))

df_true   <- data.frame(t = t_dense, true = y_true_dense)
df_points <- data.frame(t = t_points, y = y_values)


estimation.plot <- ggplot() +
  #draw the true function as a black line
  geom_line(data = df_true, aes(x = t, y = true), linewidth = 1.2, color = "black") +
  #add a seperate linetype, color for each lambda
  geom_line(data = df_pred,
            aes(x = t, y = predicted, color = lambda_lab, linetype = lambda_lab),
            linewidth = 0.9, alpha = 0.9) +
  #plot observed data points as big points, but remove aesthetics from the previous layers
  geom_point(data = df_points, aes(x = t, y = y), size = 2, alpha = 0.8, inherit.aes = FALSE) +
  #Name legend
  scale_color_discrete(name = expression(lambda)) +
  scale_linetype_discrete(name = expression(lambda)) +
  #Adapt title for various Basis functions used
  labs(title = sprintf("Penalized Fourier LS: True vs Predictions (K = %d)", K),
       x = "t", y = "f(t)") +
  #use the minimal theme with a larger text size
  theme_minimal(base_size = 13) +
  #make legend key wider.  Might help since we use different linetypes
  theme(legend.title = element_text(),
        legend.key.width = unit(1.6, "lines"))
estimation.plot

#Problem 4
#Write code to do FPCA
fpca_cust <- function(F, t, K = 3){
  
  #center the functional data
  mu <- colMeans(F)
  mu.fun.mat <- matrix(mu, nrow = nrow(F), ncol = ncol(F), byrow = TRUE)
  F.centered <- F - mu
  
  #assuming evenly spaced time points
  dt <- diff(as.vector(t))[1]
  
  #estimating covariance operator
  n <- nrow(F.centered)
  C_hat <- crossprod(F.centered)/n
  
  #eigen decomposition
  decomp <- eigen(C_hat, symmetric = TRUE)
  lambda <- decomp$values
  psi <- decomp$vectors
  
  #computing coefficients
  scores <- F.centered %*% (psi * dt)
  
  #getting into the plotting mechanisms now
  #plot, want plot of mean function
  df_mu <- data.frame(t = as.vector(t), mu = mu)
  mean_plot <- ggplot(df_mu, aes(x = t, y = mu)) +
    geom_line(linewidth = 1) +
    labs(title = "Mean Function", x = "t", y = "Mean f(t)") +
    theme_minimal()
  print(mean_plot)
  
  #getting the K dominant directions
  psiK <- psi[, seq_len(K), drop = FALSE]
  colnames(psiK) <- paste0("PC", seq_len(K))
  df_psi <- data.frame(t = as.vector(t), psiK)
  
  #pivotting longer for plotting
  psi_long <- tidyr::pivot_longer(df_psi, -t, names_to = "Component", values_to = "Value")
  eigen.vector_plot <- ggplot(psi_long, aes(x = t, y = Value, color = Component)) +
    geom_line(linewidth = 1) +
    labs(title = paste0("Leading Eigenfunctions (Top ", K, ")"),
         x = "t", y = "Eigenfunction Value") +
    theme_minimal()
  
  print(eigen.vector_plot)
  
  #scree plot
  df_lambda <- data.frame(Index = seq_along(lambda), Eigenvalue = lambda)
  
  scree_plot <- ggplot(df_lambda, aes(x = Index, y = Eigenvalue)) +
    geom_point(size = 2) +
    geom_line() +
    labs(title = "Eigenvalues (Scree Plot)", x = "Component", y = "Eigenvalue") +
    theme_minimal()
  
  print(scree_plot)
  
  #print the variance explained by the K PCs
  total_var <- sum(lambda)
  pve <- lambda / total_var 
  cum_pve <- cumsum(pve)
  
  cat("\nProportion of Variance Explained (PVE):\n")
  cat(paste(sprintf("  PC%-2d: %6.2f%%   |  Cumulative: %6.2f%%",
                    1:K, 100*pve[1:K], 100*cum_pve[1:K]),
            collapse = "\n"), "\n\n")
  
  #doing the correlation plot
  scores_df <- as.data.frame(scores[, seq_len(K), drop = FALSE])
  colnames(scores_df) <- paste0("Coefficient:", seq_len(K))
  
  top_coef <- ggpairs(scores_df, title = paste("Pairwise Scatter of Top", K, "FPCA Coefficients"))
  print(top_coef)
}


#read in data
#dataset 1
prob4.dat.1 <- readMat("Datasets/HW1/Problem 4/DataFile1_0_1.mat")

#implementing FPCA.  Need to center the functions
Fmatrix <- prob4.dat.1$f
ts <- prob4.dat.1$t

prob4.ques1 <- fpca_cust(F = Fmatrix, t = ts)


#now do this for the second dataset in the question
prob4.dat.2 <- readMat("Datasets/HW1/Problem 4/DataFile1_0_2.mat")

#implementing FPCA.  Need to center the functions
Fmatrix <- prob4.dat.2$f
ts <- prob4.dat.2$t

prob4.ques2 <- fpca_cust(F = Fmatrix, t = ts)


#problem 5
#Generate functional data as follows:  let f_i(t) = a_iPhi_1(t) + b_iPhi_(t)
#where phi_1(t) = sqrt(cos(wi*t)), phi_2(t) = sqrt(2)cos(3*pi*t)
#a_i, b_i ~ N(0,1)

prob5 <- function(n = 20, m = 50){
  #create time vector
  t <- seq(0,1,length.out=m)
  
  #phi functions
  phi1 <- sqrt(2)*cos(pi*t)
  phi2 <- sqrt(2)*cos(3*pi*t)
  
  #random coefficients
  a <- rnorm(n)
  b <- rnorm(n)
  
  #putting this together to create functional data
  F.dat <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:n) {
    F.dat[i, ] <- a[i] * phi1 + b[i] * phi2
  }
  
  #could probably create a list to return more stuff but eh
  out <- list(Fmatrix = F.dat, t = t)
  return(out)
}

#create functional data
F.dat <- prob5() 

#FPCA portion
prob5 <- fpca_cust(F = F.dat$Fmatrix, t = F.dat$t)



#Problem 6
prob6.dat <- readMat("Datasets/HW1/Problem 6/RegressionDataFile.mat")

#set up data
Fmat <- prob6.dat$f0
t <- prob6.dat$t
y <- prob6.dat$y0

#choose a basis.  Using B splines for this, bc we used Fourier for the previous.
#need to choose the number of basis functions
K <- 10
B <- bs(t, df=K, intercept = TRUE)

#need to discretize 
dt <- t[5] - t[4] #evenly spaced

#projecting functions onto basis space, with the dt added in
X <- Fmat %*% (B*dt)

#fit
fit <- lm(y ~ X - 1)
summary(fit)
beta_hat <- coef(fit)

beta_hat_fun <- as.vector(B %*% beta_hat)

#estimating residual variance
resid <- resid(fit)
df_res <- nrow(X) - ncol(X)
sigma2 <- sum(resid^2) / df_res

#covariance matrix of the coefficients
#sigma^2 * crossprod of X
XtX <- crossprod(X)             
Vtheta <- sigma2 * solve(XtX)
BV <- B %*% Vtheta
beta_var <- rowSums(BV*B)
se_beta <- sqrt(beta_var)

#critical values
crit <- qt(0.975, df = df_res)
beta_lo <- beta_hat_fun - crit * se_beta
beta_hi <- beta_hat_fun + crit * se_beta

df_plot <- data.frame(
  t = as.numeric(t),
  beta = beta_hat_fun,
  lo = beta_lo,
  hi = beta_hi
)

ggplot(df_plot, aes(x = t, y = beta)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "skyblue") +
  geom_line(linewidth = 1) +
  labs(
    title = expression(hat(beta)(t)),
    y = expression(hat(beta)(t)),
    x = "t"
  ) +
  theme_minimal(base_size = 14)

ridge_confidence_intervals <- function(lambda){
  #extract the coefficients
  theta_ridge <- as.numeric(coef(fit_glmnet, s = lambda))[-1]
  
  #get the beta curve
  beta_ridge <- as.vector(B %*% theta_ridge)
  
  #need to compute residuals manually
  yhat <- as.vector(X %*% theta_ridge)
  rss  <- sum((y - yhat)^2)
}

#comparing penalized regression (equivalent to ridge) to ridge
lambda_grid <- c(0.01, 0.1, 1, 10)
fit_glmnet <- glmnet(
  x = X, y = y,
  alpha = 0,                
  lambda = lambda_grid,     
  intercept = FALSE,
  standardize = FALSE
)

ridge_df <- do.call(rbind, lapply(lambda_grid, ridge_confidence_intervals))



#effective df for ridge
df_lambda <- function(lambda) sum(d / (d + lambda))
