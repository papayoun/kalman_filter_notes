rm(list = ls())
source("integrated_ou.R")
rho <- diag(1, 2); rho[1, 2] <- 0.7; rho[2, 1] <- -0.7
sigma <- diag(1, 2)
mu <- rnorm(2)

Delta <- runif(1, 0, 1)
modele <- IOU$new(rho, mu, sigma)
mat <- modele$compute_model_matrices()
eigen(mat$hidden_cov_matrix)
