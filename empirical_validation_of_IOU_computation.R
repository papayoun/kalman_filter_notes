rm(list = ls())
library(expm)
# Paramètres ---------------------------------------------------------------

set.seed(123)

sigma <- rWishart(1, 2, diag(1, 2))[,, 1]
diag_rho <- runif(2, 0, 0.5)
anti_diagrho <- c(1, -1) * runif(1, 0, min(diag_rho))
rho <- matrix(anti_diagrho, nrow = 2, ncol = 2)
diag(rho) <- diag_rho
# rho <- diag(runif(2))
mu <- rnorm(2)
# Calculs utiles 
I_2 <- diag(1, 2)
I_kron_rho <- kronecker(I_2, rho)
rho_kron_I <- kronecker(rho, I_2)
rho_t_kron_I <- kronecker(t(rho), I_2)
rho_kron_rho_inv <- solve(I_kron_rho  + rho_kron_I)
sst_vec <- as.numeric(sigma %*% t(sigma))
tF <- 0.3
exp_rho_delta <- expm::expm(-rho * tF)
exp_rho_t_delta <- expm::expm(-t(rho) * tF)

# Mes  calculs ------------------------------------------------------------


source("integrated_ou.R")
modele <- IOU$new(rho, mu, sigma, 1)
mat <- modele$simulate(tF, 2)
my_mat <- purrr::map(modele$get_model_matrices(), 1)$hid_cov

# Covariance OU
# set.seed(123)
t1 <- 0.1
(cov_Pierre <- (expm::expm(-rho_kron_I * (tF - t1)) -
     expm::expm(-(rho_kron_I * tF + I_kron_rho * t1))) %*%
    rho_kron_rho_inv %*%
  sst_vec %>% 
  matrix(nrow = 2, ncol = 2))

# Calculs Paul ------------------------------------------------------------

S <- matrix(rho_kron_rho_inv %*% sst_vec, 
            nrow = 2, ncol = 2)
K_Y <- S - exp_rho_delta  %*% S %*% exp_rho_t_delta
(cov_Paul <- S %*% expm::expm(-t(rho)*(tF - t1)) - 
  expm::expm(-rho * t1) %*% S %*% expm::expm(-t(rho) * tF))
K_XY <- S %*% t(solve(rho)) %*% (I_2 - t(exp_rho_delta )) - 
  (I_2 - exp_rho_delta) %*% solve(rho) %*% S %*%exp_rho_t_delta

var_Paul <- function(t){
  rho_inv <- solve(rho)
  delta_M <- rho_inv %*% S + S %*% t(rho_inv)
  e_rho_t <-  expm::expm(-rho * t)
  term1 <- delta_M * t
  term2 <- (I_2 - e_rho_t) %*% rho_inv %*% delta_M
  term3 <- delta_M %*% t(rho_inv) %*% (I_2 - t(e_rho_t))
  term4 <- rho_inv %*% S %*% t(rho_inv)
  term5 <- e_rho_t %*% term4 %*% t(e_rho_t)
  
}


my_simu <- function(delta = 1e-5){
  out <- rep(NA, 4)
  temps <- 0
  x0 <- c(10, 10)
  while(temps < tF){
    xt <- x0 -rho %*% (x0 - mu) * delta + sigma %*% rnorm(2, 0, sqrt(delta))
    temps <- temps + delta
    if(all.equal(temps, t1) == T){
      out[1:2] <- as.numeric(xt)
    }
    if(all.equal(temps, tF) == T){
      out[3:4] <- as.numeric(xt)
    }
    x0 = xt
  }
  return(out)
}

my_simu_exact <- function(){
  out <- rep(NA, 4)
  temps <- 0
  x0 <- c(10, 10)
  exp_rho_1 <- expm::expm(-rho * t1)
  x1 <- mixtools::rmvnorm(1, as.numeric(mu + exp_rho_1 %*% (x0 - mu)),
                          S - exp_rho_1 %*% S %*% t(exp_rho_1)) %>% 
    as.numeric()
  exp_rho_2 <- expm::expm(-rho * (tF - t1))
  x2 <- mixtools::rmvnorm(1, as.numeric(mu + exp_rho_2 %*% (x1 - mu)),
                          S - exp_rho_2 %*% S %*% t(exp_rho_2)) %>% 
    as.numeric()
  out <- c(x1, x2)
  return(out)
}

res <- parallel::mclapply(1:0000, function(i) my_simu_exact(),
                mc.cores = 11)
mat_res <- 
  res %>% do.call(what = rbind)
cov(mat_res)
