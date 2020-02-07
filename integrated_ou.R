library(R6)
library(tidyverse)
IOU <- R6Class(classname = "IOU",
               # Public --------------------------------------------------------
               public = list(
                 # Attributes --------------------------------------------------
                 observations = NULL,
                 hidden_states = NULL,
                 times = NULL,
                 # Methods -----------------------------------------------------
                 initialize = function(rho_,
                                       mu_, # A vector of 
                                       sigma_,
                                       obs_var_,
                                       m0_ = NULL,
                                       var0_ = NULL) {
                   private$validate(rho_, mu_, sigma_, obs_var_)
                   private$rho <- rho_
                   private$mu <- mu_
                   private$sigma <- sigma_
                   private$obs_var <- obs_var_
                   private$dimension <- length(mu)
                   if(is.null(m0_)){
                     private$m0 <- rep(0, private$dimension)
                   }
                   if(is.null(var0_)){
                     private$var0 <- diag(1, private$dimension)
                   }
                   else{
                     private$check_and_set_m0_var0(m0_, var0_)
                   }
                   private$compute_sigma_square()
                   private$compute_kronecker_terms()
                 }, # End initialize method
                 print = function(){
                   cat("mu: \n")
                   base::print(private$mu)
                   cat("rho: \n")
                   base::print(private$rho)
                   cat("sigma: \n")
                   base::print(private$sigma)
                   cat("sigma_square: \n")
                   base::print(private$sigma_square)
                 }, # End print method
                 simulate = function(times_ = NULL, n = NULL){
                   private$check_and_set_times(times_, n)
                   n <- length(self$times) - 1
                   private$compute_model_matrices()
                   invisible(self)
                 }, # End simulate method
                 get_model_matrices = function(){
                   list(hid_inter = private$hidden_dyn_intercepts,
                        hid_dyn = private$hidden_dyn_matrices,
                        hid_cov = private$hidden_cov_matrices,
                        obs_dyn = private$obs_dyn_matrices,
                        obs_cov = private$obs_cov_matrices)
                 }
               ),
               # Private -------------------------------------------------------
               private = list(
                 # Attributes --------------------------------------------------
                 dimension = NULL,
                 rho = NULL,
                 rho_inv = NULL,
                 mu = NULL,
                 sigma = NULL,
                 sigma_square = NULL,
                 sigma_square_vec = NULL,
                 kron_sum_rho = NULL,
                 kron_sum_rho_inv = NULL,
                 regular_time_stamp = NULL,
                 hidden_dyn_intercepts = NULL,
                 hidden_dyn_matrices = NULL,
                 hidden_cov_matrices = NULL,
                 obs_dyn_matrices = NULL,
                 obs_cov_matrices = NULL,
                 # Methods -----------------------------------------------------
                 ## Methods to compute attributes
                 compute_sigma_square = function(){
                   private$sigma_square = private$sigma %*% t(private$sigma)
                   private$sigma_square_vec <- as.numeric(private$sigma_square)
                 },
                 compute_kronecker_terms = function(){
                   Id <- diag(1, private$dimension)
                   private$kron_sum_rho <- kronecker(private$rho, 
                                                     Id) + 
                     kronecker(Id, private$rho)
                   private$kron_sum_rho_inv <- solve(private$kron_sum_rho)
                 },
                 compute_model_matrices = function(n){
                   private$rho_inv = solve(private$rho)
                   I_dim = diag(1, private$dimension)
                   I_2dim = diag(1, 2 * private$dimension)
                   I_kron_rho <- kronecker(I_dim,
                                           private$rho)
                   I_kron_rho_inv <- solve(I_kron_rho)
                   I_kron_rho_inv_sq <- I_kron_rho_inv %*% I_kron_rho_inv
                   rho_kron_I <- kronecker(private$rho,
                                           I_dim)
                   rho_kron_I_inv <- solve(rho_kron_I)
                   rho_kron_I_inv_sq <- rho_kron_I_inv %*% rho_kron_I_inv
                   compute_all_hidden_terms <- function(delta_t){
                     e_rho_t <- expm::expm(-private$rho * delta_t)
                     e_2rho_t <-  expm::expm(-private$kron_sum_rho * delta_t)
                     e_rho_kron_I_t <- expm::expm(-rho_kron_I * delta_t)
                     e_I_kron_rho_t <- expm::expm(-I_kron_rho * delta_t)
                     # Dynamics terms
                     # First, the intercept
                     intercept_vel <- (I_dim - e_rho_t) %*% private$mu
                     intercept_pos <- (delta_t * I_dim +
                                         private$rho_inv * e_rho_t) %*% mu
                     hid_dyn_intercept <- c(intercept_vel,
                                            intercept_pos)
                     # Then, the matrix, that will multiply the vector of
                     # hidden states
                     dyn_matrix_vel <- cbind(e_rho_t,
                                             matrix(0,
                                                    nrow = private$dimension,
                                                    ncol = private$dimension))
                     dyn_matrix_pos <- cbind(- private$rho_inv %*% e_rho_t,
                                             I_dim)
                     hid_dyn_matrix <- rbind(dyn_matrix_vel,
                                             dyn_matrix_pos)
                     # Covariance matrix
                     var_vel <- (private$kron_sum_rho_inv %*% 
                                   (I_2dim - e_2rho_t) %*% 
                                   private$sigma_square_vec) %>% 
                       matrix(ncol = private$dimension, 
                              nrow = private$dimension)
                     var_pos <- (private$kron_sum_rho_inv %*% 
                       (
                         (I_kron_rho_inv + rho_kron_I_inv) * delta_t + # First term
                           rho_kron_I_inv_sq %*% (e_rho_kron_I_t - I_2dim) + # Second
                           I_kron_rho_inv_sq %*% (e_I_kron_rho_t - I_2dim) + # Third
                           I_kron_rho_inv %*% rho_kron_I_inv %*%  # Fourth
                           (e_rho_kron_I_t + e_I_kron_rho_t - # Still Fourth
                              e_2rho_t - I_2dim) # Still fourth
                       ) %*% 
                       private$sigma_square_vec) %>% 
                       matrix(ncol = private$dimension,
                              nrow = private$dimension)
                     cov_pos_vel <- (private$kron_sum_rho_inv %*% 
                                       (
                                         rho_kron_I_inv + 
                                           I_kron_rho_inv %*% e_2rho_t -
                                           (rho_kron_I_inv + 
                                              I_kron_rho_inv) %*% e_rho_kron_I_t
                                       ) %*% 
                                       private$sigma_square_vec) %>% 
                       matrix(ncol = private$dimension,
                              nrow = private$dimension)
                     hid_cov_matrix <- rbind(
                       cbind(var_vel, t(cov_pos_vel)),
                       cbind(cov_pos_vel, var_pos)
                     )
                     list(cst  = hid_dyn_intercept,
                          dyn = hid_dyn_matrix,
                          cov = hid_cov_matrix
                     )
                   }
                   time_lags <- diff(self$times)
                   if(regular_time_stamp){
                     dt <- time_lags[1]
                     terms <- compute_all_hidden_terms(dt)
                     private$hidden_dyn_intercepts <- replicate(n, 
                                                                terms[["cst"]],
                                                                simplify = F)
                     private$hidden_dyn_matrices <- replicate(n, 
                                                      terms[["dyn"]],
                                                      simplify = F)
                     private$hidden_cov_matrices <- replicate(n, 
                                                      terms[["cov"]],
                                                      simplify = F)
                   }
                   else{
                     private$hidden_dyn_intercepts <- rep(list(rep(NA, 
                                                                   private$dimension),
                                                               n))
                     private$hidden_dyn_intercepts <- 
                       private$hidden_dyn_intercepts <- 
                        rep(list(matrix(NA, private$dimension, private$dimension)), n)
                     for(i in 1:n){
                       dt <- time_lags[i]
                       terms <- compute_all_hidden_terms(dt)
                       private$hidden_dyn_intercepts[[i]] <- terms[["cst"]]
                       private$hidden_dyn_matrices[[i]]  <- terms[["dyn"]]
                       private$hidden_cov_matrices[[i]] <-  terms[["cov"]]
                     }
                   }
                   obs_dyn_matrix <- cbind(diag(0, 2), diag(1, 2))
                   private$obs_dyn_matrices <- replicate(n,
                                                         obs_dyn_matrix,
                                                         simplify = F)
                   private$obs_cov_matrices <- replicate(n,
                                                         private$obs_var,
                                                         simplify = F)
                   
                 },
                 ## Check methods
                 vector2matrix = function(my_vec, my_dim){
                   if(!(length(my_vec) == 1 | length(my_vec) == my_dim)){
                     stop("Length of rho and sigma should be one, or the
                          length of mu")
                   }
                   if(length(my_vec) == 1){
                     my_vec <- diag(my_vec, my_dim)
                   }
                   else{
                     my_vec <- diag(my_vec)
                   }
                   return(my_vec)
                 }, # End vector2matrix method
                 check_dim = function(dim_mu, rho_, sigma_, obs_var_){
                   if(any(dim(rho_) != dim_mu)){
                     stop(paste0("rho must be a ",
                                 dim_mu, " by ", dim_mu,
                                 "square matrix"))
                   }
                   if(any(dim(sigma_) != dim_mu)){
                     stop(paste0("sigma must be a ",
                                 dim_mu, " by ", dim_mu,
                                 " square matrix"))
                   }
                   if(any(dim(obs_var_) != dim_mu)){
                     stop(paste0("obs_var must be a ",
                                 dim_mu, " by ", dim_mu,
                                 " square matrix"))
                   }
                 }, # End checkdim method
                 check_and_set_m0_var0 = function(m0_, var0_){
                   # First m0
                   if(!is.vector(m0_)){
                     stop("m0 should be a vector")
                   }
                   else if(length(m0_) != length(mu_)){
                     stop("m0 should have the same length has mu")
                   }
                   else{
                     private$m0 <- m0_
                   }
                   # Then var0
                   if(is.vector(var0_)){
                     if(length(var0_) == private$dimension){
                       private$var0 <- diag(var0_)
                     }
                     else if(length(var0_) == 1){
                       private$var0 <- diag(var0_, private$dimension)
                     }
                     else{
                       stop("If var0 is a vector, it should be of length
                            1 or the same length as mu")
                     }
                   }
                   else if(is.matrix(var0_)){
                     if(any(dim(var0_) != private$dimension)){
                       stop("Matrix var0 must be a length(mu)*length(mu) 
                            square matrix")
                     }
                     else if(!(isSymmetric(var0_) & 
                               all(eigen(var0_)$values > 0))){
                       stop("Matrix var0 must be symmetric and positive definite")
                     }
                     else{
                       private$var0 <- var0_
                     }
                   }
                   else{
                     stop("var0 must be either a vector (of length 1 of length(mu)),
                          or a length(mu)*length(mu) square matrix ")
                   }
                 },
                 check_and_set_times = function(times_, n){
                   if(is.null(times_)){
                     if(is.null(self$times)){
                       stop("Simulation times must be specified")
                     }
                     else{
                       times_ <- self$times
                       warning("Simulation times were set to previous 
                               model times attribute")
                     }
                   }
                   # Checking whether it is a vector
                   if(!is.vector(times_)){
                     stop("times must be a vector")
                   }
                   # Creating regular vector if times_ is of length 1
                   if(length(times_) == 1){
                     if(is.null(n)){
                       stop("When times has length one, n must be specified")
                     }
                     else if(n <= 1){
                       stop("n must be at least equal to 2")
                     }
                     else{
                       times_ <- seq(from = 0, by = times_, length.out = n + 1)
                     }
                   }
                   # Checking that times are increasing
                   if(any(diff(times_) <= 0)){
                     stop("Simulation times must be a vector of strictly
                          increasing values")
                   }
                   # Setting the value of the regular_time_stamp attribute
                   if(length(unique(diff(times_))) > 1){
                     private$regular_time_stamp <- FALSE
                   }
                   else{
                     private$regular_time_stamp <- TRUE
                   }
                   self$times <- times_
                 }, # End check and set times method
                 validate = function(rho_,
                                     mu_,
                                     sigma_,
                                     obs_var_){
                   stopifnot(is.numeric(rho_),
                             is.numeric(mu_),
                             is.numeric(sigma_),
                             is.numeric(obs_var_))
                   if(!is.vector(mu_)){
                     stop("mu must be a vector")
                   }
                   dim_mu <- length(mu_)
                   if(is.vector(rho_)){
                     rho_ <- private$vector2matrix(rho_, dim_mu)
                     warning(paste("rho has been transformed to a ",
                                   dim_mu, " by ", dim_mu,
                                   "matrix"))
                   }
                   if(is.vector(sigma_)){
                     sigma_ <- private$vector2matrix(sigma_, dim_mu)
                     warning(paste("sigma has been transformed to a ",
                                   dim_mu, " by ", dim_mu,
                                   "matrix"))
                   }
                   if(is.vector(obs_var_)){
                     obs_var_ <- private$vector2matrix(obs_var_, dim_mu)
                     warning(paste("obs_var_ has been transformed to a ",
                                   dim_mu, " by ", dim_mu,
                                   "matrix"))
                   }
                   private$check_dim(dim_mu, rho_, sigma_, obs_var_)
                 }
               ))
