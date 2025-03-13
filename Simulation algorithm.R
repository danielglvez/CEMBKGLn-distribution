library(numDeriv)
library(parallel)
library(doParallel)
library(VGAM)



g <- function(x,a,b){   
  return( 6 * dkumar(x,a,b) * pkumar(x,a,b)  * (1- pkumar(x,a,b) ) )
}

G <- function(x,a,b){   
  return( 3 * pkumar(x,a,b)^2 - 2 * pkumar(x,a,b)^3 )
}

log_like <- function(parametros, data) {
  p <- parametros[1]
  alpha <- parametros[2]
  a <- parametros[3]
  b <- parametros[4]
  
  x =  data
  
  fx <-  -(alpha * G(x,a,b)^(alpha - 1) * g(x,a,b) * (1-p))/( log(p) * ( 1 - (1-p)*G(x,a,b)^(alpha) ) )
  
  llike <- log(fx)
  
  -sum(llike)
}

inverse_logit <- function(x) {
  return(1 / (1 + exp(-x)))  
}

log_like_r <- function(parametros, data) {
  p <- parametros[1]
  alpha <- parametros[2]
  a <- parametros[3]
  b <- parametros[4]
  
  x =  data
  
  fx <-  -(alpha * G(x,a,b)^(alpha - 1) * g(x,a,b) * (1-p))/( log(p) * ( 1 - (1-p)*G(x,a,b)^(alpha) ) )
  
  llike <- log(fx)
  
  -sum(llike)
}




quantile_f<- function(u,parametros){
  p = parametros[1]
  alpha =  parametros[2]
  a = parametros[3]
  b = parametros[4]
  
  
  c <- ( (1-p^u)/(1-p) )^(1/alpha)
  z <- 0.5 - sin( asin(1-2*c)/3  )
  
  valor <- ( 1 - (1-z)^(1/b)  )^(1/a)
  return(valor)
}






f_sim = function(n, param_vector ){
  
  alpha_mle <- c()
  p_mle <- c()
  a_mle <- c()
  b_mle <- c()
  
  i <- 1
  N <- 10000
  while( i <= N){
    
    while(T){  
      muestra = quantile_f( runif(n),  param_vector) 
      
      init_guess = optim( c(inverse_logit(runif(1,-1,1)),runif(1),runif(1),runif(1)), fn = log_like_r, data = muestra)
      init_guess = c( inverse_logit(init_guess$par[1]), init_guess$par[2], init_guess$par[3], init_guess$par[4] )
      
      mle_estimates <-  try( optim( init_guess, fn = log_like, data = muestra, method = "L-BFGS-B", 
                                    lower = c(0,0,0,0), upper = c(1,Inf,Inf,Inf), hessian = TRUE), silent <- T)
      if (class(mle_estimates) != "try-error") break
      
    }
    H <- solve(hessian(func = log_like, 
                       x = c(mle_estimates$par[1], mle_estimates$par[2], 
                             mle_estimates$par[3], mle_estimates$par[4]), 
                       data = muestra))
    
    if (all(diag(H) > 0) && all(!is.na(diag(H)))) {
      
      alpha_mle[i] <- mle_estimates$par[2]
      p_mle[i] <-  mle_estimates$par[1]
      a_mle[i] <- mle_estimates$par[3]
      b_mle[i] <-  mle_estimates$par[4]
      i <- i+1
    }
    cat("Iteration:", i, "\r")
  }
  info_alpha <- c( mean(alpha_mle), ( mean(alpha_mle) - param_vector[2] ), sqrt( mean( (alpha_mle  - param_vector[2] )^2 ) ) )
  info_p <- c( mean(p_mle), ( mean(p_mle) - param_vector[1] ), sqrt( mean( (p_mle - param_vector[1] )^2 ) ) )
  info_a <- c( mean(a_mle), ( mean(a_mle) - param_vector[3] ), sqrt( mean( (a_mle - param_vector[3] )^2 ) ) )
  info_b <- c( mean(b_mle), ( mean(b_mle) - param_vector[4] ), sqrt( mean( (b_mle - param_vector[4] )^2 ) ) )
  
  table_info <- t(cbind(info_p, info_alpha, info_a, info_b))
  colnames(table_info) <- c("Mean", "Bias", "RMSE")
  rownames(table_info) <- c("phi", "alpha", "a", "b")
  
  return(table_info)
}


