library(ggplot2)
library(qqplotr)
library(unitBSQuantReg)
library(unitquantreg)
library(VGAM)
library(numDeriv)


# Data --------------------------------------------------------------------


data = read.csv("RiskSurvey.csv")

x = data[, c("SIZELOG","INDCOST") ]
x_data = data[, c("SIZELOG","INDCOST") ]
x_data = as.matrix(x_data)

X = as.matrix(x)

X = cbind(1,X)

y = data$FIRMCOST/100



# Functions --------------------------------------------------------------------

g = function(x,a,b){
  return( 6 * dkumar(x,a,b) * pkumar(x,a,b)  * (1- pkumar(x,a,b) ) )
}

G = function(x,a,b){
  return( 3 * pkumar(x,a,b)^2 - 2 * pkumar(x,a,b)^3 )
}

inverse_logit <- function(x) {
  return(1 / (1 + exp(-x)))  
}

logito = function(betas){
  
  
  return( exp(X%*%c(betas[1],betas[2],betas[3]) )/(1+exp(X%*%c(betas[1],betas[2],betas[3]) )) )
  
}



fy_r = function(param){  # re parametrized density
  
  x = X
  
  q = param[1]
  alpha = param[2]
  a = param[3]
  beta0 = param[4]
  beta1 = param[5]
  beta2 = param[6]
  
  y = y
  
  
  
  
  c = ( (1-inverse_logit(q)^(0.5))/(1-inverse_logit(q)) )^(1/alpha)
  z = 0.5 - sin( asin(1-2*c)/3  )
  b = log( (1-z) )/log( 1-logito(param[4:10])^(a) ) 
  b = as.vector(b)
  
  t1 = - alpha * G(y,a,b)^(alpha - 1) * g(y,a,b) * (1-inverse_logit(q))
  t2 = log(inverse_logit(q)) * ( 1 - (1-inverse_logit(q))*G(y,a,b)^(alpha) )
  return(t1/t2)
}

Fy_r = function(param){     # re parametrized acum
  
  q = param[1]
  alpha = param[2]
  a = param[3]
  beta0 = param[4]
  beta1 = param[5]
  beta2 = param[6]
  
  
  c = ( (1-inverse_logit(q)^(0.5))/(1-inverse_logit(q)) )^(1/alpha)
  z = 0.5 - sin( asin(1-2*c)/3  )
  b = log( (1-z) )/log( 1-logito(param[4:10])^(a) ) 
  b = as.vector(b)
  
  
  t1 = log(1 - (1-inverse_logit(q))*G(y,a,b)^(alpha) )
  t2 = log(inverse_logit(q))
  return(t1/t2)
}




log_like_fy_r <- function(param) {
  -sum(log(fy_r(param)))
}




# Estimation -------------------------------------------------------

###---- CEMBKGLn ----###

parametros_reg =  optim(par = c(inverse_logit(0.6306047), 0.1776758,  0.5275131,  3.1746858, -0.6121712 , 1.3228549), fn = log_like_fy_r)

parameters_reg = c(parametros_reg$par[1],parametros_reg$par[2],parametros_reg$par[3],parametros_reg$par[4],parametros_reg$par[5],parametros_reg$par[6])



# log likelihood

-log_like_fy_r(parameters_reg) 


# bic

6*log(length(y))+2*log_like_fy_r(parameters_reg) 


# aic

6*2+2*log_like_fy_r(c( parametros_reg$par[1],parametros_reg$par[2],parametros_reg$par[3],parametros_reg$par[4],parametros_reg$par[5],parametros_reg$par[6])) 


# check hessian

solve(hessian( log_like_fy_r, parameters_reg))
      
# sd's

sqrt( diag(solve(hessian( log_like_fy_r, parameters_reg ) )))

# p - values betas

hessian_sd = sqrt( diag(solve(hessian( log_like_fy_r, parameters_reg ) )))

sd_betas = hessian_sd[4:6]

2*(1- pnorm( abs(parametros_reg$par[4:6])/sd_betas ) )

format(2*(1- pnorm( abs(parametros_reg$par[4:6])/sd_betas ) ),scientific = FALSE)

# randomized quantile residuals 

r_is = Fy_r(parameters_reg)
r_is = qnorm(r_is)

chek_ris = data.frame(values = r_is)

ggplot(chek_ris, aes(sample = values)) +
  stat_qq_band(distribution = "norm", alpha = 0, color = "black", linetype = "solid") + 
  stat_qq_line(distribution = "norm", color = "red", linetype = "dashed", size = 0.5) + 
  stat_qq_point(distribution = "norm", shape = 1, size = 2) + 
  labs(title = "CEMBKGLn",
       y = "ri", x = "Normal Quantiles") +
  theme_minimal()


###---- UBS ----###


df_data = data.frame("y" = y, "x" = X[,2:3])

taus <- 0.5
fits <- lapply(taus, function(TAU) unitBSQuantReg(y ~ x.SIZELOG + x.INDCOST, data = df_data, tau = TAU))
lapply(fits, coef)


fits[[1]]$loglik

4*log(73)-2*fits[[1]]$loglik
4*2-2*fits[[1]]$loglik

summary(fits[[1]])

r_is_UBS = data.frame(values = residuals( fits[[1]], type = c("quantile") ))

ggplot(r_is_UBS, aes(sample = values)) +
  stat_qq_band(distribution = "norm", alpha = 0, color = "black", linetype = "solid") + 
  stat_qq_line(distribution = "norm", color = "red", linetype = "dashed", size = 0.5) + 
  stat_qq_point(distribution = "norm", shape = 1, size = 2) + 
  labs(title = "UBS",
       y = "ri", x = "Normal Quantiles") +
  theme_minimal()



###---- KM ----###


fit_kum =  unitquantreg( y ~ x.SIZELOG + x.INDCOST, data = df_data, tau = 0.5, family = "kum")


fit_kum$loglik

4*log(73) - 2*fit_kum$loglik
2*4-2*fit_kum$loglik


summary(fit_kum)

r_is_KW = data.frame(values = residuals( fit_kum, type = c("quantile") ))

ggplot(r_is_KW, aes(sample = values)) +
  stat_qq_band(distribution = "norm", alpha = 0, color = "black", linetype = "solid") + 
  stat_qq_line(distribution = "norm", color = "red", linetype = "dashed", size = 0.5) + 
  stat_qq_point(distribution = "norm", shape = 1, size = 2) + 
  labs(title = "KM",
       y = "ri", x = "Normal Quantiles") +
  theme_minimal()



###---- Jhonson SB  ----###


fit_sb =  unitquantreg( y ~ x.SIZELOG + x.INDCOST, data = df_data, tau = 0.5, family = "johnsonsb")


fit_sb$loglik

4*log(73) - 2*fit_sb$loglik
2*4-2*fit_sb$loglik


summary(fit_sb)


r_is_SB = data.frame(values = residuals( fit_sb, type = c("quantile") ))

ggplot(r_is_SB, aes(sample = values)) +
  stat_qq_band(distribution = "norm", alpha = 0, color = "black", linetype = "solid") + 
  stat_qq_line(distribution = "norm", color = "red", linetype = "dashed", size = 0.5) + 
  stat_qq_point(distribution = "norm", shape = 1, size = 2) + 
  labs(title = "Johnson S_b",
       y = "ri", x = "Normal Quantiles") +
  theme_minimal()



# Vuong -------------------------------------------------------------------



values_CEMBKGLn = fy_r(parameters_reg) 
values_ubs = dubs(y, exp(X%*%(fits[[1]]$coefficients)[1:3])/(1 + exp(X%*%(fits[[1]]$coefficients)[1:3]) ) , (fits[[1]]$coefficients)[4],0.5  )
values_km = dkum(y, exp(X%*%fit_kum$coefficients$mu)/(1 + exp(X%*%fit_kum$coefficients$mu )  ) , (fit_kum$coefficients$theta),0.5  )
values_sb =djohnsonsb(y, exp(X%*%fit_sb$coefficients$mu )/(1 + exp(X%*%fit_sb$coefficients$mu )  ) , (fit_sb$coefficients$theta),0.5  )


sum_num1 = (1/sqrt( length(y)) ) * sum(  log(values_CEMBKGLn/values_ubs  )  )
w1 = mean(  log(values_CEMBKGLn/values_ubs  )^2 ) - ( mean( log(values_CEMBKGLn/values_ubs  )) )^2
t1 = sum_num1/sqrt(w1)

2*(1-pnorm(t1))



sum_num2 = (1/sqrt( length(y)) ) * sum(  log(values_CEMBKGLn/values_km  )  )
w2 = mean(  log(values_CEMBKGLn/values_km  )^2 ) - ( mean( log(values_CEMBKGLn/values_km  )) )^2
t2 = sum_num2/sqrt(w2)

2*(1-pnorm(t2))


sum_num3 = (1/sqrt( length(y)) ) * sum(  log(values_CEMBKGLn/values_sb  )  )
w3 = mean(  log(values_CEMBKGLn/values_sb  )^2 ) - ( mean( log(values_CEMBKGLn/values_sb  )) )^2
t3 = sum_num3/sqrt(w3)

2*(1-pnorm(t3))


### Example

num = (1/sqrt( length(y)) ) * sum(  log(values_km/values_sb )  )
den = (sum( ( log(values_km/values_sb  ) )^2   )/73) - (  (1/73) * sum(  log(values_km/values_sb  )  )  )^2
value_v = num/sqrt(den) 


# manual

c( value_v, 2*(1-pnorm(abs(value_v))) )

# using a function

vuong.test(fit_kum, fit_sb)
