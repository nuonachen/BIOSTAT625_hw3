
fit_linear_regression <- function(x, y){
  temp = rep(1, length(y))
  x = cbind(temp, x)
  return(as.vector(solve(t(x) %*% x) %*% t(x) %*% y))
}

get_fitted_values <- function(x, y){
  temp = rep(1, length(y))
  x = cbind(temp, x)
  weights = solve(t(x) %*% x) %*% t(x) %*% y
  temp = sweep(x, MARGIN = 2, weights, "*")
  f_values = rep(0, length(y))
  for(i in 1:nrow(temp)){
    f_values[i] = sum(temp[i, ])
  }
  return(f_values)
}

get_residual <- function(x, y){
  f_values = get_fitted_values(x, y)
  r = y - f_values
  return(r)
}

get_residual_std_error <- function(x, y){
  k = ncol(x)
  residual = get_residual(x, y)
  SSE = sum(residual^2)
  rse = sqrt(SSE/(length(residual) - (1+k)))
  return(rse)
}

get_coef_stdErr <- function(x, y){
  temp = sum((get_fitted_values(x, y) - y)^2)
  temp2 = rep(1, length(y))
  x = cbind(temp2, x)
  results = rep(0, ncol(x))
  for(i in 1:ncol(x)){
    results[i] = sqrt(temp/sum((x[,i] - mean(x[,i]))^2)) *sqrt(1/(nrow(x)-2))
  }
  return(results)

}

get_tValue <- function(x, y){
  results = fit_linear_regression(x, y) / get_coef_stdErr(x, y)
  return(results)
}

get_mult_Rsquared <- function(x, y){
  ssyy = sum((y - mean(y)) ^2)
  sse = sum(get_residual(x, y) ^2)
  return((ssyy-sse)/ssyy)
}

get_adjusted_Rsquared <- function(x, y){
  k = ncol(x)
  sse = sum(get_residual(x, y)^2)
  ssyy = sum((y-mean(y)) ^2)
  result = 1-(sse/ssyy) *(nrow(x)-1)/(nrow(x)-(k+1))
  return(result)
}

get_f_stat <- function(x, y){
  n = nrow(x)
  k = ncol(x)
  sse = sum(get_residual(x, y)^2)
  ssyy = sum((y-mean(y)) ^2)
  result = ((ssyy-sse)/k) / (sse/(n-(k+1)))
  return((result))
}

summary_table <- function(x, y){
  data = matrix(c(fit_linear_regression(x, y), get_coef_stdErr(x, y), get_tValue(x, y)), nrow = ncol(x)+1, ncol = 3)
  rownames(data) = c("Intercept", colnames(x))
  colnames(data) = c("Estimate", "Std. Error", "t value")
  print(data)
  print(paste("Residual standard error: ", get_residual_std_error(x, y)))
  print(paste("Multiple R-squared: ", get_mult_Rsquared(x, y)))
  print(paste("Adjusted R-squared: ", get_adjusted_Rsquared(x, y)))
  print(paste("F Statistics: ", get_f_stat(x, y)))

}


