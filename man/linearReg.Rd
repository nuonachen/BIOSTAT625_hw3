\name{linear regression}
\alias{linear regression}
\title{Linear Regression}
\usage{
fit_linear_regression(x, y)
get_fitted_values(x, y)
get_residual(x, y)
get_residual_std_error(x, y)
get_coef_stdErr(x, y)
get_tValue(x, y)
get_mult_Rsquared(x, y)
get_adjusted_Rsquared(x, y)
get_f_stat(x, y)
summary_table(x, y)
}
\arguments{
\item(x){independent variables/features as a matrix, n*m}
\item(y){dependent variables/outcomes as a matrix, n*1}
}
\value{
\item(fit_linear_regression){return the weights/coefficients of the linear regression model}
\item(get_fitted_values){return the fitted values of the linear regression model}
\item(get_residual){return the residuals of the linear regression model}
\item(get_residual_std_error){return the residual standard error of the linear regression model}
\item(get_coef_stdErr){return the standard errors of the coefficients}
\item(get_tValue){return the t scores of the coefficients}
\item(get_mult_Rsquared){return the multiple R squared of the linear regression model}
\item(get_adjusted_Rsquared){return the adjusted R squared of the linear regression model}
\item(get_f_stat){return the F statistic of the linear regression model}
\item(summary_table){print all of above computed values}

}
\description{
Fit a linear regression model with the given features and outcomes data and compute various attributes of the model.
}

\examples{
library(tidyverse)
data(cars)
x = as.matrix(cars["speed"])
y = as.matrix(cars["dist"])
weights = fit_linear_regression(x, y)
fitted_values = get_fitted_values(x, y)
residuals = get_residual(x, y)
residuals_std_error = get_residual_std_error(x, y)
coef_stdErr = get_coef_stdErr(x, y)
tValue = get_tValue(x, y)
mult_Rsquared = get_mult_Rsquared(x, y)
adjusted_Rsquared = get_adjusted_Rsquared(x, y)
f_stat = get_f_stat(x, y)
summary_table(x, y)
}
