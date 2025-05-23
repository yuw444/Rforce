% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpiu.R
\name{true_Y_numerical_form}
\alias{true_Y_numerical_form}
\title{Numerical form to calculate the true Y given the true hazard and other parameters}
\usage{
true_Y_numerical_form(
  t,
  constant_baseline_hazard = FALSE,
  baseline_hazard = 1,
  a_shape_weibull,
  sigma_scale_weibull,
  sigma_scale_gamma,
  lambdaZ,
  lambda
)
}
\arguments{
\item{t}{the time point to check}

\item{constant_baseline_hazard}{logical; whether constant baseline hazard is used}

\item{baseline_hazard}{a scalar, the constant baseline hazard}

\item{a_shape_weibull}{the parameter of weibull distribution when simulate the data}

\item{sigma_scale_weibull}{the parameter of weibull distribution when simulate the data}

\item{sigma_scale_gamma}{the parameter of gamma distribution when simulate the frality term}

\item{lambdaZ}{a vector, recurrent event hazard rate}

\item{lambda}{a scalar, the hazard rate of the stopping time}
}
\value{
a scalar, the mean number of event at time \code{t}
}
\description{
Numerical form to calculate the true Y given the true hazard and other parameters
}
