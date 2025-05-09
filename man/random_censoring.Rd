% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/censor.R
\name{random_censoring}
\alias{random_censoring}
\title{apply random censoring in the given dataset}
\usage{
random_censoring(data_to_convert, event_rate, seed = 926, verbose = FALSE)
}
\arguments{
\item{data_to_convert}{the dataset to be converted; including the columns of \code{Id}, \code{Status} and \code{Time}}

\item{event_rate}{the desired event rate}

\item{seed}{seed for RNG, default is 926}

\item{verbose}{print out details}
}
\value{
the censored dataset
}
\description{
censoring the given dataset based on the given event rate
}
\examples{
list_data_to_convert <- compo_sim()
df_converted <- random_censoring(
    list_data_to_convert$dataset,
    0.8
)
str(df_converted)

}
