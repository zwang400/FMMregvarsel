\name{r_sq_post}
\alias{r_sq_post}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate R-square Values for all MCMC Samples
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to calculate the R-square values for all posterior samples returned by MCMC simulations.
}
\usage{
r_sq_post(sim_res, data, scl=10001:1e5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim_res}{
%%     ~~Describe \code{x} here~~
list, MCMC simulation results returned by simulation_func, simulation_split or simulation_func_rpms

}
\item{data}{dataframe, consists of response y and covariates X}
\item{scl}{numeric vector, the index of remaining samples after burn-in}
}
\details{
Given the posterior samples of K, and c, for each iteration, within each cluster, fit a linear regression, obtain the R^2 value, then for each iteration, we obtain k_post R^2 values for k_post different clusters. Do such calculations for all posterior samples.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The R-square values for all clusters in all posterior samples are returned as a numeric vector.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Zhen Wang \email{zwangiowa@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##generate simulation data
simulation_data <- data_gen_func()

##FBMM without VS, hyper-prior for beta_bel
simulation_1 <- simulation_func(simulation_data$X, simulation_data$y, prior="Bessel", SS=F)

##calculate the R^2 values
r_sq_samples <- r_sq_post(simulation_1, simulation_data)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
