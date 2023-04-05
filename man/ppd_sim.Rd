\name{ppd_sim}
\alias{ppd_sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate the Posterior Predictive Density (ppd) for New Data Subjects
}
\description{
A function to calculate the posterior predictive density for a new subject with response y* and covariates X*.
}
\usage{
ppd_sim(sim.res, x.new, y.new, scl=10001:1e5, thin=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim.res}{
%%     ~~Describe \code{x} here~~
list, MCMC simulation results returned by simulation_func, simulation_split or simulation_func_rpms
}
\item{x.new}{numeric vector of length D, covariates for the new subject}
\item{y.new}{numeric, the new y value to estimate the ppd at}
\item{scl}{numeric vector, the index of remaining samples after burn-in}
\item{thin}{numeric, the thinning factor to speed up calculations, posterior samples for every thin-th iterations are used to calculate ppd}
}
\details{
See the Appendix of the paper for formulas of the ppd.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The ppd of x.new at y.new is returned.
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

x.new <- rnorm(6)
y.new <- 0.5
##calculate ppd
ppd.new.sample <- ppd_sim(simulation_1, x.new, y.new, scl=10001:1e5, thin=20)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
