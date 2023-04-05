\name{posterior_inf_split}
\alias{posterior_inf_split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Inference for Split Models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculate different metrics defined in the paper to evaluate performance of different split models in : clustering accuracy, parameter estimation accuracy and variable selection accuracy.
}
\usage{
post_inf_split(sim_res, scl = 10001:1e5, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim_res}{list, MCMC simulation results, return from the simulation functions defined in the package
%%     ~~Describe \code{x} here~~
}
\item{scl}{numeric vector, index of remaining samples after burn-in}
\item{data}{list, the data input used to run the simulations}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Runs inference for N-split models.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
\item{true_size}{numeric vector, the true cluster sizes in the data}
\item{measure_km_ari}{numeric, a measurement defined to measure how far away the clusters are from each other, calculated by the ARI values of the true clutering membership and the clustering result of the K-means method, a larger value generally means more distant clusters, thus easier clustering problem}
\item{auto_corr_k}{numeric vector, auto-correlations for posterior of K}
\item{geweke_k}{numeric, the Z-statistic for the geweke diagnostics for posterior of K}
\item{post_k}{numeric vector, posterior distribution of K}
\item{mse_a1, mse_a2, mse_a3}{numeric vectors, the quartiles of the mean squared errors of alpha for the 3 clusters}
\item{mse_b1, mse_b2, mse_b3}{numeric vectors, the quartiles of the mean squared errors of beta for the 3 clusters}
\item{mse_psi}{numeric vectors, the quartiles of the mean squared errors of psi}
\item{ARI}{numeric vector, the (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1) quantiles of ARI values calculated by all posterior samples of c}
\item{c_bin}{numeric vector, sizes of clusters obtained by the Binder loss}
\item{c_vi}{numeric vector, sizes of clusters obtained by the VI loss}
\item{ari_point_bin}{numeric, ARI values calculated by c_bin and the truth}
\item{ari_point_vi}{numeric, ARI values calculated by c_vi and the truth}
\item{FS}{numeric vector, the mean False Spaecity for each cluster}
\item{MS}{numeric vector, the mean Missed Spaecity for each cluster}
}
\references{
%% ~put references to the literature/web site here ~
Chris Fraley and Adrian E Raftery. Model-based clustering, discriminant analysis, and density esti- mation. Journal of the American Statistical Association, 97(458):611–631, 2002.

Riccardo Corradin, Antonio Canale, and Bernardo Nipoti. Bnpmix: An r package for bayesian non- parametric modeling via pitman-yor mixtures. Journal of Statistical Software, 100(15):1–33, 2021.
}
\author{
%%  ~~who you are~~
Zhen Wang \email{zwangiowa@gmail.com}
}
\note{
%%  ~~further notes~~
The MSE values, FS, MS values for VS accuracy are calculated with samples of K = M = 3. ARI values are calculated with all posterior samples.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##generate simulation data for split models
simulation_data <- data_gen_split()

##FBMM with VS, hyper-prior for beta_bel
simulation_1 <- simulation_split(simulation_data$W, simulation_data$Z, simulation_data$y, prior="Bessel")

##run posterior inference, burn-in the first 20k samples
post_fbmm_split_1 <- post_inf_split(simulation_1, 20001:1e5, simulation_data)


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
