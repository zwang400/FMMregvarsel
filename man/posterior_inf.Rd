\name{posterior_inf}
\alias{posterior_inf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Inference for Various Models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculate different metrics defined in the paper to evaluate performance of different models in : clustering accuracy, parameter estimation accuracy and variable selection accuracy. Besides the N-IFPP models, functions for inference of RPMS, BNPmix (P-Y mixture models, M9 and M10) and mclust (model-based clustering, M7) are also included.
}
\usage{
post_inf(sim_res, scl = 10001:1e5, data)

post_inf_rpms(sim_res, scl = 10001:1e5, data)

post_inf_bnpmix(sim_res, data)

mclust_vs(data)
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
post_inf runs inference for N-IFPP models, post_inf_rpms runs inference results for the RPMS, post_inf_bnpmix runs inference for the two P-Y mixture models (M9, M10). The current versions of the functions can only handle the case of M_true = 3 as defined in the default simulation data setups.

mclust_vs does inference for the model-based clustering model (M7) introduced by Fraley and Raftrey (2002, JASA), two different types of VS procedures were implemented to adapt to the regression setup, which is not considered in the original model. The two methods are:
1) implemented the variable selection methods by the clustvarsel package, obtain the VS and clustering results, then a common set of variables are selected for all clusters, fit a linear regression model within each cluster to obtain estimates of alpha and beta;
2) run the mclust package for clustering with both the response and covariates, based on the clustering results with fixed K=3, run the variable selection procedure by the clustvarsel package for each cluster so a cluster-specific VS can be achieved, then within each cluster, fit a linear regression model to obtain estimates of alpha and beta;
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
post_inf, post_inf_rpms and post_inf_bnpmix return:

\item{true_size}{numeric vector, the true cluster sizes in the data}
\item{measure_km_ari}{numeric, a measurement defined to measure how far away the clusters are from each other, calculated by the ARI values of the true clutering membership and the clustering result of the K-means method, a larger value generally means more distant clusters, thus easier clustering problem}
\item{auto_corr_k}{numeric vector, auto-correlations for posterior of K}
\item{geweke_k}{numeric, the Z-statistic for the geweke diagnostics for posterior of K}
\item{post_k}{numeric vector, posterior distribution of K}
\item{mse_a1, mse_a2, mse_a3}{numeric vectors, the quartiles of the mean squared errors of alpha for the 3 clusters}
\item{mse_b1, mse_b2, mse_b3}{numeric vectors, the quartiles of the mean squared errors of beta for the 3 clusters}
\item{ARI}{numeric vector, the (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1) quantiles of ARI values calculated by all posterior samples of c}
\item{c_bin}{numeric vector, sizes of clusters obtained by the Binder loss}
\item{c_vi}{numeric vector, sizes of clusters obtained by the VI loss}
\item{ari_point_bin}{numeric, ARI values calculated by c_bin and the truth}
\item{ari_point_vi}{numeric, ARI values calculated by c_vi and the truth}
\item{FS}{numeric vector, the mean False Spaecity for each cluster}
\item{MS}{numeric vector, the mean Missed Spaecity for each cluster}
\item{ARI_bnpmix}{numeric vector, the ARI values calculated from all posterior samples of c by BNPmix}

mclust_vs returns:
\item{clust_no_g}{numeric vector, the sizes of clusters by mclust when G is not pre-specified}
\item{ari_no_g}{numeric, ARI values for clustering results by mclust when G is not pre-specified}
\item{cluster_no_g_vs}{numeric, the sizes of clusters by clustvarsel when G is not pre-specified}
\item{ari_no_g_vs}{numeric, ARI values for clustering results by clustvarsel when G is not pre-specified}
\item{clust_g3}{numeric vector, the sizes of clusters by mclust when G is fixed at 3}
\item{ari_g3}{numeric, ARI values for clustering results by mclust when G is fixed at 3}
\item{clust_g3_vs}{numeric, ARI values for clustering results by clustvarsel when G is fixed at 3}
\item{ari_g3_vs}{numeric, ARI values for clustering results by clustvarsel when G is fixed at 3}
\item{se_a1.1, se_a2.1, se_a3.1}{squared-errors of alpha values for M7, method 1)}
\item{se_b1.1, se_b2.1, se_b3.1}{squared-errors of alpha values for M7, method 1)}
\item{FS1.1, FS2.1, FS3.1}{FS values for each cluster for M7, method 1)}
\item{MS1.1, MS2.1, MS3.1}{MS values for each cluster for M7, method 1)}
\item{se_a1.2, se_a2.2, se_a3.2}{squared-errors of alpha values for M7, method 2)}
\item{se_b1.2, se_b2.2, se_b3.2}{squared-errors of alpha values for M7, method 2)}
\item{FS1.2, FS2.2, FS3.2}{FS values for each cluster for M7, method 2)}
\item{MS1.2, MS2.2, MS3.2}{MS values for each cluster for M7, method 2)}
%% ...
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
##generate simulation data
simulation_data <- data_gen_func()

##FBMM with VS, hyper-prior for beta_bel
simulation_1 <- simulation_func(simulation_data$X, simulation_data$y, prior="Bessel")

##run posterior inference, burn-in the first 20k samples
post_fbmm_vs_1 <- post_inf(simulation_1, 20001:1e5, simulation_data)


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
