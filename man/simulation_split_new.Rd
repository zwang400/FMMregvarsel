\name{simulation_split_new}
\alias{simulation_split_new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC Simulations for the Bayesian Finite Mixture Regression Models with Cluster-Specific Variable Selection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A simulation function running MCMC simulations for various split models. Accepted data should be one-dim response y and multi-dim (D >= 1) covariates W and Z. Here Z are the covariates being modeled in the model, while W are the fixed covariates. In this split mode, regression coefficients $\psi$ for the fixed W are cluster-specific.
}
\usage{
simulation_split_new(W, Z, y, prior = "Dirichlet", SS = TRUE, N = 1e5, gamma_hyperprior = TRUE,
                     gamma_fixed = 1, a_gamma = 10, b_gamma = 10, a_unif = 0,
                     a_w = 1, b_w = 1, Lambda = 3, a_bessel = 2,
                     b_bessel_hyperprior = TRUE, b_bessel_fixed = 1.1, a_b_bessel = 1, b_b_bessel = 10,
                     mu = 0, a_tau = 1, b_tau = 1, a_zeta = 0, b_zeta = 1,
                     a_psi = 0, a_4 = 1, b_4 = 1,
                     a_lambda = 0.01, b_lambda = 0.005, a_alpha = 0, b_alpha = 0.01,
                     a_eta = 5, b_eta = 2, M_init = 6,
                     lambda_init = 2, alpha_init = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{W}{numeric matrix, covariates matrix that is fixed in a split model}
\item{X}{numeric matrix, covariates matrix that is modeled with Gaussian distributions in a split model}
\item{y}{numeric vector, response}
\item{prior}{char, priors on mixture weights, or different algorithm. "Dirichlet" for FDMM with conditional algorithm, "Bessel" for FBMM with conditional algorithm, "Uniform" for FUMM with conditional algorithm.}
\item{SS}{logical, if TRUE, spike and slab prior is specified on beta for cluster-specific variable selection; if FALSE, a continuous Gaussian prior is specified to beta, no variable selection in implemented}
\item{N}{numeric, number of iterations}
\item{gamma_hyperprior}{logical, if TRUE, for FDMMM with conditional and marginal algorithm, a hyper-prior is specified to the concentration parameter gamma, otherwise, gamma is fixed}
\item{gamma_fixed}{numeric, if gamma_hyperprior=FALSE, the fixed value for gamma}
\item{a_gamma, b_gamma}{numeric, if gamma_hyperprior=TRUE, a Gamma(a_gamma, b_gamma) is assigned to gamma}
\item{a_unif}{numeric, prior for S in FUMM is Uniform(a_unif, 1)}
\item{a_w, b_w}{numeric, a Beta(a_w, b_w) is assigned to the SS weights w}
\item{Lambda}{numeric, M ~ Posson_1(Lambda), a shifted Poisson distribution}
\item{a_bessel}{numeric, for FBMM, alpha_bel is fixed at this value}
\item{b_bessel_hyperprior}{logical, if TRUE, a hyperprior is specified to beta_bel}
\item{b_bessel_fixed}{numeric, if b_bessel_hyperprior=FALSE, beta_bel is fixed at this value}
\item{a_b_bessel, b_b_bessel}{numeric, if b_bessel_hyperprior=TRUE, a Gamma(a_b_bessel, b_b_bessel) hyper-prior is assigned to beta_bel-1}
\item{mu}{numeric, mean of the slab part of the SS prior on beta}
\item{a_tau, b_tau}{numeric, a Gamma(a_tau, b_tau) hyper-prior for the precision parameter tau of the slab part of the SS prior on beta}
\item{a_zeta, b_zeta}{numeric, a Normal(a_tau, b_tau) hyper-prior for the mean parameter of X, zeta}
\item{a_lambda, b_lambda}{numeric, a Gamma(a_lambda, b_lambda) hyper-prior for the precision parameter of y, lambda}
\item{a_alpha, b_alpha}{numeric, a Normal(a_alpha, b_alpha) hyper-prior for the intercept parameter alpha}
\item{a_eta, b_eta}{numeric, a Gamma(a_eta, b_eta) hyper-prior for the precision parameter of X, eta}
\item{a_psi}{numeric, the prior mean for psi, the regression coefficients corresponding to W}
\item{a_4, b_4}{numeric, a Gamma(a_4, b_4) hyper-prior is assigned to b_psi, the prior precison of psi}
\item{M_init}{numeric, initial value of M}
\item{lambda_init}{numeric, initial value of lambda}
\item{alpha_init}{numeric, initial value of alpha}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Conditional algorithms are impelemted for both FDMM and FBMM with variable selection.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{M_post}{numeric vector of length N, posterior samples of M}
\item{c_post}{numeric matrix of dim N*n, each row represents posterior samples of c for one iteration}
\item{k_post}{numeric vector of length N, posterior samples of K}
\item{U}{numeric vector of length N, posterior samples of U, the auxiliary variable for the conditional algorithm of N-IFPP models}
\item{gamma_post}{numeric vector of length N, posterior samples of the concentration parameter gamma in both FDMM models}
\item{Beta_post}{list, each element of the list represents posterior samples of all beta components in each iteration, of length M_post*d2}
\item{alpha_post}{list, each element of the list represents posterior samples of all alpha components in each iteration, of length M_post}
\item{zeta_post}{list, each element of the list represents posterior samples of all zeta components in each iteration, of length M_post*d2}
\item{lambda_post}{numeric vector, posterior samples of lambda, precision of y}
\item{eta_post}{numeric vector, posterior samples of eta, precision of X}
\item{tau_post}{numeric vector, posterior samples of tau, precision of the slab part in SS prior}
\item{w_post}{numeric matrix of dim N*D, each row represents posterior smaples of the SS weights, w, in one iteration}
\item{weight_post}{list, each element represents posterior samples of the mixture weights, omega, of length M_post}
\item{psi_post}{list, each element of the list represents posterior samples of all psi components in each iteration, of length M_post*d1}
\item{b_psi_post}{numeric vector, posterior smaples of b_psi}
\item{a_lambda}{pre-specified values}
\item{b_lambda}{pre-specified values}
\item{a_eta}{pre-specified values}
\item{b_eta}{pre-specified values}
\item{y}{response input}
\item{W}{covariates input that is fixed}
\item{Z}{covariates input that is modeled}
}
\references{
%% ~put references to the literature/web site here ~
Is infinity that far? A Bayesian nonparametric perspective of finite mixture models. Raffaele Argiento and Maria De Iorio, Ann. Statist. 50(5): 2641-2663, 2022.
}
\author{
%%  ~~who you are~~
Zhen Wang \email{zwangiowa@gmail.com}

}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##generate data for the split model
data_split <- data_gen_split()

###split model with FBMM
simulation_split_fbmm <- simulation_split_new(data_split$W, data_split$Z, data_split$y, prior = "Bessel")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
