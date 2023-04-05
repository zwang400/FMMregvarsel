\name{simulation_func_rpms}
\alias{simulation_func_rpms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC Simulations for the RPMS Model}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A simulation function running MCMC simulations for the RPMS model introduced in Barcella et al (SIM, 2016). Accepted data should be one-dim response y and multi-dim (D >= 1) covariates X.
}
\usage{
simulation_func_rpms(X, y, SS = TRUE, N = 1e5, a_w = 1, b_w = 1, mu = 0, a_tau = 1, b_tau = 1,
                     a_zeta = 0, b_zeta = 1, a_lambda = 0.01, b_lambda = 0.005, a_alpha = 0,
                     b_alpha = 0.01, a_eta = 5, b_eta = 2, a_adp = 1, b_adp = 1,
                     k_init = 6, lambda_init = 2, alpha_dp_init = 1, alpha_init = 0, m_aux = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{X}{numeric matrix, covariates matrix
%%     ~~Describe \code{x} here~~
}
\item{y}{numeric vector, response}
\item{SS}{logical, if TRUE, spike and slab prior is specified on beta for cluster-specific variable selection; if FALSE, a continuous Gaussian prior is specified to beta, no variable selection in implemented}
\item{N}{numeric, number of iterations}
\item{a_w, b_w}{numeric, a Beta(a_w, b_w) is assigned to the SS weights w}
\item{mu}{numeric, mean of the slab part of the SS prior on beta}
\item{a_tau, b_tau}{numeric, a Gamma(a_tau, b_tau) hyper-prior for the precision parameter tau of the slab part of the SS prior on beta}
\item{a_zeta, b_zeta}{numeric, a Normal(a_tau, b_tau) hyper-prior for the mean parameter of X, zeta}
\item{a_lambda, b_lambda}{numeric, a Gamma(a_lambda, b_lambda) hyper-prior for the precision parameter of y, lambda}
\item{a_alpha, b_alpha}{numeric, a Normal(a_alpha, b_alpha) hyper-prior for the intercept parameter alpha}
\item{a_eta, b_eta}{numeric, a Gamma(a_eta, b_eta) hyper-prior for the precision parameter of X, eta}
\item{a_adp, b_adp}{numeric, a Gamma(a_adp, b_adp) hyper-prior for the concentration parameter, alpha, of DP}
\item{k_init}{numeric, initial value of K}
\item{lambda_init}{numeric, initial value of lambda}
\item{alpha_dp_init}{numeric, initial value of alpha of DP}
\item{alpha_init}{numeric, initial value of alpha}
\item{m_aux}{numeric, the number of the auxiliary states, m, in Algorithm 8 of Neal (20000, JCGS)}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
A marginal algorithm is implemented, Algorithm 8 of Neal (2000, JCGS) is applied in this algorithm.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{c_post}{numeric matrix of dim N*n, each row represents posterior samples of c for one iteration}
\item{k_post}{numeric vector of length N, posterior samples of K}
\item{alpha_dp_post}{numeric vector of length N, posterior samples of concentration parameter alpha of DP}
\item{Beta_post}{list, each element of the list represents posterior samples of all beta components in each iteration, of length M_post*D}
\item{alpha_post}{list, each element of the list represents posterior samples of all alpha components in each iteration, of length M_post}
\item{zeta_post}{list, each element of the list represents posterior samples of all zeta components in each iteration, of length M_post*D}
\item{lambda_post}{numeric vector, posterior samples of lambda, precision of y}
\item{eta_post}{numeric vector, posterior samples of eta, precision of X}
\item{tau_post}{numeric vector, posterior samples of tau, precision of the slab part in SS prior}
\item{w_post}{numeric matrix of dim N*D, each row represents posterior smaples of the SS weights, w, in one iteration}
\item{a_lambda}{pre-specified values}
\item{b_lambda}{pre-specified values}
\item{a_eta}{pre-specified values}
\item{b_eta}{pre-specified values}
\item{X}{covariates input}
\item{y}{response input}
}
\references{
%% ~put references to the literature/web site here ~
William Barcella, Maria De Iorio, Gianluca Baio, and James Malone-Lee. Variable selection in co- variate dependent random partition models: an application to urinary tract infection. Statistics in Medicine, 35(8):1373–1389, 2016.

Radford M. Neal. Markov chain sampling methods for Dirichlet process mixture models. Journal of Computational and Graphical Statistics, 9(2):249–265, 2000.
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

##RPMS with VS
simulation_rpms_1 <- simulation_func_rpms(simulation_data$X, simulation_data$y)

##RPMS without VS
simulation_rpms_2 <- simulation_func_rpms(simulation_data$X, simulation_data$y, SS=F)


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
