\name{data_gen_func}
\alias{data_gen_func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate Simulation Data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generate the simulation data for various models. The data consists of single-dim response y and multi-dim covariates X. For subjects in different clusters, (alpha, beta, zeta, lambda) are different. The function can generate correlated data, by specifying the correlation matrix among X, and balanced/inbalanced data by specifying the cluster probabilities.

For the split model, the data consists of y, fixed covariates W and random covariates Z. For subjects in different clusters, (alpha, beta, psi, zeta, lambda) are different.

For covariates following t and Gamma distributions, only independent covariates are geneated.
}
\usage{
data_gen_func(n = 500, alpha_true =  c(0.1, -0.6, 0.5),
             beta_true = rbind(c(0, 0, -0.5, 0, 0.5, 0),
                              c(-0.7, 0, 0.4, 0, 0, 0),
                              c(0.6, 0, 0, 0, -0.4, 0)),
             lambda = c(2, 2, 2), zeta_sep = 1,
             eta = 1, sample_prob = c(1, 1, 1),
             cor_mtx = NULL, rho = rep(0, 3))

data_gen_split(n = 500, alpha_true = c(0.1, -0.6, 0.5),
              beta_true = rbind(c(0, 0, -0.5, 0, 0.5, 0),
                                c(-0.7, 0, 0.4, 0, 0, 0),
                                c(0.6, 0, 0, 0, -0.4, 0)),
              psi_true = rep(c(-1,0,1),3),
              W_mean = 0, lambda = c(2, 2, 2), zeta_sep = 1,
              eta = 1, sample_prob = c(1, 1, 1), rho = 0.5)

data_gen_gamma(n = 500, alpha_true =  c(0.1, -0.6, 0.5),
           beta_true = rbind(c(0, 0, -0.5, 0, 0.5, 0),
                             c(-0.7, 0, 0.4, 0, 0, 0),
                             c(0.6, 0, 0, 0, -0.4, 0)),
           lambda = c(2, 2, 2), zeta_sep = 1, sample_prob = c(1, 1, 1))

data_gen_t(n = 500, alpha_true =  c(0.1, -0.6, 0.5),
           beta_true = rbind(c(0, 0, -0.5, 0, 0.5, 0),
                             c(-0.7, 0, 0.4, 0, 0, 0),
                             c(0.6, 0, 0, 0, -0.4, 0)),
           lambda = c(2, 2, 2), zeta_sep = 1, sample_prob = c(1, 1, 1))


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size
%%     ~~Describe \code{x} here~~
}
  \item{alpha_true}{numeric vector, each element represents intercept for each cluster, length should equal the number of clusters}
  \item{beta_true}{numeric matrix, each row represents regression coefficients for each cluster, number of rows should equal the number of clusters}
  \item{lambda}{numeric vector, each element represent precision of response y in each cluster}
  \item{zeta_sep}{numeric, the difference of zeta values in different clusters, the true zeta values are set to be seq(0, true_M-1) * zeta_sep}
  \item{eta}{numeric, precision of covariates X, for split model, the precision of Z}
  \item{sample_prob}{numeric vector, ratios of sizes of clusters}
  \item{cor_mtx}{numeric matrix, correlation matrix of X, default value is for D=6, M=3, X1-X2, X3-X4, X5-X6 have blocked-wise correlation with coefficients given by rho (no specified in definition)}
  \item{rho}{numeric vector, correlation coefficients among X in the case described above. For split model, rho is a numeric value, which is the pair-wise correlation among columns of W}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
For data_gen_func:

The regression coefficients alpha, beta, the mean of X, zeta, and the precision of y, lambda are cluster-specific.

If no cor_mtx value is specified: 1) if rho is a vector of zero's, generated samples have independent X's; 2) if non-zero values are specified in rho, it's assumed to be the D=6, M=3 case, with a blocked-wise correlations of X1-X2, X3-X4 and X5-X6, and correlation coefficients given by rho, if a different setup (eg, different D or M) is desired, cor_mtx must be specified.

The true value of $\zeta$ for the three groups are 0, zeta_sep and 2*zeta_sep respectively.

For data_gen_split:

The regression coefficients alpha, beta for W, psi for Z, the mean of Z, zeta, and the precision of y, lambda, are cluster-specific.

W has a pair-wise correlation of rho.

For data_gen_gamma:

Covariates are independently following a Gamma(z^2, z) distributions, where for the M groups, z = 2+m*zeta_sep, m = 0, 1, ..., M-1.

For data_gen_t:

Covariates are independently following the t_M/sqrt(M/(M-2)) + z, where t_M is the Student’s t-distribution with a degree of freedom of M, z = m*zeta_sep, m = 0, 1, ..., M-1.

}
\value{
A list is returned:
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{y}{numeric vector, response with length n}
\item{X}{numeric matrix, covariates matrix with dimension M*D}
\item{W}{numeric matrix, fixed covariates matrix in split model}
\item{Z}{numeric matrix, random covariates matrix in split model}
\item{index}{numeric vector with length n, group membership indicator ranging from 1 to M for each subject}
\item{alpha_true}{specified true alpha values}
\item{beta_true}{specified true beta values}
\item{psi_true}{specified true psi values for split model}
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
##generate data with independent X
sample_data_1 <- data_gen_func()

##generate inbalanced data with ratio 1:2:10
sample_data_2 <- data_gen_func(sample_prob = c(1,2,10))

##generate data where D=6, M=3, and X has blocked-wise correlations of X1-X2, X3-X4 and X5-X6, and correlation coefficients (0.2,0.5,0.8)
sample_data_3 <- data_gen_func(rho = c(0.2,0.5,0.8))

##generate data where correlation matrix among X is given by:
require(Matrix)
cor_mtx = bdiag(matrix(c(1/1, 0.8,0.8, 1/1), 2, 2),
                matrix(c(1/1, 0.2,0.2, 1/1), 2, 2),
                matrix(c(1/1, 0.9,0.9, 1/1), 2, 2))
sample_data_4 <- data_gen_func(cor_mtx=cor_mtx)


##generate data for the split model, where W has pair-wise correlation of 0.3
sample_split_data_1 <- data_gen_split(rho=0.3)

##generate data with gamma covariates
sample_data_gamma <- data_gen_gamma()

##generate data with t covariates
sample_data_t <- data_gen_t()

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
