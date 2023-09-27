################################
#sim functions
#updated 23/09/24: update step 1 of marginal algorithm to fix a bug
################################
simulation_func <- function(X, y, SS = TRUE, N = 1e5, gamma_hyperprior = TRUE,
                         gamma_fixed = 1, a_gamma = 10, b_gamma = 10, a_unif = 0,
                         #for h ~ Unif(a_unif, 1)
                         algorithm = "MAR", prior = "Dirichlet",
                         #algorithm can be "MAR" or "Con",
                         #"Dynamic_FDMM" (a conditional algorithm for dynamic FDMM),
                         #"MAR_Dir_MH" for the marginal algorithm of FDMM based on Miller-Harrison (2018, JASA)
                         #prior can be "Dirichlet", "Bessel" or "Uniform" (unif not available for MAR yet)
                         a_w = 1, b_w = 1, Lambda = 3, a_bessel = 2,
                         b_bessel_hyperprior = TRUE, b_bessel_fixed = 1.1, a_b_bessel = 1, b_b_bessel = 10,
                         mu = 0, a_tau = 1, b_tau = 1, a_zeta = 0, b_zeta = 1,
                         a_lambda = 0.01, b_lambda = 0.005, a_alpha = 0, b_alpha = 0.01,
                         a_eta = 5, b_eta = 2, L_dynamic = 10, M_init = 6,
                         lambda_init = 2, alpha_init = 0){
    #X, y are data
    #gamma is parameter of Dirichlet, w is SS weights, Lambda is parameter of the Poisson prior of M
    #M_init, lambda_init are initial values of M and lambda
    n <- length(y)
    D <- ncol(X)
    if(algorithm == "Con" && prior == "Dirichlet"){
        #initializations
        U <- rep(1, N)
        M <- M_init
        S <- rep(1, M)
        Beta <- matrix(rep(-0.1, D*M), ncol = D)
        alpha <- rep(alpha_init, M)
        zeta <- matrix(rep(1, D*M), ncol = D)
        lambda <- lambda_init
        eta <- 1
        gamma <- gamma_fixed
        tau <- rep(1, D)
        w <- rep(0, D)
        c <- c(1:k, sample(1:k, n-k, replace = T)) #changed to make sure all k indies appear
        M_post <- rep(0, N)
        k_post <- rep(0, N)
        c_post <- matrix(0, ncol = n, nrow = N)
        Beta_post <- list(0)
        alpha_post <- list(0)
        zeta_post <- list(0)
        lambda_post <- rep(0, N)
        eta_post <- rep(0, N)
        gamma_post <- rep(1, N)
        weight_post <- list(0)
        tau_post <- matrix(1, ncol = D, nrow = N)
        w_post <- matrix(0.5, ncol = D, nrow = N)

        group_member <- function(m, Beta, zeta, lambda, alpha, eta){
            llkhd <- rep(0, n)
            for(i in 1:n){
                llkhd_x <- 0
                for(j in 1:D){
                    llkhd_x <- llkhd_x + dnorm(X[i,j], zeta[m, j], sqrt(1/eta), log = T)
                }
                llkhd[i] <- dnorm(y[i], alpha[m] + X[i,]%*%Beta[m,], sqrt(1/lambda), log = T) + llkhd_x
            }
            llkhd
        }

        #group_member <- function(m, Beta, zeta, lambda, alpha, eta){
        #    dnorm(y, alpha[m] + X%*%Beta[m,], sqrt(1/lambda), log = T) +
        #        apply(dnorm(t(X), zeta[m, ], sqrt(1/eta), log = T), 2, sum)
        #}

        #log of posterior pdf of \gamma
        pdf.gamma <- function(gamma, Lambda, u, a_gamma, b_gamma, par_vec){
            k <- length(par_vec)
            psi <- 1/((1+u)^gamma)
            log(Lambda*psi + k) + Lambda*psi + k*log(psi) + (a_gamma-1)*log(gamma) - b_gamma*gamma +
                sum(lgamma(gamma+par_vec) - lgamma(gamma))
        }

        for(i in 1:N){
            #step 1
            T <- sum(S)
            u <- rgamma(1, n, T)
            U[i] <- u

            #step 2
            mem_wgt <- matrix(0, nrow = n, ncol = M)
            for(m in 1:M){
                mem_wgt[ ,m] <- log(S[m]) + group_member(m, Beta, zeta, lambda, alpha, eta)
            }

            for(j in 1:n){
                wgt_shift <- max(mem_wgt[j, ]) #exp normalization trick
                c[j] <- sample(1:M, size = 1, prob = exp(mem_wgt[j, ] - wgt_shift) /
                                   sum(exp(mem_wgt[j, ] - wgt_shift)))
            }
            unique_c <- unique(c)
            k <- length(unique_c)
            k_post[i] <- k
            #relabeling
            Beta <- Beta[unique_c, ]
            zeta <- zeta[unique_c, ]
            alpha <- alpha[unique_c]
            b<-rep(0, n)
            for(j in 1:k){
                b[which(c == unique_c[j])] <- j
            }
            c <- b
            c_post[i, ] <- c

            #step 3a(need to further check)
            p_pois <- Lambda/((u+1)^gamma)
            p_coe <- Lambda/(((u+1)^gamma) * k + Lambda)
            p <- rbinom(1, 1, p_coe)
            M_na <- ifelse(p == 0, rpois(1, p_pois), rpois(1, p_pois)+1)
            M <- k + M_na
            M_post[i] <- M

            par_vec <- as.vector(table(c))
            #step 3b(need to check if k < M, i.e., if M_na == 0)
            if(k < M){
                S <- c(rgamma(k, par_vec+gamma, u+1), rgamma(M-k, gamma, u+1))
                Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = D)) #fill in the M-k gap
                zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = D))
                alpha <- c(alpha, rep(0, M-k))
                for(m in 1:k){
                    V_m <- which(c == m)
                    for(d in 1:D){
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                    }
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m, , drop=F] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                    for(d in 1:D){
                        C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d, drop=F] %*% Beta[m, -d])*
                                                          X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d, drop=F] %*%
                                                                                         Beta[m, -d])*X[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                    }
                }
                for(m in (k+1):M){
                    for(d in 1:D){
                        p_w <- rbinom(1, 1, w[d])
                        Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau)))
                    }
                    zeta[m, ] <- rnorm(D, a_zeta, 1/sqrt(b_zeta))
                }
                alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
            }
            else{
                if(k == 1){
                    Beta <- matrix(Beta, nrow = 1)
                    zeta <- matrix(zeta, nrow = 1)
                }
                S <- rgamma(k, par_vec+gamma, u+1)
                for(m in 1:k){
                    V_m <- which(c == m)
                    for(d in 1:D){
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                    }
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m, , drop=F] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                    for(d in 1:D){
                        C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                          X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                         Beta[m, -d])*X[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                    }
                }
            }

            #sample \lambda and \eta
            lambda <- rgamma(1, a_lambda + n/2, b_lambda + sum((y-alpha[c]-apply(Beta[c,]*X, 1, sum))^2)/2)
            eta <- rgamma(1, a_eta + n*D/2, b_eta + sum((X - zeta[c,])^2)/2)

            #update \tau and w
            if(SS == TRUE){
                for(d in 1:D){
                    w[d] <- rbeta(1, a_w + sum(Beta[,d] == 0), b_w + sum(Beta[,d] != 0))
                }
            }

            for(d in 1:D){
                nd_Plus <- sum(Beta[,d] != 0)
                tau[d] <- rgamma(1, a_tau + nd_Plus/2, b_tau +sum(Beta[,d]^2)/2)
            }

            #update \gamma
            if(gamma_hyperprior == TRUE){
                gamma.new <- rgamma(1, 10, 10)
                acpt.mh.gamma <- min(pdf.gamma(gamma.new, Lambda, u, a_gamma, b_gamma, par_vec) -
                                         pdf.gamma(gamma, Lambda, u, a_gamma, b_gamma, par_vec) +
                                         dgamma(gamma, 10, 10, log = T) - dgamma(gamma.new, 10, 10, log = T), 0)
                gamma <- ifelse(runif(1) <= exp(acpt.mh.gamma), gamma.new, gamma)
            }

            Beta_post[[i]] <- as.vector(Beta)
            alpha_post[[i]] <- alpha
            zeta_post[[i]] <- as.vector(zeta)
            lambda_post[i] <- lambda
            eta_post[i] <- eta
            tau_post[i,] <- tau
            gamma_post[i] <- gamma
            w_post[i,] <- w
            weight_post[[i]] <- S/sum(S)
        }
        sim_res <- list(M_post = M_post, c_post = c_post, k_post = k_post, U = U, gamma_post = gamma_post,
                        Beta_post = Beta_post, alpha_post = alpha_post, zeta_post = zeta_post,
                        lambda_post = lambda_post, eta_post = eta_post, tau_post = tau_post,
                        w_post = w_post, weight_post = weight_post, a_lambda = a_lambda, b_lambda = b_lambda,
                        a_eta = a_eta, b_eta = b_eta, X = X, y = y)
    }
    if(algorithm == "Con" && prior == "Bessel"){
        #initializations
        U <- rep(1, N)
        M <- M_init
        #this may matter as the chain becomes sticky
        S <- rep(1, M)
        Beta <- matrix(rep(-0.1, D*M), ncol = D)
        alpha <- rep(alpha_init, M)
        zeta <- matrix(rep(1, D*M), ncol = D)
        lambda <- lambda_init
        eta <- 1
        b_bessel <- b_bessel_fixed
        tau <- rep(1, D)
        w <- rep(0, D)
        c <- c(1:k, sample(1:k, n-k, replace = T)) #changed to make sure all k indies appear
        M_post <- rep(0, N)
        k_post <- rep(0, N)
        c_post <- matrix(0, ncol = n, nrow = N)
        Beta_post <- list(0)
        alpha_post <- list(0)
        zeta_post <- list(0)
        lambda_post <- rep(0, N)
        eta_post <- rep(0, N)
        b_bessel_post <- rep(1, N)
        weight_post <- list(0)
        tau_post <- matrix(1, ncol = D, nrow = N)
        w_post <- matrix(0.5, ncol = D, nrow = N)

        group_member <- function(m, Beta, zeta, lambda, alpha, eta){
            llkhd <- rep(0, n)
            for(i in 1:n){
                llkhd_x <- 0
                for(j in 1:D){
                    llkhd_x <- llkhd_x + dnorm(X[i,j], zeta[m, j], sqrt(1/eta), log = T)
                }
                llkhd[i] <- dnorm(y[i], alpha[m] + X[i,]%*%Beta[m,], sqrt(1/lambda), log = T) + llkhd_x
            }
            llkhd
        }

        #group_member <- function(m, Beta, zeta, lambda, alpha, eta){
        #    dnorm(y, alpha[m] + X%*%Beta[m,], sqrt(1/lambda), log = T) +
        #        apply(dnorm(t(X), zeta[m, ], sqrt(1/eta), log = T), 2, sum)
        #}

        #log-density, without the constant term
        bessel.exp.tilt.core <- function(x, a, b, u){
            lgamma(2*x + a) -lgamma(x + 1) - lgamma(x + a + 1) - 2*x* (log(2) + log(b + u))
        }

        bessel.gamma.tilt <- function(l, a, b, u, n_j){
            ans <- -log(hyperg_2F1((n_j+a)/2, (n_j+a+1)/2, a+1, 1/(u+b)^2)) + lgamma(a+1) +
                lgamma(2*l+a+n_j) - lgamma(a+n_j) - lgamma(l+1) -lgamma(a+l+1) - 2*l*log(2*u + 2*b)
            exp(ans)
        }

        bessel.gamma.tilt.core <- function(l, a, b, u, n_j){
            lgamma(2*l+a+n_j) - lgamma(l+1) -lgamma(a+l+1) - 2*l*log(2*u + 2*b)
        }
        #for exp-tilted Bessel
        sample.exp.bel <- function(a, b, u){
            q <- 1/(b+u)^2
            U <- 1
            tol <- 0
            while(U > tol){
                l.sam <- rgeom(1, 1-q)
                U <- log(runif(1))
                tol <- bessel.exp.tilt.core(l.sam, a, b, u) - dgeom(l.sam, 1-q, log = T) + log(1-q) + log(a)
            }
            rgamma(1, 2*l.sam + a, u + b)
        }

        #for gamma-tilted Bessel
        sample.gamma.bel <- function(a, b, u, n_j){
            e <- ifelse((u+b)^2 > 2, 2/(u+b)^2, (b/(u+b))^2)
            coef <- ifelse((u+b)^2 > 2, 2, b^2)
            coef1 <- 1
            coef2 <- a+1
            coef3 <- (a+n_j)/2
            coef4 <- (a+n_j+1)/2
            delta <- (coef*(coef1+coef2) - (coef3+coef4))^2 - 4*(coef-1)*(coef*coef1*coef2 - coef3*coef4)
            if(delta <= 0){
                ls <- 2*n_j
            }else{
                ls <- max(ceiling((coef3+coef4-coef*(coef1+coef2) + sqrt(delta))/(2*(coef-1))), 2*n_j)
            }
            p_vec <- bessel.gamma.tilt(0:ls, a, b, u, n_j)
            p <- sum(p_vec)
            uu <- runif(1)
            #when extreme cases happen, where p_vec returns NA element
            #we use max(ls, 100) as the cutoff to do approximate calculations
            if(is.na(p)){
                new_cap = max(ls, 100)
                p_vec_new <- exp(bessel.gamma.tilt.core(0:new_cap, a, b, u, n_j))
                l.sam <- sample(0:new_cap, 1, prob = p_vec_new)
            }else{
                if(uu <= p){
                    l.sam <- sample(0:ls, 1, prob = p_vec)
                }else{
                    U <- 1
                    tol <- 0
                    while(U > tol){
                        l.sam <- rgeom(1, 1-e)
                        U <- log(runif(1))
                        tol <- bessel.gamma.tilt.core(l.sam, a, b, u, n_j) - bessel.gamma.tilt.core(ls+1, a, b, u, n_j)
                        + dgeom(ls+1, 1-e, log = T) - dgeom(l.sam, 1-e, log = T)
                    }
                }
            }
            rgamma(1, 2*l.sam + a + n_j, u + b)
        }

        #log pdf of \beta - 1
        pdf.b.bel <- function(bb_bessel, a_bessel, Lambda, u, a_b_bessel, b_b_bessel, par_vec){
            k <- length(par_vec)
            b_bessel <- bb_bessel + 1
            log_psi <- a_bessel*(log(b_bessel+sqrt(b_bessel^2-1)) - log(b_bessel+u+sqrt((b_bessel+u)^2-1)))
            res <- log(Lambda*exp(log_psi) + k) + Lambda*exp(log_psi) + (a_b_bessel-1)*log(bb_bessel) -
                b_b_bessel*bb_bessel
            for(j in 1:k){
                res <- res + log(a_bessel) + a_bessel*log(b_bessel+sqrt(b_bessel^2-1)) - a_bessel*log(2) -
                    (par_vec[j]+a_bessel)*log(u+b_bessel) + lgamma(a_bessel + par_vec[j]) - lgamma(a_bessel + 1) +
                    log(hyperg_2F1((par_vec[j] + a_bessel)/2, (par_vec[j] + a_bessel + 1)/2,
                                   a_bessel + 1, (u + b_bessel)^(-2)))
            }
            res
        }

        #main program
        for(i in 1:N){
            #step 1
            T <- sum(S)
            u <- rgamma(1, n, T)
            U[i] <- u

            #step 2
            mem_wgt <- matrix(0, nrow = n, ncol = M)
            for(m in 1:M){
                mem_wgt[ ,m] <- log(S[m]) + group_member(m, Beta, zeta, lambda, alpha, eta)
            }

            for(j in 1:n){
                wgt_shift <- max(mem_wgt[j, ]) #exp normalization trick
                c[j] <- sample(1:M, size = 1, prob = exp(mem_wgt[j, ] - wgt_shift) /
                                   sum(exp(mem_wgt[j, ] - wgt_shift)))
            }

            unique_c <- unique(c)
            k <- length(unique_c)
            k_post[i] <- k
            #relabeling
            Beta <- Beta[unique_c, ]
            zeta <- zeta[unique_c, ]
            alpha <- alpha[unique_c]
            b<-rep(0, n)
            for(j in 1:k){
                b[which(c == unique_c[j])] <- j
            }
            c <- b
            c_post[i, ] <- c

            #step 3a(need to further check)
            log_psi <- a_bessel*(log(b_bessel+sqrt(b_bessel^2-1)) - log(b_bessel+u+sqrt((b_bessel+u)^2-1)))
            p_pois <- Lambda*exp(log_psi)
            p_coe <- exp(log(k) - log(k+Lambda*exp(log_psi)))
            p <- rbinom(1, 1, p_coe)
            M_na <- ifelse(p == 1, rpois(1, p_pois), rpois(1, p_pois)+1)
            M <- k + M_na
            M_post[i] <- M

            #step 3b(need to check if k < M, i.e., if M_na == 0)
            S <- rep(1, M)
            if(k < M){
                Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = D)) #fill in the M-k gap
                zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = D))
                alpha <- c(alpha, rep(0, M-k))
                for(m in 1:k){
                    V_m <- which(c == m)
                    S[m] <- sample.gamma.bel(a_bessel, b_bessel, u, length(V_m))
                    for(d in 1:D){
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                    }
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                    for(d in 1:D){
                        C_md <- sqrt(tau[d]/(tau[d] + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau[d]*mu^2)*
                            exp(((mu*tau[d] + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                             X[V_m, d]))^2)/(2*(tau[d] + lambda*sum(X[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                         Beta[m, -d])*X[V_m, d]) +
                                                                             mu * tau[d])/(tau[d] + lambda*sum(X[V_m, d]^2)),
                                                                     1/sqrt(tau[d] + lambda*sum(X[V_m, d]^2))))
                    }
                }
                for(m in (k+1):M){
                    S[m] <- sample.exp.bel(a_bessel, b_bessel, u)
                    for(d in 1:D){
                        p_w <- rbinom(1, 1, w[d])
                        Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau[d])))
                    }
                    zeta[m, ] <- rnorm(D, a_zeta, 1/sqrt(b_zeta))
                }
                alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
            }
            else{
                if(k == 1){
                    Beta <- matrix(Beta, nrow = 1)
                    zeta <- matrix(zeta, nrow = 1)
                }
                for(m in 1:k){
                    V_m <- which(c == m)
                    S[m] <- sample.gamma.bel(a_bessel, b_bessel, u, length(V_m))
                    for(d in 1:D){
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                    }
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                    for(d in 1:D){
                        C_md <- sqrt(tau[d]/(tau[d] + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau[d]*mu^2)*
                            exp(((mu*tau[d] + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                             X[V_m, d]))^2)/(2*(tau[d] + lambda*sum(X[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                         Beta[m, -d])*X[V_m, d]) +
                                                                             mu * tau[d])/(tau[d] + lambda*sum(X[V_m, d]^2)),
                                                                     1/sqrt(tau[d] + lambda*sum(X[V_m, d]^2))))
                    }
                }
            }

            #sample \lambda and \eta
            lambda <- rgamma(1, a_lambda + n/2, b_lambda + sum((y-alpha[c]-apply(Beta[c,]*X, 1, sum))^2)/2)
            eta <- rgamma(1, a_eta + n*D/2, b_eta + sum((X - zeta[c,])^2)/2)

            #update \tau and w
            if(SS == TRUE){
                for(d in 1:D){
                    w[d] <- rbeta(1, a_w + sum(Beta[,d] == 0), b_w + sum(Beta[,d] != 0))
                }
            }

            for(d in 1:D){
                nd_Plus <- sum(Beta[,d] != 0)
                tau[d] <- rgamma(1, a_tau + nd_Plus/2, b_tau +sum(Beta[,d]^2)/2)
            }

            #update \beta of bessel
            if(b_bessel_hyperprior == TRUE){
                par_vec <- as.vector(table(c))
                bb_bessel <- b_bessel - 1
                bb_bessel.new <- rgamma(1, 2, 10)
                acpt.mh.bbel <- min(pdf.b.bel(bb_bessel.new, a_bessel, Lambda, u, a_b_bessel, b_b_bessel, par_vec) -
                                        pdf.b.bel(bb_bessel, a_bessel, Lambda, u, a_b_bessel, b_b_bessel, par_vec) +
                                        dgamma(bb_bessel, 2, 10, log = T) -
                                        dgamma(bb_bessel.new, 2, 10, log = T), 0)
                b_bessel <-ifelse(runif(1) <= exp(acpt.mh.bbel), bb_bessel.new, bb_bessel) + 1
            }


            Beta_post[[i]] <- as.vector(Beta)
            alpha_post[[i]] <- alpha
            zeta_post[[i]] <- as.vector(zeta)
            lambda_post[i] <- lambda
            eta_post[i] <- eta
            tau_post[i,] <- tau
            w_post[i,] <- w
            b_bessel_post[i] <- b_bessel
            weight_post[[i]] <- S/sum(S)
        }
        sim_res <- list(M_post = M_post, c_post = c_post, k_post = k_post, U = U, b_bessel_post = b_bessel_post,
                        Beta_post = Beta_post, alpha_post = alpha_post, zeta_post = zeta_post,
                        lambda_post = lambda_post, eta_post = eta_post, tau_post = tau_post,
                        w_post = w_post, weight_post = weight_post, a_lambda = a_lambda, b_lambda = b_lambda,
                        a_eta = a_eta, b_eta = b_eta, X = X, y = y)
    }
    if(algorithm == "Con" && prior == "Uniform"){
        #initializations
        U <- rep(1, N)
        M <- M_init
        S <- rep(1, M)
        Beta <- matrix(rep(-0.1, D*M), ncol = D)
        alpha <- rep(alpha_init, M)
        zeta <- matrix(rep(1, D*M), ncol = D)
        lambda <- lambda_init
        eta <- 1
        tau <- rep(1, D)
        w <- rep(0, D)
        c <- c(1:k, sample(1:k, n-k, replace = T)) #changed to make sure all k indies appear
        M_post <- rep(0, N)
        k_post <- rep(0, N)
        c_post <- matrix(0, ncol = n, nrow = N)
        Beta_post <- list(0)
        alpha_post <- list(0)
        zeta_post <- list(0)
        lambda_post <- rep(0, N)
        eta_post <- rep(0, N)
        weight_post <- list(0)
        tau_post <- matrix(1, ncol = D, nrow = N)
        w_post <- matrix(0.5, ncol = D, nrow = N)

        group_member <- function(m, Beta, zeta, lambda, alpha, eta){
            llkhd <- rep(0, n)
            for(i in 1:n){
                llkhd_x <- 0
                for(j in 1:D){
                    llkhd_x <- llkhd_x + dnorm(X[i,j], zeta[m, j], sqrt(1/eta), log = T)
                }
                llkhd[i] <- dnorm(y[i], alpha[m] + X[i,]%*%Beta[m,], sqrt(1/lambda), log = T) + llkhd_x
            }
            llkhd
        }

        #group_member <- function(m, Beta, zeta, lambda, alpha, eta){
        #    dnorm(y, alpha[m] + X%*%Beta[m,], sqrt(1/lambda), log = T) +
        #        apply(dnorm(t(X), zeta[m, ], sqrt(1/eta), log = T), 2, sum)
        #}


        for(i in 1:N){
            #step 1
            T <- sum(S)
            u <- rgamma(1, n, T)
            U[i] <- u

            #step 2
            mem_wgt <- matrix(0, nrow = n, ncol = M)
            for(m in 1:M){
                mem_wgt[ ,m] <- log(S[m]) + group_member(m, Beta, zeta, lambda, alpha, eta)
            }

            for(j in 1:n){
                wgt_shift <- max(mem_wgt[j, ]) #exp normalization trick
                c[j] <- sample(1:M, size = 1, prob = exp(mem_wgt[j, ] - wgt_shift) /
                                   sum(exp(mem_wgt[j, ] - wgt_shift)))
            }
            unique_c <- unique(c)
            k <- length(unique_c)
            k_post[i] <- k
            #relabeling
            Beta <- Beta[unique_c, ]
            zeta <- zeta[unique_c, ]
            alpha <- alpha[unique_c]
            b<-rep(0, n)
            for(j in 1:k){
                b[which(c == unique_c[j])] <- j
            }
            c <- b
            c_post[i, ] <- c

            #step 3a(need to further check)
            psi_u <- (exp(-a_unif*u)-exp(-u))/u
            p_pois <- Lambda*psi_u
            p_coe <- Lambda*psi_u/(k+Lambda*psi_u)
            p <- rbinom(1, 1, p_coe)
            M_na <- ifelse(p == 0, rpois(1, p_pois), rpois(1, p_pois)+1)
            M <- k + M_na
            M_post[i] <- M

            #step 3b(need to check if k < M, i.e., if M_na == 0)
            S <- rep(1, M)
            if(k == 1){
                if(k < M){
                    Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = D)) #fill in the M-k gap
                    zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = D))
                    alpha <- c(alpha, rep(0, M-k))
                    for(m in 1:k){
                        V_m <- which(c == m)
                        S[m] <- rtrunc(1, spec="gamma", a=0, b=1, shape=length(V_m)+1, rate=u)
                        for(d in 1:D){
                            zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                    (eta * length(V_m) + b_zeta),
                                                1/sqrt(eta * length(V_m) + b_zeta))
                        }
                        alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                              (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                        for(d in 1:D){
                            C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                                exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                              X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                            theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                            Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                             Beta[m, -d])*X[V_m, d]) +
                                                                                 mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                         1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                        }
                    }
                    for(m in (k+1):M){
                        S[m] <- rtrunc(1, spec="exp", a=0, b=1, rate=u)
                        for(d in 1:D){
                            p_w <- rbinom(1, 1, w[d])
                            Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau)))
                        }
                        zeta[m, ] <- rnorm(D, a_zeta, 1/sqrt(b_zeta))
                    }
                    alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
                }
                else{
                    Beta <- matrix(Beta, nrow = 1)
                    zeta <- matrix(zeta, nrow = 1)
                    for(m in 1:k){
                        V_m <- which(c == m)
                        S[m] <- rtrunc(1, spec="gamma", a=0, b=1, shape=length(V_m)+1, rate=u)
                        for(d in 1:D){
                            zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                    (eta * length(V_m) + b_zeta),
                                                1/sqrt(eta * length(V_m) + b_zeta))
                        }
                        alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                              (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                        for(d in 1:D){
                            C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                                exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                              X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                            theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                            Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                             Beta[m, -d])*X[V_m, d]) +
                                                                                 mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                         1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                        }
                    }
                }
            }
            else if(k < M){
                Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = D)) #fill in the M-k gap
                zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = D))
                alpha <- c(alpha, rep(0, M-k))
                for(m in 1:k){
                    V_m <- which(c == m)
                    S[m] <- rtrunc(1, spec="gamma", a=0, b=1, shape=length(V_m)+1, rate=u)
                    for(d in 1:D){
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                    }
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                    for(d in 1:D){
                        C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                          X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                         Beta[m, -d])*X[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                    }
                }
                for(m in (k+1):M){
                    S[m] <- rtrunc(1, spec="exp", a=0, b=1, rate=u)
                    for(d in 1:D){
                        p_w <- rbinom(1, 1, w[d])
                        Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau)))
                    }
                    zeta[m, ] <- rnorm(D, a_zeta, 1/sqrt(b_zeta))
                }
                alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
            }
            else{
                for(m in 1:k){
                    V_m <- which(c == m)
                    S[m] <- rtrunc(1, spec="gamma", a=0, b=1, shape=length(V_m)+1, rate=u)
                    for(d in 1:D){
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                    }
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                    for(d in 1:D){
                        C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                          X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                         Beta[m, -d])*X[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                    }
                }
            }

            #sample \lambda and \eta
            lambda <- rgamma(1, a_lambda + n/2, b_lambda + sum((y-alpha[c]-apply(Beta[c,]*X, 1, sum))^2)/2)
            eta <- rgamma(1, a_eta + n*D/2, b_eta + sum((X - zeta[c,])^2)/2)

            #update \tau and w
            if(SS == TRUE){
                for(d in 1:D){
                    w[d] <- rbeta(1, a_w + sum(Beta[,d] == 0), b_w + sum(Beta[,d] != 0))
                }
            }

            for(d in 1:D){
                nd_Plus <- sum(Beta[,d] != 0)
                tau[d] <- rgamma(1, a_tau + nd_Plus/2, b_tau +sum(Beta[,d]^2)/2)
            }

            Beta_post[[i]] <- as.vector(Beta)
            alpha_post[[i]] <- alpha
            zeta_post[[i]] <- as.vector(zeta)
            lambda_post[i] <- lambda
            eta_post[i] <- eta
            tau_post[i,] <- tau
            w_post[i,] <- w
            weight_post[[i]] <- S/sum(S)
        }
        sim_res <- list(M_post = M_post, c_post = c_post, k_post = k_post, U = U,
                        Beta_post = Beta_post, alpha_post = alpha_post, zeta_post = zeta_post,
                        lambda_post = lambda_post, eta_post = eta_post, tau_post = tau_post,
                        w_post = w_post, weight_post = weight_post, a_lambda = a_lambda, b_lambda = b_lambda,
                        a_eta = a_eta, b_eta = b_eta, X = X, y = y)
    }
    if(algorithm == "Dynamic_FDMM"){
        #initializations
        U <- rep(1, N)
        M <- M_init
        S <- rep(1/M, M)
        Beta <- matrix(rep(-0.1, D*M), ncol = D)
        alpha <- rep(alpha_init, M)
        zeta <- matrix(rep(1, D*M), ncol = D)
        lambda <- lambda_init
        eta <- 1
        gamma <- gamma_fixed
        tau <- rep(1, D)
        w <- rep(0, D)
        c <- c(1:k, sample(1:k, n-k, replace = T)) #changed to make sure all k indies appear
        M_post <- rep(0, N)
        k_post <- rep(0, N)
        c_post <- matrix(0, ncol = n, nrow = N)
        Beta_post <- list(0)
        alpha_post <- list(0)
        zeta_post <- list(0)
        lambda_post <- rep(0, N)
        eta_post <- rep(0, N)
        gamma_post <- rep(1, N)
        weight_post <- list(0)
        tau_post <- matrix(1, ncol = D, nrow = N)
        w_post <- matrix(0.5, ncol = D, nrow = N)

        group_member <- function(m, Beta, zeta, lambda, alpha, eta){
            llkhd <- rep(0, n)
            for(i in 1:n){
                llkhd_x <- 0
                for(j in 1:D){
                    llkhd_x <- llkhd_x + dnorm(X[i,j], zeta[m, j], sqrt(1/eta), log = T)
                }
                llkhd[i] <- dnorm(y[i], alpha[m] + X[i,]%*%Beta[m,], sqrt(1/lambda), log = T) + llkhd_x
            }
            llkhd
        }

        #group_member <- function(m, Beta, zeta, lambda, alpha, eta){
        #    dnorm(y, alpha[m] + X%*%Beta[m,], sqrt(1/lambda), log = T) +
        #        apply(dnorm(t(X), zeta[m, ], sqrt(1/eta), log = T), 2, sum)
        #}


        for(i in 1:N){
            #step 1, sample membership vector c, here S stands for mixture weights (after normalized)
            mem_wgt <- matrix(0, nrow = n, ncol = M)
            for(m in 1:M){
                mem_wgt[ ,m] <- log(S[m]) + group_member(m, Beta, zeta, lambda, alpha, eta)
            }

            for(j in 1:n){
                wgt_shift <- max(mem_wgt[j, ]) #exp normalization trick
                c[j] <- sample(1:M, size = 1, prob = exp(mem_wgt[j, ] - wgt_shift) /
                                   sum(exp(mem_wgt[j, ] - wgt_shift)))
            }
            unique_c <- unique(c)
            k <- length(unique_c)
            k_post[i] <- k
            #relabeling
            Beta <- Beta[unique_c, ]
            zeta <- zeta[unique_c, ]
            alpha <- alpha[unique_c]
            b<-rep(0, n)
            for(j in 1:k){
                b[which(c == unique_c[j])] <- j
            }
            c <- b
            c_post[i, ] <- c

            #step3a, update M
            #function to calculate the posterior prob mass of M: (5.3) of Fru ̈hwirth-Schnatter(2021, BA)
            M_post_mass_log <- function(Lambda, gamma, M, k, par_vec){
                (M-1)*log(Lambda) - (k-1)*log(M) - lgamma(M-k+1) + sum(lgamma(par_vec+gamma/M) - lgamma(1+gamma/M))
            }

            #function to estimate the prob mass vector of M
            #the truncation is done when the prob is decreasing & the new prob mass is less then 1e-5
            #of the sum of the previous ones
            M_post_mass_est <- function(Lambda, gamma, k, par_vec){
                post_mass <- numeric(0)
                M <- k
                mass <- 0
                new_mass <- 1
                new_mass_norm <- 1
                post_mass_normalized <- 2
                while(mass < new_mass ||  (new_mass_norm / (sum(post_mass_normalized) - new_mass_norm) > 1e-5)){
                    mass <- M_post_mass_log(Lambda, gamma, M, k, par_vec)
                    post_mass <- c(post_mass, mass)
                    M <- M + 1
                    new_mass <- M_post_mass_log(Lambda, gamma, M, k, par_vec)
                    new_post_mass <- c(post_mass, new_mass)
                    post_mass_normalized <- exp(new_post_mass - max(new_post_mass)) / sum(exp(new_post_mass - max(new_post_mass)))
                    new_mass_norm <- post_mass_normalized[length(post_mass_normalized)]
                }
                post_mass_normalized
            }

            par_vec <- as.vector(table(c))
            M_post_prob <- M_post_mass_est(Lambda, gamma, k, par_vec)
            length_M_est <- length(M_post_prob)
            M <- sample((k:(k+length_M_est-1)), 1, replace = F, prob = M_post_prob)
            M_post[i] <- M

            #step 2a, sample \theta^*
            if(k == 1){
                if(k < M){
                    Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = D)) #fill in the M-k gap
                    zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = D))
                    alpha <- c(alpha, rep(0, M-k))
                    for(m in 1:k){
                        V_m <- which(c == m)
                        for(d in 1:D){
                            zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                    (eta * length(V_m) + b_zeta),
                                                1/sqrt(eta * length(V_m) + b_zeta))
                        }
                        alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                              (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                        for(d in 1:D){
                            C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                                exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                              X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                            theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                            Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                             Beta[m, -d])*X[V_m, d]) +
                                                                                 mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                         1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                        }
                    }
                    for(m in (k+1):M){
                        for(d in 1:D){
                            p_w <- rbinom(1, 1, w[d])
                            Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau)))
                        }
                        zeta[m, ] <- rnorm(D, a_zeta, 1/sqrt(b_zeta))
                    }
                    alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
                }
                else{
                    Beta <- matrix(Beta, nrow = 1)
                    zeta <- matrix(zeta, nrow = 1)
                    for(m in 1:k){
                        V_m <- which(c == m)
                        for(d in 1:D){
                            zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                    (eta * length(V_m) + b_zeta),
                                                1/sqrt(eta * length(V_m) + b_zeta))
                        }
                        alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                              (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                        for(d in 1:D){
                            C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                                exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                              X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                            theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                            Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                             Beta[m, -d])*X[V_m, d]) +
                                                                                 mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                         1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                        }
                    }
                }
            }
            else if(k < M){
                Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = D)) #fill in the M-k gap
                zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = D))
                alpha <- c(alpha, rep(0, M-k))
                for(m in 1:k){
                    V_m <- which(c == m)
                    for(d in 1:D){
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                    }
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                    for(d in 1:D){
                        C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                          X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                         Beta[m, -d])*X[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                    }
                }
                for(m in (k+1):M){
                    for(d in 1:D){
                        p_w <- rbinom(1, 1, w[d])
                        Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau)))
                    }
                    zeta[m, ] <- rnorm(D, a_zeta, 1/sqrt(b_zeta))
                }
                alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
            }
            else{
                for(m in 1:k){
                    V_m <- which(c == m)
                    for(d in 1:D){
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                    }
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                    for(d in 1:D){
                        C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                          X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                         Beta[m, -d])*X[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                    }
                }
            }

            #step 3b, update gamma

            #log pdf of log(gamma), note the last term is the Jacobian of the transformation
            log_gamma_post <- function(log_gamma, a_gamma, b_gamma, M, k, par_vec){
                -b_gamma*exp(log_gamma) + log_gamma*(a_gamma-1+k) + lgamma(exp(log_gamma)) - lgamma(n+exp(log_gamma)) +
                    sum(lgamma(par_vec+exp(log_gamma)/M) - lgamma(1+exp(log_gamma)/M)) + log_gamma
            }

            #a random walk MH for log(gamma)
            if(gamma_hyperprior == TRUE){
                add_gamma <- rnorm(1, 0, 1)
                log_gamma_new <- log(gamma)+add_gamma
                accept.random.walk.gamma <- min((log_gamma_post(log_gamma_new, a_gamma, b_gamma, M, k, par_vec) -
                                                     log_gamma_post(log(gamma), a_gamma, b_gamma, M, k, par_vec)), 0)
                gamma <- exp(ifelse(runif(1) <= exp(accept.random.walk.gamma), log_gamma_new, log(gamma)))
            }

            #sample \lambda and \eta
            lambda <- rgamma(1, a_lambda + n/2, b_lambda + sum((y-alpha[c]-apply(Beta[c,]*X, 1, sum))^2)/2)
            eta <- rgamma(1, a_eta + n*D/2, b_eta + sum((X - zeta[c,])^2)/2)

            #update \tau and w
            if(SS == TRUE){
                for(d in 1:D){
                    w[d] <- rbeta(1, a_w + sum(Beta[,d] == 0), b_w + sum(Beta[,d] != 0))
                }
            }

            for(d in 1:D){
                nd_Plus <- sum(Beta[,d] != 0)
                tau[d] <- rgamma(1, a_tau + nd_Plus/2, b_tau +sum(Beta[,d]^2)/2)
            }

            #update S, mixture weights
            if(k < M){
                par_vec_complete <- c(par_vec, rep(0, M-k))
            }else{
                par_vec_complete <- par_vec
            }
            S <- rdirichlet(1, rep(gamma/M, M) + par_vec_complete)


            Beta_post[[i]] <- as.vector(Beta)
            alpha_post[[i]] <- alpha
            zeta_post[[i]] <- as.vector(zeta)
            lambda_post[i] <- lambda
            eta_post[i] <- eta
            tau_post[i,] <- tau
            gamma_post[i] <- gamma
            w_post[i,] <- w
            weight_post[[i]] <- S
        }
        sim_res <- list(M_post = M_post, c_post = c_post, k_post = k_post, U = U, gamma_post = gamma_post,
                        Beta_post = Beta_post, alpha_post = alpha_post, zeta_post = zeta_post,
                        lambda_post = lambda_post, eta_post = eta_post, tau_post = tau_post,
                        w_post = w_post, weight_post = weight_post, a_lambda = a_lambda, b_lambda = b_lambda,
                        a_eta = a_eta, b_eta = b_eta, X = X, y = y)
    }
    if(algorithm == "MAR_Dir_MH"){
        #update step i of updating c, 2023-04-20
        #initializations
        k <- M_init
        Beta <- matrix(rep(-0.1, D*k), ncol = D)
        alpha <- rep(alpha_init, k)
        zeta <- matrix(rep(1, D*k), ncol = D)
        lambda <- lambda_init
        eta <- 1
        gamma <- gamma_fixed
        tau <- rep(1, D)
        w <- rep(0, D)
        c <- c(1:k, sample(1:k, n-k, replace = T)) #changed to make sure all k indies appear
        k_post <- rep(0, N)
        c_post <- matrix(0, ncol = n, nrow = N)
        Beta_post <- list(0)
        alpha_post <- list(0)
        zeta_post <- list(0)
        lambda_post <- rep(0, N)
        eta_post <- rep(0, N)
        u_post <- rep(0, N)
        gamma_post <- rep(1, N)
        tau_post <- matrix(1, ncol = D, nrow = N)
        w_post <- matrix(0.5, ncol = D, nrow = N)
        u <- 1

        group_member <- function(i, m, Beta, zeta, lambda, alpha, eta){
            sum(dnorm(X[i,], zeta[m, ], sqrt(1/eta), log = T)) +
                dnorm(y[i], alpha[m] + X[i,]%*%Beta[m,], sqrt(1/lambda), log = T)
        }

        #log of posterior pdf of \gamma
        pdf.gamma <- function(gamma, Lambda, u, a_gamma, b_gamma, par_vec){
            k <- length(par_vec)
            psi <- 1/((1+u)^gamma)
            log(Lambda*psi + k) + Lambda*psi + k*log(psi) + (a_gamma-1)*log(gamma) - b_gamma*gamma +
                sum(lgamma(gamma+par_vec) - lgamma(gamma))
        }

        #log of posterior pdf of u
        pdf.u <- function(n, gamma, Lambda, u, par_vec){
            k <- length(par_vec)
            psi <- 1/((1+u)^gamma)
            log(Lambda*psi + k) + Lambda*psi - (k*gamma+n)*log(1+u) + (n-1)*log(u)
        }

        #a function to pre-compute V_n(k)
        #log of the sum term for k
        V_k <- function(n, t, k, dg, Lambda){
            lgamma(k+1) - lgamma(k-t+1) + dpois(k-1, Lambda, log = T) - lgamma(k*dg+n) + lgamma(k*dg)
        }

        #calculate V_n(t+1)/V_n(t), use k=1e5 as upper bound
        V_t_tp1 <- function(n, t, dg, Lambda){
            k <- 100
            v_t <- V_k(n, t, t:k, dg, Lambda)
            v_tp1 <- V_k(n, t+1, (t+1):k, dg, Lambda)
            v_t_max <- max(v_t)
            v_tp1_max <- max(v_tp1)
            v_tp1_max + log(sum(exp(v_tp1 - v_tp1_max))) - v_t_max - log(sum(exp(v_t - v_t_max)))
        }

        if(gamma_hyperprior == FALSE){
            pre_comp_terms <- 10
            v_n_t_ratio <- rep(0, pre_comp_terms)
            for(t in 1:pre_comp_terms){
                v_n_t_ratio[t] <- V_t_tp1(n, t, gamma, Lambda)
            }
        }

        for(i in 1:N){
            if(gamma_hyperprior == TRUE){
                #pre-compute V_n(t+1)/V_n(t), pre-compute 50 terms
                #note if we fix \gamma, this pre-compute step can be done outside of this for loop
                pre_comp_terms <- 10
                v_n_t_ratio <- rep(0, pre_comp_terms)
                for(t in 1:pre_comp_terms){
                    v_n_t_ratio[t] <- V_t_tp1(n, t, gamma, Lambda)
                }
            }

            #step 1, update c
            c_size <- as.vector(table(c))
            for(j in 1:n){
                k_c <- length(c_size)
                k_minus <- ifelse(c_size[c[j]] == 1, k_c-1, k_c)
                h <- k_minus + L_dynamic
                if(c_size[c[j]] > 1){
                    #extend Beta, zeta and alpha
                    Beta_ext <- matrix(0, ncol = D, nrow = L_dynamic)
                    for(mm in 1:L_dynamic){
                        pp_w <- rbinom(D, 1, w)
                        Beta_ext[mm,] <- ifelse(pp_w == 1, rep(0, D), rnorm(D, mu, 1/sqrt(tau)))
                    }
                    zeta_ext <- matrix(rnorm(D*L_dynamic, a_zeta, 1/sqrt(b_zeta)), ncol = D)
                    alpha_ext <- rnorm(L_dynamic, a_alpha, 1/sqrt(b_alpha))
                    Beta <- rbind(Beta, Beta_ext)
                    zeta <- rbind(zeta, zeta_ext)
                    alpha <- c(alpha, alpha_ext)

                    #update cluster sizes after removing jthe subject
                    n_size <- c_size
                    n_size[c[j]] <- c_size[c[j]]-1
                }else{
                    Beta_ext <- matrix(0, ncol = D, nrow = (L_dynamic-1))
                    for(mm in 1:(L_dynamic-1)){
                        pp_w <- rbinom(D, 1, w)
                        Beta_ext[mm,] <- ifelse(pp_w == 1, rep(0, D), rnorm(D, mu, 1/sqrt(tau)))
                    }
                    zeta_ext <- matrix(rnorm(D*(L_dynamic-1), a_zeta, 1/sqrt(b_zeta)), ncol = D)
                    alpha_ext <- rnorm(L_dynamic-1, a_alpha, 1/sqrt(b_alpha))

                    #update cluster sizes
                    if(c[j] < max(c)){
                        c_size[c[j]: (k_c-1)] <- c_size[(c[j]+1) : k_c]
                        c_size[k_c] <- 1
                    }

                    #re-labeling
                    c[which(c > c[j])] <- c[which(c > c[j])] - 1
                    c[j] <- k_minus+1
                    cj <- c[j]

                    Beta_ord <- rbind(Beta[-cj, ], Beta[cj, ])
                    zeta_ord <- rbind(zeta[-cj, ], zeta[cj, ])
                    alpha_ord <- c(alpha[-cj], alpha[cj])
                    Beta <- rbind(Beta_ord, Beta_ext)
                    zeta <- rbind(zeta_ord, zeta_ext)
                    alpha <- c(alpha_ord, alpha_ext)

                    #update cluster sizes after removing jthe subject
                    n_size <- c_size[1:k_minus]
                }

                #calculating joint mixture density
                mix_den <- rep(0, h)
                for(hh in 1:h){
                    mix_den[hh] <- group_member(j, hh, Beta, zeta, lambda, alpha, eta)
                }

                #calculate llkhood for ith subject
                mem_wgt <- rep(0, h)
                mem_wgt[1:k_minus] <- log(n_size + gamma) + mix_den[1:k_minus]
                ratio_term <- ifelse(k_minus > pre_comp_terms, V_t_tp1(n, k_minus, gamma, Lambda), v_n_t_ratio[k_minus])
                mem_wgt[(k_minus+1):h] <- ratio_term + log(gamma) + mix_den[(k_minus+1):h] - log(L_dynamic)
                mem_shift <- max(mem_wgt)

                #update c[j], using normalization trick
                c[j] <- sample(1:h, size = 1, prob =
                                   exp(mem_wgt - mem_shift)/ sum(exp(mem_wgt - mem_shift)))
                if(c[j] > k_minus){
                    c[j] <- k_minus + 1

                    #update cluster sizes
                    c_size <- c(n_size, 1)
                }else{
                    #update cluster sizes
                    c_size <- n_size
                    c_size[c[j]] <- n_size[c[j]]+1
                }
            }

            unique_c <- unique(c)
            k <- length(unique_c)
            k_post[i] <- k
            #relabeling
            Beta <- Beta[unique_c, ]
            zeta <- zeta[unique_c, ]
            alpha <- alpha[unique_c]
            b<-rep(0, n)
            for(j in 1:k){
                b[which(c == unique_c[j])] <- j
            }
            c <- b
            c_post[i, ] <- c

            #update \theta
            if(k == 1){
                Beta <- matrix(Beta, nrow = 1)
                zeta <- matrix(zeta, nrow = 1)
            }
            for(m in 1:k){
                V_m <- which(c == m)
                for(d in 1:D){
                    zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                            (eta * length(V_m) + b_zeta),
                                        1/sqrt(eta * length(V_m) + b_zeta))
                }
                alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                      (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                for(d in 1:D){
                    C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                        exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                      X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                    theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                    Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                     Beta[m, -d])*X[V_m, d]) +
                                                                         mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                 1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                }
            }


            #sample \lambda and \eta
            lambda <- rgamma(1, a_lambda + n/2, b_lambda + sum((y-alpha[c]-apply(Beta[c,]*X, 1, sum))^2)/2)
            eta <- rgamma(1, a_eta + n*D/2, b_eta + sum((X - zeta[c,])^2)/2)

            #update \tau and w
            if(SS == TRUE){
                for(d in 1:D){
                    w[d] <- rbeta(1, a_w + sum(Beta[,d] == 0), b_w + sum(Beta[,d] != 0))
                }
            }

            for(d in 1:D){
                nd_Plus <- sum(Beta[,d] != 0)
                tau[d] <- rgamma(1, a_tau + nd_Plus/2, b_tau +sum(Beta[,d]^2)/2)
            }

            #update \gamma
            #note that in this case we'll need to update u!
            if(gamma_hyperprior == TRUE){
                par_vec <- as.vector(table(c))
                #update u
                u.new <- rgamma(1, 50, 1) #need to think about this proposal carefully
                acpt.mh.u <- min(pdf.u(n, gamma, Lambda, u.new, par_vec) -
                                     pdf.u(n, gamma, Lambda, u, par_vec) +
                                     dgamma(u, 50, 1, log = T) - dgamma(u.new, 50, 1, log = T), 0)
                u <- ifelse(runif(1) <= exp(acpt.mh.u), u.new, u)

                #update \gamma
                gamma.new <- rgamma(1, 1, 2)
                acpt.mh.gamma <- min(pdf.gamma(gamma.new, Lambda, u, a_gamma, b_gamma, par_vec) -
                                         pdf.gamma(gamma, Lambda, u, a_gamma, b_gamma, par_vec) +
                                         dgamma(gamma, 1, 2, log = T) - dgamma(gamma.new, 1, 2, log = T), 0)
                gamma <- ifelse(runif(1) <= exp(acpt.mh.gamma), gamma.new, gamma)
            }

            Beta_post[[i]] <- as.vector(Beta)
            alpha_post[[i]] <- alpha
            zeta_post[[i]] <- as.vector(zeta)
            lambda_post[i] <- lambda
            eta_post[i] <- eta
            tau_post[i,] <- tau
            gamma_post[i] <- gamma
            w_post[i,] <- w
            u_post[i] <- u
        }
        #here, M_post = k_post is purely for ease of use of post process function,
        #we do not really update M.
        sim_res <- list(M_post = k_post, c_post = c_post, k_post = k_post, gamma_post = gamma_post,
                        Beta_post = Beta_post, alpha_post = alpha_post, zeta_post = zeta_post,
                        lambda_post = lambda_post, eta_post = eta_post, tau_post = tau_post,
                        w_post = w_post, u_post = u_post, a_lambda = a_lambda, b_lambda = b_lambda,
                        a_eta = a_eta, b_eta = b_eta, X = X, y = y)
    }
     if(algorithm == "MAR"){
        #initializations
        k <- M_init
        Beta <- matrix(rep(-0.1, D*k), ncol = D)
        alpha <- rep(alpha_init, k)
        zeta <- matrix(rep(1, D*k), ncol = D)
        lambda <- lambda_init
        eta <- 1
        u <- 1
        gamma <- gamma_fixed
        b_bessel <- b_bessel_fixed
        tau <- rep(1, D)
        w <- rep(0, D)
        c <- c(1:k, sample(1:k, n-k, replace = T)) #changed to make sure all k indies appear
        k_post <- rep(0, N)
        c_post <- matrix(0, ncol = n, nrow = N)
        Beta_post <- list(0)
        alpha_post <- list(0)
        zeta_post <- list(0)
        lambda_post <- rep(0, N)
        eta_post <- rep(0, N)
        u_post <- rep(0, N)
        gamma_post <- rep(1, N)
        tau_post <- matrix(1, ncol = D, nrow = N)
        w_post <- matrix(0.5, ncol = D, nrow = N)
        b_bessel_post <- rep(1, N)

        #lilelihood function for the ith subject in the mth component
        group_member <- function(i, m, Beta, zeta, lambda, alpha, eta){
            sum(dnorm(X[i,], zeta[m, ], sqrt(1/eta), log = T)) +
                dnorm(y[i], alpha[m] + X[i,]%*%Beta[m,], sqrt(1/lambda), log = T)
        }

        #log of posterior pdf of \gamma of FDMM
        log.pdf.gamma <- function(gamma, Lambda, u, a_gamma, b_gamma, cluster_size){
            k <- length(cluster_size)
            psi <- 1/((1+u)^gamma)
            log(Lambda*psi + k) + Lambda*psi + k*log(psi) + (a_gamma-1)*log(gamma) - b_gamma*gamma +
                sum(lgamma(gamma+cluster_size) - lgamma(gamma))
        }


        #log pdf of \beta - 1 of FBMM
        pdf.b.bel <- function(bb_bessel, a_bessel, Lambda, u, a_b_bessel, b_b_bessel, cluster_size){
            k <- length(cluster_size)
            b_bessel <- bb_bessel + 1
            log_psi <- a_bessel*(log(b_bessel+sqrt(b_bessel^2-1)) - log(b_bessel+u+sqrt((b_bessel+u)^2-1)))
            res <- log(Lambda*exp(log_psi) + k) + Lambda*exp(log_psi) + (a_b_bessel-1)*log(bb_bessel) -
                b_b_bessel*bb_bessel
            for(j in 1:k){
                res <- res + log(a_bessel) + a_bessel*log(b_bessel+sqrt(b_bessel^2-1)) - a_bessel*log(2) -
                    (cluster_size[j]+a_bessel)*log(u+b_bessel) + lgamma(a_bessel + cluster_size[j]) - lgamma(a_bessel + 1) +
                    log(hyperg_2F1((cluster_size[j] + a_bessel)/2, (cluster_size[j] + a_bessel + 1)/2,
                                   a_bessel + 1, (u + b_bessel)^(-2)))
            }
            res
        }

        #log of psi
        log_small_psi <- function(u, prior){
            if(prior == "Dirichlet"){
                -gamma*log(u+1)
            }else if(prior == "Bessel"){
                a_bessel*(log(b_bessel+sqrt(b_bessel^2-1)) - log(b_bessel+u+sqrt((b_bessel+u)^2-1)))
            }else{
                NA
            }
        }

        #log of the Psi
        log_big_psi <- function(Lambda, small_psi, k){
            (k-1)*log(Lambda) + log(Lambda*small_psi + k) + Lambda*(small_psi - 1)
        }

        #log of kappa
        log_kappa <- function(n_par, u, prior){
            if(prior == "Dirichlet"){
                lgamma(gamma+n_par) - lgamma(gamma) - (gamma+n_par)*log(1+u)
            }else if(prior == "Bessel"){
                log(a_bessel) + a_bessel*log(b_bessel+sqrt(b_bessel^2-1)) - a_bessel*log(2) -
                    (n_par+a_bessel)*log(u+b_bessel) + lgamma(a_bessel + n_par) - lgamma(a_bessel + 1) +
                    log(hyperg_2F1((n_par + a_bessel)/2, (n_par + a_bessel + 1)/2,
                                   a_bessel + 1, (u + b_bessel)^(-2)))
            }else{
                NA
            }
        }

        #log of ratio of kappa in updating c
        log_ratio_kappa <- function(n_par, u, prior){
            if(prior == "Dirichlet"){
                log(n_par + gamma) - log(1+u)
            }else if(prior == "Bessel"){
                log(a_bessel + n_par) - log(u + b_bessel) +
                    log(hyperg_2F1((n_par + 1 + a_bessel)/2, (n_par + 1 + a_bessel + 1)/2,
                                   a_bessel + 1, (u + b_bessel)^(-2))) -
                    log(hyperg_2F1((n_par + a_bessel)/2, (n_par + a_bessel + 1)/2,
                                   a_bessel + 1, (u + b_bessel)^(-2)))
            }else{
                NA
            }
        }

        #log pdf of u
        log.pdf.u <- function(n, Lambda, u, cluster_size, prior){
            k <- length(cluster_size)
            small_psi <- exp(log_small_psi(u, prior))
            log_big_psi_value <- log_big_psi(Lambda, small_psi, k)
            (n-1)*log(u) + log_big_psi_value + sum(log_kappa(cluster_size, u, prior))
        }

        for(i in 1:N){
            c_size <- as.vector(table(c))
            #step 1, update c
            for(j in 1:n){
                c_j_temp <- c[j]
                k_c <- length(c_size)
                k_minus <- ifelse(c_size[c_j_temp] == 1, k_c-1, k_c)
                h <- k_minus + L_dynamic

                if(c_size[c_j_temp] > 1){
                    #extend Beta, zeta and alpha
                    Beta_ext <- matrix(0, ncol = D, nrow = L_dynamic)
                    for(mm in 1:L_dynamic){
                        pp_w <- rbinom(D, 1, w)
                        Beta_ext[mm,] <- ifelse(pp_w == 1, rep(0, D), rnorm(D, mu, 1/sqrt(tau)))
                    }
                    zeta_ext <- matrix(rnorm(D*L_dynamic, a_zeta, 1/sqrt(b_zeta)), ncol = D)
                    alpha_ext <- rnorm(L_dynamic, a_alpha, 1/sqrt(b_alpha))
                    Beta <- rbind(Beta, Beta_ext)
                    zeta <- rbind(zeta, zeta_ext)
                    alpha <- c(alpha, alpha_ext)

                    #update cluster sizes after removing jth subject
                    n_size <- c_size
                    n_size[c_j_temp] <- c_size[c_j_temp]-1
                }else{
                    Beta_ext <- matrix(0, ncol = D, nrow = (L_dynamic-1))
                    for(mm in 1:(L_dynamic-1)){
                        pp_w <- rbinom(D, 1, w)
                        Beta_ext[mm,] <- ifelse(pp_w == 1, rep(0, D), rnorm(D, mu, 1/sqrt(tau)))
                    }
                    zeta_ext <- matrix(rnorm(D*(L_dynamic-1), a_zeta, 1/sqrt(b_zeta)), ncol = D)
                    alpha_ext <- rnorm(L_dynamic-1, a_alpha, 1/sqrt(b_alpha))

                    #update cluster sizes
                    if(c_j_temp < max(c)){
                        c_size[c_j_temp: (k_c-1)] <- c_size[(c_j_temp+1) : k_c]
                        c_size[k_c] <- 1

                        #re-labeling
                        c[which(c > c_j_temp)] <- c[which(c > c_j_temp)] - 1
                        c[j] <- k_minus+1
                    }

                    Beta_ord <- rbind(Beta[-c_j_temp, ], Beta[c_j_temp, ])
                    zeta_ord <- rbind(zeta[-c_j_temp, ], zeta[c_j_temp, ])
                    alpha_ord <- c(alpha[-c_j_temp], alpha[c_j_temp])
                    Beta <- rbind(Beta_ord, Beta_ext)
                    zeta <- rbind(zeta_ord, zeta_ext)
                    alpha <- c(alpha_ord, alpha_ext)

                    #update cluster sizes after removing jth subject
                    n_size <- c_size[1:k_minus]
                }

                #calculating joint mixture density
                mix_den <- rep(0, h)
                for(hh in 1:h){
                    mix_den[hh] <- group_member(j, hh, Beta, zeta, lambda, alpha, eta)
                }

                #calculate llkhood for ith subject
                mem_wgt <- rep(0, h)
                small_psi_u <- exp(log_small_psi(u, prior))
                log_ratio_Psi <- log(Lambda) + log(Lambda*small_psi_u + k_minus + 1) - log(Lambda*small_psi_u + k_minus)
                mem_wgt[1:k_minus] <- log_ratio_kappa(n_size, u, prior) + mix_den[1:k_minus]
                mem_wgt[(k_minus+1):h] <- log_ratio_Psi + log_kappa(1, u, prior) + mix_den[(k_minus+1):h] - log(L_dynamic)
                mem_shift <- max(mem_wgt)

                #update c[j], using normalization trick
                c[j] <- sample(1:h, size = 1, prob =
                                   exp(mem_wgt - mem_shift)/ sum(exp(mem_wgt - mem_shift)))
                if(c[j] > k_minus){
                    c[j] <- k_minus + 1

                    #update cluster sizes
                    c_size <- c(n_size, 1)
                }else{
                    #update cluster sizes
                    c_size <- n_size
                    c_size[c[j]] <- n_size[c[j]]+1
                }
            }

            unique_c <- unique(c)
            k <- length(unique_c)
            k_post[i] <- k
            #relabeling
            Beta <- Beta[unique_c, ]
            zeta <- zeta[unique_c, ]
            alpha <- alpha[unique_c]
            b<-rep(0, n)
            for(j in 1:k){
                b[which(c == unique_c[j])] <- j
            }
            c <- b
            c_post[i, ] <- c

            #update \theta
            if(k == 1){
                Beta <- matrix(Beta, nrow = 1)
                zeta <- matrix(zeta, nrow = 1)
            }
            for(m in 1:k){
                V_m <- which(c == m)
                for(d in 1:D){
                    zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(X[V_m, d]))/
                                            (eta * length(V_m) + b_zeta),
                                        1/sqrt(eta * length(V_m) + b_zeta))
                }
                alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - X[V_m,  , drop=F] %*% Beta[m, ]))/
                                      (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))
                for(d in 1:D){
                    C_md <- sqrt(tau/(tau + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                        exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*% Beta[m, -d])*
                                                      X[V_m, d]))^2)/(2*(tau + lambda*sum(X[V_m, d]^2))))
                    theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                    Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d , drop=F] %*%
                                                                                     Beta[m, -d])*X[V_m, d]) +
                                                                         mu * tau)/(tau + lambda*sum(X[V_m, d]^2)),
                                                                 1/sqrt(tau + lambda*sum(X[V_m, d]^2))))
                }
            }

            #sample \lambda and \eta
            lambda <- rgamma(1, a_lambda + n/2, b_lambda + sum((y-alpha[c]-apply(Beta[c,]*X, 1, sum))^2)/2)
            eta <- rgamma(1, a_eta + n*D/2, b_eta + sum((X - zeta[c,])^2)/2)

            #update \tau and w
            if(SS == TRUE){
                for(d in 1:D){
                    w[d] <- rbeta(1, a_w + sum(Beta[,d] == 0), b_w + sum(Beta[,d] != 0))
                }
            }

            for(d in 1:D){
                nd_Plus <- sum(Beta[,d] != 0)
                tau[d] <- rgamma(1, a_tau + nd_Plus/2, b_tau +sum(Beta[,d]^2)/2)
            }

            #update u, independent M-H with a Gamma proposal
            u.new <- rgamma(1, 50, 1) #need to think about this proposal carefully
            acpt.mh.u <- min(log.pdf.u(n, Lambda, u.new, c_size, prior) -
                                 log.pdf.u(n, Lambda, u, c_size, prior) +
                                 dgamma(u, 50, 1, log = T) - dgamma(u.new, 50, 1, log = T), 0)
            u <- ifelse(runif(1) <= exp(acpt.mh.u), u.new, u)

            #update \gamma of FDMM and \beta of FBMM
            if(prior == "Dirichlet"){
                #update \gamma
                if(gamma_hyperprior == TRUE){
                    gamma.new <- rgamma(1, 1, 2)
                    acpt.mh.gamma <- min(log.pdf.gamma(gamma.new, Lambda, u, a_gamma, b_gamma, c_size) -
                                             log.pdf.gamma(gamma, Lambda, u, a_gamma, b_gamma, c_size) +
                                             dgamma(gamma, 1, 2, log = T) - dgamma(gamma.new, 1, 2, log = T), 0)
                    gamma <- ifelse(runif(1) <= exp(acpt.mh.gamma), gamma.new, gamma)
                }
            }

            if(prior == "Bessel"){
                if(b_bessel_hyperprior == TRUE){
                    bb_bessel <- b_bessel - 1
                    bb_bessel.new <- rgamma(1, 2, 10)
                    acpt.mh.bbel <- min(pdf.b.bel(bb_bessel.new, a_bessel, Lambda, u, a_b_bessel, b_b_bessel, c_size) -
                                            pdf.b.bel(bb_bessel, a_bessel, Lambda, u, a_b_bessel, b_b_bessel, c_size) +
                                            dgamma(bb_bessel, 2, 10, log = T) -
                                            dgamma(bb_bessel.new, 2, 10, log = T), 0)
                    b_bessel <-ifelse(runif(1) <= exp(acpt.mh.bbel), bb_bessel.new, bb_bessel) + 1
                }
            }

            Beta_post[[i]] <- as.vector(Beta)
            alpha_post[[i]] <- alpha
            zeta_post[[i]] <- as.vector(zeta)
            lambda_post[i] <- lambda
            eta_post[i] <- eta
            tau_post[i,] <- tau
            gamma_post[i] <- gamma
            b_bessel_post[i] <- b_bessel
            w_post[i,] <- w
            u_post[i] <- u

        }
        #here, M_post = k_post is purely for ease of use of post process function,
        #we do not really update M.
        sim_res <- list(M_post = k_post, c_post = c_post, k_post = k_post, gamma_post = gamma_post,
                        b_bessel_post = b_bessel_post, Beta_post = Beta_post, alpha_post = alpha_post,
                        zeta_post = zeta_post, lambda_post = lambda_post, eta_post = eta_post,
                        tau_post = tau_post, w_post = w_post, u_post = u_post, a_lambda = a_lambda,
                        b_lambda = b_lambda, a_eta = a_eta, b_eta = b_eta, X = X, y = y)
    }
    sim_res
}












