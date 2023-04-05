simulation_split <- function(W, Z, y, prior = "Dirichlet", SS = TRUE, N = 1e5, gamma_hyperprior = TRUE,
                             gamma_fixed = 1, a_gamma = 1, b_gamma = 1,
                             a_w = 1, b_w = 1, Lambda = 3, a_bessel = 2,
                             b_bessel_hyperprior = TRUE, b_bessel_fixed = 1.1, a_b_bessel = 1, b_b_bessel = 10,
                             mu = 0, a_tau = 1, b_tau = 1, a_zeta = 0, b_zeta = 1,
                             a_psi = 0, a_4 = 1, b_4 = 1,
                             a_lambda = 0.01, b_lambda = 0.005, a_alpha = 0, b_alpha = 0.01,
                             a_eta = 5, b_eta = 2, M_init = 6,
                             lambda_init = 2, alpha_init = 0){
    #X, y are data
    #gamma is parameter of Dirichlet, w is SS weights, Lambda is parameter of the Poisson prior of M
    #M_init, lambda_init are initial values of M and lambda
    n <- length(y)
    d1 <- ncol(W)
    d2 <- ncol(Z)
    if(prior == "Dirichlet"){
        #initializations
        U <- rep(1, N)
        M <- M_init
        S <- rep(1, M)
        Beta <- matrix(rep(-0.1, d2*M), ncol = d2)
        alpha <- rep(alpha_init, M)
        zeta <- matrix(rep(1, d2*M), ncol = d2)
        lambda <- lambda_init
        eta <- 1
        gamma <- gamma_fixed
        tau <- rep(1, d2)
        w <- rep(0, d2)
        c <- rep(0, n)
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
        tau_post <- matrix(1, ncol = d2, nrow = N)
        w_post <- matrix(0.5, ncol = d2, nrow = N)
        psi <- rep(-0.1, d1)
        psi_post <- matrix(1, ncol = d1, nrow = N)
        b_psi <- 1
        b_psi_post <- rep(1, N)

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
        #        apply(dnorm(X, zeta[m, ], sqrt(1/eta), log = T), 1, sum)
        #}



        #log of posterior pdf of \gamma
        pdf.gamma <- function(gamma, Lambda, u, a_gamma, b_gamma, par_vec){
            k <- length(par_vec)
            psi_u <- 1/((1+u)^gamma)
            log(Lambda*psi_u + k) + Lambda*psi_u + k*log(psi_u) + (a_gamma-1)*log(gamma) - b_gamma*gamma +
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
                mem_wgt[ ,m] <- log(S[m]) + group_member(m, Beta, zeta, psi, lambda, alpha, eta)
            }

            for(j in 1:n){
                wgt_shift <- max(mem_wgt[j, ]) #exp normalization trick
                c[j] <- sample(1:M, size = 1, prob = exp(mem_wgt[j, ] - wgt_shift) /
                                   sum(exp(mem_wgt[j, ] - wgt_shift)))
            }
            k <- length(unique(c))
            k_post[i] <- k
            #relabeling
            Beta <- Beta[unique(c), ]
            zeta <- zeta[unique(c), ]
            alpha <- alpha[unique(c)]
            b<-rep(0, n)
            for(j in 1:k){
                b[which(c == unique(c)[j])] <- j
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

            #step 3b(need to check if k < M, i.e., if M_na == 0)
            S <- rep(1, M)
            if(k == 1){
                if(k < M){
                    Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = d2)) #fill in the M-k gap
                    zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = d2))
                    alpha <- c(alpha, rep(0, M-k))
                    for(m in 1:k){
                        V_m <- which(c == m)
                        S[m] <- rgamma(1, length(V_m)+gamma, u+1)
                        #update zeta and beta in this way cause this seems to give better results, tho both methods should work okay.

                        #update \alpha
                        alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - W[V_m, ] %*% psi - Z[V_m, ] %*% Beta[m, ]))/
                                              (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))

                        for(d in 1:d2){
                            #update \zeta
                            zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(Z[V_m, d]))/
                                                    (eta * length(V_m) + b_zeta),
                                                1/sqrt(eta * length(V_m) + b_zeta))
                            #update \beta
                            C_md <- sqrt(tau/(tau + lambda*sum(Z[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                                exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - W[V_m, ] %*% psi - Z[V_m, -d] %*% Beta[m, -d])*
                                                              Z[V_m, d]))^2)/(2*(tau + lambda*sum(Z[V_m, d]^2))))
                            theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                            Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - W[V_m, ] %*% psi
                                                                                         - Z[V_m, -d] %*%
                                                                                             Beta[m, -d])*Z[V_m, d]) +
                                                                                 mu * tau)/(tau + lambda*sum(Z[V_m, d]^2)),
                                                                         1/sqrt(tau + lambda*sum(Z[V_m, d]^2))))
                        }

                    }
                    for(m in (k+1):M){
                        S[m] <- rgamma(1, gamma, u+1)
                        for(d in 1:d2){
                            p_w <- rbinom(1, 1, w[d])
                            Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau)))
                        }
                        zeta[m, ] <- rnorm(d2, a_zeta, 1/sqrt(b_zeta))
                    }
                    alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
                }
                else{
                    Beta <- matrix(Beta, nrow = 1)
                    zeta <- matrix(zeta, nrow = 1)
                    for(m in 1:k){
                        V_m <- which(c == m)
                        S[m] <- rgamma(1, length(V_m)+gamma, u+1)

                        #update \alpha
                        alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - W[V_m, ] %*% psi - Z[V_m, ] %*% Beta[m, ]))/
                                              (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))

                        for(d in 1:d2){
                            #update \zeta
                            zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(Z[V_m, d]))/
                                                    (eta * length(V_m) + b_zeta),
                                                1/sqrt(eta * length(V_m) + b_zeta))
                            #update beta
                            C_md <- sqrt(tau/(tau + lambda*sum(Z[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                                exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] -  W[V_m, ] %*% psi - Z[V_m, -d] %*% Beta[m, -d])*
                                                              Z[V_m, d]))^2)/(2*(tau + lambda*sum(Z[V_m, d]^2))))
                            theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                            Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] -  W[V_m, ] %*% psi
                                                                                         - Z[V_m, -d] %*%
                                                                                             Beta[m, -d])*Z[V_m, d]) +
                                                                                 mu * tau)/(tau + lambda*sum(Z[V_m, d]^2)),
                                                                         1/sqrt(tau + lambda*sum(Z[V_m, d]^2))))
                        }
                    }
                }
            }
            else if(k < M){
                Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = d2)) #fill in the M-k gap
                zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = d2))
                alpha <- c(alpha, rep(0, M-k))
                for(m in 1:k){
                    V_m <- which(c == m)
                    S[m] <- rgamma(1, length(V_m)+gamma, u+1)
                    #update zeta and beta in this way cause this seems to give better results, tho both methods should work okay.

                    #update \alpha
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - W[V_m, ] %*% psi - Z[V_m, ] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))

                    for(d in 1:d2){
                        #update \zeta
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(Z[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                        #update \beta
                        C_md <- sqrt(tau/(tau + lambda*sum(Z[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - W[V_m, ] %*% psi - Z[V_m, -d] %*% Beta[m, -d])*
                                                          Z[V_m, d]))^2)/(2*(tau + lambda*sum(Z[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - W[V_m, ] %*% psi
                                                                                     - Z[V_m, -d] %*%
                                                                                         Beta[m, -d])*Z[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(Z[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(Z[V_m, d]^2))))
                    }

                }
                for(m in (k+1):M){
                    S[m] <- rgamma(1, gamma, u+1)
                    for(d in 1:d2){
                        p_w <- rbinom(1, 1, w[d])
                        Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau)))
                    }
                    zeta[m, ] <- rnorm(d2, a_zeta, 1/sqrt(b_zeta))
                }
                alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
            }
            else{
                for(m in 1:k){
                    V_m <- which(c == m)
                    S[m] <- rgamma(1, length(V_m)+gamma, u+1)

                    #update \alpha
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - W[V_m, ] %*% psi - Z[V_m, ] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))

                    for(d in 1:d2){
                        #update \zeta
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(Z[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                        #update \beta
                        C_md <- sqrt(tau/(tau + lambda*sum(Z[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] -  W[V_m, ] %*% psi - Z[V_m, -d] %*% Beta[m, -d])*
                                                          Z[V_m, d]))^2)/(2*(tau + lambda*sum(Z[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] -  W[V_m, ] %*% psi
                                                                                     - Z[V_m, -d] %*%
                                                                                         Beta[m, -d])*Z[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(Z[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(Z[V_m, d]^2))))
                    }
                }
            }

            #update \psi
            var_mtx <- solve(lambda*(t(W) %*% W) + diag(b_psi, d1))
            psi <- mvrnorm(1, var_mtx %*% (colSums(sweep(W, MARGIN = 1, (y - alpha[c] - apply(Beta[c,]*Z, 1, sum)), "*")) *
                                               lambda + rep(a_psi*b_psi,d1)), var_mtx)


            #sample \lambda and \eta
            lambda <- rgamma(1, a_lambda + n/2, b_lambda + sum((y-alpha[c]-W %*%psi-apply(Beta[c,]*Z, 1, sum))^2)/2)
            eta <- rgamma(1, a_eta + n*d2/2, b_eta + sum((Z - zeta[c,])^2)/2)

            #update \tau and w
            if(SS == TRUE){
                for(d in 1:d2){
                    w[d] <- rbeta(1, a_w + sum(Beta[,d] == 0), b_w + sum(Beta[,d] != 0))
                }
            }

            for(d in 1:d2){
                nd_Plus <- sum(Beta[,d] != 0)
                tau[d] <- rgamma(1, a_tau + nd_Plus/2, b_tau +sum(Beta[,d]^2)/2)
            }

            #update \b_psi
            b_psi <- rgamma(1, a_4+d1/2, b_4+sum((psi - a_psi)^2)/2)

            #update \gamma
            if(gamma_hyperprior == TRUE){
                par_vec <- as.vector(table(c))
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
            weight_post[[i]] <- S/sum(S)
            psi_post[i, ] <- psi
            b_psi_post[i] <- b_psi
        }
        sim_res <- list(M_post = M_post, c_post = c_post, k_post = k_post, U = U, gamma_post = gamma_post,
                        Beta_post = Beta_post, alpha_post = alpha_post, zeta_post = zeta_post, psi_post = psi_post,
                        lambda_post = lambda_post, eta_post = eta_post, tau_post = tau_post, b_psi_post = b_psi_post,
                        w_post = w_post, weight_post = weight_post, a_lambda = a_lambda, b_lambda = b_lambda,
                        a_eta = a_eta, b_eta = b_eta, W = W, Z = Z, y = y)
    }
    if(prior == "Bessel"){
        #initializations
        U <- rep(1, N)
        M <- M_init
        #this may matter as the chain becomes sticky
        S <- rep(1, M)
        Beta <- matrix(rep(-0.1, d2*M), ncol = d2)
        alpha <- rep(alpha_init, M)
        zeta <- matrix(rep(1, d2*M), ncol = d2)
        lambda <- lambda_init
        eta <- 1
        b_bessel <- b_bessel_fixed
        tau <- rep(1, d2)
        w <- rep(0, d2)
        c <- rep(0, n)
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
        tau_post <- matrix(1, ncol = d2, nrow = N)
        w_post <- matrix(0.5, ncol = d2, nrow = N)
        psi <- rep(-0.1, d1)
        psi_post <- matrix(1, ncol = d1, nrow = N)
        b_psi <- 1
        b_psi_post <- rep(1, N)

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
        #        apply(dnorm(X, zeta[m, ], sqrt(1/eta), log = T), 1, sum)
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
                mem_wgt[ ,m] <- log(S[m]) + group_member(m, Beta, zeta, psi, lambda, alpha, eta)
            }

            for(j in 1:n){
                wgt_shift <- max(mem_wgt[j, ]) #exp normalization trick
                c[j] <- sample(1:M, size = 1, prob = exp(mem_wgt[j, ] - wgt_shift) /
                                   sum(exp(mem_wgt[j, ] - wgt_shift)))
            }
            k <- length(unique(c))
            k_post[i] <- k
            #relabeling
            Beta <- Beta[unique(c), ]
            zeta <- zeta[unique(c), ]
            alpha <- alpha[unique(c)]
            b<-rep(0, n)
            for(j in 1:k){
                b[which(c == unique(c)[j])] <- j
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
            if(k == 1){
                if(k < M){
                    Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = d2)) #fill in the M-k gap
                    zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = d2))
                    alpha <- c(alpha, rep(0, M-k))
                    for(m in 1:k){
                        V_m <- which(c == m)
                        S[m] <- sample.gamma.bel(a_bessel, b_bessel, u, length(V_m))
                        #update zeta and beta in this way cause this seems to give better results, tho both methods should work okay.

                        #update \alpha
                        alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - W[V_m, ] %*% psi - Z[V_m, ] %*% Beta[m, ]))/
                                              (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))

                        for(d in 1:d2){
                            #update \zeta
                            zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(Z[V_m, d]))/
                                                    (eta * length(V_m) + b_zeta),
                                                1/sqrt(eta * length(V_m) + b_zeta))
                            #update \beta
                            C_md <- sqrt(tau/(tau + lambda*sum(Z[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                                exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - W[V_m, ] %*% psi - Z[V_m, -d] %*% Beta[m, -d])*
                                                              Z[V_m, d]))^2)/(2*(tau + lambda*sum(Z[V_m, d]^2))))
                            theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                            Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - W[V_m, ] %*% psi
                                                                                         - Z[V_m, -d] %*%
                                                                                             Beta[m, -d])*Z[V_m, d]) +
                                                                                 mu * tau)/(tau + lambda*sum(Z[V_m, d]^2)),
                                                                         1/sqrt(tau + lambda*sum(Z[V_m, d]^2))))
                        }

                    }
                    for(m in (k+1):M){
                        S[m] <- sample.exp.bel(a_bessel, b_bessel, u)
                        for(d in 1:d2){
                            p_w <- rbinom(1, 1, w[d])
                            Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau)))
                        }
                        zeta[m, ] <- rnorm(d2, a_zeta, 1/sqrt(b_zeta))
                    }
                    alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
                }
                else{
                    Beta <- matrix(Beta, nrow = 1)
                    zeta <- matrix(zeta, nrow = 1)
                    for(m in 1:k){
                        V_m <- which(c == m)
                        S[m] <- sample.gamma.bel(a_bessel, b_bessel, u, length(V_m))

                        #update \alpha
                        alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - W[V_m, ] %*% psi - Z[V_m, ] %*% Beta[m, ]))/
                                              (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))

                        for(d in 1:d2){
                            #update \zeta
                            zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(Z[V_m, d]))/
                                                    (eta * length(V_m) + b_zeta),
                                                1/sqrt(eta * length(V_m) + b_zeta))

                            C_md <- sqrt(tau/(tau + lambda*sum(Z[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                                exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] -  W[V_m, ] %*% psi - Z[V_m, -d] %*% Beta[m, -d])*
                                                              Z[V_m, d]))^2)/(2*(tau + lambda*sum(Z[V_m, d]^2))))
                            theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                            Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] -  W[V_m, ] %*% psi
                                                                                         - Z[V_m, -d] %*%
                                                                                             Beta[m, -d])*Z[V_m, d]) +
                                                                                 mu * tau)/(tau + lambda*sum(Z[V_m, d]^2)),
                                                                         1/sqrt(tau + lambda*sum(Z[V_m, d]^2))))
                        }
                    }
                }
            }
            else if(k < M){
                Beta <- rbind(Beta, matrix(0, nrow = M-k, ncol = d2)) #fill in the M-k gap
                zeta <- rbind(zeta, matrix(0, nrow = M-k, ncol = d2))
                alpha <- c(alpha, rep(0, M-k))
                for(m in 1:k){
                    V_m <- which(c == m)
                    S[m] <- sample.gamma.bel(a_bessel, b_bessel, u, length(V_m))
                    #update zeta and beta in this way cause this seems to give better results, tho both methods should work okay.

                    #update \alpha
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - W[V_m, ] %*% psi - Z[V_m, ] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))

                    for(d in 1:d2){
                        #update \zeta
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(Z[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))
                        #update \beta
                        C_md <- sqrt(tau/(tau + lambda*sum(Z[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] - W[V_m, ] %*% psi - Z[V_m, -d] %*% Beta[m, -d])*
                                                          Z[V_m, d]))^2)/(2*(tau + lambda*sum(Z[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - W[V_m, ] %*% psi
                                                                                     - Z[V_m, -d] %*%
                                                                                         Beta[m, -d])*Z[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(Z[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(Z[V_m, d]^2))))
                    }

                }
                for(m in (k+1):M){
                    S[m] <- sample.exp.bel(a_bessel, b_bessel, u)
                    for(d in 1:d2){
                        p_w <- rbinom(1, 1, w[d])
                        Beta[m, d] <- ifelse(p_w == 1, 0, rnorm(1, mu, 1/sqrt(tau)))
                    }
                    zeta[m, ] <- rnorm(d2, a_zeta, 1/sqrt(b_zeta))
                }
                alpha[(k+1):M] <- rnorm(M-k, a_alpha, 1/sqrt(b_alpha))
            }
            else{
                for(m in 1:k){
                    V_m <- which(c == m)
                    S[m] <- sample.gamma.bel(a_bessel, b_bessel, u, length(V_m))

                    #update \alpha
                    alpha[m] <- rnorm(1, (a_alpha*b_alpha + lambda*sum(y[V_m] - W[V_m, ] %*% psi - Z[V_m, ] %*% Beta[m, ]))/
                                          (length(V_m) * lambda + b_alpha), 1/sqrt(length(V_m) * lambda + b_alpha))

                    for(d in 1:d2){
                        #update \zeta
                        zeta[m, d] <- rnorm(1, (a_zeta * b_zeta + eta * sum(Z[V_m, d]))/
                                                (eta * length(V_m) + b_zeta),
                                            1/sqrt(eta * length(V_m) + b_zeta))

                        C_md <- sqrt(tau/(tau + lambda*sum(Z[V_m, d]^2)))*exp(-0.5*tau*mu^2)*
                            exp(((mu*tau + lambda*sum((y[V_m] - alpha[m] -  W[V_m, ] %*% psi - Z[V_m, -d] %*% Beta[m, -d])*
                                                          Z[V_m, d]))^2)/(2*(tau + lambda*sum(Z[V_m, d]^2))))
                        theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                        Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] -  W[V_m, ] %*% psi
                                                                                     - Z[V_m, -d] %*%
                                                                                         Beta[m, -d])*Z[V_m, d]) +
                                                                             mu * tau)/(tau + lambda*sum(Z[V_m, d]^2)),
                                                                     1/sqrt(tau + lambda*sum(Z[V_m, d]^2))))
                    }
                }
            }

            #update \psi
            var_mtx <- solve(lambda*(t(W) %*% W) + diag(b_psi, d1))
            psi <- mvrnorm(1, var_mtx %*% (colSums(sweep(W, MARGIN = 1, (y - alpha[c] - apply(Beta[c,]*Z, 1, sum)), "*")) *
                                               lambda + rep(a_psi*b_psi,d1)), var_mtx)


            #sample \lambda and \eta
            lambda <- rgamma(1, a_lambda + n/2, b_lambda + sum((y-alpha[c]-W %*%psi-apply(Beta[c,]*Z, 1, sum))^2)/2)
            eta <- rgamma(1, a_eta + n*d2/2, b_eta + sum((Z - zeta[c,])^2)/2)

            #update \tau and w
            if(SS == TRUE){
                for(d in 1:d2){
                    w[d] <- rbeta(1, a_w + sum(Beta[,d] == 0), b_w + sum(Beta[,d] != 0))
                }
            }

            for(d in 1:d2){
                nd_Plus <- sum(Beta[,d] != 0)
                tau[d] <- rgamma(1, a_tau + nd_Plus/2, b_tau +sum(Beta[,d]^2)/2)
            }

            #update \b_psi
            b_psi <- rgamma(1, a_4+d1/2, b_4+sum((psi - a_psi)^2)/2)


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
            psi_post[i, ] <- psi
            b_psi_post[i] <- b_psi
        }
        sim_res <- list(M_post = M_post, c_post = c_post, k_post = k_post, U = U, b_bessel_post = b_bessel_post,
                        Beta_post = Beta_post, alpha_post = alpha_post, zeta_post = zeta_post, psi_post = psi_post,
                        lambda_post = lambda_post, eta_post = eta_post, tau_post = tau_post, b_psi_post = b_psi_post,
                        w_post = w_post, weight_post = weight_post, a_lambda = a_lambda, b_lambda = b_lambda,
                        a_eta = a_eta, b_eta = b_eta, W = W, Z = Z, y = y)
    }
    sim_res
}


