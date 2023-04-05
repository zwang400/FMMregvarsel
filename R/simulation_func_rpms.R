#simulation function for rpms
simulation_func_rpms <- function(X, y, SS = TRUE, N = 1e5, a_w = 1, b_w = 1, mu = 0, a_tau = 1, b_tau = 1,
                     a_zeta = 0, b_zeta = 1, a_lambda = 0.01, b_lambda = 0.005, a_alpha = 0,
                     b_alpha = 0.01, a_eta = 5, b_eta = 2, a_adp = 1, b_adp = 1,
                     k_init = 6, lambda_init = 2, alpha_dp_init = 1, alpha_init = 0, m_aux = 10){
    n <- length(y)
    D <- ncol(X)
    k <- k_init
    Beta <- matrix(rep(-0.1, D*k), ncol = D)
    alpha_dp <- alpha_dp_init
    alpha <- rep(alpha_init, k)
    zeta <- matrix(rep(1, D*k), ncol = D)
    lambda <- lambda_init
    eta <- 1
    tau <- rep(1, D)
    w <- rep(0, D)
    c <- sample(1:k, n, replace = T) #use a random start
    k_post <- rep(0, N)
    c_post <- matrix(0, ncol = n, nrow = N)
    Beta_post <- list(0)
    alpha_post <- list(0)
    zeta_post <- list(0)
    lambda_post <- rep(0, N)
    eta_post <- rep(0, N)
    tau_post <- matrix(1, ncol = D, nrow = N)
    w_post <- matrix(0.5, ncol = D, nrow = N)
    alpha_dp_post <- rep(0, N)

    group_member <- function(i, m, Beta, zeta, lambda, alpha, eta){
        llkhd_x <- 0
        for(j in 1:D){
            llkhd_x <- llkhd_x + dnorm(X[i,j], zeta[m, j], sqrt(1/eta), log = T)
        }
        llkhd <- dnorm(y[i], alpha[m] + X[i,]%*%Beta[m,], sqrt(1/lambda), log = T) + llkhd_x
        llkhd
    }


    for(i in 1:N){
        #updating group membership c
        for(j in 1:n){
            c_minus <- c[-j]
            cj <- c[j]
            k_minus <- length(unique(c_minus))
            h <- k_minus + m_aux
            if(any(c_minus == c[j])){
                #extend Beta, zeta and alpha
                Beta_ext <- matrix(0, ncol = D, nrow = m_aux)
                for(mm in 1:m_aux){
                    pp_w <- rbinom(D, 1, w)
                    Beta_ext[mm,] <- ifelse(pp_w == 1, rep(0, D), rnorm(D, mu, 1/sqrt(tau)))
                }
                zeta_ext <- matrix(rnorm(D*m_aux, a_zeta, 1/sqrt(b_zeta)), ncol = D)
                alpha_ext <- rnorm(m_aux, a_alpha, 1/sqrt(b_alpha))
                Beta <- rbind(Beta, Beta_ext)
                zeta <- rbind(zeta, zeta_ext)
                alpha <- c(alpha, alpha_ext)
            }else{
                c[j] <- k_minus + 1
                Beta_ext <- matrix(0, ncol = D, nrow = (m_aux-1))
                for(mm in 1:(m_aux-1)){
                    pp_w <- rbinom(D, 1, w)
                    Beta_ext[mm,] <- ifelse(pp_w == 1, rep(0, D), rnorm(D, mu, 1/sqrt(tau)))
                }
                zeta_ext <- matrix(rnorm(D*(m_aux-1), a_zeta, 1/sqrt(b_zeta)), ncol = D)
                alpha_ext <- rnorm(m_aux-1, a_alpha, 1/sqrt(b_alpha))

                #need to relabel first
                b<-rep(0, n-1)
                for(kk in 1:k_minus){
                    b[which(c_minus == unique(c_minus)[kk])] <- kk
                }
                c_minus <- b
                c[-j] <- c_minus

                Beta_ord <- rbind(Beta[-cj, ], Beta[cj, ])
                zeta_ord <- rbind(zeta[-cj, ], zeta[cj, ])
                alpha_ord <- c(alpha[-cj], alpha[cj])
                Beta <- rbind(Beta_ord, Beta_ext)
                zeta <- rbind(zeta_ord, zeta_ext)
                alpha <- c(alpha_ord, alpha_ext)
            }

            #size for existing clusters
            n_size <- rep(0, k_minus)
            for(cc in 1:k_minus){
                n_size[cc] <- sum(c_minus == cc)
            }

            #calculating joint mixture density
            mix_den <- rep(0, h)
            for(hh in 1:h){
                mix_den[hh] <- group_member(j, hh, Beta, zeta, lambda, alpha, eta)
            }
            mem_wgt <- rep(0, h)
            mem_wgt[1:k_minus] <- log(n_size) + mix_den[1:k_minus]
            mem_wgt[(k_minus+1):h] <- log(alpha_dp) - log(m_aux) + mix_den[(k_minus+1):h]
            mem_shift <- max(mem_wgt)
            c[j] <- sample(1:h, size = 1, prob =
                               exp(mem_wgt - mem_shift)/ sum(exp(mem_wgt - mem_shift)))
            if(c[j] > k_minus){
                c[j] <- k_minus + 1
            }
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

        #updating \alpha_dp
        u <- rbeta(1, alpha_dp+1, n)
        p_adp <- rbinom(1, 1, (a_adp + k - 1)/ (a_adp + k - 1 + n*(b_adp - log(u))))
        alpha_dp <- ifelse(p_adp == 1, rgamma(1, a_adp+k, b_adp - log(u)),
                           rgamma(1, a_adp+k-1, b_adp - log(u)))
        alpha_dp_post[i] <- alpha_dp

        #updating clustering parameters
        if(k == 1){
            Beta <- matrix(Beta, nrow = 1)
            zeta <- matrix(zeta, nrow = 1)
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
                    C_md <- sqrt(tau[d]/(tau[d] + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau[d]*mu^2)*
                        exp(((mu*tau[d] + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d, drop=F] %*% Beta[m, -d])*
                                                         X[V_m, d]))^2)/(2*(tau[d] + lambda*sum(X[V_m, d]^2))))
                    theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                    Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d, drop=F] %*%
                                                                                     Beta[m, -d])*X[V_m, d]) +
                                                                         mu * tau[d])/(tau[d] + lambda*sum(X[V_m, d]^2)),
                                                                 1/sqrt(tau[d] + lambda*sum(X[V_m, d]^2))))
                }
            }
        }
        else{
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
                    C_md <- sqrt(tau[d]/(tau[d] + lambda*sum(X[V_m, d]^2)))*exp(-0.5*tau[d]*mu^2)*
                        exp(((mu*tau[d] + lambda*sum((y[V_m] - alpha[m] - X[V_m, -d, drop=F] %*% Beta[m, -d])*
                                                         X[V_m, d]))^2)/(2*(tau[d] + lambda*sum(X[V_m, d]^2))))
                    theta_md <- rbinom(1, 1, w[d]/(w[d] + (1 - w[d])*C_md))
                    Beta[m, d] <- ifelse(theta_md == 1, 0, rnorm(1, (lambda*sum((y[V_m] - alpha[m] - X[V_m, -d, drop=F] %*%
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

        Beta_post[[i]] <- as.vector(Beta)
        alpha_post[[i]] <- alpha
        zeta_post[[i]] <- as.vector(zeta)
        lambda_post[i] <- lambda
        eta_post[i] <- eta
        tau_post[i,] <- tau
        w_post[i,] <- w
    }
    sim_res <- list(c_post = c_post, k_post = k_post, alpha_dp_post = alpha_dp_post,
                    Beta_post = Beta_post, alpha_post = alpha_post, zeta_post = zeta_post,
                    lambda_post = lambda_post, eta_post = eta_post, tau_post = tau_post,
                    w_post = w_post, a_lambda = a_lambda, b_lambda = b_lambda,
                    a_eta = a_eta, b_eta = b_eta, X = X, y = y)
    sim_res
}
