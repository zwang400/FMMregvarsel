#generate simulation data
data_gen_func <- function(n = 500, alpha_true =  c(0.1, -0.6, 0.5),
                         beta_true = rbind(c(0, 0, -0.5, 0, 0.5, 0),
                                           c(-0.7, 0, 0.4, 0, 0, 0),
                                           c(0.6, 0, 0, 0, -0.4, 0)),
                         lambda = c(2, 2, 2), zeta_sep = 1,
                         eta = 1, sample_prob = c(1, 1, 1),
                         cor_mtx = NULL, rho = rep(0, 3)){
    #correlation are blocked-wise, X1-X2, X3-X4, X5-X6 are correlated
    #correlation are the 3 elements in rho
    true_M <- nrow(beta_true)
    zeta <- seq(0, true_M-1) * zeta_sep + 2
    D <- ncol(beta_true)
    sd <- sqrt(1/lambda)
    index <- sample(1:true_M, n, prob = sample_prob, replace = T)
    X <- matrix(0, nrow = n, ncol = D)
    y <- rep(0, n)
    if(is.null(cor_mtx) && all(rho == 0)){
        Beta_init <- matrix(0, nrow = n, ncol = D)
        for(i in 1:n){
            X[i, ] <- rnorm(D, zeta[index[i]], 1/sqrt(eta))
            Beta_init[i, ] <- beta_true[index[i], ]
            y[i] <- alpha_true[index[i]] + rnorm(1, X[i, ] %*% Beta_init[i, ], sd[index[i]])
        }
    }else{
        if(is.null(cor_mtx)){
            cor_mtx <- bdiag(matrix(c(1/eta, rho[1], rho[1], 1/eta), 2, 2),
                             matrix(c(1/eta, rho[2], rho[2], 1/eta), 2, 2),
                             matrix(c(1/eta, rho[3], rho[3], 1/eta), 2, 2))
        }
        for(m in 1:true_M){
            X[which(index==m), ] <- mvrnorm(n = sum(index==m), mu = rep(zeta[m], D), Sigma = cor_mtx)
            y[which(index==m)] <- alpha_true[m] + rnorm(n = sum(index==m), X[which(index==m), ] %*% beta_true[m, ], sd[m])
        }
    }
    list(X = X, y = y, index= index, alpha_true = alpha_true, beta_true = beta_true)
}



data_gen_split <- function(n = 500, alpha_true = c(0.1, -0.6, 0.5),
                          beta_true = rbind(c(0, 0, -0.5, 0, 0.5, 0),
                                            c(-0.7, 0, 0.4, 0, 0, 0),
                                            c(0.6, 0, 0, 0, -0.4, 0)),
                          psi_true = rep(c(0.3,0,8),3),
                          W_mean = 0, lambda = c(2, 2, 2), zeta_sep = 1,
                          eta = 1, sample_prob = c(1, 1, 1), rho = 0.5){
    zeta <- c(0, zeta_sep, 2*zeta_sep)
    d1 <- length(psi_true)
    d2 <- ncol(beta_true)
    sd <- sqrt(1/lambda)
    index <- sample(1:3, n, prob = sample_prob, replace = T)
    W <- mvrnorm(n, rep(W_mean, d1), matrix(rep(rho, d1*d1), ncol = d1) + diag(1-rho, d1))
    Z <- t(mvrnorm(d2, zeta[index], diag(1/eta, n)))
    y <- rep(0, n)

    for(i in 1:n){
        y[i] <- alpha_true[index[i]] + rnorm(1, Z[i, ] %*% beta_true[index[i], ] +
                                                 W[i, ] %*% psi_true, sd[index[i]])
    }

    list(Z = Z, y = y, W = W, index= index, alpha_true = alpha_true, beta_true = beta_true, psi_true = psi_true)
}

data_gen_split_new <- function(n = 500, alpha_true = c(0.1, -0.6, 0.5),
                           beta_true = rbind(c(0, 0, -0.5, 0, 0.5, 0),
                                             c(-0.7, 0, 0.4, 0, 0, 0),
                                             c(0.6, 0, 0, 0, -0.4, 0)),
                           psi_true = rbind(c(0.1, 0.5, -0.6),
                                            c(-0.7, -0.2, 0.5),
                                            c(0.9, -0.4, 0.1)),
                           W_mean = 0, lambda = c(2, 2, 2), zeta_sep = 1,
                           eta = 1, sample_prob = c(1, 1, 1), rho = 0.5){
    zeta <- c(0, zeta_sep, 2*zeta_sep)
    d1 <- ncol(psi_true)
    d2 <- ncol(beta_true)
    sd <- sqrt(1/lambda)
    index <- sample(1:3, n, prob = sample_prob, replace = T)
    W <- mvrnorm(n, rep(W_mean, d1), matrix(rep(rho, d1*d1), ncol = d1) + diag(1-rho, d1))
    Z <- t(mvrnorm(d2, zeta[index], diag(1/eta, n)))
    y <- rep(0, n)

    for(i in 1:n){
        y[i] <- alpha_true[index[i]] + rnorm(1, Z[i, ] %*% beta_true[index[i], ] +
                                                 W[i, ] %*% psi_true[index[i], ], sd[index[i]])
    }

    list(Z = Z, y = y, W = W, index= index, alpha_true = alpha_true, beta_true = beta_true, psi_true = psi_true)
}

data_gen_gamma <- function(n = 500, alpha_true =  c(0.1, -0.6, 0.5),
                           beta_true = rbind(c(0, 0, -0.5, 0, 0.5, 0),
                                             c(-0.7, 0, 0.4, 0, 0, 0),
                                             c(0.6, 0, 0, 0, -0.4, 0)),
                           lambda = c(2, 2, 2), zeta_sep = 1,
                           sample_prob = c(1, 1, 1)){
    true_M <- nrow(beta_true)
    eta <- seq(0, true_M-1) * zeta_sep + 2
    zeta <- eta^2
    D <- ncol(beta_true)
    sd <- sqrt(1/lambda)
    index <- sample(1:true_M, n, prob = sample_prob, replace = T)
    X <- matrix(0, nrow = n, ncol = D)
    y <- rep(0, n)
    Beta_init <- matrix(0, nrow = n, ncol = D)
    for(i in 1:n){
        X[i, ] <- rgamma(D, zeta[index[i]], eta[index[i]])
        Beta_init[i, ] <- beta_true[index[i], ]
        y[i] <- alpha_true[index[i]] + rnorm(1, X[i, ] %*% Beta_init[i, ], sd[index[i]])
    }
    list(X = X, y = y, index= index, alpha_true = alpha_true, beta_true = beta_true)
}

data_gen_t <- function(n = 500, alpha_true =  c(0.1, -0.6, 0.5),
                           beta_true = rbind(c(0, 0, -0.5, 0, 0.5, 0),
                                             c(-0.7, 0, 0.4, 0, 0, 0),
                                             c(0.6, 0, 0, 0, -0.4, 0)),
                           lambda = c(2, 2, 2), zeta_sep = 1,
                           sample_prob = c(1, 1, 1)){
    true_M <- nrow(beta_true)
    zeta <- seq(0, true_M-1) * zeta_sep + 2
    D <- ncol(beta_true)
    sd <- sqrt(1/lambda)
    index <- sample(1:true_M, n, prob = sample_prob, replace = T)
    X <- matrix(0, nrow = n, ncol = D)
    y <- rep(0, n)
    Beta_init <- matrix(0, nrow = n, ncol = D)
    for(i in 1:n){
        X[i, ] <- rt(D, true_M)/sqrt(true_M/(true_M-2)) + zeta[index[i]]
        Beta_init[i, ] <- beta_true[index[i], ]
        y[i] <- alpha_true[index[i]] + rnorm(1, X[i, ] %*% Beta_init[i, ], sd[index[i]])
    }
    list(X = X, y = y, index= index, alpha_true = alpha_true, beta_true = beta_true)
}








