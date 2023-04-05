#measurements functions
get_measure_i <- function(index, k, i, d, X, y){
    prop_i <- matrix(0, nrow = sum(index == i), ncol = k)
    for(j in 1:sum(index==i)){
        not_i <- which(index != i)
        data_not_i <- cbind(y[not_i], X[not_i,])
        point_j <- c(y[which(index == i)[j]], X[which(index == i)[j],])
        ind_d <- not_i[which(distance(rbind(point_j, data_not_i))[1,][-1] < d)]
        #indeces in M_i^d
        for(l in (1:k)){
            prop_i[j, l] <- sum(index[ind_d]==l)/sum(index==l)
        }
    }
    prop <- apply(prop_i, 2, mean)
    prop
}


measure_def <- function(index, d, X, y){
    k <- length(unique(index))
    prop_measure <- matrix(0, nrow=k, ncol=k)
    for(i in 1:k){
        prop_measure[i, ] <- get_measure_i(index, k, i, d, X, y)
    }
    prop_measure
}

#k-means measure, compared with ari
measure_km <- function(index, X, y){
    n <- length(y)
    km_data <- kmeans(cbind(y,X), 3, nstart = 25)
    index_km <- km_data$cluster
    ARI(index, index_km)
}

#calculate the metrics given the posterior samples
post_inf_split <- function(sim_res, scl = 10001:1e5, data){
    index <- data$index
    X <- data$X
    y <- data$y
    k_post <- sim_res$k_post
    M_post <- sim_res$M_post
    post_k <- tabulate(k_post[scl])/length(scl)
    D <- length(sim_res$Beta_post[[1]])/M_post[1]
    d1 <- length(sim_res$psi_post[[1]])/M_post[1]
    beta_true <- data$beta_true
    psi_true <- data$psi_true
    alpha_true <- data$alpha_true
    psi_true <- data$psi_true
    #MSE
    indx_3 <- which(k_post[scl] == 3 & M_post[scl] == 3)+scl[1]-1
    if(length(indx_3) > 1){
        c_mtx <- sim_res$c_post[indx_3, ]
        ls <- label.switching(method = "ECR-ITERATIVE-1", z = c_mtx, K = 3)

        alpha_mtx <- alpha.ordered <- matrix(0, length(indx_3), 3)
        beta_mtx <- beta.ordered <- matrix(0, length(indx_3), 3*D)
        psi_mtx <- psi.ordered <- matrix(0, length(indx_3), 3*d1)
        beta_1 <- beta_2 <- beta_3 <- matrix(0, length(indx_3), D)
        psi_1 <- psi_2 <- psi_3 <- matrix(0, length(indx_3), d1)
        for(i in 1:length(indx_3)){
            beta_mtx[i, ] <- (sim_res$Beta_post)[[indx_3[i]]]
            psi_mtx[i, ] <- (sim_res$psi_post)[[indx_3[i]]]
            alpha_mtx[i, ] <- (sim_res$alpha_post)[[indx_3[i]]]
            beta.ordered[i, ] <- as.vector(matrix(beta_mtx[i,], ncol = D)[ls$permutations[[1]][i,], ])
            psi.ordered[i, ] <- as.vector(matrix(psi_mtx[i,], ncol = d1)[ls$permutations[[1]][i,], ])
            alpha.ordered[i, ] <- alpha_mtx[i, ][ls$permutations[[1]][i,]]
            beta_1[i, ] <- matrix(beta.ordered[i,], nrow = 3)[1, ]
            beta_2[i, ] <- matrix(beta.ordered[i,], nrow = 3)[2, ]
            beta_3[i, ] <- matrix(beta.ordered[i,], nrow = 3)[3, ]
            psi_1[i, ] <- matrix(psi.ordered[i,], nrow = 3)[1, ]
            psi_2[i, ] <- matrix(psi.ordered[i,], nrow = 3)[2, ]
            psi_3[i, ] <- matrix(psi.ordered[i,], nrow = 3)[3, ]
        }
        #to calculate MSE, will need to figure out exactly which one is b1, b2, b3...
        #decide the order by the sum of the first column, which actually compares the
        #first element of each \beta vector
        order_beta <-order(c(sum(beta_1[,1]), sum(beta_2[,1]), sum(beta_3[,1])))
        order_psi <-order(c(sum(psi_1[,1]), sum(psi_2[,1]), sum(psi_3[,1])))
        order_true <- order(beta_true[,1])
        order_true_psi <- order(psi_true[,1])
        beta_list <- list(beta_1, beta_2, beta_3)
        psi_list <- list(psi_1, psi_2, psi_3)

        order_a_true <- order(alpha_true)
        order_a <- order(apply(alpha.ordered, 2, sum))
        mse_a1 <- (alpha.ordered[, order_a[1]] - alpha_true[order_a_true[1]])^2
        mse_a2 <- (alpha.ordered[, order_a[2]] - alpha_true[order_a_true[2]])^2
        mse_a3 <- (alpha.ordered[, order_a[3]] - alpha_true[order_a_true[3]])^2

        mse_b1 <- apply((sweep(beta_list[[order_beta[1]]], 2, beta_true[order_true[1],]))^2, 1, sum)
        mse_b2 <- apply((sweep(beta_list[[order_beta[2]]], 2, beta_true[order_true[2],]))^2, 1, sum)
        mse_b3 <- apply((sweep(beta_list[[order_beta[3]]], 2, beta_true[order_true[3],]))^2, 1, sum)

        mse_p1 <- apply((sweep(psi_list[[order_psi[1]]], 2, psi_true[order_true_psi[1],]))^2, 1, sum)
        mse_p2 <- apply((sweep(psi_list[[order_psi[2]]], 2, psi_true[order_true_psi[2],]))^2, 1, sum)
        mse_p3 <- apply((sweep(psi_list[[order_psi[3]]], 2, psi_true[order_true_psi[3],]))^2, 1, sum)


        fdr1 <- fdr2 <- fdr3 <- error1 <- error2 <- error3 <- rep(0, length(indx_3))
        #FDR rate
        for(i in 1:length(indx_3)){
            fdr1[i] <- sum(beta_true[order_true[1],][which(beta_list[[order_beta[1]]][i, ] == 0)] != 0)/ sum(beta_list[[order_beta[1]]][i, ] == 0)
            fdr2[i] <- sum(beta_true[order_true[2],][which(beta_list[[order_beta[2]]][i, ] == 0)] != 0)/ sum(beta_list[[order_beta[2]]][i, ] == 0)
            fdr3[i] <- sum(beta_true[order_true[3],][which(beta_list[[order_beta[3]]][i, ] == 0)] != 0)/ sum(beta_list[[order_beta[3]]][i, ] == 0)
            error1[i] <- sum(beta_list[[order_beta[1]]][i, which(beta_true[order_true[1],] == 0)] != 0)/ sum(beta_true[order_true[1],] == 0)
            error2[i] <- sum(beta_list[[order_beta[2]]][i, which(beta_true[order_true[2],] == 0)] != 0)/ sum(beta_true[order_true[2],] == 0)
            error3[i] <- sum(beta_list[[order_beta[3]]][i, which(beta_true[order_true[3],] == 0)] != 0)/ sum(beta_true[order_true[3],] == 0)
        }

    }else{
        mse_a1 <- mse_a2 <- mse_a3 <- mse_b1 <- mse_b2 <- mse_b3 <- mse_p1 <- mse_p2 <- mse_p3 <- 0
        fdr1 <- fdr2 <- fdr3 <- error1 <- error2 <- error3 <- 1
    }

    c_post <- sim_res$c_post
    ari.res <- apply(c_post[scl, ], 1, ARI, index)
    ############
    #"mcclust" package to compute point estimate
    ############
    c_burn <- c_post[scl,]
    psm <- comp.psm(c_burn)
    est_bin <- minbinder.ext(psm, c_burn, method = "all", include.greedy = T)
    est_vi <- minVI(psm, c_burn, method = "all", include.greedy = T)

    c_bin <- est_bin$cl[1,]
    c_vi <- est_vi$cl[1,]

    ari_point_bin <- ARI(index, c_bin)
    ari_point_vi <- ARI(index, c_vi)

    #convergence diagnostics
    auto_corr_k <- autocorr(mcmc(k_post[scl]))[,,1]
    geweke_k <- geweke.diag(mcmc(k_post[scl]))$z
    #measurements
    true_size <- tabulate(index)
    #measure_nei <- measure_def(index, 20, X, y)
    #measure_nei <- c(measure_nei[1,2],measure_nei[1,3], measure_nei[2,3])
    measure_km_ari <- measure_km(index, X, y)
    ##what we need:
    #true_size, measure_nei, measure_km_ari
    #post_k, mse_b1, mse_b2, mse_b3, ARI
    #point estimate: c_bin, c_vi, ari_point_bin, ari_point_vi
    results <- list(true_size=unlist(true_size), measure_km_ari=measure_km_ari,
                    auto_corr_k = as.vector(auto_corr_k), geweke_k = geweke_k,
                    post_k=post_k, mse_a1=quantile(mse_a1), mse_a2=quantile(mse_a2),
                    mse_a3=quantile(mse_a3), mse_b1=quantile(mse_b1), mse_b2=quantile(mse_b2),
                    mse_b3=quantile(mse_b3), mse_p1 = quantile(mse_p1),
                    mse_p2 = quantile(mse_p2), mse_p3 = quantile(mse_p3),
                    ARI=quantile(ari.res, probs = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1), na.rm = T),
                    c_bin=tabulate(c_bin), c_vi=tabulate(c_vi), ari_point_bin=ari_point_bin,
                    ari_point_vi=ari_point_vi,
                    FS = c(mean(fdr1, na.rm = T), mean(fdr2, na.rm = T), mean(fdr3, na.rm = T)),
                    MS = c(mean(error1), mean(error2), mean(error3)))
    results
}
