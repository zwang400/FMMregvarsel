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
post_inf <- function(sim_res, scl = 10001:1e5, data){
    index <- data$index
    X <- data$X
    y <- data$y
    k_post <- sim_res$k_post
    M_post <- sim_res$M_post
    post_k <- tabulate(k_post[scl])/length(scl)
    D <- length(sim_res$Beta_post[[1]])/M_post[1]
    beta_true <- data$beta_true
    alpha_true = data$alpha_true
    #MSE
    indx_3 <- which(k_post[scl] == 3 & M_post[scl] == 3)+scl[1]-1
    if(length(indx_3) > 1){
        c_mtx <- sim_res$c_post[indx_3, ]
        ls <- label.switching(method = "ECR-ITERATIVE-1", z = c_mtx, K = 3)

        alpha_mtx <- alpha.ordered <- matrix(0, length(indx_3), 3)
        beta_mtx <- beta.ordered <- matrix(0, length(indx_3), 3*D)
        beta_1 <- beta_2 <- beta_3<- matrix(0, length(indx_3), D)
        for(i in 1:length(indx_3)){
            beta_mtx[i, ] <- (sim_res$Beta_post)[[indx_3[i]]]
            alpha_mtx[i, ] <- (sim_res$alpha_post)[[indx_3[i]]]
            beta.ordered[i, ] <- as.vector(matrix(beta_mtx[i,], ncol = D)[ls$permutations[[1]][i,], ])
            alpha.ordered[i, ] <- alpha_mtx[i, ][ls$permutations[[1]][i,]]
            beta_1[i, ] <- matrix(beta.ordered[i,], nrow = 3)[1, ]
            beta_2[i, ] <- matrix(beta.ordered[i,], nrow = 3)[2, ]
            beta_3[i, ] <- matrix(beta.ordered[i,], nrow = 3)[3, ]
        }
        #to calculate MSE, will need to figure out exactly which one is b1, b2, b3...
        #decide the order by the sum of the first column, which actually compares the
        #first element of each \beta vector
        order_beta <-order(c(sum(beta_1[,1]), sum(beta_2[,1]), sum(beta_3[,1])))
        order_true <- order(beta_true[,1])
        beta_list <- list(beta_1, beta_2, beta_3)

        order_a_true <- order(alpha_true)
        order_a <- order(apply(alpha.ordered, 2, sum))
        mse_a1 <- (alpha.ordered[, order_a[1]] - alpha_true[order_a_true[1]])^2
        mse_a2 <- (alpha.ordered[, order_a[2]] - alpha_true[order_a_true[2]])^2
        mse_a3 <- (alpha.ordered[, order_a[3]] - alpha_true[order_a_true[3]])^2

        mse_b1 <- apply((sweep(beta_list[[order_beta[1]]], 2, beta_true[order_true[1],]))^2, 1, sum)
        mse_b2 <- apply((sweep(beta_list[[order_beta[2]]], 2, beta_true[order_true[2],]))^2, 1, sum)
        mse_b3 <- apply((sweep(beta_list[[order_beta[3]]], 2, beta_true[order_true[3],]))^2, 1, sum)

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
        mse_a1 <- mse_a2 <- mse_a3 <- mse_b1 <- mse_b2 <- mse_b3 <- 0
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
                    mse_b3=quantile(mse_b3), ARI=quantile(ari.res, probs = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1), na.rm = T),
                    c_bin=tabulate(c_bin), c_vi=tabulate(c_vi), ari_point_bin=ari_point_bin,
                    ari_point_vi=ari_point_vi,
                    FS = c(mean(fdr1, na.rm = T), mean(fdr2, na.rm = T), mean(fdr3, na.rm = T)),
                    MS = c(mean(error1), mean(error2), mean(error3)))
    results
}

#for RPMS
post_inf_rpms <- function(sim_res, scl = 10001:1e5, data){
    index <- data$index
    X <- data$X
    y <- data$y
    k_post <- sim_res$k_post
    post_k <- tabulate(k_post[scl])/length(scl)
    D <- ncol(X)
    beta_true <- data$beta_true
    alpha_true = data$alpha_true
    #MSE
    indx_3 <- which(k_post[scl] == 3)+scl[1]-1
    if(length(indx_3) > 1){
        c_mtx <- sim_res$c_post[indx_3, ]
        ls <- label.switching(method = "ECR-ITERATIVE-1", z = c_mtx, K = 3)

        alpha_mtx <- alpha.ordered <- matrix(0, length(indx_3), 3)
        beta_mtx <- beta.ordered <- matrix(0, length(indx_3), 3*D)
        beta_1 <- beta_2 <- beta_3<- matrix(0, length(indx_3), D)
        for(i in 1:length(indx_3)){
            beta_mtx[i, ] <- (sim_res$Beta_post)[[indx_3[i]]]
            alpha_mtx[i, ] <- (sim_res$alpha_post)[[indx_3[i]]]
            beta.ordered[i, ] <- as.vector(matrix(beta_mtx[i,], ncol = D)[ls$permutations[[1]][i,], ])
            alpha.ordered[i, ] <- alpha_mtx[i, ][ls$permutations[[1]][i,]]
            beta_1[i, ] <- matrix(beta.ordered[i,], nrow = 3)[1, ]
            beta_2[i, ] <- matrix(beta.ordered[i,], nrow = 3)[2, ]
            beta_3[i, ] <- matrix(beta.ordered[i,], nrow = 3)[3, ]
        }
        #to calculate MSE, will need to figure out exactly which one is b1, b2, b3...
        #decide the order by the sum of the first column, which actually compares the
        #first element of each \beta vector
        order_beta <-order(c(sum(beta_1[,1]), sum(beta_2[,1]), sum(beta_3[,1])))
        order_true <- order(beta_true[,1])
        beta_list <- list(beta_1, beta_2, beta_3)

        order_a_true <- order(alpha_true)
        order_a <- order(apply(alpha.ordered, 2, sum))
        mse_a1 <- (alpha.ordered[, order_a[1]] - alpha_true[order_a_true[1]])^2
        mse_a2 <- (alpha.ordered[, order_a[2]] - alpha_true[order_a_true[2]])^2
        mse_a3 <- (alpha.ordered[, order_a[3]] - alpha_true[order_a_true[3]])^2

        mse_b1 <- apply((sweep(beta_list[[order_beta[1]]], 2, beta_true[order_true[1],]))^2, 1, sum)
        mse_b2 <- apply((sweep(beta_list[[order_beta[2]]], 2, beta_true[order_true[2],]))^2, 1, sum)
        mse_b3 <- apply((sweep(beta_list[[order_beta[3]]], 2, beta_true[order_true[3],]))^2, 1, sum)

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
        mse_a1 <- mse_a2 <- mse_a3 <- mse_b1 <- mse_b2 <- mse_b3 <- 0
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
                    mse_b3=quantile(mse_b3), ARI=quantile(ari.res, probs = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1), na.rm = T),
                    c_bin=tabulate(c_bin), c_vi=tabulate(c_vi), ari_point_bin=ari_point_bin,
                    ari_point_vi=ari_point_vi,
                    FS = c(mean(fdr1, na.rm = T), mean(fdr2, na.rm = T), mean(fdr3, na.rm = T)),
                    MS = c(mean(error1), mean(error2), mean(error3)))
    results
}

#for BNPmix
post_inf_bnpmix <- function(sim_res, data){
    index <- data$index
    X <- data$X
    y <- data$y
    k_post <- apply(sim_res$clust, 1, function(x) length(unique(x)))
    post_k <- tabulate(k_post)/length(k_post)
    D <- length(sim_res$beta[[1]])/k_post[1]
    beta_true <- data$beta_true
    alpha_true = data$alpha_true
    indx_3 <- which(k_post == 3)
    if(length(indx_3) > 1){
        c_mtx <- sim_res$clust[indx_3, ]+1
        ls <- label.switching(method = "ECR-ITERATIVE-1", z = c_mtx, K = 3)

        beta_mtx <- beta.ordered <- matrix(0, length(indx_3), 3*D)
        beta_1 <- beta_2 <- beta_3<- matrix(0, length(indx_3), D)
        for(i in 1:length(indx_3)){
            beta_mtx[i, ] <- as.vector((sim_res$beta)[[indx_3[i]]])
            beta.ordered[i, ] <- as.vector(matrix(beta_mtx[i,], ncol = D)[ls$permutations[[1]][i,], ])
            beta_1[i, ] <- matrix(beta.ordered[i,], nrow = 3)[1, ]
            beta_2[i, ] <- matrix(beta.ordered[i,], nrow = 3)[2, ]
            beta_3[i, ] <- matrix(beta.ordered[i,], nrow = 3)[3, ]
        }
        #to calculate MSE, will need to figure out exactly which one is b1, b2, b3...
        #decide the order by the sum of the first column, which actually compares the
        #first element of each \beta vector
        order_beta <-order(c(sum(beta_1[,2]), sum(beta_2[,2]), sum(beta_3[,2])))
        order_true <- order(beta_true[,1])
        beta_list <- list(beta_1[,-1], beta_2[,-1], beta_3[,-1])

        alpha.ordered <- cbind(beta_1[,1], beta_2[,1], beta_3[,1])
        order_a_true <- order(alpha_true)
        order_a <- order(apply(alpha.ordered, 2, sum))
        mse_a1 <- (alpha.ordered[, order_a[1]] - alpha_true[order_a_true[1]])^2
        mse_a2 <- (alpha.ordered[, order_a[2]] - alpha_true[order_a_true[2]])^2
        mse_a3 <- (alpha.ordered[, order_a[3]] - alpha_true[order_a_true[3]])^2

        mse_b1 <- apply((sweep(beta_list[[order_beta[1]]], 2, beta_true[order_true[1],]))^2, 1, sum)
        mse_b2 <- apply((sweep(beta_list[[order_beta[2]]], 2, beta_true[order_true[2],]))^2, 1, sum)
        mse_b3 <- apply((sweep(beta_list[[order_beta[3]]], 2, beta_true[order_true[3],]))^2, 1, sum)
    }else{
        mse_a1 <- mse_a2 <- mse_a3 <- mse_b1 <- mse_b2 <- mse_b3 <- 0
    }

    ARI_post <- apply(sim_res$clust, 1, ARI, index)

    results <- list(post_k=post_k, ARI_bnpmix = quantile(ARI_post,
                                                         probs = c(0, 0.05, 0.1,
                                                                   0.15, 0.2, 0.25,
                                                                   0.5, 0.75, 1),
                                                         na.rm = T),
                    mse_a1=quantile(mse_a1), mse_a2=quantile(mse_a2),
                    mse_a3=quantile(mse_a3), mse_b1=quantile(mse_b1), mse_b2=quantile(mse_b2),
                    mse_b3=quantile(mse_b3))
    results
}

#for mclust
mclust_vs <- function(data){
    data.index <- data$index
    beta_true <- data$beta_true
    alpha_true <- data$alpha_true
    data1 <- as.data.frame(cbind(data$X, data$y))

    #find best clustering with no fixed k
    mod1 <- Mclust(data1)
    clust_no_g <- tabulate(mod1$classification)
    ari_no_g <- ARI(data.index, mod1$classification)

    mod2 <- clustvarsel(data1, direction = "backward", verbose = F)
    data2 <- data1[, mod2$subset, drop = F]
    mod22 <- Mclust(data2)
    cluster_no_g_vs <- tabulate(mod22$classification)
    ari_no_g_vs <- ARI(data.index, mod22$classification)

    #fix k = 3, and run clustvarsel separatly for 3 clusters
    mode.mcl <- Mclust(data1, G=3)
    clust_g3 <- tabulate(mode.mcl$classification)
    ari_g3 <- ARI(data.index, mode.mcl$classification)

    #best clustering after VS
    mod.vs <- clustvarsel(data1, G = 3, direction = "backward",  verbose = F)
    data1.v <- data1[, mod.vs$subset, drop = F]
    mod.ff <- Mclust(data1.v, G=3)
    clust_g3_vs <- tabulate(mod.ff$classification)
    ari_g3_vs <- ARI(data.index, mod.ff$classification)


    #doing VS for 3 clusters respectively
    ind1 <- which(mod.ff$classification == 1)
    ind2 <- which(mod.ff$classification == 2)
    ind3 <- which(mod.ff$classification == 3)

    #regression for 3 clusters respectively, based on 1 single VS
    data2.v <- data1[, mod.vs$subset[mod.vs$subset != 7], drop = F]
    data2.v$y <- data$y

    alpha1 <- alpha2 <- alpha3 <- 0
    beta1 <- beta2 <- beta3 <- rep(0, ncol(data1)-1)

    reg1.v <- lm(y~., data = data2.v[ind1,])
    alpha1 <- reg1.v$coefficients[1]
    beta1[mod.vs$subset[mod.vs$subset != 7]] <- reg1.v$coefficients[-1]

    reg2.v <- lm(y~., data = data2.v[ind2,])
    alpha2 <- reg2.v$coefficients[1]
    beta2[mod.vs$subset[mod.vs$subset != 7]] <- reg2.v$coefficients[-1]

    reg3.v <- lm(y~., data = data2.v[ind3,])
    alpha3 <- reg3.v$coefficients[1]
    beta3[mod.vs$subset[mod.vs$subset != 7]] <- reg3.v$coefficients[-1]

    #alpha vector
    alpha_seq <- as.numeric(c(alpha1, alpha2, alpha3))

    order_b <- order(c(beta1[1], beta2[1], beta3[1]))
    beta.mtx <- rbind(beta1, beta2, beta3)
    order_true <- order(beta_true[,1])

    #se of alpha
    se_a1.1 <- (alpha_seq[order_b[1]] - alpha_true[order_true[1]])^2
    se_a2.1 <- (alpha_seq[order_b[2]] - alpha_true[order_true[2]])^2
    se_a3.1 <- (alpha_seq[order_b[3]] - alpha_true[order_true[3]])^2


    #se of beta
    se_b1 <- sum((beta.mtx[order_b[1],] - beta_true[order_true[1],])^2)
    se_b2 <- sum((beta.mtx[order_b[2],] - beta_true[order_true[2],])^2)
    se_b3 <- sum((beta.mtx[order_b[3],] - beta_true[order_true[3],])^2)

    #fdr
    fdr1 <- sum(beta_true[order_true[1],][which(beta.mtx[order_b[1],] == 0)] != 0)/ sum(beta.mtx[order_b[1],] == 0)
    fdr2 <- sum(beta_true[order_true[2],][which(beta.mtx[order_b[2],] == 0)] != 0)/ sum(beta.mtx[order_b[2],] == 0)
    fdr3 <- sum(beta_true[order_true[3],][which(beta.mtx[order_b[3],] == 0)] != 0)/ sum(beta.mtx[order_b[3],] == 0)

    #error
    err1 <- sum(beta.mtx[order_b[1],][which(beta_true[order_true[1],] == 0)] != 0 )/ sum(beta_true[order_true[1],] == 0)
    err2 <- sum(beta.mtx[order_b[2],][which(beta_true[order_true[2],] == 0)] != 0 )/ sum(beta_true[order_true[2],] == 0)
    err3 <- sum(beta.mtx[order_b[3],][which(beta_true[order_true[3],] == 0)] != 0 )/ sum(beta_true[order_true[3],] == 0)

    result_1 <- list(se_b1.1 = se_b1, se_b2.1 = se_b2,
                     se_b3.1 = se_b3, se_a1.1 = se_a1.1, se_a2.1 = se_a2.1,
                     se_a3.1 = se_a3.1, FS1.1 = fdr1, FS2.1 = fdr2, FS3.1 = fdr3,
                     MS1.1 = err1, MS2.1 = err2, MS3.1 = err3)

    #regression for 3 clusters respectively, based on 3 VS respectively
    ind1 <- which(mode.mcl$classification == 1)
    ind2 <- which(mode.mcl$classification == 2)
    ind3 <- which(mode.mcl$classification == 3)

    data.f1 <- data1[ind1,]
    data.f2 <- data1[ind2,]
    data.f3 <- data1[ind3,]

    if(nrow(data.f1) <= ncol(data.f1)){
        mod.v1 <- clustvarsel(data.f1, G = 3, direction = "backward",  verbose = F, emModels2 = "EII")
    } else{
        mod.v1 <- clustvarsel(data.f1, G = 3, direction = "backward",  verbose = F)
    }

    if(nrow(data.f2) <= ncol(data.f2)){
        mod.v2 <- clustvarsel(data.f2, G = 3, direction = "backward",  verbose = F, emModels2 = "EII")
    } else{
        mod.v2 <- clustvarsel(data.f2, G = 3, direction = "backward",  verbose = F)
    }

    if(nrow(data.f3) <= ncol(data.f3)){
        mod.v3 <- clustvarsel(data.f3, G = 3, direction = "backward",  verbose = F, emModels2 = "EII")
    } else{
        mod.v3 <- clustvarsel(data.f3, G = 3, direction = "backward",  verbose = F)
    }

    alpha1 <- alpha2 <- alpha3 <- 0
    beta1 <- beta2 <- beta3 <- rep(0, ncol(data1)-1)

    data.v1 <- data.f1[, mod.v1$subset[mod.v1$subset != 7], drop = F]
    data.v1$y <- data.f1[,7]

    reg.v1 <- lm(y~., data = data.v1)
    alpha1 <- reg.v1$coefficients[1]
    beta1[mod.v1$subset[mod.v1$subset != 7]] <- reg.v1$coefficients[-1]

    data.v2 <- data.f2[, mod.v2$subset[mod.v2$subset != 7], drop = F]
    data.v2$y <- data.f2[,7]

    reg.v2 <- lm(y~., data = data.v2)
    alpha2 <- reg.v2$coefficients[1]
    beta2[mod.v2$subset[mod.v2$subset != 7]] <- reg.v2$coefficients[-1]


    data.v3 <- data.f3[, mod.v3$subset[mod.v3$subset != 7], drop = F]
    data.v3$y <- data.f3[,7]

    reg.v3 <- lm(y~., data = data.v3)
    alpha3 <- reg.v3$coefficients[1]
    beta3[mod.v3$subset[mod.v3$subset != 7]] <- reg.v3$coefficients[-1]

    #alpha vector
    alpha_vs <- as.numeric(c(alpha1, alpha2, alpha3))

    order_b <- order(c(beta1[1], beta2[1], beta3[1]))
    beta.mtx <- rbind(beta1, beta2, beta3)
    order_true <- order(beta_true[,1])

    #se of alpha
    se_a1.2 <- (alpha_vs[order_b[1]] - alpha_true[order_true[1]])^2
    se_a2.2 <- (alpha_vs[order_b[2]] - alpha_true[order_true[2]])^2
    se_a3.2 <- (alpha_vs[order_b[3]] - alpha_true[order_true[3]])^2

    #se of beta
    se_b1 <- sum((beta.mtx[order_b[1],] - beta_true[order_true[1],])^2)
    se_b2 <- sum((beta.mtx[order_b[2],] - beta_true[order_true[2],])^2)
    se_b3 <- sum((beta.mtx[order_b[3],] - beta_true[order_true[3],])^2)

    #fdr
    fdr1 <- sum(beta_true[order_true[1],][which(beta.mtx[order_b[1],] == 0)] != 0)/
        sum(beta.mtx[order_b[1],] == 0)
    fdr2 <- sum(beta_true[order_true[2],][which(beta.mtx[order_b[2],] == 0)] != 0)/
        sum(beta.mtx[order_b[2],] == 0)
    fdr3 <- sum(beta_true[order_true[3],][which(beta.mtx[order_b[3],] == 0)] != 0)/
        sum(beta.mtx[order_b[3],] == 0)

    #error
    err1 <- sum(beta.mtx[order_b[1],][which(beta_true[order_true[1],] == 0)] != 0 )/
        sum(beta_true[order_true[1],] == 0)
    err2 <- sum(beta.mtx[order_b[2],][which(beta_true[order_true[2],] == 0)] != 0 )/
        sum(beta_true[order_true[2],] == 0)
    err3 <- sum(beta.mtx[order_b[3],][which(beta_true[order_true[3],] == 0)] != 0 )/
        sum(beta_true[order_true[3],] == 0)

    result_2 <- list(se_b1.2 = se_b1, se_b2.2 = se_b2,
                     se_b3.2 = se_b3, se_a1.2 = se_a1.2, se_a2.2 = se_a2.2, se_a3.2 = se_a3.2,
                     FS1.2 = fdr1, FS2.2 = fdr2, FS3.2 = fdr3,
                     MS1.2 = err1, MS2.2 = err2, MS3.2 = err3)
    result0 <- list(clust_no_g = clust_no_g, ari_no_g = ari_no_g, cluster_no_g_vs =
                        cluster_no_g_vs, ari_no_g_vs = ari_no_g_vs, clust_g3 = clust_g3,
                    ari_g3 = ari_g3, clust_g3_vs = clust_g3_vs, ari_g3_vs = ari_g3_vs)
    results = c(result0, result_1, result_2)
    results
}
