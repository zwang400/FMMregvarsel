r_sq_post <- function(sim_res, data, scl=10001:1e5){
    #####################
    #data: data frame containing all covariates and response, note response has to have name "y"
    #####################
    D <- ncol(data)-1
    y <- data$y
    k_post <- sim_res$k_post[scl]
    #M_post <- sim_res$M_post[scl]
    c_post <- sim_res$c_post[scl,]
    r_sqr <- numeric(0)

    for(i in 1:nrow(c_post)){
        for(j in 1:k_post[i]){
            n_j <- which(c_post[i, ] == j)
            if(length(n_j) > D){
                r_sqr <- c(r_sqr, summary(lm(y~., data[n_j, ]))$r.square)
            }
        }
    }
    r_sqr
}
