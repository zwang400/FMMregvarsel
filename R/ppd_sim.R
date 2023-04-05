ppd_sim <- function(sim.res, x.new, y.new, scl=10001:1e5, thin=20){
    #some constants
    a_lambda <- sim.res$a_lambda
    b_lambda = sim.res$b_lambda
    a_eta <- sim.res$a_eta
    b_eta <- sim.res$b_eta
    X <- sim.res$X
    y <- sim.res$y
    n <- length(y)
    D <- ncol(X)
    a_star <- a_lambda + n/2
    a_jing <- n*D/2 + a_eta

    Beta_post <- sim.res$Beta_post
    zeta_post <- sim.res$zeta_post
    alpha_post <- sim.res$alpha_post
    M_post <- sim.res$M_post
    c_post <- sim.res$c_post
    weight_post <- sim.res$weight_post

    sample.size <- length(scl)
    sample.index <- seq(scl[1], scl[sample.size], thin)
    answ1 <- answ2 <- rep(0, length(sample.index))
    k <- 1
    const <- -0.5*log(2*pi) + lgamma(a_star+0.5) - lgamma(a_star)
    #calculate each iteration
    for(i in sample.index){
        beta.sample <- matrix(Beta_post[[i]], ncol = D)
        zeta.sample <- matrix(zeta_post[[i]], ncol = D)
        alpha.sample <- alpha_post[[i]]
        b_star <- b_lambda + sum((y - apply((beta.sample[c_post[i, ],])*X, 1, sum) - alpha.sample[c_post[i, ]])^2)/2
        b_jing <- b_eta + sum((X - zeta.sample[c_post[i, ],])^2)/2
        M.sample <- M_post[i]
        #log-scale, for each term from 1:M
        #numerator
        term.1 <- log(weight_post[[i]]) - (a_star + 1/2) *
            log(((y.new - x.new%*%t(beta.sample) - alpha.sample)^2)/2 + b_star) -
            (a_jing + D/2) * log(apply((apply(zeta.sample, 1, function(s) x.new-s))^2, 2, sum)/2 + b_jing) +
            a_star*log(b_star) + a_jing*log(b_jing)

        #denominator
        term.2 <- log(weight_post[[i]]) - (a_jing + D/2) *
            log(apply((apply(zeta.sample, 1, function(s) x.new-s))^2, 2, sum)/2 + b_jing) + a_jing*log(b_jing)

        #n.const.anw <- mean(c(max(term.1), max(term.1)))
        answ1[k] <- log(sum(exp(term.1)))
        answ2[k] <- log(sum(exp(term.2)))

        k <- k + 1
    }
    n.const <- mean(c(max(answ1), max(answ2)))
    exp(const + log(mean(exp(answ1 - n.const))) - log(mean(exp(answ2 - n.const))))
}
