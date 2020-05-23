##
#
createContainer <- function(n_arms, alpha, beta) {

    return list(
        # priors
        alpha = alpha,
        beta = beta,
        # outcomes
        n_trials = rep(0, n_arms)
        n_success = rep(0, n_arms)
    )
}

updateContainer <- function(mab_container, arm_index, is_success) {

    mab_container$n_trials[arm_index] <- mab_container$n_trials + 1
    mab_container$n_success[arm_index] <- mab_container$n_success + is_success

    return(mab_container)

}

sampler <- function(mab_container) {

    a <- mab_container$n_success + mab_container$alpha
    b <- mab_container$n_trials - mab_container$n_success + mab_container$beta

    theta <- mapply(
        FUN = rbeta,
        shape1 = a,
        shape2 = b,
        MoreArgs = list(n = 1)
    )

    return(theta)
}



