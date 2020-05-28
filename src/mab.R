# Create tracker list for the multi-arm bandit
# n_arms: the number of arms to track
# alpha: prior beta distribution parameter
# beta: prior beta distribution parameter
mabCreateContainer <- function(n_arms, alpha, beta) {

    return list(
        # priors
        alpha = alpha,
        beta = beta,
        # outcomes
        n_trials = rep(0, n_arms)
        n_success = rep(0, n_arms)
    )
}

# Update the tracker list for the MAB
# mab_container: mab tracker list
# arm_index: arm that was recommended
# is_success: bool indicating the result of the arm
mabUpdateContainer <- function(mab_container, arm_index, is_success) {

    mab_container$n_trials[arm_index] <- mab_container$n_trials + 1
    mab_container$n_success[arm_index] <- mab_container$n_success + is_success

    return(mab_container)

}

# Make a recommendation from the current MAB
# mab_container: mab tracker list
mabRecommend <- function(mab_container) {

    a <- mab_container$n_success + mab_container$alpha
    b <- mab_container$n_trials - mab_container$n_success + mab_container$beta

    theta <- mapply(
        FUN = rbeta,
        shape1 = a,
        shape2 = b,
        MoreArgs = list(n = 1)
    )

    return(which.max(theta))

}




