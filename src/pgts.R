required(pgdraw) # for Polya-Gamma distribution
source('./utility.R')

## Create tracker list for PG-TS
# mu: p \times 1 matrix for the prior mean.
# sigma: p \times p matrix for the prior covariance. Must be positive definite.
# n_sample: The number of samples iterations for the Gibbs sampler.
# window_size: The size of the rolling window to keep for the data history.
pgtsCreateContainer <- function(
    mu,
    sigma,
    n_sample,
    window_size
) {
    omega = solve(sigma)

    return (
        list(
            # priors
            mu = mu,
            sigma = sigma,
            omega = omega,
            omega_x_mu = omega %*% mu, # needed for each update

            # gibbs sampler param
            n_sample = n_sample,
            window_size = window_size

            # data
            t = -1,
            x = matrix(
                data = NA, 
                nrow = window_size,
                ncol = dim(mu)[1]
            ),
            y = matrix(
                data = NA,
                nrow = window_size,
                ncol = 1
            ),

            # coefficients
            theta = matrix(
                data = NA,
                nrow = dim(mu)[1],
                ncol = 1
            )
        )
    )
}

## Update the pg-ts tracker list
# pgts_container: The pg-ts tracker list.
# arm_index: The index for the arm suggested.
# instance_features: Row vector of instance features.
# arm_features: k \timex p matrix of arm features, each row for one arm.
# is_success: Bool indicating if it was a success or not.
pgtsUpdateValues <- function(
    pgts_container,
    arm_index,
    instance_features,
    arms,
    is_success
) {
    # update time
    pgts_container$t = pgts_container$t + 1
    i = (pgts_container$t %% pgts_container$window_size) + 1

    # update observations
    pgts_container$x[i,] = c(instance_features, arm_features[arm_index,])
    pgts_container$y[i,] = is_success
}

## Sampler required for generating coefficients
# pgts_container: The pg-ts tracker list.
pgtsSampler <- function(pgts_container) {

    kappa <- pgts_container$y - 1/2

    if (pgts_container$t == 0) {
        # initialize theta to be the specified prior mean
        curr_theta <- pgts_container$mu
    } else {
        curr_theta <- pgts_container$theta
    }

    if (pgts_container$t+1 >= window_size) {
        x <- pgts_container$x
        y <- pgts_container$y
    } else {
        # take the subset of the data matrix that currently exists
        # constructor because I don't trust R's magic type converting
        x <- matrix(
            data = pgts_container$x[1:(pgts_container$t+1),], 
            nrow = pgts_container$t + 1
        )
        y <- matrix(
            data = pgts_container$y[1,(pgts_container$t+1)], 
            nrow = pgts_container$t + 1
        )
    } 

    # Gibbs sampler
    for (i in 1:pgts_container$n_sample) {
        z = pgdraw::pgdraw(1, c(x %*% curr_theta))   

        v <- solve((t(x) %*% diag(z) %*% x) + pgts_container$omega)
        m <- v %*% (t(x) %*% kappa + pgts_container$omega_x_mu)

        curr_theta <- mvtnorm::rmvnorm(
            n = 1,
            mean = m,
            sigma = v
        )
    }

    return(curr_theta)
}

## Get arm recommendation from pg-ts for available data
# pgts_container: pg-ts tracker list.
# instance: vector of features for the instance
# arms: arm tracker list.
pgtsRecommender <- function(
    pgts_container, 
    instance,
    arms
) {
    # run gibbs sampler
    pgts_container$theta[,1] <- pgtsSampler(pgts_container)
    # calculate prob of success for each arm with drawn coefficients
    arm_probs = calcChoiceProbability(
        instance = instance, 
        arms = arms, 
        coefficients = pgts_container$theta
    )

    return which.max(arm_probs)
}