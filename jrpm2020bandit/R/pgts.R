require(pgdraw) # for Polya-Gamma distribution
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
    window_size,
    diagonal = False
) {
    if (diagonal) {
        sigma = diag(sigma)
        omega = 1 / sigma
        omega_x_mu = mu * omega
    } else {
        omega = solve(sigma)
        omega_x_mu = omega %*% mu
    }

    return (
        list(
            # priors
            mu = mu,
            sigma = sigma,
            omega = omega,
            omega_x_mu = omega_x_mu, # needed for each update

            # gibbs sampler param
            n_sample = n_sample,
            window_size = window_size,
            diagonal = diagonal,

            # data
            t = 0,
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
            theta = mu
        )
    )
}

## Sampler required for generating coefficients
# pgts_container: The pg-ts tracker list.
pgtsSampler <- function(pgts_container) {

    if (pgts_container$t >= pgts_container$window_size) {
        x <- pgts_container$x
        y <- pgts_container$y
    } else {
        # take the subset of the data matrix that currently exists
        x <- pgts_container$x[1:pgts_container$t, , drop = F]
        y <- pgts_container$y[1:pgts_container$t, 1, drop = F]
    } 

    kappa <- y - 1/2
    curr_theta <- pgts_container$theta

    # Gibbs sampler
    for (i in 1:pgts_container$n_sample) {
        z = pgdraw::pgdraw(1, c(x %*% curr_theta))

        if (pgts_container$diagonal) {
            v_diag <- 1 / colSums(x^2 * z)
            v <- diag(v_diag)
            m <- v_diag * (t(x) %*% kappa + pgts_container$omega_x_mu)
        } else {
            v <- solve(
                (t(x) %*% diag(z, nrow = length(z)) %*% x) 
                + pgts_container$omega
            )
            m <- v %*% (t(x) %*% kappa + pgts_container$omega_x_mu)
        }

        curr_theta <- matrix(
            data = MASS::mvrnorm(
                n = 1,
                mu = m,
                Sigma = v
            ),
            ncol = 1
        )
    }

    return(curr_theta)
}

## Update the pg-ts tracker list
# pgts_container: The pg-ts tracker list.
# arm_index: The index for the arm suggested.
# instance_features: Row vector of instance features.
# arms: arm list
# is_success: Bool indicating if it was a success or not.
pgtsUpdateContainer <- function(
    pgts_container,
    arm_index,
    instance_features,
    arms,
    is_success
) {
    # update time
    pgts_container$t = pgts_container$t + 1
    i = pgts_container$t %% pgts_container$window_size
    # why must R indexing start at 1... don't want to refactor
    if (i == 0) i <- pgts_container$window_size

    # update observations
    pgts_container$x[i,] = c(
        arms$features[arm_index,],
        makeInstanceArmInteractions(
            instance = instance_features,
            arms = arms
        )[arm_index,]
    )
    pgts_container$y[i,1] = is_success

    # run gibbs sampler
    pgts_container$theta <- pgtsSampler(pgts_container)

    return(pgts_container)
}

## Get arm recommendation from pg-ts for available data
# pgts_container: pg-ts tracker list.
# instance: vector of features for the instance
# arms: arm tracker list.
pgtsRecommend <- function(
    pgts_container, 
    instance,
    arms
) {
    # calculate prob of success for each arm with drawn coefficients
    arm_probs = calcChoiceProbability(
        instance = instance, 
        arms = arms, 
        coefficients = pgts_container$theta
    )

    return(which.max(arm_probs))
}