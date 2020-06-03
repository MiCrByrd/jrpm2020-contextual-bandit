banditSimulation <- function(
    horizon,
    n_arms,
    n_instance_features,
    n_arm_features,
    alpha,
    beta,
    n_sample,
    window_size
) {
    if (length(n_sample) != length(window_size)) {
        stop('n_sample and window_size must be the same length')
    }
    n_pgts_settings <- length(n_sample)

    # arms and underlying linear relationship
    arms <- createArms(
        n_arms = n_arms,
        n_features = n_arm_features
    )
    true_coef <- createCoefficients(
        n_arm_features = n_arm_features,
        n_instance_features = n_instance_features
    )

    # performance evaluated by regret at each time point for each method
    random_regret <- rep(NA, horizon)
    mab_regret <- rep(NA, horizon)
    pgts_regret <- matrix(
        data = rep(NA, horizon*n_pgts_settings),
        ncol = n_pgts_settings
    )

    # tracker for mab and pgts
    mab_tracker <- mabCreateContainer(
        n_arms = n_arms,
        alpha = alpha,
        beta = beta
    )
    pgts_trackers <- lapply(
        X = 1:n_pgts_settings,
        FUN = function(j) {
                pgtsCreateContainer(
                mu = matrix(
                    data = rep(0, true_coef$p),
                    nrow = true_coef$p,
                    ncol = 1
                ),
                sigma = diag(1, true_coef$p),
                n_sample = n_sample[j],
                window_size = window_size[j]
            )
        }
    )

    for (i in 1:horizon) { 
        instance <- createInstance(n_instance_features)

        # make recommendations
        random_suggest <- sample(
            x = 1:n_arms, 
            size = 1
        )
        mab_suggest <- mabRecommend(mab_tracker)
        pgts_suggests <- lapply(
            X = pgts_trackers,
            FUN = function(tracker) {
                pgtsRecommend(
                    pgts_container = tracker,
                    instance = instance,
                    arms = arms
                )
            }
        )

        # calculate truth, identify success, and calculate regret
        succ_prob <- calcChoiceProbability(
            instance = instance,
            arms = arms,
            coefficients = true_coef$coefficients
        )

        random_success <- rbinom(
            n = 1, 
            size = 1, 
            prob = succ_prob[random_suggest]
        )
        mab_success <- rbinom(
            n = 1,
            size = 1,
            prob = succ_prob[mab_suggest]
        )
        pgts_successes <- lapply(
            X = pgts_suggests,
            FUN = function(suggest) {
                rbinom(
                    n = 1,
                    size = 1,
                    prob = succ_prob[suggest]
                )
            }
        )

        ideal_prob <- max(succ_prob)
        random_regret[i] <- ideal_prob - succ_prob[random_suggest]
        mab_regret[i] <- ideal_prob - succ_prob[mab_suggest]
        pgts_regret[i,] <- sapply(
            X = pgts_suggests,
            FUN = function(suggest) {
                ideal_prob - succ_prob[suggest]
            }
        )

        # update
        mab_tracker <- mabUpdateContainer(
            mab_container = mab_tracker,
            arm_index = mab_suggest,
            is_success = mab_success
        )
        pgts_trackers <- lapply(
            X = 1:n_pgts_settings,
            FUN = function(j) {
                pgtsUpdateContainer(
                    pgts_container = pgts_trackers[[j]],
                    arm_index = pgts_suggests[[j]],
                    instance_features = instance,
                    arms = arms,
                    is_success = pgts_successes[[j]]
                )
            }
        )
    }

    return(
        list(
            random_regret = random_regret,
            mab_regret = mab_regret,
            pgts_regret = pgts_regret
        )
    )
}