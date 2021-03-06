## Create arms for the simulation study.
# n_arms: number of arms for the experiment
# n_features: number of features for each arm.
createArms <- function(n_arms, n_features) {
    return(
        list(
            n_arms = n_arms,
            n_features = n_features,
            features = matrix(
                data = rnorm(n = n_arms * n_features),
                nrow = n_arms,
                ncol = n_features
            )
        )
    )
}

## Create a new instance to be given a recommendation.
# n_features: number of features for the instance.
createInstance <- function(n_features) {
    return(rnorm(n = n_features))
}

## Create the true coefficients for a simulation study.
# n_arm_features: Number of arm features.
# n_instance_features: Number of instance features.
createCoefficients <- function(
    n_arm_features, 
    n_instance_features
) {
    # the total number of features
    # don't directly use instance as it doesn't distinguish info
    # includes interactions between instance and arm
    p = n_arm_features + n_arm_features*n_instance_features

    return(
        list(
            p = p,
            coefficients = matrix(
                data = rnorm(n = p),
                nrow = p,
                ncol = 1
            )
        )
    )
}

## For two feature matrices of same number of rows, create feature interactions.
# m1: Matrix#1 
# m2: Matrix#2 of the same number of rows as m1
makeInstanceArmInteractions <- function(instance, arms) {

    broadcast_instance <- t(matrix(
            data = instance,
            nrow = length(instance),
            ncol = arms$n_arms
    ))

    return (
        do.call(
            what = cbind,
            args = lapply(1:ncol(broadcast_instance), function(i) {
                broadcast_instance[,i] * arms$features
            })
        )
    )
}

## Calculate the probability of success for each arm
# instance: instance features
# arms: arm tracker list
# coefficients: p \times 1 matrix of model coefficients
calcChoiceProbability <- function(
    instance, 
    arms, 
    coefficients
) {
    features <- cbind(
        # arm features
        arms$features,
        # interaction of instance and arm features - hacky but idk how else
        makeInstanceArmInteractions(
            instance = instance,
            arms = arms
        )
    )

    xb <- features %*% coefficients
    return(c(1 / (1 + exp(-xb))))
}

