createArms <- function(n_arms, n_features) {
    return(
        list(
            n_arms = n_arms,
            n_features = n_features,
            features = matrix(
                data = rnorm(n = n_arms * n_cts_features),
                nrow = n_arms,
                ncol = n_cts_features
            )
        )
    )
}

createInstance <- function(n_features) {
    return(rnorm(n = n_features))
}

createCoefficients <- function(
    n_arms, 
    n_arm_features, 
    n_instance_features
) {
    # the total number of features
    # includes interactions between instance and arm
    p = n_arm_features + n_instance_features + n_arm_features*n_instance_features

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

calcChoiceProbability <- function(
    instance, 
    arms, 
    coefficients
) {
    broadcast_instance <- t(matrix(
            data = instance,
            nrow = length(instance)
            ncol = arms$n_arms
    ))

    features <- cbind(
        # instance features
        broadcast_instance,
        # arm features
        arms$features,
        # interaction of instance and arm features - hacky but idk how else
        do.call(
            what = cbind,
            args = lapply(1:length(instance), function(i) {
                broadcast_instance[,i] * arms$features
            })
        )
    )

    xb <- features %*% coefficients$coefficients
    return(1 / (1 + exp(-xb)))
}

