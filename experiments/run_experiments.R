library(parallel)

n_replicates <- 25
# adjust for the system
c1 <- makeCluster(6)

## Small #arms Small #features
exp1 <- parLapply(c1, 1:n_replicates, function(i) {
    library(jrpm2020bandit)
    banditSimulation(
        horizon = 1000,
        n_arms = 10,
        n_instance_features = 2,
        n_arm_features = 2,
        alpha = 1, 
        beta = 1,
        n_sample = c(25, 25, 1, 1),
        window_size = c(1000, 100, 1000, 100)
    )
})
save(exp1, file = 'experiment1-results.Rdata')

## Small #arms Large #features
exp2 <- parLapply(c1, 1:n_replicates, function(i) {
    library(jrpm2020bandit)
    banditSimulation(
        horizon = 1000,
        n_arms = 10,
        n_instance_features = 2,
        n_arm_features = 9,
        alpha = 1, 
        beta = 1,
        n_sample = c(25, 25, 1, 1),
        window_size = c(1000, 100, 1000, 100)
    )
})
save(exp2, file = 'experiment2-results.Rdata')

## Large #arms Small #features
exp3 <- parLapply(c1, 1:n_replicates, function(i) {
    library(jrpm2020bandit)
    banditSimulation(
        horizon = 1000,
        n_arms = 100,
        n_instance_features = 2,
        n_arm_features = 2,
        alpha = 1, 
        beta = 1,
        n_sample = c(25, 25, 1, 1),
        window_size = c(1000, 100, 1000, 100)
    )
})
save(exp3, file = 'experiment3-results.Rdata')

## Large #arms Large #features
exp4 <- parLapply(c1, 1:n_replicates, function(i) {
    library(jrpm2020bandit)
    banditSimulation(
        horizon = 1000,
        n_arms = 100,
        n_instance_features = 2,
        n_arm_features = 9,
        alpha = 1, 
        beta = 1,
        n_sample = c(25, 25, 1, 1),
        window_size = c(1000, 100, 1000, 100)
    )
})
save(exp4, file = 'experiment4-results.Rdata')
