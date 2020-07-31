library(ggplot2)

# each result is call exp{#}, where # is the experiment number
load('./experiment1-results.Rdata') 
load('./experiment2-results.Rdata')
load('./experiment3-results.Rdata')
load('./experiment4-results.Rdata')

#' @description Convience function used to enumerate summary results for each
#' instance of pg-ts.
#' 
#' @param n_pgts The number of pg-ts settings being summarized.
makePgtsNames <- function(n_pgts) {
    output <- rep(NA, 3*n_pgts)
    for (i in 1:n_pgts) {
        start <- 3*(i-1) + 1
        end <- start + 2
        output[start:end] <- c(
            paste('lcb_pgts', i, sep = ""),
            paste('mid_pgts', i, sep = ""),
            paste('ucb_pgts', i, sep = "") 
        )
    }
    return(output)
}

#' @description Convience for summarizing method performance. Calculates a measure
#'  of center with some percentile bands for the method.
#' 
#' @param cumregret_t A vector of the cumulative regret at time t for replicates 
#'  of the same method.
#' @param center_fn A callable, which takes a vector and gives a measure of center.
#' @param bounds A vector of length 2 that is between 0 and 1 giving the percentiles
#'  to report for the method.
summarizeRegretTime <- function(cumregret_t, center_fn, bounds) {
    bands <- quantile(cumregret_t, probs = bounds)
    center <- center_fn(cumregret_t)
    return(c(bands[1], center, bands[2]))
}

#' @description Calculates the summary of the cumulative regret at each time
#'  for the random, mab, and pg-ts methods. Note: we do not report the bands
#'  for the random, as it is only used for a baseline.
#' 
#' @param result_list list of results from each iteration of the experiment.
#' @param summary_type either 'average' or 'quantile'. The 'average' option gives
#'  average regret for each time point for random, mab, and pgts. The 'quantile' 
#'  option gives the median regret for each time point for random, mab, and pgts.
#' @param bounds vector of length 2, giving the percentile bounds to present for
#'  the mab and pgts methods. Each element should be in (0,1) and the first 
#'  element is less than the second.
#' @param pgts_index vector for the index of different pgts settings to use.
summarizeExperiment <- function(
    result_list, 
    summary_type,
    bounds, 
    pgts_index
) {
    if (bounds[1] > bounds[2]) {
        stop('Lower bound must be less than the upper bound.')
    }

    n_rep <- length(result_list)
    n_time <- length(result_list[[1]]$random_regret)
    n_pgts <- length(pgts_index)

    random_cumregret <- matrix(NA, nrow = n_time, ncol = n_rep)
    mab_cumregret <- matrix(NA, nrow = n_time, ncol = n_rep)
    pgts_cumregret <- rep(list(matrix(NA, nrow = n_time, ncol = n_rep)), n_pgts)

    for (i in 1:n_rep) {
        random_cumregret[,i] <- cumsum(result_list[[i]]$random_regret)
        mab_cumregret[,i] <- cumsum(result_list[[i]]$mab_regret) 
        for (j in 1:n_pgts) {
            pgts_cumregret[[j]][,i] <- cumsum(result_list[[i]]$pgts_regret[,pgts_index[j]])
        }
    }

    if (summary_type == 'average') {
        getCenter <- mean
    } else {
        getCenter <- median
    }

    # time, random cumregret, (lcb, middle, ucb) for mab, (lcb, middle, ucb) for each pgts
    plot_info <- matrix(NA, nrow = n_time, ncol = 1+1+3+(3*n_pgts))
    plot_info[,1] <- 1:n_time
    plot_info[,2] <- apply(random_cumregret, 1, getCenter)
    plot_info[,3:5] <- t(apply(
        mab_cumregret, 
        1, 
        summarizeRegretTime, 
        center_fn = getCenter, 
        bounds = bounds
    ))
    for (i in 1:n_pgts) {
        start <- 5 + (3*(i-1))+1
        end <- start + 2
        plot_info[,start:end] <- t(apply(
            pgts_cumregret[[i]], 
            1, 
            summarizeRegretTime, 
            center_fn = getCenter, 
            bounds = bounds
        ))
    }

    plot_df <- as.data.frame(plot_info)
    colnames(plot_df) <- c(
        c('Time', 'random', 'lcb_mab', 'mid_mab', 'ucb_mab'),
        makePgtsNames(n_pgts)
    )

    return(plot_df)
}

plotRegretSummary <- function(summary_df, title) {
    ggplot(summary_df, aes_string(x = 'Time')) +
        geom_line(aes_string(y = 'random')) +
        geom_line(aes_string(y = 'mid_mab'), col = 'blue') + 
        geom_ribbon(aes_string(ymin = 'lcb_mab', ymax = 'ucb_mab'), col = 'blue', alpha = .2) +
        geom_line(aes_string(y = 'mid_pgts1'), col = 'darkgreen') + 
        geom_ribbon(aes_string(ymin = 'lcb_pgts1', ymax = 'ucb_pgts1'), col = 'darkgreen', alpha = .2) +
        geom_line(aes_string(y = 'mid_pgts2'), col = 'orange') + 
        geom_ribbon(aes_string(ymin = 'lcb_pgts2', ymax = 'ucb_pgts2'), col = 'orange', alpha = .2) +
        ggtitle(title) +
        xlab('Time') +
        ylab('Cumulative Regret') +
        theme(plot.title = element_text(hjust = 0.5))
}

exp1_summary <- summarizeExperiment(exp1, 'average', c(.05, .95), c(1,2))
exp2_summary <- summarizeExperiment(exp2, 'average', c(.05, .95), c(1,2))
exp3_summary <- summarizeExperiment(exp3, 'average', c(.05, .95), c(1,2))
exp4_summary <- summarizeExperiment(exp4, 'average', c(.05, .95), c(1,2))

# this is ran inside the experiments folder - change path accordingly
ggsave(
    filename = '../figures/experiment1-plot.pdf', 
    plot = plotRegretSummary(exp1_summary, 'Setting#1') # Small #Arms & Large #Features
)
# ggsave(
#     filename = '../figures/experiment2-plot.pdf', 
#     plot = plotRegretSummary(exp2_summary, 'Experiment#2 - Small #Arms & Large #Features')
# )
ggsave(
    filename = '../figures/experiment3-plot.pdf', 
    plot = plotRegretSummary(exp3_summary, 'Setting#2') # Large #Arms & Small #Features
)
# ggsave(
#     filename = '../figures/experiment4-plot.pdf', 
#     plot = plotRegretSummary(exp4_summary, 'Experiment#4 - Large #Arms & Large #Features')
# )