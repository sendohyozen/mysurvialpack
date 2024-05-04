#' Find the best cutoff threhold for a numeric var; if median grouping could not give a significant P value
#'
#' @param surv.data a dataframe with multiple varnames: column: time, vital status, var names; row: samples
#' @param time colnames of survial time
#' @param vital colnames of vital status
#' @param numeric.varName colnames  var names as colnames (must be numeric vectors)
#' @param positive.lab positive label defines event occurs; such as 'Dead'
#' @param minprop minimal proption to prevent unbanlanced grouping
#'
#' @return a best cutoff surv fit cutoff value
#' @export
#'
#' @examples surv_bestcutoff(surv.data = sur_data, time = "days_to_last_follow_up", vital = "vital_status", numeric.varName = "total_perMB", positive.lab = 'Dead')
surv_bestcutoff <- function(surv.data, time, vital, numeric.varName, positive.lab, minprop=0.3){

    # survival data
    df = data.frame(time = surv.data[[time]],
                    status = surv.data[[vital]],
                    var = surv.data[[numeric.varName]])

    # data type transforming
    df$status = ifelse(df$status == positive.lab, 1, 0)  # vital statues: event positive 1; negative, 0
    df$time = as.numeric(df$time)  # time as numeric
    df$var = as.numeric(df$var)  # var must be a numeric vector

    colnames(df) = c("time", "status", numeric.varName)

    # find the best cutoff
    best_threshold_surv <- survminer::surv_cutpoint(df,
                                                    time = "time",
                                                    event = "status",
                                                    variables = numeric.varName,  # variable to find the cutoff
                                                    minprop = minprop,     # minimal proption to prevent unbanlance
                                                    progressbar = TRUE)  # show progressbar

    return(best_threshold_surv$cutpoint$cutpoint)
}
