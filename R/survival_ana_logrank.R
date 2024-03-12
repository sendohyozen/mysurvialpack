#' Geting P value using logrank from a dataframe with multiple varnames
#'
#' @param surv.data a dataframe with multiple varnames: column: time, vital status, var names; row: samples
#' @param time  colnames of survial time
#' @param vital colnames of vital status
#' @param group  colnames  var names as colnames
#' @param positive.lab positive label defines event occurs; such as 'Dead'
#' @param factor.Levels define factor levels of the var name (group)
#'
#' @return a vector of var name and pvalue
#' @export
#'
#' @examples logrank_pv(surv.data = sur_data, time = 'days_to_last_follow_up', vital = 'vital_status', group = snp,
#'                   positive.lab = 'Dead', factor.Levels = c('wild', 'mutation'))
logrank_pv = function(surv.data, time, vital, group, positive.lab, factor.Levels){

    # var name for survial analysis
    varName = group

    # survial data
    df = data.frame(time = surv.data[[time]],
                    status = surv.data[[vital]],
                    group = surv.data[[group]])

    # data type transforming
    df$status = ifelse(df$status == positive.lab, 1, 0)  # vital statues: event positive 1; negative, 0
    df$time = as.numeric(df$time)  # time as numeric
    df$group = factor(df$group, levels = factor.Levels)  # var as factor, and levels defined


    # survial analysis- geting p value with logrank method
    x = survival::survdiff(Surv(time, status) ~ group, data = df)
    pValue= round(1-pchisq(x$chisq,df=1), 3)

    # return var name and pvalue
    res = c(varName, pValue)
}
