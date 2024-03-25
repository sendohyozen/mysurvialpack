#' Geting P value using logrank from a dataframe with multiple varnames
#'
#' @param surv.data a dataframe with multiple varnames: column: time, vital status, var names; row: samples
#' @param time  colnames of survial time
#' @param vital colnames of vital status
#' @param group  colnames  var names as colnames
#' @param positive.lab positive label defines event occurs; such as 'Dead'
#' @param factor.Levels define factor levels of the var name (group)
#' @param cutoff group dividing for numeric var value; if a give optimized cutoff value is given, default median
#'
#' @return a vector of var name and pvalue
#' @export
#'
#' @examples logrank_pv(surv.data = sur_data, time = 'days_to_last_follow_up', vital = 'vital_status', group = snp,
#'                   positive.lab = 'Dead', factor.Levels = c('wild', 'mutation'))
logrank_pv = function(surv.data, time, vital, group, cutoff,
                      positive.lab, factor.Levels){

    # var name for survial analysis
    varName = group

    # survial data
    df = data.frame(time = surv.data[[time]],
                    status = surv.data[[vital]],
                    group = surv.data[[group]])

    # if group is continous numeric value, transform into discrete groups
    if(is.numeric(df$group)){
        cat('group contains continous numeric value, transform into discrete groups, default using median, else used defined cutoff')
        if(missing(cutoff)){
             # default using median for numeric group
            df$group = ifelse(df$group > median(df$group), 'high', 'low')
            df$group = factor(df$group, levels = c('low', 'high'))
        }else{
            # given a defined cutoff value
            df$group = ifelse(df$group > cutoff, 'high', 'low')
            df$group = factor(df$group, levels = c('low', 'high'))
        }

    }else{
        cat('group is better be factors, or give a previously cutted values!')
        if(missing(factor.Levels)){
            warning("is group a factor or character data type? better give a factor level!")
        }else{
            df$group = factor(df$group, levels = factor.Levels)  # var as factor, and levels defined
        }
    }


    # data type transforming
    df$status = ifelse(df$status == positive.lab, 1, 0)  # vital statues: event positive 1; negative, 0
    df$time = as.numeric(df$time)  # time as numeric

    # survial analysis- geting p value with logrank method
    x = survival::survdiff(Surv(time, status) ~ group, data = df)
    pValue= round(1-pchisq(x$chisq,df=1), 3)

    # return var name and pvalue
    res = c(varName, pValue)
}
