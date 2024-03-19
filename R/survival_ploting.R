#' Kaplan Meier plot using survminer package
#'
#' @param surv.data a dataframe with multiple varnames: column: time, vital status, var names; row: samples
#' @param time colnames of survial time
#' @param vital colnames of vital status
#' @param group colnames  var names as colnames
#' @param positive.lab positive label defines event occurs; such as 'Dead'
#' @param factor.Levels define factor levels of the var name (group)
#' @param font.x font size of x axis title, default 15
#' @param font.y font size of y axis title,default 15
#' @param font.tickslab font size of labs in x and y,default 15
#' @param font.legend font size of legend title and labs, default 12
#' @param risk.table if show risk table, default False
#' @param title title of the plot
#' @param x.lab x lab
#' @param y.lab y lab
#' @param conf.int  if add confidence ,default False
#' @param col palette ,defalut lancet
#'
#' @return a KM plot
#' @export
#'
#' @examples surv_plot1(surv.data = sur_data, time = 'days_to_last_follow_up', vital = 'vital_status', group = snp,
#'                      positive.lab = 'Dead', factor.Levels = c('wild', 'mutation'))
#'
surv_plot1 = function(surv.data, time, vital, group, positive.lab, factor.Levels,
                      title = NULL, x.lab = 'Time', y.lab = 'Survival probability',
                      font.x = 15, font.y = 15, font.tickslab = 15, font.legend = 12,
                      risk.table = F, conf.int = F,
                      col = ggsci::pal_lancet(alpha = 0.8)(9)
                      ){

    # https://zhuanlan.zhihu.com/p/596831604?utm_id=0

    # var name for survial analysis
    varName = group

    # survival data
    df = data.frame(time = surv.data[[time]],
                    status = surv.data[[vital]],
                    group = surv.data[[group]])

    # data type transforming
    df$status = ifelse(df$status == positive.lab, 1, 0)  # vital statues: event positive 1; negative, 0
    df$time = as.numeric(df$time)  # time as numeric
    df$group = factor(df$group, levels = factor.Levels)  # var as factor, and levels defined


    # get a surv object
    fit <- survival::survfit(Surv(time,status)~group, data = df)

    # KM plot
    p = survminer::ggsurvplot(fit, data = df,
                              pval = T, # show P value
                              title = title, xlab = x.lab, ylab = y.lab,
                              legend.title = varName,
                              legend.labs = levels(df$group),
                              font.x = font.x, font.y = font.y,
                              font.tickslab = font.tickslab,  font.legend = font.legend,
                              risk.table = risk.table, risk.table.col = "strata",
                              conf.int = conf.int,
                              palette = col[1:nlevels(df$group)]
                              )

    return(p)
}








#' Custom version of Kaplan Meier plot from IM package
#' @description
#'  Custom display of survival data, including plotting lines for median survival and adding number at risk
#'
#' @param surv.data a dataframe with multiple varnames: column: time, vital status, var names; row: samples
#' @param time colnames of survial time
#' @param vital colnames of vital status
#' @param group colnames  var names as colnames
#' @param positive.lab  positive label defines event occurs; such as 'Dead'
#' @param factor.Levels define factor levels of the var name (group)
#' @param main main character; main title for plot
#' @param xLab character; x-axis label for plot
#' @param yLab character; y-axis label for plot
#' @param cols character vector; colors, default to darkgreen, darkmagenta, cyan4, darkorange, darkred
#' @param ltypes numeric vector; line types (lty) for plot
#' @param plotMedian logical; show the median survival for each group be plotted, default False
#' @param pval type of p-value computed, either "logrank" or "coxph" (the latter can only be computed on two groups)
#'             or "none" don't show
#' @param legPos legend position for legend , default "topright"
#' @param mar plotting margins, default c(12,9,3,2)
#' @param plot.nrisk logical; show number of samples for each group and time point,
#'                    be plotted underneath graph .default False
#' @param nrisk.interval numeric; spacing of time intervals for which number of samples
#'                        is given; defaults to 2
#' @param timemark  logical; if set to TRUE (default), censoring marks are added to survival curves
#' @param lwd  line width for survival curves, default 2
#' @param cexMedian numeric; relative size of font for median survival times , default 0.8
#' @param cexLegend numeric; relative size of legend font, default 0.9
#' @param ... other arguments passed to plot.survfit
#'
#' @return  a KM plot
#' @export
#'
#' @examples mySurvial::surv_plot2(surv.data = sur_data, time = 'days_to_last_follow_up', vital = 'vital_status', group = snp,
#'             positive.lab = 'Dead', factor.Levels = c('wild', 'mutation'),
surv_plot2 <- function(surv.data, time, vital, group, positive.lab, factor.Levels,
                       main='',  xLab="OS (months)", yLab="Probability of survival",
                       cols,
                       ltypes=1,
                       plotMedian=FALSE,
                       pval=c("logrank", "coxph", "none"),
                       legPos="topright",
                       mar=c(12,9,3,2),
                       plot.nrisk=FALSE,
                       nrisk.interval=2,
                       timemark=TRUE,
                       lwd=2,
                       cexMedian=0.8,
                       cexLegend=0.9,
                       ...) {

    # survival data
    df = data.frame(time = surv.data[[time]],
                    status = surv.data[[vital]],
                    group = surv.data[[group]])

    # data type transforming
    df$status = ifelse(df$status == positive.lab, 1, 0)  # vital statues: event positive 1; negative, 0
    df$time = as.numeric(df$time)  # time as numeric, month unit
    df$group = factor(df$group, levels = factor.Levels)  # var as factor, and levels defined

    # survfit object; survfit(Surv(TTE, Cens) ~ group, data=df)
    survFit <- survival::survfit(Surv(time,status)~group, data = df)

    # surfdiff object; survdiff(Surv(TTE, Cens) ~ group, data=df)
    survDiff <- survival::survdiff(Surv(time,status)~group, data = df)

    # surf group as factor
    diff.factor <- df$group

    stopifnot(is(survFit, "survfit"))
    stopifnot(is(survDiff, "survdiff"))
    stopifnot(is(diff.factor, "factor"))

    ## simple checks to decrease chance that fitted data used indeed
    fitCall <- as.character(survFit$call)
    diffCall <- as.character(survDiff$call)
    stopifnot(fitCall[2] == diffCall[2])
    stopifnot(fitCall[3] == diffCall[3])

    group.labels <- levels(diff.factor)  # group of factor levels

    ## setting colors
    if (missing(cols)) {
        cols <- c("darkgreen", "darkmagenta", "cyan4", "darkorange", "darkred")[1:nlevels(diff.factor)]
    }else{
        cols <- cols[1:nlevels(diff.factor)]
    }
    # stopifnot(length(cols) == length(group.labels))   # stop if number of color palette is not consistent with number of group levels

    # line types
    ltypes = rep(ltypes, nlevels(diff.factor))

    ## computing number of samples at each time point
    if (plot.nrisk) {
        time.pt <- seq(0, max(survFit$time), nrisk.interval)
        ix = 0
        n.risk = c()
        for (kk in 1:(length(survFit$strata))) {
            fit.n.risk = survFit$n.risk[(ix+1) : (ix + survFit$strata[kk])]
            fit.time = survFit$time[(ix+1) : (ix + survFit$strata[kk])]
            tmp = findInterval(time.pt, fit.time)
            n.risk=rbind(n.risk,
                         ifelse(tmp < length(fit.time), fit.n.risk[tmp+1], 0)
            )
            ix = ix + survFit$strata[kk]
        }
        dimnames(n.risk)[[2]] = time.pt

        if (mar[1]<4+length(group.labels)) {
            mar[1] <- 4+length(group.labels)
        }
        org.mar <- par()$mar
        par(mar=mar)
    }

    ## plot curves and axes
    plot(survFit,
         xlab=xLab,
         lwd=lwd,
         lty=ltypes,
         col=cols,
         ylab=yLab,
         main=main,
         axes=FALSE,
         mark.time=timemark,
         ...)
    box()
    axis(1,
         at=seq(0, max(survFit$time), nrisk.interval),
         seq(0, max(survFit$time), nrisk.interval)
    )
    axis(2,
         at=seq(0, max(survFit$time), 0.1),
         seq(0, max(survFit$time), 0.1),
         las=2)
    abline(h=0,
           col="darkgrey")

    ## add the survival medians
    if (plotMedian) {
        median.surv <- summary(survFit)$table[,"median"]
        wMed <- which(!is.na(median.surv))
        median.surv <- median.surv[!is.na(median.surv)]
        median.surv <- round(median.surv, digits=2)
        for (i in 1:length(median.surv)) {
            lines(x=rep(median.surv[i], 6),
                  y=seq(0,0.5,0.1),
                  lty=2,
                  col="darkgrey")
            lines(x=seq(0, median.surv[i], length.out=6),
                  y=rep(0.5,6),
                  lty=2,
                  col="darkgrey")
            text(x=median.surv[i],
                 y=-0.02,
                 labels=median.surv[i],
                 cex=cexMedian,
                 col=cols[wMed][i])
        }
    }

    ## add the number of samples underneath plot
    if (plot.nrisk) {
        for (i in 1:length(group.labels)) {
            mtext(side=1,
                  at=-1.7,
                  line=i+3.5,
                  text=group.labels[i],
                  col=cols[i],
                  adj=1,
                  cex=0.6)
            mtext(side=1,
                  at=time.pt,
                  line=i+3.5,
                  text=n.risk[i,],
                  col=cols[i],
                  cex=0.8)
        }
    }

    ##  compute statistics of pvalue
    if (pval == "logrank") {
        pval <- round(1 - pchisq(survDiff$chisq, length(survDiff$n) - 1),
                      digits=3)
        group.labels <- c(group.labels,
                          paste("log rank pval:", pval))
    }
    if (pval == "coxph") {
        # cox object; if needed
        coxphData <- coxph(Surv(time,status)~group, data = df)

        stopifnot(is(coxphData, "coxph"))
        coxCall <- as.vector(coxphData$formula)
        stopifnot(fitCall[2] == coxCall)
        stopifnot(sum(survFit$n) == coxphData$n)
        stopifnot(length(group.labels) == 2)

        tmp <- round(getHRandCIfromCoxph(coxphData),
                     digits=4)
        group.labels <- c(group.labels,
                          paste("pval:", tmp[1,"P"]),
                          paste("HR:",
                                paste0(tmp[1,"HR"], " (", tmp[1,"CIl0.95"], ";",
                                       tmp[1,"CIu0.95"], ")")))
    }

    ## plot legend
    legend(legPos,
           group.labels,
           lwd=2,
           col=c(cols, "white", "white"),
           lty=ltypes,
           legend=group.labels,
           bty="n",
           cex=cexLegend)

    ## reset mar that was changed to allow adding numbers underneath plot
    if (plot.nrisk) {
        par(mar=org.mar)
    }
}

