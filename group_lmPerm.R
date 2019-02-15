## group lmPerm

require(lmPerm)
require(AICcmodavg)


group.lmp <- function(group.vars, # should be character string
                      var.table.name, # should be character string
                      dependent.var, #should be character string
                      control.var = NULL,
                      iter = 999999,
                      Ca = .01

) {


    lmp.results <- data.frame(row.names = group.vars,
                                  estimate = rep(NA, times = length(group.vars)),
                                  var.explnd = rep(NA),
                                  # avg.var.explnd = rep(NA),
                                  F.stat.lmp = rep(NA),
                                  prob.pval = rep(NA),
                                  prob.iter = rep(NA),
                                  adj.pval = rep(NA),
                                  AICc.lm = rep(NA)
    )

    if (!is.null(control.var)) {control.var <- paste0(control.var, "+")}


    for (i in 1:length(group.vars)) {

        lmp.temp <- lmp(as.formula(paste0(dependent.var, " ~ ",
                                          control.var,
                                          var.table.name, "$", group.vars[i])),
                        model = TRUE, seqs = TRUE, x = TRUE, y = TRUE,
                        center = FALSE, maxIter = iter, Ca = Ca)
        summary.lmp.temp <- lmPerm::summary.lmp(lmp.temp)
        lm.temp <- lm(as.formula(paste0(dependent.var, " ~ ", var.table.name, "$",
                                               group.vars[i])))


        # add calculated values to table:
        lmp.results$estimate[i] <- ifelse(length(summary.lmp.temp$coefficients[,1]) > 2, NA,
                                              summary.lmp.temp$coefficients[2,1])
        lmp.results$var.explnd[i] <- summary.lmp.temp$adj.r.squared


        lmp.results$F.stat.lmp[i] <- summary.lmp.temp$fstatistic["value"]
        lmp.results$prob.pval[i] <- summary.lmp.temp$coefficients[2,3]
        lmp.results$prob.iter[i] <- summary.lmp.temp$coefficients[1,2]
        lmp.results$AICc.lm[i] <- AICc(lm.temp)[1]


        }

lmp.results$adj.pval <- p.adjust(p = lmp.results$prob.pval, method = "holm")



return(lmp.results)


}

