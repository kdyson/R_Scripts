## group lmPerm

# based on http://www.utstat.toronto.edu/~brunner/workshops/mixed/NormalWithR.pdf#:~:text=Repeated%20measures%20analysis%20with%20R%20Summary%20for%20experienced,to%20the%20model%20for%20the%20random%20subject%20effect.
# with information from https://stats.stackexchange.com/questions/7240/proportion-of-explained-variance-in-a-mixed-effects-model?rq=1



require(lme4) #
require(car)  # car also has vif function for multivariate models...
require(MuMIn)


group.lmer <- function(var.names, # should be character string
                      in.data, 
                      dependent.var, #should be character string
                      control.var = NULL,
                      random.subject.effect = NULL # should be character string in the format 
                                                   # (1 | site)

) {

    # Want to capture MSE, F, Pr(>F), estimate and standard error, AICc
    
    lmer.results <- data.frame(row.names = var.names,
                               var.names = var.names,
                                  slope.estimate = rep(NA, times = length(var.names)),
                                  std.error = rep(NA),
                                  Mean.Sq = rep(NA),
                                  F.stat = rep(NA),
                                  prob.pval = rep(NA),
                                  adj.pval = rep(NA),
                                  AICc = rep(NA),
                                  AIC.REML = rep(NA),
                                  logLik = rep(NA)
    )

    if (!is.null(control.var)) {control.var <- paste0(control.var, "+")}
    if (!is.null(random.subject.effect)) {random.subject.effect <- 
                                                paste0("+ ", random.subject.effect)}
    
i=1
    for (i in 1:length(var.names)) {

        lmer.temp <- lmer(formula = as.formula(paste0(dependent.var, " ~ ",
                                          control.var,
                                          var.names[i],
                                          random.subject.effect
                                          )),
                        data = in.data)
        cT.temp <- MuMIn::coefTable(lmer.temp)
        indx <- (2 + length(control.var))

        # add calculated values to table:
        lmer.results$slope.estimate[i] <- cT.temp[indx, 1]
        lmer.results$std.error[i] <- cT.temp[indx, 2]
        
        lmer.results$Mean.Sq[i] <- anova(lmer.temp)[var.names[i],3]

        lmer.results$F.stat[i] <- car::Anova(lmer.temp, test = "F")[var.names[i],"F"]
        lmer.results$prob.pval[i] <- car::Anova(lmer.temp, test = "F")[var.names[i],"Pr(>F)"]
        
        
        lmer.results$AICc[i] <- AICc(lmer.temp)[1] # this is for the model...
        lmer.results$AIC.REML[i] <- lme4::llikAIC(lmer.temp)[[2]]
        lmer.results$logLik[i] <- lme4::llikAIC(lmer.temp)[[1]]

        }


lmer.results$adj.pval <- p.adjust(p = lmer.results$prob.pval, method = "holm")

return(lmer.results)


}

