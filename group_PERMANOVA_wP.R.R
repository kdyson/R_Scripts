
## Function for sequentially testing multiple variables in PERMANOVA, then
## adjusting pseudo-F values. Also calculates multivariate dispersion for
## categorical variables.

# First, a quick function to calculate omega squared based on:
    # https://academic.oup.com/bioinformatics/article/31/15/2461/188732
    # http://www.real-statistics.com/multiple-regression/other-measures-effect-size-anova/ 
    # http://psych.colorado.edu/~willcutt/pdfs/Olejnik_2003.pdf if you want to make this generalized...
    # https://gist.github.com/arnoud999/e677516ed45e9a11817e for some r code of generalization

PERMANOVA.omega2 <- function(adonis2.object, num.control.vars) {
    

    SSE <- adonis2.object[ 1 + num.control.vars, "SumOfSqs"] 
    dfE <- adonis2.object[ 1 + num.control.vars, "Df"]
    MSE <- adonis2.object["Residual", "SumOfSqs"]/adonis2.object["Residual", "Df"]      # Mean Square for Error
    SST <- sum(adonis2.object$SumOfSqs)
    
    omega2 <- (SSE - dfE*MSE) / (SST + MSE)
    
    return(omega2)
    
}



## Note that I have some defaults for adonis2: two processors (parallel = 2),
## 99999 permutations, and using the bray-curtis distance matrix (method =
## "bray")

## Note also that this does not test for interaction effects.

## Variable definitions:
     ## var.names should be a vector of column names, one for each variable to be tested
          # Might consider adding functionality to take column names or position numbers, but for now...
     ## var.table is the table where the independent variables are.
     ## var.table.c is a string of the table name: eventually generate this automatically.
     ## species.table is the species/site matrix (e.g. created by matrify)
     ## species.table.c is a string of the table name: eventually generate this automatically.
     ## control.vars are any variables that need to be included before testing variables--character string, e.g. "var1 + var2"
     ## num.control.vars specifies the number of control variables used. Replace this later with something deriving it from control.vars text...
     ## by.adonsi2 is to pass through for the by = argument in adonis2


group.PERMANOVA <- function(var.names, var.table, var.table.c, control.vars = "", species.table, species.table.c, num.control.vars = 0, by.adonis2 = "terms", perms = 99999) {

     ## need to order columns otherwise the order of the columns and the order of col.numbers is incorrect.

     var.table <- var.table[ , order(colnames(var.table))]
     var.names <- var.names[order(var.names)]



     output <- data.frame(row.names = var.names,
                                         var.explnd = rep(NA, times = length(var.names)),
                                         avg.var.explnd = rep(NA),
                                         pseudo.F = rep(NA),
                                         F.pval = rep(NA),
                                         adj.F.pval = rep(NA),
                                         disp.F.pval = rep(NA),
                                         adj.disp.F.pval = rep(NA)
                              )
     col.numbers <- which(colnames(var.table) %in% var.names)


     ## Attempting to get the names from the tables
     #var.table.t <- as.character(quote(var.table, env = ))
     #species.table.t <- as.character(var.table)


     ## Populate the output table.

     for (i in 1:length(var.names)) {


          temp <- adonis2(formula = as.formula(paste(species.table.c, "~", control.vars, "+", var.names[i])),
                          parallel = 2, permutations = perms, method = "bray", by = by.adonis2, data = var.table)

          output$var.explnd[i] <- temp$SumOfSqs[1 + num.control.vars] / sum(temp$SumOfSqs)
          output$avg.var.explnd[i] <- output$var.explnd[i] / temp$Df[1 + num.control.vars]
          output$omega2[i] <- PERMANOVA.omega2(temp, num.control.vars)
          output$pseudo.F[i] <- temp$F[num.control.vars + 1]
          output$F.pval[i] <- temp$`Pr(>F)`[num.control.vars + 1]
          


               # the ordering is so that col.names[i] works!

               if (is.character(var.table[[col.numbers[i]]]) == TRUE | is.factor(var.table[[col.numbers[i]]]) == TRUE) {

                    temp.disp <- anova(betadisper(vegdist(species.table, method = "bray"),
                                                  as.factor(var.table[[col.numbers[i]]])))

                    output$disp.F.pval[i] <- temp.disp$`Pr(>F)`[1]

               }

     }

     output$adj.F.pval <- p.adjust(p = output$F.pval, method = "holm")
     output$adj.disp.F.pval <- p.adjust(p = output$disp.F.pval, method = "holm")

return(output)
     
     
     
}


## with time, test:
    #             paste0("\"", sample.covariates, "\"" )



## Test Code...

# var.names     <- non.management.vars.f <- c("DissolvedArea", "Proportion_ForeignBorn_B99051e5", "DistrictName",
#                                           "MedianHouseholdIncome_B19013e1", "dissolved_parcel_500m_buffer_impervious_500m_mean")
# var.table     <- sample.data.nf
# var.table.c   <- "sample.data.nf"
# species.table <- sqrt.matrify.OTU.genus.nf
# species.table.c <- "sqrt.matrify.OTU.genus.nf"
# control.vars  <- "vegetation_class + gc.simple"
# num.control.vars <- 2
# by.adonis2    <- "terms"




# test <- group.PERMANOVA(var.names = non.management.vars.f, var.table =  sample.data, var.table.c = "sample.data",
#                 species.table = sqrt.matrify.OTU.genus, species.table.c = "sqrt.matrify.OTU.genus",
#                 control.vars = "sample.data$vegetation_class + sample.data$gc.simple", num.control.vars = 2,
#                 by.adonis2 = "terms")
#
# adonis2(formula = sqrt.matrify.OTU.genus ~ sample.data$vegetation_class + sample.data$gc.simple +
#              sample.data$Proportion_ForeignBorn_B99051e5,
#         parallel = 2, permutations = 99999, method = "bray", by = "terms")
#
# adonis2(formula = as.formula(paste(species.table.c, "~", control.vars, "+ sample.data$", var.names[i])),
#         parallel = 2, permutations = 99999, method = "bray", by = "terms")
