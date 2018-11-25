## Create AICc tables for PERMANOVA output.
require(vegan)
require(tibble)

## -- For two variables: ------------------------------------

AICc.table.2var <- function(sig.vars, control.var.char = NULL, c.var = 0, matrix.char, perm = 999, type = "AICc", method = "bray") {

    varcomb.2.AICc <- tibble(variables = rep("var.name", choose(length(sig.vars),2)),
                             AICc.values = rep(0),
                             `Pseudo-_F_` = rep(0),
                             `p-value` = rep(0),
                             `Var Explnd` = rep(0),
                             Model = rep("model"))

    if (is.character(control.var.char) == TRUE & c.var == 0) {c.var = 1} 
        
    combo.list <- combn(x = sig.vars, m = 2, simplify = FALSE)

    if (!is.null(control.var.char)) {
        control.var.char <- paste0(control.var.char, " +")
        }

    for (r in 1:choose(length(sig.vars),2)) {
        # label the row with variable names
        varcomb.2.AICc[r,1] <- paste(control.var.char, paste(combo.list[[r]], collapse = " + "))
        
        # create a temporary PERMANOVA to take info from
        
        temp <- adonis2(
            as.formula(paste0(
                matrix.char,
                " ~ ",
                control.var.char,
                paste0(combo.list[[r]], collapse =  "+")
            )),
            permutations = perm,
            method = method,
            by = NULL
        )
        
        
        varcomb.2.AICc[r, 2] <- AICc.PERMANOVA2(temp)[type]
        varcomb.2.AICc[r, 3] <- temp$F[1]
        varcomb.2.AICc[r, 4] <- temp$`Pr(>F)`[1]
        varcomb.2.AICc$`Var Explnd`[r] <- temp$SumOfSqs[1] / temp$SumOfSqs[3]
        

        r <- r + 1

    }

    varcomb.2.AICc$`Delta AICc` <- varcomb.2.AICc$AICc.values -
        min(varcomb.2.AICc$AICc.values)
    varcomb.2.AICc$`Relative Likelihood` <- exp((min(varcomb.2.AICc$AICc.values) -
                                                     varcomb.2.AICc$AICc.values)/2)

    # Relative likelihood compared with best model; see
    # https://en.wikipedia.org/wiki/Likelihood_function
    
    
return(varcomb.2.AICc)

}



## -- For N variables: ---------------------------------------------------


AICc.table.Nvar <- function(sig.vars, control.var.char = NULL, c.var = 0, matrix.char, perm = 999, n.var = 1, composite = FALSE, type = "AICc", method = "bray") {

    if (n.var > length(sig.vars)) { stop("n.var greater than number of significant variables")}

    if (!is.null(control.var.char)) {
        control.var.char <- paste0(control.var.char, " + ")
    }
    if (is.character(control.var.char) == TRUE & c.var == 0) {c.var = 1} 
    

    varcomb.N.AICc <- tibble(variables = rep("var.name", choose(length(sig.vars), n.var)),
                             AICc.values = rep(0),
                             `Pseudo-_F_` = rep(0),
                             `p-value` = rep(0),
                             `Var Explnd` = rep(0),
                             Model = rep("model"))

    combo.list <- combn(x = sig.vars, m = n.var, simplify = FALSE)


    for (r in 1:choose(length(sig.vars), n.var)) {
        # label the row with variable names
        varcomb.N.AICc[r,1] <- paste(control.var.char, paste(combo.list[[r]], collapse = " and "))

        # create a temporary PERMANOVA to take info from
        temp <- adonis2(
            as.formula(paste0(
                matrix.char,
                " ~ ",
                control.var.char,
                paste0(combo.list[[r]], collapse =  "+")
            )),
            permutations = perm,
            method = method,
            by = NULL
        )
        
        
        varcomb.N.AICc[r,2] <- AICc.PERMANOVA2(temp)[type]

        varcomb.N.AICc[r,3] <- temp$`F`[1]
        varcomb.N.AICc[r,4] <- temp$`Pr(>F)`[1]
        varcomb.N.AICc$`Var Explnd`[r] <- temp$SumOfSqs[1] / temp$SumOfSqs[3]
        

        r <- r + 1

    }


    # Calculate diagnostic variables:

if (composite == FALSE) {
    varcomb.N.AICc$`Delta AICc` <- varcomb.N.AICc$AICc.values -
        min(varcomb.N.AICc$AICc.values)
    varcomb.N.AICc$`Relative Likelihood` <- exp((min(varcomb.N.AICc$AICc.values) -
                                                     varcomb.N.AICc$AICc.values)/2)
    }


    return(varcomb.N.AICc)

}


## -- Wrapper function: ---------------------------------------------------

# comb.incl is which # of variables for the combinations you want to include.
# e.g. all combinations of 2 variables, 3 variables, etc.

AICc.table.all <- function(sig.vars, control.var.char = NULL, matrix.char, perm = 999, comb.incl = 1, extra.var = FALSE, extra.var.char = NULL, type = "AICc", method = "bray") {

    varcomb.all <- data.frame()

    # If there is a control variable, create a one variable model with only
    # control variable, for comparison with rest of proposed models
    # 
    if (!is.null(control.var.char)) {

        temp <- AICc.table.Nvar(sig.vars = control.var.char, control.var.char = NULL,
                                matrix.char = matrix.char, n.var = 1, composite = TRUE,
                                type = type, method = method)

        varcomb.all <- rbind(varcomb.all, temp)


    }

    # Iterate through the comb.incl. e.g. all 1 var models, then 2 var models, then 3 var models...
    for (i in comb.incl) {
    
        temp <- AICc.table.Nvar(sig.vars = sig.vars, control.var.char = control.var.char,
                                matrix.char = matrix.char, n.var = i, composite = TRUE,
                                type = type, method = method)
    
        varcomb.all <- rbind(varcomb.all, temp)
    
    
    }

    
    
    
    # If you want to include a non-significant variable for comparison...
    if (extra.var == TRUE) {
        
        for (i in 1:length(extra.var.char)) {
        temp <- AICc.table.Nvar(sig.vars = extra.var.char[i], control.var.char = control.var.char,
                                matrix.char = matrix.char, n.var = 1, composite = TRUE, 
                                type = type, method = method)
        
        varcomb.all <- rbind(varcomb.all, temp)
        
        }
    }
    
    
    varcomb.all$`Delta AICc` <- varcomb.all$AICc.values -
        min(varcomb.all$AICc.values)
    varcomb.all$`Relative Likelihood` <- exp((min(varcomb.all$AICc.values) -
                                                  varcomb.all$AICc.values)/2)



    return(varcomb.all)

}







## Testing
# sig.vars <- sigvars.symbionfNGS[-1]
# matrix.char <- "sqrt(symbio.transpose.nf)"
# control.var.char <- sigvars.symbionfNGS[1]
#
#
# AICc.table.2var(sig.vars = sigvars.symbionfNGS[-1], control.var.char = sigvars.symbionfNGS[1], matrix.char = "sqrt(symbio.transpose.nf)")[3]
#
# AICc.table.Nvar(sig.vars = sigvars.symbionfNGS[-1], control.var.char = sigvars.symbionfNGS[1], matrix.char = "sqrt(symbio.transpose.nf)", n.var = 1)[3]
#
# test <- AICc.table.all(sig.vars = sigvars.symbionfNGS[-1], control.var.char = sigvars.symbionfNGS[1],
#                        matrix.char = "sqrt(symbio.transpose.nf)", comb.incl = 1)
#

