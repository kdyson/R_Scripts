# Function to calculate AICc for PERMANOVA. Requires input from adonis {vegan}


AICc.PERMANOVA <- function(adonis.model) {
    
    # check to see if object is an adonis model...
    
    if (!(adonis.model$aov.tab[1,1] >= 1))
        stop("object not output of adonis {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model
    # Calculating AICc using residual sum of squares (RSS) since I don't think that adonis returns something I can use as a liklihood function...
    
    RSS <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "SumsOfSqs"]
    MSE <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "MeanSqs"]
    
    k <- ncol(adonis.model$model.matrix)# + 1 # add one for error variance
    
    nn <- nrow(adonis.model$model.matrix)
    
    # AIC : 2*k + n*ln(RSS)
    # AICc: AIC + [2k(k+1)]/(n-k-1)

    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
    
    # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
    # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
    
    AIC <- 2*k + nn*log(RSS)
    AIC.g <- 2*k + nn * (1 + log( 2 * pi * RSS / nn))
    AIC.MSE <- 2*k + nn * log(MSE)
    AIC.pi <- k + nn*(1 + log( 2*pi*RSS/(nn-k) )   )
    AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
    AICc.MSE <- AIC.MSE + (2*k*(k + 1))/(nn - k - 1)
    AICc.pi <- AIC.pi + (2*k*(k + 1))/(nn - k - 1)
    
    output <- list("AIC" = AIC, "AIC.g" = AIC.g, "AICc" = AICc,
                   "AIC.MSE" = AIC.MSE, "AICc.MSE" = AICc.MSE,
                   "AIC.pi" = AIC.pi, "AICc.pi" = AICc.pi, "k" = k)
    
    return(output)   

}

AICc.PERMANOVA2 <- function(adonis2.model) {
    
    # check to see if object is an adonis2 model...
    
    if (is.na(adonis2.model$SumOfSqs[1]))
        stop("object not output of adonis2 {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model Calculating AICc
    # using residual sum of squares (RSS or SSE) since I don't think that adonis
    # returns something I can use as a liklihood function... o wait maximum likelihood may be MSE.
    
    RSS <- adonis2.model$SumOfSqs[ length(adonis2.model$SumOfSqs) - 1 ]
    MSE <- RSS / adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    nn <- adonis2.model$Df[ length(adonis2.model$Df) ] + 1
    
    k <- nn - adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    
    # AIC : 2*k + n*ln(RSS)
    # AICc: AIC + [2k(k+1)]/(n-k-1)
    
    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
    
    # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
    # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
    
    AIC <- 2*k + nn*log(RSS)
    AIC.g <- 2*k + nn * (1 + log( 2 * pi * RSS / nn))
    AIC.MSE <- 2*k + nn * log(MSE)
    AIC.pi <- k + nn*(1 + log( 2*pi*RSS/(nn-k) )   )
    AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
    AICc.MSE <- AIC.MSE + (2*k*(k + 1))/(nn - k - 1)
    AICc.pi <- AIC.pi + (2*k*(k + 1))/(nn - k - 1)
    
    output <- list("AIC" = AIC, "AIC.g" = AIC.g, "AICc" = AICc,
                   "AIC.MSE" = AIC.MSE, "AICc.MSE" = AICc.MSE,
                   "AIC.pi" = AIC.pi, "AICc.pi" = AICc.pi, "k" = k)
    
    return(output)   
    
}


## Test code:

# AICc.PERMANOVA(size.age.PERMANOVA)
# 
# adonis.model <- size.qual.PERMANOVA
# 
# size.qual.PERMANOVA <- adonis(formula = matrify.shrub.bysite ~ as.factor(sample.covariates$BldgQuality_lookup),
#                               permutations = 99999, method = "bray")
# size.qual.PERMANOVA$aov.tab[2,3]
# size.qual.PERMANOVA$model.matrix

