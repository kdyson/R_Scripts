# Function to calculate AICc for PERMANOVA. Requires input from adonis {vegan}


AICc.PERMANOVA <- function(adonis.model) {
    
    # check to see if object is an adonis model...
    
    if (!(adonis.model$aov.tab[1,1] >= 1))
        stop("object not output of adonis {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model
    # Calculating AICc using residual sum of squares (RSS) since I don't think that adonis returns something I can use as a liklihood function...
    
    RSS <- adonis.model$aov.tab[which(rownames(adonis.model$aov.tab) == "Residuals"), 2]
    
    k <- ncol(adonis.model$model.matrix)
    
    nn <- nrow(adonis.model$model.matrix)
    
    # AIC : 2*k + n*ln(RSS/n)
    # AICc: AIC + [2k(k+1)]/(n-k-1)

    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
    
    AIC <- 2*k + (nn*log(RSS/nn))
    AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
    
    output <- list("AIC" = AIC, "AICc" = AICc)
    
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

