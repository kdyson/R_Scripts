
## ----- Data manip. funs. | transpose ---------------------------

speciesSite <- function(inputOTU, rare = 0) {
    #input OTU should only be data, no text/ID fields.
    # inputs Species as columns, outputs Species as rows.
    temp <- inputOTU[ , colSums(inputOTU) > rare]
    temp <- as.data.frame(inputOTU)
    temp.col.names <- colnames(temp)
    temp <- as.data.frame(t(temp))
    rownames(temp) <- temp.col.names
    
    return(temp)
}

siteSpecies <- function(inputOTU, rowNames, rare = 0) {
    #input OTU should only be data, no text/descriptions.
    #inputs Species as rows, outputs Species as columns
    temp <- as.data.frame(inputOTU, row.names = rowNames)
    temp <- as.data.frame(t(temp))
    colnames(temp) <- rowNames
    temp <- temp[rowSums(temp) > rare,]
    
    return(temp)
    
}


## ----- Data manip. funs. | Compositional Matrices ----------------

compMatrix <- function(inputMatrix, z.warning = 0.8){
    # square-root Bayesian-multiplicative replacement of zeros with the cmultRepl()
    # function (Ladin et al., 2021) 
    
    # Most other code in GitHub uses CZM. e.g. See
    # https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-Biplot
    # and https://raw.githubusercontent.com/ggloor/CoDaSeq/6ff864aade46cd3c8b0eff3bb54d5460775f92cd/CoDaSeq/vignettes/CoDaSeq_vignette.Rnw
    # This latter contends that this is the most principled method.
    
    zeroRepl <- cmultRepl(inputMatrix,
                          label = 0,
                          method = "CZM",
                          z.warning = z.warning)
    
    output <- cdt.acomp(x = zeroRepl) %>% 
        as.data.frame()
    #rownames(output) <- rownames(inputMatrix)
    
    ## see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7811025/
    
    return(output)
}


## ----- Data Summary | Reads ------------------------------------------
# Creates a summary table for number of eDNA reads found; total number of species variants, etc. 
# note this is customized for Eco Mol output tables. Column names will differ for other providers.


ecoMolSummary <- function(ecoMolInput) {
    output <- tibble(
        sumReads = dplyr::select(ecoMolInput, asvAbsoluteAbundance) %>% sum(),
        ASVCount = dplyr::select(ecoMolInput, ASVHeader) %>% unique() %>% nrow(),
        genusCount = dplyr::select(ecoMolInput, genusBLASTn) %>% unique() %>% nrow()
    )
    
    return (output)
}



## ----- Data Summary | Alpha Diversity --------------------------

# Note these two are based on the VEGAN function, but with modifications to preserve data order. Per this issue I created, changing the data order is the desired behavior and they won't fix it, so this is my workaround: https://github.com/vegandevs/vegan/issues/552

`specnumberMOD` <-
    function(x, groups, MARGIN = 1)
    {
        if (!missing(groups)) {
            if (length(groups) == 1)
                groups <- rep(groups, nrow(x))
            groups <- factor(groups, levels = unique(groups)) # this preserves the data order!
            x <- aggregate(x, list(groups), max) # max is used because the actual number doesn't matter, just that it is over 0
            rownames(x) <- x[,1] 
            x <- x[,-1]
        }
        if (length(dim(x)) > 1)
            apply(x > 0, MARGIN, sum)
        else
            sum(x > 0)
    }

`diversityMOD` <- function (x, index = "shannon", groups, equalize.groups = FALSE, 
                            MARGIN = 1, base = exp(1)) 
{
    x <- drop(as.matrix(x))
    if (!is.numeric(x)) 
        stop("input data must be numeric")
    if (any(x < 0, na.rm = TRUE)) 
        stop("input data must be non-negative")
    if (!missing(groups)) {
        if (MARGIN == 2) 
            x <- t(x)
        if (length(groups) == 1) 
            groups <- rep(groups, NROW(x))
        if (equalize.groups) 
            x <- decostand(x, "total")
        groups <- factor(groups, levels = unique(groups)) # this preserves the data order!
        x <- aggregate(x, list(groups), sum)
        rownames(x) <- x[, 1]
        x <- x[, -1, drop = FALSE]
        if (MARGIN == 2) 
            x <- t(x)
    }
    INDICES <- c("shannon", "simpson", "invsimpson")
    index <- match.arg(index, INDICES)
    if (length(dim(x)) > 1) {
        total <- apply(x, MARGIN, sum)
        x <- sweep(x, MARGIN, total, "/")
    }
    else {
        x <- x/(total <- sum(x))
    }
    if (index == "shannon") 
        x <- -x * log(x, base)
    else x <- x * x
    if (length(dim(x)) > 1) 
        H <- apply(x, MARGIN, sum, na.rm = TRUE)
    else H <- sum(x, na.rm = TRUE)
    if (index == "simpson") 
        H <- 1 - H
    else if (index == "invsimpson") 
        H <- 1/H
    if (any(NAS <- is.na(total))) 
        H[NAS] <- NA
    H
}




alphaMetrics <- function(inputOTUSiSp, groupNames, replNames){
    
    alphaTable <- tibble(
        siteNames = rownames(inputOTUSiSp),
        siteType = groupNames,
        replicate = replNames,
        speciesRichness = specnumberMOD(inputOTUSiSp , MARGIN = 1),
        shannonRichness = vegan::diversity(inputOTUSiSp, index = "shannon"),
        effectiveSR = exp(shannonRichness),
        invSimpson = vegan::diversity(inputOTUSiSp, index = "invsimpson")
    )
    
    temp <- alphaTable %>% group_by(siteType) %>% summarise(siteType_n = n())
    
    alphaTable <- left_join(x = alphaTable, y = temp)
    
    return(alphaTable)
    
}

alphaGroupMetrics <- function(inputOTUSiSp, groupNames) {
    alphaTable <- tibble(
        siteType = groupNames %>% unique(), # unique preserves order of groupNames, no alph reordering
        speciesRichness = specnumberMOD(inputOTUSiSp, MARGIN = 1, groups = groupNames),
        shannonRichness = diversityMOD(inputOTUSiSp, index = "shannon", groups = groupNames),
        effectiveSR = exp(shannonRichness),
        invSimpson = diversityMOD(inputOTUSiSp, index = "invsimpson", groups = groupNames)
    )
    
    temp <- as_tibble(groupNames) %>% group_by(value) %>% summarise(siteType_n = n())
    
    alphaTable <- left_join(x = alphaTable, y = temp, join_by(siteType == value))
    
    return(alphaTable)
    
}


## ----- Plotting functions | Aitchison -------------------------------
# helper function--note this is taken from other sources. I couldn't identify the original source of the code but e.g. it is used here: http://sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

get_lower_tri <- function(inpMatrix){
    inpMatrix[upper.tri(inpMatrix, diag = T)]<- NA
    return(inpMatrix)
}



aitHeatmap <- function(inputDist, fillColor1 = "blue", fillColor2 = "orange", textPlease = FALSE){
    graph <-
        inputDist %>%
        as.matrix() %>%
        get_lower_tri() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("plot1") %>%
        pivot_longer(-c(plot1),
                     names_to = "plot2",
                     values_to = "distance",
                     values_drop_na = T) %>%
        ggplot(aes(x = plot1, y = plot2, fill = distance)) + 
        geom_raster() +
        
        scale_fill_gradient(low = fillColor1, high = fillColor2,  
                            name="Aitchison\nDistance") +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            #panel.border = element_blank(),
            panel.background = element_blank(),
            #axis.ticks = element_blank(),
            legend.justification = c(1, 0),
            legend.position = c(0.5, 0.7),
            legend.direction = "horizontal",
            axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1),
            axis.text.y = element_text(size = 12))+
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                     title.position = "top", title.hjust = 0.5))
    
    if(textPlease == TRUE){
        graph <- graph +
            geom_text(aes(label = round(distance)))
        
        
    }
    
    
    return(graph)   
}


# Create function that takes distance matrix, switches to long form, then
# creates relevant columns. Then can feed it into a box plot.

aitComparison <- function(inputDist, remap = NULL, repeatSamples = FALSE, fillColor = NULL, levelsPlot = NULL, plotPlease = TRUE) {
    # remap should have three columns: the original site names, pretty site names,
    # and the type of site (for grouping)
    temp <-
        inputDist %>%
        as.matrix() %>%
        get_lower_tri() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("plot1") %>%
        pivot_longer(
            -c(plot1),
            names_to = "plot2",
            values_to = "distance",
            values_drop_na = T
        )
    
    if(!is.null(remap)) {
        # make this a remap so that the names aren't so awful. need to add else to handle if remap is NULL
        old <- remap[[1]]
        new <- remap[[2]]
        temp$plot1[ temp$plot1 %in% old] <- new[base::match(temp$plot1, old)]
        temp$plot2[ temp$plot2 %in% old] <- new[match(temp$plot2, old)]
        
        if(repeatSamples == TRUE) {
            type <- remap[[3]]
            temp$type1 <- type[match(temp$plot1, new)]
            temp$type2 <- type[match(temp$plot2, new)]
            
        }
        
    }
    
    if(repeatSamples == FALSE){
        temp <- temp %>%
            mutate(pair = paste0(plot1, "-", plot2)) 
    }
    
    if (repeatSamples == TRUE) {
        temp <- temp %>%
            mutate(pair = ifelse(
                type1 < type2,
                paste(type1, type2, sep = "-"),
                paste(type2, type1, sep = "-")
            ))
    } 
    
    print(temp$pair)
    
    if(is.null(fillColor) == TRUE) {
        library(RColorBrewer)
        fillColor <- brewer.pal(length(unique(temp$pair)),"Set1")
    }
    
    if(is.null(levelsPlot) == TRUE){
        levelsPlot <- sort(unique(temp$pair))
    }
    
    
    #make the graph here.
    if(repeatSamples == TRUE & plotPlease == TRUE){
        temp <- temp %>%
            mutate(pair = factor(pair, levels = levelsPlot)) %>%
            ggplot(aes(y = distance, x = pair)) +
            stat_boxplot(geom = "errorbar",
                         width = 0.25) +
            geom_boxplot() +
            geom_jitter(aes(color = pair), width = 0.05, size = 3) +
            scale_x_discrete(limits = levelsPlot) +
            theme(legend.position = "none",
                  plot.margin = margin(t = .5,  # Top margin
                                       r = .5,  # Right margin
                                       b = .5,  # Bottom margin
                                       l = 1.5,  # Left margin
                                       unit = "cm"
                  )) +
            theme(axis.text.x = element_text(
                angle = 45,
                vjust = 1,
                size = 10,
                hjust = 1
            )) +
            scale_color_manual(values = fillColor) +
            labs(x = element_blank(),
                 y = "Aitchison Distance") 
    }
    
    
    
    if(repeatSamples == FALSE & plotPlease == TRUE){
        temp <- temp %>%
            ggplot(aes(y = distance, x = pair)) +
            geom_bar(aes(y = distance, fill = pair), stat = "identity") +
            geom_text(aes(label = round(distance,1)), vjust = -0.5) +
            scale_x_discrete(limits = levelsPlot) +
            theme(legend.position = "none") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                             size = 10, hjust = 1)) +
            scale_fill_manual(values = fillColor) +
            labs(x = element_blank(),
                 y = "Aitchison Distance"
            )
    }
    
    
    return(temp)
}




