## This is a function that runs multipatt a specified number of times, then outputs two key things: 

        # 1. A summary table of the species over 100 runs and
        # 2. A plot showing the different indicator species made in ggplot (OPTIONAL)

# Reminder to user: check that matrix site order and cluster order are the same

# repeat.multipatt(matrix.name = matrify.shrub.bysite, cluster.name = gc.cluster$gc.group.3, plot.please = FALSE, freq.cutoff = .1)
# repeat.multipatt(matrix.name = matrify.tree.species, cluster.name = tree.abs.cluster$tree.abs, func.name = "r.g")
# repeat.multipatt(matrix.name = matrify.foraging,
#                  cluster.name = shrub.abs.cluster$shrub.abs,
#                  xlab.input = "Bird Incidence Indicator Species \nof shrub clusters")




repeat.multipatt <- function(repeats = 100, matrix.name, cluster.name, p.cutoff = .05, func.name = "IndVal.g", phi = FALSE, plot.please = TRUE, plot.colors = c("deeppink2", "cyan3", "yellowgreen", "gray 50", "violet"), freq.cutoff = .5, xlab.input = "Indicator Species \n(freq > 50%)", ylab.input = "Mean Indicator Value over 100 Runs", quiet = TRUE, stat.cutoff = .75) {

    require(ggplot2, dplyr, indicspecies)
    
    
if (freq.cutoff != .5) {
    
    print("You probably want to change the default Y axis label!")
    
}
    
    if (length(unique(cluster.name)) > 5) {
        
        stop("Please use fewer than 5 clusters (or edit the function).")
    }
    
    if (func.name == "r.g" & phi == TRUE & quiet == FALSE) {
        
            matrix.name <- as.data.frame(ifelse(matrix.name > 0, 1, 0))
            print("you are using the r.g function to calculate Pearson's phi coefficient of association.")
    }
        
    if (func.name == "r.g" & phi == FALSE & quiet == FALSE) {
        
        print("you are using the r.g function to calculate the point biserial correlation coefficient.")
    }


mp.sign.dump <- data.frame()
mp.AB.dump <- data.frame()
i <- 1

for (i in 1:repeats) {
    multipatt.out <- multipatt(x = matrix.name,
                              cluster = cluster.name, func = func.name)
    
    mp.sign <- multipatt.out$sign[complete.cases(multipatt.out$sign),]
    mp.sign <- mp.sign[mp.sign$p.value <= p.cutoff, ]
    mp.sign$species <- rownames(mp.sign)
    
        if ( length(multipatt.out$A) > 0 ) {
            mp.A <- as.data.frame(multipatt.out$A[complete.cases(multipatt.out$A), ] , make.row.names = FALSE)
            mp.A <- mp.A[which(rownames(mp.A) %in% mp.sign$species), ]
                colnames(mp.A) <- paste("A.", colnames(mp.A), sep = "")
            
                
            mp.B <- as.data.frame(multipatt.out$B[complete.cases(multipatt.out$B), ], make.row.names = FALSE)
            mp.B <- mp.B[which(rownames(mp.B) %in% mp.sign$species), ]
                colnames(mp.B) <- paste("B.", colnames(mp.B), sep = "")
            
            mp.AB <- cbind.data.frame(mp.A, mp.B)
            mp.AB$species <- rownames(mp.AB)
            
                if (nrow(mp.AB) > 0) {
                    mp.AB$i <- i
                    }
        
            }
    
        if (nrow(mp.sign) > 0) {
            mp.sign$i <- i
            }
    
    mp.sign.dump <- rbind(mp.sign.dump, mp.sign, make.row.names = FALSE)
    
    
    if (func.name %in% c("IndVal", "IndVal.g")) {
            mp.AB.dump <- rbind(mp.AB.dump, mp.AB, make.row.names = FALSE)
        }    
    
    i <- i + 1
}



if (!(any(colnames(mp.sign.dump) == "s.3"))) {mp.sign.dump$s.3 = rep(0, times = length(mp.sign.dump$s.1))}
if (!(any(colnames(mp.sign.dump) == "s.4"))) {mp.sign.dump$s.4 = rep(0, times = length(mp.sign.dump$s.1))}
if (!(any(colnames(mp.sign.dump) == "s.5"))) {mp.sign.dump$s.5 = rep(0, times = length(mp.sign.dump$s.1))}


# Make summary tables

mp.summary <- group_by(mp.sign.dump, species) %>%
                summarize(
                           count.sp = n(),
                           frequency.sp = round(count.sp/repeats, digits = 3),
                           mean.stat = round(mean(stat), digits = 3),
                           group = paste(
                                          (if (any(s.1 == 1)) {"1."}), 
                                          (if (any(s.2 == 1)) {".2."}),
                                          (if (any(s.3 == 1)) {".3."}),
                                          (if (any(s.4 == 1)) {".4."}),
                                          (if (any(s.5 == 1)) {".5"}), 
                                                 sep = ""),
                           min.p.val = min(p.value),
                           max.p.val = max(p.value),
                           all.p.vals = paste(p.value, collapse = ", ")
                           )

mp.summary$group <- ifelse(grepl("[1-9]..[1-9]", mp.summary$group),
                                gsub(pattern = "\\.\\.", replacement = " & ", mp.summary$group), mp.summary$group)
mp.summary$group <- paste("Cluster", gsub(pattern = "\\.", replacement = "", mp.summary$group), sep = " ")




if (func.name %in% c("IndVal", "IndVal.g")) {
        mp.AB.summary  <- group_by(mp.AB.dump, species) %>%
                            summarize_each(funs(mean))
            mp.AB.summary$i <- NULL
            colnames(mp.AB.summary)[-1] <- paste(colnames(mp.AB.summary[-1]), "mean", sep = ".")
        
        mp.AB.summary <- cbind(mp.AB.summary, 
                               group_by(mp.AB.dump, species) %>%
                                summarize(
                                    count.sp = n(),
                                    frequency.sp = round(count.sp/max(i), digits = 3)
                                )
                            )
}       
                        

## Output graph
                           
if (plot.please == TRUE) {

    if (nrow(subset(mp.summary, frequency.sp > freq.cutoff)) > 0) {
        
        mp.plot <- ggplot(subset(mp.summary, frequency.sp > freq.cutoff & mean.stat > stat.cutoff), 
                          aes(x = reorder(species, mean.stat))) +
                geom_bar(aes(y = mean.stat, fill = group), alpha = .75, stat = "identity") +
                geom_text(aes(y = mean.stat, label = mean.stat), hjust = 1.2) +
                facet_grid(group ~ ., scales = "free_x") + coord_flip() + guides(fill = "none") +
                xlab(xlab.input) + ylab(ylab.input) +
                scale_fill_manual(values = c(plot.colors))

        print(mp.plot)

    }
    
    else print("can't make graph; no indicator species for any clusters!")
}




# Output data tibble
if (quiet == FALSE) {
if (func.name %in% c("IndVal", "IndVal.g")) {
    return(list(mp.summary, mp.AB.summary))
    }
else return(mp.summary)
}

}
