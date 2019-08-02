require(dplyr)
require(tibble)
require(data.table)

## --- Make some nice tables for each TITAN ------------------------------

## For use with basic site x species/genus table, e.g. with mushrooms collected in the field.

extract.titan.taxa <-
    function(titan.out,
             purity.cutoff = 0.9,
             reliability.cutoff = 0.9,
             taxa.label.name = "Genus") {
        titan.out <- titan.out
        
        titan.out.filtered <- titan.out$sppmax %>%
            as.data.frame() %>%
            rownames_to_column()
        
        titan.out.filtered <-
            titan.out.filtered[titan.out.filtered$purity >= purity.cutoff &
                                   titan.out.filtered$reliability >= reliability.cutoff, ]
        titan.out.filtered <-
            titan.out.filtered[order(titan.out.filtered$zscore, decreasing = T), ] %>%
            mutate(zgrp = ifelse(maxgrp == 1, "z-", "z+"))
        
        
        titan.out.taxonomy.summary <-
            group_by(titan.out.filtered, zgrp, rowname) %>%
            summarise(
                mean.zscore = mean(zscore),
                mean.purity = mean(purity),
                mean.reliability = mean(reliability),
                mean.zenv.cp = mean(zenv.cp),
                mean.5pct.cp = mean(`5%`),
                mean.95pct.cp = mean(`95%`)
            )
        
        
        colnames(titan.out.taxonomy.summary)[1:ncol(titan.out.taxonomy.summary)] <-
            c(
                "Decreasing/Increasing Z Taxa",
                taxa.label.name,
                "Mean Z Score",
                "Mean Purity",
                "Mean Reliability",
                "Mean Env CP z-max",
                "Mean 5% CP",
                "Mean 95% CP"
            )
        
        return(ungroup(titan.out.taxonomy.summary))
        
    }




## For use with NGS data, could be made to work with phyloseq etc.

# titan.out = I5CL.titan.pctN
# taxonomy.table = raw_dada2_taxa


extract.titan.taxa.NGS <-
    function(titan.out,
             taxonomy.table,
             taxonomy.table.merge = "taxa.label.unique",
             taxa.level = Genus,
             purity.cutoff = 0.90,
             reliability.cutoff = 0.90,
             taxa.label.name = "taxa.label.unique",
             label = label) {
        titan.out <- titan.out
        taxa.level <- enquo(taxa.level)
        label <- enquo(label)
        taxonomy.table <- taxonomy.table
        
        titan.out.filtered <- titan.out$sppmax %>%
            as.data.frame()
        titan.out.filtered[, taxa.label.name] <-
            rownames(titan.out.filtered)
        titan.out.filtered <-
            titan.out.filtered[titan.out.filtered$purity >= purity.cutoff &
                                   titan.out.filtered$reliability >= reliability.cutoff, ]
        titan.out.filtered <-
            titan.out.filtered[order(titan.out.filtered$zscore, decreasing = T), ]
        titan.out.filtered <-
            merge(
                titan.out.filtered,
                taxonomy.table,
                by.x = taxa.label.name,
                by.y = taxonomy.table.merge,
                sort = F
            ) %>%
            mutate(zgrp = ifelse(maxgrp == 1, "z-", "z+"))
        
        
        titan.out.taxonomy.summary <-
            group_by(titan.out.filtered, zgrp, !!taxa.level) %>%
            summarise(
                mean.zscore = mean(zscore),
                mean.purity = mean(purity),
                mean.reliability = mean(reliability),
                mean.zenv.cp = mean(zenv.cp),
                mean.5pct.cp = mean(`5%`),
                mean.95pct.cp = mean(`95%`),
                count = n(),
                svs = paste(unique(!!label), collapse = ";")
            )
        
        colnames(titan.out.taxonomy.summary)[1] <-
            "Decreasing/Increasing Z Taxa"
        colnames(titan.out.taxonomy.summary)[3:ncol(titan.out.taxonomy.summary)] <-
            c(
                "Mean Z Score",
                "Mean Purity",
                "Mean Reliability",
                "Mean Env CP z-max",
                "Mean 5% CP",
                "Mean 95% CP",
                "Count",
                "All SVs or OTUs"
            )
        
        return(ungroup(titan.out.taxonomy.summary))
        
    }





## -- For cell filler for table assembler ---------------------

# this needs to fill cells for EITHER z1 or z2 (z- or z+ taxa)

# taxa.col.name <- quo(Genus)
# titan.taxa.tables <- list(NGS.TITAN.cnd.Genus, NGS.TITAN.sesr.Genus)
# table.col.names <- c("cnd", "sesr")
# z.score.pct.cutoff <- .5
# round.val <- 3
# which.z = "z-"
# taxon.list = all.z1.taxon
# j = 3


cell.filler <-
    function(titan.taxa.tables,
             which.z,
             taxon.list,
             taxa.col.name,
             table.col.names,
             z.score.pct.cutoff,
             round.val) {
        # create the table framework...
        
        taxa.comparison.table <- data.table(`All Taxon` = taxon.list,
                                            `Z type` = rep(which.z, length(taxon.list)))
        
        
        # ...and create the columns for the data to go into.
        
        for (i in 1:length(titan.taxa.tables)) {
            taxa.comparison.table[, table.col.names[i]] <-
                rep("", length(taxon.list))
            
        }
        
        
        
        # Now fill everything.
        
        for (i in 1:length(titan.taxa.tables)) {
            temp.table <- titan.taxa.tables[[i]]
            temp.table <-
                temp.table[temp.table$`Decreasing/Increasing Z Taxa` == which.z ,]
            
            # create a cutoff based on user input.
            z.score.percentile <-
                quantile(temp.table$`Mean Z Score`, z.score.pct.cutoff)
            
            
            
            # iterate through each taxon in your taxon.list...
            for (j in 1:length(taxon.list)) {
                check.list <- unlist(select(temp.table, !!taxa.col.name))
                
                # is this taxon in the titan output of the ith table?
                if (taxon.list[j] %in% check.list) {
                    # if yes, then give it some information like if it is a z- or z+ taxa, z score, change point
                    
                    # isolate the info for the right genus
                    temp.row <-
                        temp.table[select(temp.table, !!taxa.col.name) == taxon.list[j] &
                                       !(is.na(select(
                                           temp.table, !!taxa.col.name
                                       ))), ]
                    
                    # bold the z score if it's above the cutoff
                    z.score.text <-
                        ifelse(
                            temp.row[, "Mean Z Score"] >= z.score.percentile ,
                            yes = paste0(
                                "**",
                                temp.row[, "Mean Z Score"] %>% round(digits = round.val),
                                "**"
                            ),
                            no = paste0(temp.row[, "Mean Z Score"])
                        )
                    
                    text <-
                        paste0("Z-score: ",
                               z.score.text,
                               "; CP: ",
                               temp.row[, "Mean Env CP z-max"] %>% round(round.val))
                    
                    taxa.comparison.table[j , i + 2] <- text
                    
                    
                }
                
                # if not R will just leave it blank.
                
                
            }
        }
        
        return(taxa.comparison.table)
        
    }

## --- Multivariable summary table assembly ----------------------

combine.titan.results <-
    function(titan.taxa.tables,
             taxa.col.name,
             table.col.names,
             z.score.pct.cutoff = .5,
             round.val = 3) {
        taxa.col.name <- enquo(taxa.col.name)
        
        all.z1.taxon <- c()
        all.z2.taxon <- c()
        
        
        # create z- and z+ separately b/c some taxa will have both (esp. for higher taxonomic levels.)
        
        # create a list of all taxonomic groups for z- taxon
        for (i in 1:length(titan.taxa.tables)) {
            temp.z1.table <-
                titan.taxa.tables[[i]] %>% subset(`Decreasing/Increasing Z Taxa` == "z-")
            
            all.z1.taxon <-
                c(all.z1.taxon,
                  select(temp.z1.table, !!taxa.col.name) %>% unlist())
            
        }
        
        # create a list of all taxonomic groups for z+ taxon
        for (i in 1:length(titan.taxa.tables)) {
            temp.z2.table <-
                titan.taxa.tables[[i]] %>% subset(`Decreasing/Increasing Z Taxa` == "z+")
            
            all.z2.taxon <-
                c(all.z2.taxon,
                  select(temp.z2.table, !!taxa.col.name) %>% unlist())
            
        }
        
        
        # make that only uniques, and remove NA
        
        all.z1.taxon <-
            all.z1.taxon[!(is.na(all.z1.taxon))] %>% unique()
        all.z1.taxon <- all.z1.taxon[order(all.z1.taxon)]
        
        
        all.z2.taxon <-
            all.z2.taxon[!(is.na(all.z2.taxon))] %>% unique()
        all.z2.taxon <- all.z2.taxon[order(all.z2.taxon)]
        
        
        
        
        # and now iterate through to put text in each slot.
        
        summary.table.z1 <-
            cell.filler(
                titan.taxa.tables = titan.taxa.tables,
                which.z = "z-",
                taxon.list = all.z1.taxon,
                taxa.col.name = taxa.col.name,
                table.col.names = table.col.names,
                z.score.pct.cutoff = z.score.pct.cutoff,
                round.val = round.val
            )
        
        summary.table.z2 <-
            cell.filler(
                titan.taxa.tables = titan.taxa.tables,
                which.z = "z+",
                taxon.list = all.z2.taxon,
                taxa.col.name = taxa.col.name,
                table.col.names = table.col.names,
                z.score.pct.cutoff = z.score.pct.cutoff,
                round.val = round.val
            )
        
        
        # smash the two resulting tables together
        
        summary.table <- rbind(summary.table.z1, summary.table.z2)
        
        # alphebetize
        
        summary.table <-
            summary.table[order(summary.table$`All Taxon`),]
        
        
        return(summary.table)
        
        
    }
