#
# Make data anonymous! This script takes an input data table (csv or in memory)
# and a conversion table (lookup table) and then produces an output table (saved
# to disk or memory).
# 
# Based in part on information found at: https://stackoverflow.com/questions/35636315/replace-values-in-a-dataframe-based-on-lookup-table
#
# ... passes arguments to write.csv.
# 


anon.data <- function(input.table, lookup.table, site.column, match.column, anon.column, output.path = NA, write.to.disk = FALSE, ...){

# -- Input table -------

if(is.character(input.table) == TRUE) {
    
    input.table <- read.csv(input.table, stringsAsFactors = FALSE)
    
}

exists("input.table")

# -- Lookup table ------------

if(is.character(lookup.table) == TRUE) {
    
    lookup.table <- read.csv(lookup.table, stringsAsFactors = FALSE)
    
}

exists("lookup.table")

# -- Anon the table -----------

replace.me <- input.table[ , site.column]
lookup.match <- lookup.table[ , match.column]
replace.with <- lookup.table[ , anon.column]

anon.table <- input.table
anon.table[ , which(colnames(anon.table) == site.column)] <- replace.with[match(
    replace.me, lookup.match)]


# -- Export the output table ---------

if(is.character(output.path) == TRUE & write.to.disk == TRUE) {
    
    write.csv(x = anon.table, file = output.path, ...)
    return(anon.table)
    
}

if(is.character(output.path) == TRUE & write.to.disk == FALSE) {
    
    warning("warning: table is chr and write to disk is false, creating a default output table in Global Environment")
    
    return (anon.table)
    
}

if(is.character(output.path) == FALSE & write.to.disk == FALSE) {
    
    return (anon.table)
    
}


if(is.character(output.path) == FALSE & write.to.disk == TRUE) {
    
    warning("warning: no output path specified, csv not saved to disk.")
    return (anon.table)
    
}

}
