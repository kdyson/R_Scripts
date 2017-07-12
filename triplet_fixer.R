# This script creates a sourceable function to turn any matrix into a triplet using the labdsv matrify()
# Written 4/21/2016

# requires the labdsv package

# The function inputs should be:
    # filename or data table
    # three columns that make up the triplet.

# The function name is ez.matrify

ez.matrify <- function(filename, species.name, site.name, abundance) {
    
    library(labdsv)
    
    #check if filename is a .csv or already imported
    if (is.character(filename))
        data.for.analysis <- read.csv(file=filename, header = T) else
        data.for.analysis <- filename
    
    proper.columns <- as.data.frame(data.for.analysis[ ,c(site.name, species.name, abundance)])
    ez.matrify.output <- matrify(proper.columns)
    return(ez.matrify.output)
    
}


#test
#head(ez.matrify(filename='../../DataRepository/VegetationData/ShrubsCSV.csv', species.name = 'SpeciesTaxonomic', site.name = 'SiteStandardGroup', abundance = 'RandomTest'))

#head(ez.matrify(filename=tree.data, species.name = 'tree.species', site.name = 'site', abundance = 'tree.number'))

