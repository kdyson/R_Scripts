# likertTable for all three functions should be formatted with question columns
# and demographic (or grouping) columns, each row should be a different
# respondents answer. Cells should be e.g. likert responses (More, etc.) or
# demographic data (e.g. age groups)

## Results table for the addKruskal and addChiSquared functions should be formulated like this:

#
# activityGardensResults <-
#   tibble(
#     demographicVariable = c(
#       "Gender",
#       "Marital Status",
#       "Age",
#       "Housing",
#       "Town Type",
#       "Governorate",
#       "Employment",
#       "Work Location",
#       "Education",
#       "Before Income",
#       "After Income"
#     ),
#     afterGardensRelaxing_CHI = rep(0.0, length(demographicVariable)),
#     afterGardensRelaxing_PVAL = rep(0.0, length(demographicVariable)),
#     afterGardensRelaxing_CORPVAL = rep(0.0, length(demographicVariable)),
#     afterGardensExercise_CHI = rep(0.0, length(demographicVariable)),
#     afterGardensExercise_PVAL = rep(0.0, length(demographicVariable)),
#     afterGardensExercise_CORPVAL = rep(0.0, length(demographicVariable)),
#     afterGardensBirdPhotography_CHI = rep(0.0, length(demographicVariable)),
#     afterGardensBirdPhotography_PVAL = rep(0.0, length(demographicVariable)),
#     afterGardensBirdPhotography_CORPVAL = rep(0.0, length(demographicVariable))
#
#   )



addKruskal <-
  function(likertTable,
           resultsTable,
           questionColumns,
           demographicColumns) {
    i = min(demographicColumns) # row
    j = min(questionColumns) # column
    
    a = 1 # the first row in the resultsTable that should have results assigned.
    c = 4 # first column that should have corrected pvalues
    
    # Iterate through each demographic variable
    for (i in min(demographicColumns):max(demographicColumns)) {
      b = 2 # the first column in the resultsTable that should have results assigned to it
      
      for (j in min(questionColumns):max(questionColumns)) {
        # Create a temp table and filter out blanks.
        tempTable <- likertTable[, c(j, i)]
        tempTable <- tempTable[tempTable[2] != "",]
        
        # perform the kruskal wallis test
        tempKruskal <- kruskal.test(tempTable[, 1] ~ tempTable[, 2])
        
        # add statistic value and pvalue to the table.
        resultsTable[a, b] <- tempKruskal$statistic
        resultsTable[a, b + 1] <- tempKruskal$p.value
        
        # add 2 so next results go in the proper place
        b = b + 3
      } # end question loop
      
      a = a + 1 # go to next row for next set of demographic info
      
    } # end demographic loop
    
    # add corrected pvalue
    
    z = 1
    for (z in 1:length(questionColumns)) {
      tempCorrected <-
        p.adjust(as.vector(unlist(resultsTable[, c - 1])), method = "holm")
      
      resultsTable[, c] <- tempCorrected
      
      c = c + 3
      
    } # end addition of corrected pvalues
    
    return(resultsTable)
    
  } # end function



addChiSquared <-
  function(likertTable,
           resultsTable,
           questionColumns,
           demographicColumns,
           minVal = 25) {
    i = min(demographicColumns) # row
    j = min(questionColumns) # column
    
    a = 1 # the first row in the resultsTable that should have results assigned.
    c = 4 # first column that should have corrected pvalues
    
    # Iterate through each demographic variable
    for (i in min(demographicColumns):max(demographicColumns)) {
      b = 2 # the first column in the resultsTable that should have results assigned to it
      
      for (j in min(questionColumns):max(questionColumns)) {
        # Create a temp table and filter out blanks.
        tempTable <- likertTable[, c(j, i)]
        tempTable <-
          tempTable[tempTable[2] != "" & !is.na(tempTable[1]) ,]
        
        # Check to see if any group has less than min responses.
        testSize <-
          tempTable %>% count(tempTable[1:2]) %>% group_by(across(.cols = 2)) %>% summarise(n = sum(n))
        
        tooSmall <- testSize[testSize[2] < minVal, 1]
        print(unlist(tooSmall))
        # print(colnames(tempTable)) # use to debug
        
        # Based on this, either assign a "can't be tested"/NA indicator or the values for the chi-squared test.
        if (length(setdiff(unique(tempTable[, 2]), tooSmall)) < 2) {
          resultsTable[a, b] <- NA
          resultsTable[a, b + 1] <- NA
          
        } else {
          tempTable <- tempTable[!(tempTable[, 2] %in% tooSmall),]
          
          # perform the chi squared test
          tempChi <- chisq.test(tempTable[, 1], tempTable[, 2])
          
          # add statistic value and pvalue to the table.
          resultsTable[a, b] <- tempChi$statistic
          resultsTable[a, b + 1] <- tempChi$p.value
          
        } # end of else
        
        # Use this if you want to export a Pivot table
        # pivot <- pivot_wider(testSize,
        #                         names_from = colnames(testSize[1]),
        #                         values_from = n)
        
        # Clean up
        remove(testSize, tooSmall)
        
        # add 2 so next results go in the proper place
        b = b + 3
      } # end question loop
      
      a = a + 1 # go to next row for next set of demographic info
      
    } # end demographic loop
    
    # add corrected pvalue
    
    z = 1
    for (z in 1:length(questionColumns)) {
      tempCorrected <-
        p.adjust(as.vector(unlist(resultsTable[, c - 1])), method = "holm")
      
      resultsTable[, c] <- tempCorrected
      
      c = c + 3
      
    } # end addition of corrected pvalues
    
    return(resultsTable)
    
  } # end function



# This function needs to compare between categories within a question to see
# which categories are significantly different.

posthocChiSquared <-
  function(likertTable,
           questionColumn,
           demographicColumn,
           minVal = 25,
           correction = TRUE) {
    tempTable <- likertTable[, c(questionColumn, demographicColumn)]
    tempTable <-
      tempTable[tempTable[2] != "" & !is.na(tempTable[1]) , ]
    
    # Check to see if any group has less than min allowed responses.
    testSize <-
      tempTable %>% count(tempTable[1:2]) %>% group_by(across(.cols = 2)) %>% summarise(n = sum(n))
    
    tooSmall <- testSize[testSize[2] < minVal, 1]
    print(unlist(tooSmall))
    
    # Filter out too small categories
    tempTable <- tempTable[!(tempTable[, 2] %in% tooSmall), ]
    # and make sure demographic is a factor (this gets lost sometimes)
    if (is.factor(tempTable[, 2]) == FALSE) {
      tempTable[, 2] <- as.factor(tempTable[, 2])
    }
    
    # Make sure that the number of unique is more than two
    if (nlevels(tempTable[, 2]) < 3) {
      stop(print("Too few levels"))
    }
    
    numLevels <- nlevels(tempTable[, 2])
    
    # create the results table
    #create matrix with correct number of columns
    resultsTable <- matrix(rep(999, times = numLevels ^ 2),
                           ncol = numLevels,
                           byrow = TRUE)
    
    #define column names and row names of matrix
    tempLevels <- levels(tempTable[, 2])
    colnames(resultsTable) <- tempLevels
    rownames(resultsTable) <- tempLevels
    
    # for each [i,j] pair of factors add the pvalue for chisquared test
    i = 1 # row
    j = 1 # column
    for (i in 1:numLevels) {
      for (j in 1:numLevels) {
        if (i != j) {
          # subset for i and j levels
          testTable <-
            tempTable[tempTable[, 2] %in% tempLevels[c(i, j)] ,]
          
          # run test and assign pvalue to i,j spot
          resultsTable[i, j] <-
            chisq.test(testTable[, 1], testTable[, 2])$p.value
          
        } else {
          resultsTable[i, j] <- NA
        }
        
      } # end column loop
      
    } # end row loop
    
    # remove the lower triangle--can probably make this a filter but this is easy
    resultsTable[lower.tri(resultsTable, diag = FALSE)] <- NA
    
    # correct pvalues
    if (correction == TRUE) {
      resultsTable <-
        matrix(
          p.adjust(as.vector(resultsTable), method = 'holm'),
          ncol = numLevels,
          dimnames = list(tempLevels, tempLevels)
        )
    }
    #convert matrix to a tibble
    resultsTable <- as_tibble(resultsTable, rownames = "levels")
    
    return(resultsTable)
    
  } # end function
