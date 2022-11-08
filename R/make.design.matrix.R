# Function to make design matrix
"make.design.matrix" <- function(edesign, # Design Matrix supplied as data.frame
                                 degree = 2, # Order of polynomial; 2 means a quadratic fit
                                 time.col = 1, # Column having time-series information
                                 repl.col = 2, # Column having replicate information
                                 group.cols = c(3:ncol(edesign)) # Dummy Variable Information
                                 # Only dimensions are stored as vector
                                 # Always start at 3 #---Fix Required---
) {
  # Selection of column to be treated as control
  # The first column after the replicate is the control column
  control.label <- colnames(edesign)[group.cols][1] #---Fix Required---

  # Enter the condition if the edesign has more than 3 columns
  # Columns are checked after matrix conversion (might not be necessary) # Can be removed
  # dim(as.matrix(edesign))[2] returns column counts, ncols(edesign) is more suitable
  if (dim(as.matrix(edesign))[2] > 3) {

    # This will select the column indexes from the vector
    # Index start from 2 to avoid dummy variable trap
    # Dummy columns will always be 1 less than group.cols in length
    # There are more efficient ways to do this #---Fix Required
    dummy.cols <- group.cols[2:length(group.cols)]

    # Adding label information
    # Sep = "_" for easier interpretation  #---A fix might be required---#
    # This will create a single label vector of type (experiment'vs'control)
    treatm.label <- paste(colnames(edesign)[dummy.cols], "vs", control.label, sep = "")

    # Making group label information
    # This will make a vector with control (i.e. 3rd column always) and exp
    # i.e. other columns 4...
    groups.label <- c(control.label, treatm.label)

    # This will create a one column matrix as columns 3 is dropped due to
    # dummy variable trap
    # Ideally it will create a matrix having with columns == ncols(group.cols)-1
    matrix.dummy <- as.matrix(edesign[, dummy.cols])

    ## Shared origin
    # Variable initialization
    dummy <- NULL
    j <- 0

    # Select the value having first time-point
    # This might need a fix for scRNA-trajectory
    # Multiple cells can have 0 time points
    # Only 1 cell can have 0 time points
    # Cluster/Clusters of cells can have 0 time points
    # This has to be accounted for #---A fix might be required---#
    origen <- min(edesign[, time.col])

    #---Bug-Fix-Required---#
    # It will always select column 1 as the time-series column
    # Either make is default by removing the argument from the function
    # Or replace ```edesign[, 1]``` as ```edesign[, time.col]```
    #---
    # This selects all the cells/samples having time-point 0 or minimum
    # Select happens on the design file, therefore all the columns are selected
    origen <- edesign[edesign[, 1] == origen, ]

    # Enter a for loop
    # Run for the length of dummy variables i.e. ncols(group.cols)-1
    # i is initialized with 0
    for (i in 1:length(dummy.cols)) {

      # Create a new columns 'share', to store the sum of dummny variable
      # Eg. df[2,1] + df[2,2] + df[2,3] = SUM
      # If the dummy variable share is unique then the sum for all will be 1
      share <- apply(origen[, c(3, dummy.cols[i])], 1, sum)

      # Check condition if nothing is shared among the groups
      # This looks for the shared sum, if all of them are less than 1
      # Then that means nothing is shared among the groups
      # This returns true negatively
      if (!is.element(TRUE, share > 1)) {

        # Incrementing by 1
        j <- j + 1

        # Actual dummy variable column creation happens here
        # It runs in for-loop for i number of columns
        # i is always -1 than the groups to avoid dummy trap
        # In case of the two lineage it will run only once
        # Cbind is used to add columns one by one
        # This method only happens when the share is 1 for all
        dummy <- cbind(dummy, matrix.dummy[, i])

        # Rename the column using the experiment'vs'control label
        # Renming also happens iteratively with for loop
        colnames(dummy)[j] <- treatm.label[i]
      }

      # for-loop close
    }

    # Make a vector of all the time values
    # User-specified column of time-series index is used for the purpose
    # `time` is a one column matrix always (what if- two time series exist?)
    # Slingshot often results in two time series
    time <- as.matrix(edesign[, time.col])

    # Renaming the column of the `time` based on the edesign file supplied by
    # user. This can be something custom as it is handled internally
    #---Might need a fix---#
    colnames(time) <- colnames(edesign)[time.col]

    # Using cbind to combine the time column with the dummy column
    # In a two lineage case dummy is a single column data.frame
    #---This is calling for trouble if the data is not sorted---#
    #---Fix: A custom function will be helful to sort the data---#
    dis <- cbind(dummy, time)

    # Transfer the rownames from the user supplied edesign file to this
    # internal design file
    #---This is calling for trouble if the data is not sorted---#
    #---Fix: A custom function will be helful to sort the data---#
    rownames(dis) <- rownames(edesign)

    # Create a vector of comparison and control
    # group vector will have experiment'vs'control and control
    groups.vector <- c(colnames(dummy), control.label)

    # Extract the column names of the internal design file
    # It should have the comparison groups and time columns
    colnames.dis <- colnames(dis)

    # Combine the dummy variable design with the time column
    # All the 1 will have time variable multiplied
    # A column having unique time values per experiment
    #---Needs more explanation---#
    dis <- cbind(dis, dis[, ncol(dis)] * matrix.dummy)

    # Renaming he columns
    colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col], "x", colnames(edesign)[dummy.cols], sep = ""))

    # Added extra column for the experiment'vs'control
    groups.vector <- c(groups.vector, treatm.label)


    # Polynomial Modelling with degree two or more; 1 will be a straight line always
    # Time is given the orders
    if (degree >= 2) {

      # For loop start with quadratic fits and can go further
      for (i in 2:degree) {

        # Save colnames for the internal design file
        colnames.dis <- colnames(dis)

        # Modeled order is -1 than the supplied order
        # if 3 is supplied max-order is 2
        dis <- cbind(
          dis, # Previous
          edesign[, time.col]^i,
          edesign[, time.col]^i * edesign[, dummy.cols]
        )

        # Column names are renamed based on model formula
        #---Fix needed for '^' to make formula interpretable
        colnames(dis) <- c(
          colnames.dis, # Previous
          paste(colnames(edesign)[time.col], i, sep = ""),
          # paste(colnames(edesign)[time.col], i, sep = "^"), # Changes
          paste(colnames(edesign)[time.col], "",
            i, "x", colnames(edesign)[dummy.cols],
            sep = ""
          )
        )

        # Seven terms for no reason --- need fix or explanation
        # In this case the group vector will be equal to the product of order
        # and groups (if more than 1) + 1
        # In this case group.cols == 2 so 2*order of polynomial+1
        groups.vector <- c(groups.vector, groups.label)
      }
      # close for-loop
    } # close condition
    #---Fix needed, no check available for linear models (order==1)---#


    # Else for the condition if the edesign has less than or equal to 3 columns
    # No interaction term used in this case with the successive groups
  } else {
    print("Entered in the else nest")

    # Make the internal design matrix with time column
    dis <- as.matrix(edesign[, time.col])

    # Transfer column name 'Time', can be something arbitary supplied by user
    #---Imprrovemnet, make is order or pseudotime for trajectory inference
    colnames(dis) <- colnames(edesign)[time.col]

    # Transfer the rownames from the user supplied edesign file to this
    # internal design file
    #---This is calling for trouble if the data is not sorted---#
    #---Fix: A custom function will be helful to sort the data---#
    rownames(dis) <- rownames(edesign)

    # Design polynomial model
    #---Syntax improvemnet needed here, Different conditions used
    # Either set  (degree > 1) or  (degree >= 2) (like in previous)
    if (degree > 1) {

      # For loop start with quadratic fits and can go further
      for (i in 2:degree) {

        # Transfer time-sries colnames to dataframe
        colnames.dis <- colnames(dis)

        # Add higher order columns
        dis <- cbind(dis, edesign[, time.col]^i)

        # Rename the columns with higher order column informations
        # colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col], i, sep = ""))
        colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col], i, sep = "^")) # change
      }
      # close for-loop
    }
    # Return group vector
    # In this case the group vector will be equal to the degree
    # Group vector will only depend on the group cols
    # In this case group.cols == 1 so 1*order of polynomial
    groups.vector <- rep(colnames(edesign)[group.cols], degree)
  } # close condition

  # Return a named list of objects
  output <- list(
    dis, # internal design matrix
    groups.vector, # Unused for now
    edesign # User supplied design file
  )

  # Name list objets
  names(output) <- c("dis", "groups.vector", "edesign")

  #--Possible fix: Can be done by single command
  # output <- list(dis = dis, groups.vector = groups.vector, edesign = edesign)

  #---Syntax-Fix-needed use proper return()
  output
}
