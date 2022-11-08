p.vector <- function(
        data, # Expected Single Cell Count Table
        design, # Internal Design File made previously # Accepts list,
        # data.frame, and matrix #--Fix required: Restrict input to maSigPro internal design file
        Q = 0.05, # significance level
        MT.adjust = "BH", # Depends on library(stats) p.adjust
        min.obs = 6,
        counts = FALSE, # Logical to indicate if the data are counts
        family = NULL,
        theta = 10,
        epsilon = 1e-05,
        item = "gene") {
    
    # Condition-1: Check the datatype for the design file
    # If the supplied file is a dataframe or a matrix, 'groups.vector' and 'edesign'
    # will have null invoke
  if (is.data.frame(design) || # if it is a dataframe
      is.matrix(design)) {# if it is a matrix
      
      # Internal Design file == User supplied Matrix or Data.frame
    dis <- design
    
    # Invoke NULL
    groups.vector <- NULL
    
    # Invoke NULL
    edesign <- NULL
    
    # Else if a list from make.matrix.design()
    # Check list structure #---Fix: Better check has to be performed
  } else if (is.list(design)) {
      
      # Assign internal design file
    dis <- as.data.frame(design$dis)
    
    # Assign group vectors
    groups.vector <- design$groups.vector
    
    # Add original experimental.design file
    edesign <- design$edesign
  }

    # Conditio-2; Check for Residual Family
    # For single cell we will need negative binomial
    # Selects negative binomial if counts == T
  if (is.null(family)) { # Default gaussian
    if (!counts) { # if not equal to null
      family <- gaussian() # Gaussian family
    }
    if (counts) { # if counts ==T
      family <- negative.binomial(theta) # NB with input/default theta
    }
  }

    # Data to matrix conversion
    #-Fix: Check it before hand
  dat <- as.matrix(data)
  
  # Re-ordering via rownames
  # This will also drop the cell/ samples which donot exist in metadata
  #---Bug--: Will not work if the internal design file will have extra cells
  dat <- dat[, as.character(rownames(dis))]
  
  # Number of genes in the count matrix
  G <- nrow(dat)

  # Removing rows with many missings:
  # This will count the number of missing values 
  #---Warning: Bux fix required for some cases: Only NA is counted
  # ScRNA-Seq has many 0s, which is not counted here
  count.na <- function(x) (length(x) - length(x[is.na(x)]))
  
  # Subset the counts, aginst the min.obs, default == 6
  dat <- dat[apply(dat, 1, count.na) >= min.obs, ]

  # Removing rows with all ceros:
  # Sepcific to single cell
  #---Better-Syntax: dat <- dat[rowSums(dat) !=0,]
  sumatot <- apply(dat, 1, sum)
  counts0 <- which(sumatot == 0)
  if (length(counts0) > 0) {
    dat <- dat[-counts0, ]
  }
  
  # Save attributes of the final dataframe/ counts
  g <- dim(dat)[1] # Number of Genes/ Rows #--- better syntax:  g <- nrow(dat)
  n <- dim(dat)[2] # Number of cells/ Columns  #--- better syntax:  n <- ncol(dat)
  p <- dim(dis)[2] # Number of terms of the model or group vectors
  
  # Invoke vector equivalent to length of genes
  # Not necessary to invoke a defined length
  p.vector <- vector(mode = "numeric", length = g)

  #----------------------------------------------------------------------
  # Model fitting
  # Run for loop on each gene one by one
  # g is the total number of genes
  for (i in 1:g) {
      
    # Subset one gene at a time
      # Store in Variable y for modelling
    y <- as.numeric(dat[i, ])

    #---Bad syntax to print progress**************************************
    # Replace with progress()*********************************************
    div <- c(1:round(g / 100)) * 100
    if (is.element(i, div)) {
      print(paste(c("fitting ", item, i, "out of", g), collapse = " "))
    } # ******************************************************************
    # ********************************************************************

    # Actual model fitting
    model.glm <- glm(y ~ ., data = dis, family = family, epsilon = epsilon)
    
    print(summary(model.glm))
    stop()
    
    # Condition: 
    if (model.glm$null.deviance == 0) {
      p.vector[i] <- 1
    } else {
      model.glm.0 <- glm(y ~ 1, family = family, epsilon = epsilon)

      if (family$family == "gaussian") {
        test <- anova(model.glm.0, model.glm, test = "F")
        if (is.na(test[6][2, 1])) {
          p.vector[i] <- 1
        } else {
          p.vector[i] <- test[6][2, 1]
        }
      } else {
        test <- anova(model.glm.0, model.glm, test = "Chisq")
        if (is.na(test[5][2, 1])) {
          p.vector[i] <- 1
        } else {
          p.vector[i] <- test[5][2, 1]
        }
      }
    }
    
    # Model fitting outer loop close
  }

  #----------------------------------------------------------------------
  # P-value adjustment with stats package
  p.adjusted <- p.adjust(p.vector, method = MT.adjust, n = length(p.vector))
  genes.selected <- rownames(dat)[which(p.adjusted <= Q)]

  FDR <- sort(p.vector)[length(genes.selected)]

  SELEC <- as.matrix(as.data.frame(dat)[genes.selected, ])
  if (nrow(SELEC) == 0) {
    print("no significant genes")
  }

  p.vector <- as.matrix(p.vector)
  rownames(p.vector) <- rownames(dat)
  colnames(p.vector) <- c("p.value")

  #-------------------------------------------------------------------------
  #---Better Syntax: list.name <- list(obj.name = obj, obj2.name = obj2)
  # return named list of objects
  # Make list of object made
  output <- list(SELEC, p.vector, p.adjusted, G, g, FDR, nrow(SELEC), dis, dat, min.obs, Q, groups.vector, edesign, family)
  
  # Assign name
  names(output) <- c(
    "SELEC", "p.vector", "p.adjusted", "G", "g", "FDR", "i", "dis", "dat", "min.obs", "Q", "groups.vector", "edesign",
    "family"
  )
  
  #---Syntax-fix needed: Use proper retun statement, return(output)
  output
}
