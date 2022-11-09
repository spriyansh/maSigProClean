## Testing
source("stepback.R")
source("position.R")

## Variable selection method
T.fit <- function(data, # fit model list
                  design = data$dis, # Design file
                  step.method = "backward", # Method of variable selection
                  min.obs = data$min.obs, # Minimum number of observations to consider
                  alfa = data$Q, # Significance threshold
                  nvar.correction = FALSE,
                  family = gaussian(), #--Bug---Need.fix Default function set to gaussian(), set to Ngeative binomial
                  epsilon = 1e-05, #---Syntax-fix needed Epsilon Needs same formal
                  item = "gene" # Item gene
                  ) {
    
    # Condition to check data type
    # Why accept two formats? ---Better syntax needed to select standard class
    
  if (is.list(data)) {
    
      #---Syntax-fix needed---Pass as matrix  
    dat <- as.matrix(data$SELEC)

    # Add row at index-1 having values 1
    # Increase gene counts by 1
    dat <- rbind(c(rep(1, ncol(dat))), dat)
    
    # Group.vector transfer
    #---Syntax-fix: no need to tranfer 
    groups.vector <- data$groups.vector
    
    # This will select for the single independent variable not the 'vs'
    # Then add that in the group vector in the begining
    #---Better syntax needed---
    groups.vector <- c(groups.vector[nchar(groups.vector) == min(nchar(groups.vector))][1], groups.vector)
    
    # edesign transfer
    #---Syntax-fix: no need to tranfer 
    edesign <- data$edesign
    
    # number of gene transfer
    #---Syntax-fix: no need to tranfer
    G <- data$g
    
    # family transfer
    #---Syntax-fix: no need to tranfer 
    family <- data$family
    
    # Cruft-code// No need to accept different dataframe
    #----------------------------------------
  } else {
    G <- nrow(data)
    data <- rbind(c(rep(1, ncol(data))), data)
    dat <- as.matrix(data)
    count.na <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[apply(dat, 1, count.na) >= min.obs, ]
    groups.vector <- NULL
    edesign <- NULL
  }
    #----------------------------------------------------
    
    #---Better syntax needed
    #---No need to dataframe conversion
  dis <- as.data.frame(design)
  
  #---Cruft-code
  #---Better syntax needed
  #---Accept counts with matched data
  #---Use of S4 classes is encouraged
  dat <- dat[, as.character(rownames(dis))]
  
  # Number of Genes
  # 1 is subtracted because it was added before
  g <- (dim(dat)[1] - 1)
  
  # Number of samples/cells
  n <- dim(dat)[2]
  
  # P something (No explanation provided)
  # Equal to number of cells/ samples
  p <- dim(dis)[2]
  
  # This was equal to group.vector previously
  #---Better Synatx, used group vector instead
  vars.in <- colnames(dis)
  
  # Null invoke decleration
  #---User better syntax
  sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
  
  # Outlier Genes matrix
  # Equal to the number of cells/groups
  influ.info <- matrix(NA, nrow = nrow(dis), ncol = 1)
  
  # Transfer rownames
  # Cell/sample names tranfer
  rownames(influ.info) <- rownames(dis)
  
  # Ask Ana about this-----
  #------------------------
  #----Imapct of nvar------
  if (nvar.correction) {
    alfa <- alfa / ncol(dis)
  }
  #------------------------
  
  # For loop: fits model on each gene
  # Since first row is full of 1, loop is started at +1
  for (i in 2:(g + 1)) {
      
      # Vector of counts for gene_i
    y <- as.numeric(dat[i, ])
    
    # Transfer rownames of the genes to vector
    name <- rownames(dat)[i]
    
    # Backward selection methods
    if (step.method == "backward") {
        
        # return object
        # Stepback is custom function
      reg <- stepback(y = y, # Gene Counts for gene_i
                      d = dis, # Design file
                      alfa = alfa, # Significance level
                      family = family, # Family
                      epsilon = epsilon # Epsilon
                      )
      
      stop()
      
      ## Foward-Selection-Procedure: Check-later
    } else if (step.method == "forward") {
        
        #-- Can we use 'apply()` here?
      reg <- stepfor(y = y, d = dis, alfa = alfa, family = family, epsilon = epsilon)
      
      ## two-ways-backward-Procedure: Check-later
    } else if (step.method == "two.ways.backward") {
      reg <- two.ways.stepback(y = y, d = dis, alfa = alfa, family = family, epsilon = epsilon)
      
      ## two-ways-forward-Procedure: Check-later
    } else if (step.method == "two.ways.forward") {
      reg <- two.ways.stepfor(y = y, d = dis, alfa = alfa, family = family, epsilon = epsilon)
      
      ## If something else is supplied 
      ##---Syntax fix--could be checked before 
    } else {
      stop("stepwise method must be one of backward, forward, two.ways.backward, two.ways.forward")
    }
    
    
    
    div <- c(1:round(g / 100)) * 100
    if (is.element(i, div)) {
      print(paste(c("fitting ", item, i, "out of", g), collapse = " "))
    }
    lmf <- glm(y ~ ., data = as.data.frame(dis), family = family, epsilon = epsilon)
    result <- summary(lmf)
    novar <- vars.in[!is.element(vars.in, names(result$coefficients[, 4]))]
    influ <- influence.measures(reg)$is.inf
    influ <- influ[, c(ncol(influ) - 3, ncol(influ) - 1)]
    influ1 <- which(apply(influ, 1, all))
    if (length(influ1) != 0) {
      paste.names <- function(a) {
        paste(names(a)[a], collapse = "/")
      }
      match <- match(rownames(dis), rownames(influ))
      influ <- as.data.frame(apply(influ, 1, paste.names))
      influ.info <- cbind(influ.info, influ[match, ])
      colnames(influ.info)[ncol(influ.info)] <- name
    }
    result <- summary(reg)
    if ((!(result$aic == -Inf) & !is.na(result$aic) & family$family == "gaussian") | family$family != "gaussian") {
      k <- i

      # Computing p-values

      model.glm.0 <- glm(y ~ 1, family = family, epsilon = epsilon)

      if (family$family == "gaussian") {
        test <- anova(model.glm.0, reg, test = "F")
        p.value <- test[6][2, 1]
      } else {
        test <- anova(model.glm.0, reg, test = "Chisq")
        p.value <- test[5][2, 1]
      }
      # Computing goodness of fitting:

      bondad <- (reg$null.deviance - reg$deviance) / reg$null.deviance
      if (bondad < 0) {
        bondad <- 0
      }
      beta.coeff <- result$coefficients[, 1]
      beta.p.valor <- result$coefficients[, 4]
      coeff <- rep(0, (length(vars.in) + 1))
      if (length(novar) != 0) {
        for (m in 1:length(novar)) {
          coeff[position(dis, novar[m]) + 1] <- NA
        }
      }
      p.valor <- t <- as.numeric(rep(NA, (length(vars.in) + 1)))

      if (result$coefficients[, 4][rownames(result$coefficients) == "(Intercept)"] < alfa) {
        coeff[1] <- result$coefficients[, 1][rownames(result$coefficients) == "(Intercept)"]
        p.valor[1] <- result$coefficients[, 4][rownames(result$coefficients) == "(Intercept)"]
        t[1] <- result$coefficients[, 3][rownames(result$coefficients) == "(Intercept)"]
      }
      for (j in 2:length(coeff)) {
        if (is.element(vars.in[j - 1], rownames(result$coefficients))) {
          coeff[j] <- result$coefficients[, 1][rownames(result$coefficients) == vars.in[j - 1]]
          p.valor[j] <- result$coefficients[, 4][rownames(result$coefficients) == vars.in[j - 1]]
          t[j] <- result$coefficients[, 3][rownames(result$coefficients) == vars.in[j - 1]]
        }
      }
      if (!all(is.na(p.valor))) {
        sol <- rbind(sol, as.numeric(c(p.value, bondad, p.valor)))
        coefficients <- rbind(coefficients, coeff)
        t.score <- rbind(t.score, t)
        sig.profiles <- rbind(sig.profiles, y)
        h <- nrow(sol)
        rownames(sol)[h] <- name
        rownames(coefficients)[h] <- name
        rownames(t.score)[h] <- name
        rownames(sig.profiles)[h] <- name
      }
    }
  }
  if (!is.null(sol)) {
    sol <- as.data.frame(sol)
    coefficients <- as.data.frame(coefficients)
    coeffic <- coefficients
    t.score <- as.data.frame(t.score)
    sig.profiles <- as.data.frame(sig.profiles)
    colnames(sol) <- c("p-value", "R-squared", "p.valor_beta0", paste("p.valor_", vars.in, sep = ""))
    colnames(coefficients) <- c("beta0", paste("beta", vars.in, sep = ""))
    colnames(t.score) <- c("t.score_beta0", paste("t.score_", vars.in, sep = ""))
    colnames(sig.profiles) <- colnames(dat)
    if (!is.null(groups.vector) & !is.null(edesign)) {
      groups <- colnames(edesign)[3:ncol(edesign)]
      degree <- (length(groups.vector) / length(groups)) - 1
      for (w in 1:nrow(coefficients)) {
        A <- NULL
        col.names <- NULL
        for (l in 1:length(groups)) {
          B <- reg.coeffs(coefficients = coefficients[w, ], groups.vector = groups.vector, group = groups[l])
          cols <- paste(rep(groups[l], each = length(B)), paste("beta", c(0:(length(B) - 1)), sep = ""), sep = "_")
          A <- c(A, B)
          col.names <- c(col.names, cols)
        }
        group.coeffs <- (rbind(group.coeffs, A))
      }
      colnames(group.coeffs) <- col.names
      rownames(group.coeffs) <- rownames(coefficients)
    }
  }
  if (ncol(influ.info) > 2) {
    print(paste("Influence:", ncol(influ.info) - 1, "genes with influential data at slot influ.info. Model validation for these genes is recommended"))
  }
  influ.info <- influ.info[, -1]
  output <- list(
    sol, sig.profiles, coefficients, as.data.frame(group.coeffs), t.score, vars.in, G, g, dat, dis, step.method, groups.vector,
    edesign, influ.info
  )
  names(output) <- c(
    "sol", "sig.profiles", "coefficients", "group.coeffs", "t.score", "variables", "G", "g", "dat", "dis", "step.method",
    "groups.vector", "edesign", "influ.info"
  )
  output
}
