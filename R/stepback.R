# Integrated function

stepback <- function(y = y, # Single gene vector with rowname
                     d = d, # Design file
                     alfa = 0.05, # Significant level #---Bug--fix needed
                     family = gaussian(), # Family #---Bug--fix needed
                     epsilon = 1e-05 # Epsilon
) {

  # Fill model with all cells (gene_i across all cells/samples)
  lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon)
  
  # Save summary
  result <- summary(lm1)
  
  # Select Max Value from Pr(>|t|)
  #[-1] Donot take the b_0 even if it is significant
  #---Good Synatx-Readapt in codes: If missing remove 
  max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
  
  
 
  # Check later -----------------------------------
  if (length(result$coefficients[, 4][-1]) == 1) {
    if (max > alfa) {
      max <- 0
      lm1 <- glm(y ~ 1, family = family, epsilon = epsilon)
      stop()
    }
  }
  #------------------------------------------------
  print("skip")
  
 # Run the model for backward selection
  while (max > alfa) {
    varout <- names(result$coefficients[, 4])[result$coefficients[, 4] == max][1]
    pos <- position(matrix = d, vari = varout)
    d <- d[, -pos]
    if (length(result$coefficients[, 4][-1]) == 2) {
      min <- min(result$coefficients[, 4][-1], na.rm = TRUE)
      lastname <- names(result$coefficients[, 4][-1])[result$coefficients[, 4][-1] == min]
    }
    if (is.null(dim(d))) {
      d <- as.data.frame(d)
      colnames(d) <- lastname
    }
    lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon)
    result <- summary(lm1)
    max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
    if (length(result$coefficients[, 4][-1]) == 1) {
      max <- result$coefficients[, 4][-1]
      if (max > alfa) {
        max <- 0
        lm1 <- glm(y ~ 1, family = family, epsilon = epsilon)
      }
    }
  }
  return(lm1)
}
