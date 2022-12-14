# Integrated function

"stepfor" <- function(y = y, # Single gene vector with rowname
                      d = d, # Design file
                      alfa = 0.05, # Significant level #---Bug--fix needed
                      family = gaussian(), # Family #---Bug--fix needed
                      epsilon = 1e-05 # Epsilon
) {

  # Null variable Invoke
  # p-value
  pval <- NULL

  # Design file
  design <- NULL
  
  # Invoke j as 1 
  j <- 1
  
  # Fill model with all cells (gene_i across all cells/samples)
  # Extract summary and save it
  resul0 <- summary(glm(y ~ ., data = d, family = family, epsilon = epsilon))$coefficients[, 4]
  
  print(resul0)
  
  stop()
  
  d <- as.data.frame(d[, names(resul0)[-1]])
  for (i in 1:ncol(d)) {
    sub <- cbind(design, d[, i])
    sub <- as.data.frame(sub)
    lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon)
    result <- summary(lm2)
    pval[i] <- result$coefficients[, 4][j + 1]
  }
  min <- min(pval)
  while (min < alfa) {
    b <- pval == min
    c <- c(1:length(pval))
    pos <- c[b]
    pos <- pos[!is.na(pos)][1]
    design <- cbind(design, d[, pos])
    design <- as.data.frame(design)
    colnames(design)[j] <- colnames(d)[pos]
    j <- j + 1

    if (ncol(d) == 2) {
      lastname <- colnames(d)[!b]
    }
    d <- as.data.frame(d[, -pos])
    if (ncol(d) == 1) {
      colnames(d) <- lastname
    }

    pval <- NULL
    if (ncol(d) != 0) {
      for (i in 1:ncol(d)) {
        sub <- cbind(design, d[, i])
        sub <- as.data.frame(sub)
        lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon)
        result <- summary(lm2)
        pval[i] <- result$coefficients[, 4][j + 1]
      }
      min <- min(pval, na.rm = TRUE)
    } else {
      min <- 1
    }
  }
  if (is.null(design)) {
    lm1 <- glm(y ~ 1, family = family, epsilon = epsilon)
  } else {
    lm1 <- glm(y ~ ., data = design, family = family, epsilon = epsilon)
  }
  return(lm1)
}
