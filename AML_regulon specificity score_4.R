## R version 4.0.3
rm(list = ls())

#calcucate JS score for specificity score
#the following is the true code for specificity score (JS score)
#input is a data.frame
#the colnames represent sample names, the rownames represent gene names
JSD <- function(data) {
  data.cp <- data
  data <- log10(data + 1)
  data <- data / rowSums(data) #the sum of each row is 1
  data <- data + 0.000001
  KLD  <-  function(x, y) sum(x * log(x / y))
  JSD  <-  function(x, y) sqrt(0.5 * KLD(x, (x + y) / 2) +
    0.5 * KLD(y, (x + y) / 2))
  refcase <- diag(ncol(data)) + 0.000001
  results <- matrix(0, nrow(data), ncol(data)) #specificity score
  maxJSspec <- NULL
  sampleID <- NULL
  maxValue <- NULL
  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      results[i, j] <- 1 - JSD(data[i, ], refcase[j, ]) #specificity score, not the JS distance
    }
    maxJSspec <- c(maxJSspec, max(results[i, ]))
    idx <- which(results[i, ] == max(results[i, ]))[1] #if multiple sample have the same max JS score, choose the first sample
    sampleID <- c(sampleID, colnames(data)[idx])
    maxValue <- c(maxValue, data.cp[i, idx]) #corresponding sampleID
  }
  specificity <- data.frame(rownames(data), maxJSspec, sampleID, maxValue)
}

#res <- JSD(data)  