library(clue)
# This function helps to compute the true accuracy
# Function to align cluster labels and compute accuracy
true_accuracy <- function(output, labels) {
  # Compute the confusion matrix
  conf_matrix <- table(output, labels)
  
  # Find the best alignment using the Hungarian algorithm
  alignment <- solve_LSAP(conf_matrix, maximum = TRUE)
  
  # Align the cluster labels
  aligned_cinx <- output
  for (i in seq_along(alignment)) {
    aligned_cinx[output == i] <- alignment[i]
  }
  
  # Compute the aligned confusion matrix
  aligned_conf_matrix <- table(aligned_cinx, labels)
  
  # Calculate the clustering accuracy
  accuracy <- sum(diag(aligned_conf_matrix)) / sum(aligned_conf_matrix)
  
  # Return the accuracy
  return(accuracy)
}








library(clue)  # Ensure the 'clue' package is loaded for solve_LSAP

# Modified function to compute the true accuracy
true_accuracy1 <- function(output, labels) {
  # Compute the confusion matrix
  conf_matrix <- table(output, labels)
  
  # Check if all labels are present in both output and labels
  unique_output <- sort(unique(output))
  unique_labels <- sort(unique(labels))
  
  # Expand the confusion matrix if some labels are missing
  all_labels <- union(unique_output, unique_labels)
  full_conf_matrix <- matrix(0, nrow = length(all_labels), ncol = length(all_labels))
  
  # Map the actual confusion matrix to the full_conf_matrix
  rownames(conf_matrix) <- as.character(unique_output)
  colnames(conf_matrix) <- as.character(unique_labels)
  rownames(full_conf_matrix) <- colnames(full_conf_matrix) <- as.character(all_labels)
  
  # Populate the full confusion matrix with the actual values
  full_conf_matrix[rownames(conf_matrix), colnames(conf_matrix)] <- conf_matrix
  
  # Find the best alignment using the Hungarian algorithm
  alignment <- solve_LSAP(full_conf_matrix, maximum = TRUE)
  
  # Align the cluster labels
  aligned_output <- output
  for (i in seq_along(alignment)) {
    aligned_output[output == all_labels[i]] <- all_labels[alignment[i]]
  }
  
  # Compute the aligned confusion matrix
  aligned_conf_matrix <- table(aligned_output, labels)
  
  # Calculate the clustering accuracy
  accuracy <- sum(diag(aligned_conf_matrix)) / sum(aligned_conf_matrix)
  
  # Return the accuracy
  return(accuracy)
}
