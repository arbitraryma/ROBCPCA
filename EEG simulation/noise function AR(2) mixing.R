add_noise_to_group <- function(group_data_list, noise_df = 3) {
  # Randomly select 2 processes
  selected_processes <- sample(1:length(group_data_list), 2)
  
  # Randomly select 15% of the rows
  num_rows <- nrow(group_data_list[[1]])
  num_noisy_rows <- ceiling(0.15 * num_rows)
  selected_rows <- sample(1:num_rows, num_noisy_rows)
  
  # Add noise to the selected rows in the selected processes
  for (process_idx in selected_processes) {
    # Generate noise from t-distribution
    noise <- rt(num_noisy_rows * ncol(group_data_list[[process_idx]]), df = noise_df)
    noise_matrix <- matrix(noise, nrow = num_noisy_rows, ncol = ncol(group_data_list[[process_idx]]))
    
    # Add the noise to the selected rows
    group_data_list[[process_idx]][selected_rows, ] <- group_data_list[[process_idx]][selected_rows, ] + noise_matrix
  }
  
  return(group_data_list)
}
