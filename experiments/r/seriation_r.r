# Load required libraries
library(seriation)
library(data.table)

# Function to get neighboring indices in a matrix
get_neighbors <- function(i, j, nrows, ncols, eps_neigh) {
  neighbors <- list()
  
  for (i1 in 1:nrows) {
    for (j1 in 1:ncols) {
      dist <- sqrt((i - i1)^2 + (j - j1)^2)
      if (dist >= 0.1 && dist <= eps_neigh) {
        neighbors <- append(neighbors, list(c(i1, j1)))
      }
    }
  }
  
  return(neighbors)
}

# Function to compute the Homogeneity index
compute_homogeneity <- function(matrix, eps_neigh) {
  nrows <- nrow(matrix)
  ncols <- ncol(matrix)
  HOM <- matrix(0, nrow = nrows, ncol = ncols)

  for (i in 1:nrows) {
    for (j in 1:ncols) {
      neighbors <- get_neighbors(i, j, nrows, ncols, eps_neigh)
      if (length(neighbors) > 0) {
        HOM[i, j] <- (1 / length(neighbors)) * sum(sapply(neighbors, function(idx) {
          abs(matrix[i, j] - matrix[idx[1], idx[2]])
        }))
      }
    }
  }

  return(sum(HOM) / (nrows * ncols))
}

# Function to run seriation experiments and save results incrementally
run_seriation_experiments <- function(file, matrix_instance, output_file = "seriation_results.csv") {
  
  # Get available seriation methods
  methods <- unlist(seriation::list_seriation_methods("dist"))
  
  # Convert matrix to distance matrix
  dist_matrix <- dist(matrix_instance)
  
  # Create a new file and write the header if it doesn't exist
  file_exists <- file.exists(output_file)
  if (!file_exists) {
    fwrite(data.table(Method = character(), CPU_Time = numeric(), Homogeneity = numeric(), Order = character()), 
           output_file, row.names = FALSE)
  }
  
  # Iterate through each method and append results
  for (method in methods) {
    cat("Running seriation with method:", method, "\n")
    # Measure CPU time
    start_time <- Sys.time()
    # Apply seriation method
    seriation_result <- tryCatch({
      seriate(dist_matrix, method = method)
    }, error = function(e) {
      cat("Error with method", method, ":", e$message, "\n")
      return(NULL)
    })
    end_time <- Sys.time()
    cpu_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    # If method was successful
    if (!is.null(seriation_result)) {
      order <- get_order(seriation_result)
      # Reorder the matrix using the obtained order
      reordered_matrix <- matrix_instance[order, order]
      # Compute homogeneity
      homogeneity <- compute_homogeneity(reordered_matrix, 1)
      # Save reordered matrix
    #write.csv(reordered_matrix, paste0("reordered_matrix_", method, ".csv"), row.names = FALSE)
      # Append results to CSV
      result_row <- data.table(File=file, Method = method, CPU_Time = cpu_time, Homogeneity = homogeneity, Order = paste(order, collapse = " "))
      fwrite(result_row, output_file, append = TRUE, row.names = FALSE)
      cat("Saved results for method:", method, "\n")
    }
  }
  cat("All results saved to", output_file, "\n")
}

# Example Usage: Create a random matrix and run seriation
set.seed(42)
matrix_instance <- matrix(sample(1:10, 10 * 10, replace = TRUE), nrow = 10, ncol = 10)

# Run the function and store the results incrementally
txt_files <- list.files(path = "instances/", pattern = "\\.txt$", full.names = TRUE)

read_matrix_from_txt <- function(file_path) {
  matrix_data <- as.matrix(read.table(file_path, header = FALSE))
  return(matrix_data)
}

for (file in txt_files) {
    matrix_data <- as.matrix(read.table(file, header = FALSE))
    run_seriation_experiments(file, matrix_instance)
}