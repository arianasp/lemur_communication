library(dplyr)
library(tidyr)
library(lubridate)
library(abind)
library(igraph)

#specify time window for pooling (in seconds) - L15!


#Import and prepare data-------------------------
        data <- readRDS("output/clean_proxit.rds")
        colnames(data)[c(5,6,7)] <- c("TX_Power", "Start_Time", "Seconds_Duration")

    #time window
       window_size <- 5

    #convert Start_Time to numeric seconds
        data <- data %>%
          mutate(
            second = as.numeric(floor_date(Start_Time, "second"))
          )
    
    #get unique animal IDs
        all_animals <- unique(c(data$animal_ID1, data$animal_ID2))
    
    #determine overall time range
        min_second <- min(data$second)
        max_second <- max(data$second)
        
        # start_time <- data[data$second==min_second, "Start_Time"]
        # end_time <- data[data$second==max_second, "Start_Time"]
        

#Assign values to time window-----------------
      data <- data %>%
        mutate(
          window_index = floor((second - min_second) / window_size),
          window_start = min_second + window_index * window_size
        )
      
      #create full sequence of window start times
          time_windows <- seq(min_second, max_second, by = window_size)
      
      #average RSSI per dyad per window
          data_clean <- data %>%
            group_by(window_start, animal_ID1, animal_ID2) %>%
            summarise(RSSI = mean(RSSI, na.rm = TRUE), .groups = "drop")
      
      #create a list of matrices, one per time window 
          matrices <- list()
      
      
      #time_windows <- time_windows[1:10000] #to subset
      
      
      for (w in time_windows) {
        print(w)
        
        #empty matrix for this window
            mat <- matrix(0, nrow = length(all_animals), ncol = length(all_animals))
            rownames(mat) <- all_animals
            colnames(mat) <- all_animals
        
        #extract records for this window
            records <- data_clean %>% filter(window_start == w)
        
        if (nrow(records) == 0) {
          #no records -> matrix stays all zeros
            } else {
              for (i in 1:nrow(records)) {
                a1 <- as.character(records$animal_ID1[i])
                a2 <- as.character(records$animal_ID2[i])
                rssi <- records$RSSI[i]
                
                mat[a1, a2] <- rssi
                mat[a2, a1] <- rssi
              }
            }
        
        matrices[[as.character(w)]] <- mat
      }
      
      #combine into 3D array 
        array_3d <- abind::abind(matrices, along = 3)
        dimnames(array_3d)[[3]] <- as.character(time_windows)
        

# 
# #Smoothing-----------------------------------------------------
# #takes the neighbouring time steps (amount defined by window) and fills out NAs with the closest value (within window)
# #quite slow, and now has been replaced with window ealier in the code
# 
#   smooth_array_with_zero <- function(array_3d, window = 2) { # window is +- seconds. so window=2, means -2 to +2
# 
#       #get dimensions
#           n_animals <- dim(array_3d)[1]
#           n_seconds <- dim(array_3d)[3]
# 
#       #get time labels (seconds)
#          seconds <- as.numeric(dimnames(array_3d)[[3]])
# 
#       #create output array (copy of input)
#          result <- array_3d  # We'll modify this in place
# 
#       #loop over each dyad (i,j) and each time point t
#           for (i in 1:n_animals) {
#             print(i)
#             for (j in 1:n_animals) {
#               print(j)
#               #extract the time series for this dyad (i,j)
#                  values <- result[i, j, ] 
# 
#               #loop over each time point
#                   for (t in 1:n_seconds) {
#                     if (values[t] != 0) next 
# 
#                 #find all non-zero values within time window
#                     t_current <- seconds[t]
#                     t_min <- t_current - window
#                     t_max <- t_current + window
# 
#                 #find indices within time window
#                    valid_indices <- which(seconds >= t_min & seconds <= t_max)
# 
#                 #filter to non-zero values
#                     valid_values <- values[valid_indices]
#                     valid_times <- seconds[valid_indices]
# 
#                 #keep only non-zero values
#                     valid_indices <- valid_indices[valid_values != 0]
#                     valid_times <- valid_times[valid_values != 0]
# 
#                     if (length(valid_indices) == 0) next  # No data in window
# 
#                         #calculate time distances
#                         distances <- abs(valid_times - t_current)
# 
#                         #find minimum distance
#                         min_dist <- min(distances)
# 
#                         #get indices with minimum distance
#                         candidates <- valid_indices[distances == min_dist]
# 
#                         #randomly pick one if multiple
#                         if (length(candidates) > 1) {
#                           chosen_idx <- sample(candidates, 1)
#                 } else {
#                   chosen_idx <- candidates
#                 }
# 
#                 # Get the value from the chosen time
#                 chosen_value <- values[chosen_idx]
# 
#                 # Replace 0 with chosen value
#                 result[i, j, t] <- chosen_value
#           }
#         }
#       }
# 
#       return(result)
#     }
# 
# 
# # Example: smooth with window = 2 seconds
#     smoothed_array <- smooth_array_with_zero(array_3d, window = 2)
# 
# # #check
# #     dyad_values <- smoothed_array["62601", "76927", ]
# #     dyad_values

      
      
#ave (window size in name)
      
    #convert back to posix time
      current_seconds <- as.numeric(dimnames(array_3d)[[3]])
      posix_times <- as.POSIXct(current_seconds, origin = "1970-01-01", tz = "UTC")
      dimnames(array_3d)[[3]] <- as.character(posix_times)
    
    #save
    save(array_3d, file = paste0("proxit_array_ws", window_size, ".RData"))
      

    
#Calculating transitivity (triadic closure)--------------------------------

      #change to binary
          binary_array <- array(0, dim = dim(array_3d), dimnames = dimnames(array_3d))
          binary_array[array_3d != 0] <- 1
      
      #remove self-loops (diagonal), just in case
          for (t in 1:dim(binary_array)[3]) {
            diag(binary_array[, , t]) <- 0
          }
      
      #summarize probability of connection by number of shared neighbors
          n_animals <- dim(binary_array)[1]
          n_seconds <- dim(binary_array)[3]
          animal_names <- dimnames(binary_array)[[1]]
          
      #store results
          results_list <- vector("list", n_seconds)
      
          for (t in 1:n_seconds) {
            A <- binary_array[, , t]
            
            #common neighbors matrix: CN[i,j] = number of shared connections
              CN <- A %*% A
              diag(CN) <- NA  # exclude self
            
            #get upper triangle only (avoid double-counting dyads)
                idx <- which(upper.tri(CN), arr.ind = TRUE)
                
                df <- data.frame(
                  animal1 = animal_names[idx[, 1]],
                  animal2 = animal_names[idx[, 2]],
                  common_neighbors = CN[idx],
                  connected = A[idx],
                  second = t
                )
                
                results_list[[t]] <- df
              }
          

      #combine all time points
         all_results <- bind_rows(results_list)
         

         
#Summarise probability-----------------------------
       summary_df <- all_results %>%
         group_by(common_neighbors) %>%
         summarise(
           n_pairs = n(),
           n_connected = sum(connected),
           prop_connected = mean(connected),
           se = sqrt(prop_connected * (1 - prop_connected) / n_pairs)
         )
       
       saveRDS(all_results, file = paste0("output/proxit_commonNeighbors_ws", window_size, ".rds"))
       
       
#Change array names into posix---------------------------------------
         current_seconds <- as.numeric(dimnames(binary_array)[[3]])
         posix_times <- as.POSIXct(current_seconds, origin = "1970-01-01", tz = "UTC")
         dimnames(binary_array)[[3]] <- as.character(posix_times)

      #save
         save(binary_array, file = paste0("proxit_array_ws", window_size, "_bin.RData"))
         


#Plot —--------------------
#Probability of Connection vs. Shared Neighbors
#double check the plot - Claud based! apparently it's running a glm internally (geom_smooth)?
#limit adjusted to be the max od common neighbours


    png("output/sharedNeighbors_vs_ProbDirConnection.png")
    
        ggplot(summary_df, aes(x = common_neighbors, y = prop_connected)) +
          geom_point(size = 2) +
          geom_errorbar(aes(ymin = prop_connected - 1.96 * se,
                            ymax = prop_connected + 1.96 * se), width = 0.2) +
          geom_smooth(method = "glm", method.args = list(family = "binomial"),
                      aes(weight = n_pairs), se = TRUE, color = "slategrey", fill = "slategrey") +
          labs(
            x = "Number of Shared Neighbors (common connections)",
            y = "Probability of Direct Connection"
          ) +
          scale_x_continuous(limits = c(0, max(summary_df$common_neighbors)), expand = expansion(mult = c(0, 0.05))) +
          scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +  
          theme(
            panel.grid = element_blank(),
            axis.title.x = element_text(size = 16, margin = margin(t = 15)),
            axis.title.y = element_text(size = 16, margin = margin(r = 15)),
            axis.text = element_text(size = 14),
            axis.line = element_line(color = "black"),
            plot.title = element_blank()
          )

    dev.off()
    
    # #checks
    # #range and distribution of common_neighbors
    #   table(all_results$common_neighbors)
    # 
    # #connected
    #   table(all_results$connected)


