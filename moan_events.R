#This R script processes animal communication and proximity data obtained by Vesper tags and produces a event-based dataset linking 
# vocal initiation and response events and spatial proximity.


#Data Loading & Preprocessing
##Loads call data (vocalizations) and individual metadata (sex, age)
##Filters calls by type (male/female), orders by time, and computes inter-call intervals

#Event Detection
##Groups calls into "events" using a 4-second time threshold
##Identifies the first caller as the initiator of each event
##Tracks participants, call counts, and response patterns

#Proximity Data Integration
##Loads a binary proximity array (5-second intervals) indicating dyad connectivity
##For each call, finds the nearest proximity timestamp and checks if the initiator-responder dyad was connected
##Adds self-responses (animals responding to themselves) with connected_at_call = 1

#Analysis & Visualization
##Creates event_pairs_full: a comprehensive dataframe with call timing, connectivity status, and individual attributes
##Generates a bar plot showing response rates by dyad connectivity (connected vs. not connected) and initiator sex



# Last modified: 23.10.2026 Lily Broos



Sys.setenv(TZ='UTC')

library(fields)


#PARAMS
call_types <- c('m','w')
call_cols <- c('darkblue','darkred')
chain_time <- 4 #time for chaining together events - in seconds
cols_sex <- c('darkblue','darkred')


#DIRECTORY
#setwd('~/Dropbox/lemur_communication')

#LOAD DATA
ids <- read.csv("C:/Users/Lily/Dropbox/lemur_communication/rawdata_csv/ids_matrix.csv")
calls <- read.csv("C:/Users/Lily/Dropbox/lemur_communication/Matrixes_Updated/calls_df_20261021.csv")

#PROCESS
calls$t <- as.POSIXct(calls$time, tz = 'UTC')
calls <- calls[which(calls$call_type %in% call_types),] #get only specified call types
calls$col <- call_cols[match(calls$call_type, call_types)]
calls <- calls[order(calls$t),]

#make male/female all lower case
ids$sex <- tolower(ids$sex)

#add colors column to ids data frame
ids$color <- NA
male_idxs <- which(ids$sex=='male' & !is.na(ids$age))
ids$color[male_idxs] <- sapply(ids$age[male_idxs]/max(ids$age,na.rm=T), function(x){return(adjustcolor(cols_sex[1], alpha.f = x))})
female_idxs <- which(ids$sex=='female' & !is.na(ids$age))
ids$color[female_idxs] <- sapply(ids$age[female_idxs]/max(ids$age,na.rm=T), function(x){return(adjustcolor(cols_sex[2], alpha.f = x))})
ids$color[which(is.na(ids$color))] <- 'gray'

#inter-call-interval
calls$time_to_last <- NA
calls$time_to_last[2:nrow(calls)] <- calls$t[2:nrow(calls)] - calls$t[1:(nrow(calls)-1)]

#PLOT VISUALIZATION
date_plot <- '2026-04-15'
tmp <- calls[which(as.Date(calls$t)==date_plot),]
t0 <- min(tmp$t)
dt <- (tmp$t - t0) / 60 / 60
plot(NULL, xlim = c(min(dt),max(dt)), ylim = c(0, max(tmp$ind_idx)))
arrows(dt, tmp$ind_idx-1, dt, tmp$ind_idx, lwd = 1, length = 0, col = tmp$col)

#Find events
calls$event_id <- NA
calls$event_id[1] <- 1
curr_id <- 1
for(i in 2:nrow(calls)){
  if(calls$time_to_last[i] > chain_time){
    curr_id <- curr_id + 1
  }
  calls$event_id[i] <- curr_id
}

#Number of calls per event
events_table <- table(calls$event_id)
histo <- hist(events_table, breaks=seq(1,max(events_table),1), plot = F)
plot(histo$mids, cumsum(histo$counts)/sum(histo$counts), xlab = '# calls', ylab = 'Cumulative probability', log = 'xy', type = 'l', lwd = 4)
plot(histo$mids, (1-cumsum(histo$counts)/sum(histo$counts)), xlab = 'N calls', ylab = 'P(> N calls)', log = 'x', type = 'l', lwd = 4)

#Make a table of events
events <- data.frame(event_id = seq(1,length(events_table)))
events$tf <- events$t0 <- as.POSIXct(NA, tz = 'UTC')
events$initiator <- events$n_calls <- events$n_callers <- NA
events$participants <- rep(list(c(NA)))
for(i in 1:nrow(events)){
  tmp <- calls[which(calls$event_id == events$event_id[i]),]
  events$t0[i] <- min(tmp$t) #event start time
  events$tf[i] <- max(tmp$t) #event_end_time
  events$initiator[i] <- tmp$animal_ID[1] #event initiator id
  events$n_calls[i] <- nrow(tmp) #number of calls
  events$n_callers[i] <- length(unique(tmp$animal_ID)) #number of unique callers
  events$participants[i] <- list(unique(tmp$animal_ID))
}

#events matrix (rows = inds, cols = events, T = participated, F = did not)
events_mat <- response_mat <- matrix(0, nrow = nrow(ids), ncol = nrow(events))
for(i in 1:nrow(events)){
  particip <- events$participants[i][[1]]
  init <- which(ids$id_code == events$initiator[i])
  idxs <- which(ids$id_code %in% particip)
  events_mat[idxs,i] <- 1
  idxs_resp <- setdiff(idxs, init)
  response_mat[idxs_resp,i] <- 1
}

#summary statistics
ids$events_participated <- ids$events_initiated <- ids$tot_calls <- NA
n_calls <- table(calls$animal_ID)
ids$tot_calls <- as.numeric(n_calls[match(as.character(ids$id_code), names(n_calls))])
n_inits <- table(events$initiator)
ids$events_initiated <- as.numeric(n_inits[match(as.character(ids$id_code), names(n_inits))])
ids$events_participated <- rowSums(events_mat)
ids$events_responded <- rowSums(response_mat)

#probability of response given a particular initiator
#row = initiator, col = responder
p_resp <- matrix(NA, nrow = nrow(ids), ncol = nrow(ids))
for(i in 1:nrow(ids)){
  events_ind <- events_mat[,which(events$initiator==ids$id_code[i])]
  events_ind[i,] <- NA #don't count responding to yourself
  if(ncol(events_ind)>0){
    p_resp[i,] <- rowMeans(events_ind)
  }
}

#Overall probability for each initiator was responded to (averaged across all responders)
ids$p_responded_to <- rowMeans(p_resp, na.rm=T)

#Overall probability that each individual responded to someone else (averaged across all initiators)
ids$p_responds <- colMeans(p_resp, na.rm=T)

#sort by sex and then age
ids_sort <- ids[order(ids$sex, ids$age, decreasing = T),]
p_resp_sort <- p_resp[order(ids$sex, ids$age, decreasing = T),order(ids$sex, ids$age, decreasing = T)]



#PLOTS

#make barplots for total number of calls, event participation, and event initiation
dev.off()
#quartz()
par(mfrow=c(1,3))
barplot(ids_sort$events_initiated, horiz=T, names.arg = ids_sort$id_code, las = 2, col = ids_sort$color, cex.names = .8, xlab = '# events initiated', cex.lab=1.5)
barplot(ids_sort$events_participated, horiz=T, names.arg = ids_sort$id_code, las = 2, col = ids_sort$color, cex.names = .8, xlab = '# events participated', cex.lab=1.5)
barplot(ids_sort$tot_calls, horiz=T, names.arg = ids_sort$id_code, las = 2, col = ids_sort$color, cex.names = .8, xlab = 'total # calls', cex.lab=1.5)

dev.off()
#quartz()
image.plot(p_resp_sort, xaxt = 'n', yaxt = 'n', col = viridis(256))
Map(axis, side=2, at=seq(0,1,length.out = nrow(p_resp_sort)), col.axis=ids_sort$color, labels=ids_sort$id, lwd=0, las=1)
Map(axis, side=1, at=seq(0,1,length.out = nrow(p_resp_sort)), col.axis=ids_sort$color, labels=ids_sort$id_code, lwd=0, las=2)






# Create event based dataframe ----
## 1. Base: one row per call, with event info -----
event_pairs <- calls[, c("event_id", "animal_ID", "t", "time_to_last")]

# Rename for consistency
names(event_pairs)[names(event_pairs) == "animal_ID"] <- "id_code"
names(event_pairs)[names(event_pairs) == "t"] <- "time_call"

# Add date
event_pairs$date <- as.Date(event_pairs$time_call)

# Initialize columns
event_pairs$initiator <- 0
event_pairs$responder <- 1
event_pairs$initial_responder <- FALSE
event_pairs$dyadic_rssi <- NA_real_
event_pairs$dominance_rank <- NA_real_
event_pairs$sex <- NA_character_
event_pairs$age <- NA_character_


## 2. Mark the initiator for each event----
# (first call in event, by time)
for (i in 1:nrow(events)) {
  ev_id <- events$event_id[i]
  init_id <- events$initiator[i]
  
  # Subset calls for this event, sorted by time_call
  ev_calls <- event_pairs[event_pairs$event_id == ev_id, ]
  ev_calls <- ev_calls[order(ev_calls$time_call), ]  # ensure sorted
  
  # First call = initiator
  first_idx <- 1  # because sorted
  idx <- which(
    event_pairs$event_id == ev_id &
      event_pairs$id_code == init_id &
      event_pairs$time_call == ev_calls$time_call[first_idx]
  )
  
  if (length(idx) == 1) {
    event_pairs$initiator[idx] <- 1
    event_pairs$responder[idx] <- 0
  }
}


## 3. Mark the *first responder* ----
## (first non-initiator call *after* initiator)
for (i in 1:nrow(events)) {
  ev_id <- events$event_id[i]
  init_id <- events$initiator[i]
  
  # Subset event
  ev_calls <- event_pairs[event_pairs$event_id == ev_id, ]
  ev_calls <- ev_calls[order(ev_calls$time_call), ]
  
  # Get initiator's call time
  init_row <- ev_calls[ev_calls$initiator == 1, ]
  if (nrow(init_row) == 0) next  # no initiator found — skip event
  
  time_i <- init_row$time_call[1]  # should be only one
  
  # Find first non-initiator call *after* time_i
  non_init <- ev_calls$initiator == 0
  resp_calls <- ev_calls[non_init, ]
  
  if (nrow(resp_calls) == 0) next  # no responders — skip
  
  valid_times <- resp_calls$time_call[resp_calls$time_call > time_i]
  if (length(valid_times) == 0) next  # no calls after initiator — skip
  
  first_resp_time <- min(valid_times)
  
  # Mark that row as initial_responder
  idx <- which(
    event_pairs$event_id == ev_id &
      event_pairs$time_call == first_resp_time &
      event_pairs$initiator == 0  # ensure it's not the initiator
  )
  
  
  # Handle ties (multiple calls at same time) — mark all
  event_pairs$initial_responder[idx] <- TRUE
}



## 4. Merge individual attributes ----
## (ids uses 'animal_ID' as ID column)
event_pairs <- merge(
  event_pairs,
  ids[, c("id_code", "sex", "age")],  # ← use 'animal_ID' here
  by = "id_code",                      # ← left side uses 'id_code'
  #by.y = "animal_ID",                    # ← right side uses 'animal_ID'
  all.x = TRUE
)

# Rename age.y → age, sex.y → sex, and drop age.x, sex.x
names(event_pairs)[names(event_pairs) == "age.y"] <- "age"
names(event_pairs)[names(event_pairs) == "sex.y"] <- "sex"

# Remove .x columns if they exist
event_pairs <- event_pairs[, !names(event_pairs) %in% c("age.x", "sex.x", "dyadic_rssi")]



## 6. Add proxtit data to df ----
library(dplyr)
library(data.table)
setDT(event_pairs)


#proxtit <- readRDS(file.path("C:/Users/Lily/Dropbox/lemur_communication/rawdata_csv/proxtit/cleaned_proxit.rds"))
prox_array_bin <- file.path("C:/Users/Lily/Dropbox/lemur_communication/processed/proxit_array_ws5_bin.RData")
load(prox_array_bin)

# Check data 
head(binary_array)
str(binary_array)

## Start time of binary array 13.04.2026 09:53:35 - 16.04.2026 16:27:17 with increments of 5 seconds 

### a. Calculate timeframe of proximity data----
# Define start and end times
start_time <- as.POSIXct("13.04.2026 09:53:35", format = "%d.%m.%Y %H:%M:%S", tz = "Europe/Zurich")
end_time <- as.POSIXct("16.04.2026 16:27:17", format = "%d.%m.%Y %H:%M:%S", tz = "Europe/Zurich")
interval_seconds <- 5

# Calculate duration and expected number of points
duration_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
expected_n <- floor(duration_seconds / interval_seconds) + 1
last_expected_time <- start_time + (expected_n - 1) * interval_seconds

#########


# Create dyads_per_event
initiators <- event_pairs[initiator == 1, .(event_id, initiator_id = id_code)]
responder_pairs <- event_pairs[responder == 1, .(event_id, responder_id = id_code)]
dyads_per_event <- initiators[responder_pairs, on = "event_id"]
dyads_per_event[, dyad_key := paste(pmin(initiator_id, responder_id), 
                                    pmax(initiator_id, responder_id), 
                                    sep = "_")]

##########
# Get unique (event, initiator, responder) dyads, excluding self-responses
dyads_per_event_unique <- dyads_per_event[
  initiator_id != responder_id,  # exclude self-responses
  .(event_id, initiator_id, responder_id)
]

# Explode into one row per responder (already done, but ensure uniqueness)
dyads_per_event_all <- unique(dyads_per_event_unique[
  , .(responder_id = responder_id),
  by = .(event_id, initiator_id)
])

# Now join to responder rows
event_pairs_test <- event_pairs[responder == 1]

event_pairs_dyads_test <- event_pairs_test[
  dyads_per_event_all, 
  on = c("event_id", "id_code" = "responder_id"), 
  allow.cartesian = TRUE
]

### b. Get connected dyads at call initiation-----
time_names <- as.numeric(dimnames(binary_array)[[3]])
id_to_idx <- setNames(seq_along(dimnames(binary_array)[[1]]), dimnames(binary_array)[[1]])

# Convert time_call to numeric timestamp
event_pairs_dyads_test[, time_call_ts := as.numeric(
  #as.POSIXct(
    time_call, tz = "UTC")
 # )
]

# Find nearest timestamp index (vectorized safely)
event_pairs_dyads_test[, time_idx := sapply(time_call_ts, function(x) which.min(abs(time_names - x)))]

# Get dyad indices (sorted order matches binary_array)
event_pairs_dyads_test[, dyad1_idx := sapply(pmin(initiator_id, id_code), function(x) id_to_idx[[as.character(x)]])]
event_pairs_dyads_test[, dyad2_idx := sapply(pmax(initiator_id, id_code), function(x) id_to_idx[[as.character(x)]])]

# Extract binary value
event_pairs_dyads_test[, connected_at_call := 
                         mapply(function(d1, d2, t) binary_array[d1, d2, t],
                                dyad1_idx, dyad2_idx, time_idx)]

# Handle missing dyads (NA in binary_array) → set to 0 for now
event_pairs_dyads_test[is.na(connected_at_call), connected_at_call := 0]

# Add the proxit start time corresponding to each time index
event_pairs_dyads_test[,  prox_time_start := as.POSIXct(
  as.numeric(start_time) + (as.integer(time_idx) - 1) * 5,
  origin = "1970-01-01", tz = "Europe/Zurich")]


### c. Add self-responses ------
# Get self-responder rows (initiator == responder)
self_responders <- event_pairs[
 # event_id == test_event_id & 
    responder == 1 & 
    id_code == event_pairs[#event_id == test_event_id & 
      initiator == 1, id_code][1]
]


# Add the new columns (copy from event_pairs_dyads_test structure)
self_responders[, `:=`(
  initiator_id = id_code,  # self-responder: initiator_id = id_code
  connected_at_call = 1,
  time_call_ts = as.numeric(as.POSIXct(time_call, tz = "UTC")),
  time_idx = NA_integer_,
  dyad1_idx = NA_integer_,
  dyad2_idx = NA_integer_
)]


# Work on a copy of the original dataset
event_pairs_full <- copy(event_pairs)

# Create the new columns if they do not already exist
new_cols <- c("initiator_id", "time_call_ts", "time_idx", "dyad1_idx", "dyad2_idx", "connected_at_call")

for (col in new_cols) {
  if (!col %in% names(event_pairs_full)) {
    event_pairs_full[, (col) := NA_real_]
  }
}

# Ensure the join table contains only the required columns
update_data <- unique(
  event_pairs_dyads_test[,
    .(
      event_id,
      id_code,
      time_call,
      initiator_id,
      time_call_ts,
      time_idx,
      dyad1_idx,
      dyad2_idx,
      connected_at_call,
      prox_time_start
    )
  ],
  by = c("event_id", "id_code", "time_call")
)

# LEFT JOIN: update matching rows in the original dataset
event_pairs_full[
  update_data,
  on = .(event_id, id_code, time_call),
  `:=`(
    initiator_id = i.initiator_id,
    time_call_ts = i.time_call_ts,
    time_idx = i.time_idx,
    dyad1_idx = i.dyad1_idx,
    dyad2_idx = i.dyad2_idx,
    connected_at_call = i.connected_at_call,
    prox_time_start = i.prox_time_start
  )
]

# Identify initiators for the test event
initiators_in_event <- unique(event_pairs[initiator == 1, id_code])


# Update existing self-responder rows directly in event_pairs_full
event_pairs_full[  #event_id == test_event_id &
    responder == 1 &
    id_code %in% initiators_in_event,
  `:=`(
    initiator_id = id_code,
    time_call_ts = as.numeric(as.POSIXct(time_call, tz = "UTC")),
    time_idx = NA_integer_,
    dyad1_idx = NA_integer_,
    dyad2_idx = NA_integer_,
    connected_at_call = 1L,
    prox_time_start = as.POSIXct(NA, tz = "Europe/Zurich")
    
  )
]

# Identify self-responses within each event
self_dyads <- unique(
  dyads_per_event[
    initiator_id == responder_id,
    .(event_id, id_code = responder_id)
  ]
)

# Set initiator_id and connected_at_call for self-responders
event_pairs_full[
  responder == 1 &
    id_code %in% initiators_in_event,  # self-responders are those who initiated
  `:=`(
    initiator_id = id_code,  # self: initiator_id = id_code
    connected_at_call = 1L
  )
]


