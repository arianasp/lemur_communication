#Playing around with moan events
Sys.setenv(TZ='UTC')

library(fields)

#PARAMS
call_types <- c('m','w')
call_cols <- c('darkblue','darkred')
chain_time <- 4 #time for chaining together events - in seconds
cols_sex <- c('darkblue','darkred')


#DIRECTORY
setwd('~/Dropbox/lemur_communication/rawdata_csv/')

#LOAD DATA
ids <- read.csv('ids_matrix.csv')
calls <- read.csv('calls_df_new.csv')

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
date_plot <- '2026-04-14'
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
events_mat <- matrix(0, nrow = nrow(ids), ncol = nrow(events))
for(i in 1:nrow(events)){
  particip <- events$participants[i][[1]]
  idxs <- which(ids$id_code %in% particip)
  events_mat[idxs,i] <- 1
}

#summary statistics
ids$events_participated <- ids$events_initiated <- ids$tot_calls <- NA
n_calls <- table(calls$animal_ID)
ids$tot_calls <- as.numeric(n_calls[match(as.character(ids$id_code), names(n_calls))])
n_inits <- table(events$initiator)
ids$events_initiated <- as.numeric(n_inits[match(as.character(ids$id_code), names(n_inits))])
ids$events_participated <- rowSums(events_mat)

#sort by sex and then age
ids_sort <- ids[order(ids$sex, ids$age, decreasing = T),]

#make barplots for total number of calls, event participation, and event initiation
quartz()
par(mfrow=c(1,3))
barplot(ids_sort$events_initiated, horiz=T, names.arg = ids_sort$id_code, las = 2, col = ids_sort$color, cex.names = .8, xlab = '# events initiated', cex.lab=1.5)
barplot(ids_sort$events_participated, horiz=T, names.arg = ids_sort$id_code, las = 2, col = ids_sort$color, cex.names = .8, xlab = '# events participated', cex.lab=1.5)
barplot(ids_sort$tot_calls, horiz=T, names.arg = ids_sort$id_code, las = 2, col = ids_sort$color, cex.names = .8, xlab = 'total # calls', cex.lab=1.5)

#probability of response given a particular initiator
#row = initiator, col = responder
p_resp <- matrix(NA, nrow = nrow(ids), ncol = nrow(ids))
for(i in 1:nrow(ids)){
  events_ind <- events_mat[,which(events$initiator==ids$id_code[i])]
  events_ind[i,] <- NA #don't respond to yourself
  if(ncol(events_ind)>0){
    p_resp[i,] <- rowMeans(events_ind)
  }
}

p_responded_to <- rowMeans(p_resp, na.rm=T)
p_responds <- colMeans(p_resp, na.rm=T)

par(mfrow=c(1,1))
image.plot(p_resp)
     