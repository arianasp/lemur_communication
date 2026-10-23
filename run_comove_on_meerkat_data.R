source('get_comoves.R')

load('~/Dropbox/lemur_communication/Mapping vocal interactions in space and time differentiates signal broadcast versus signal exchange in meerkat groups/data/HM2017_COORDINATES_all_sessions_anon.RData')

#---Generate simulated magnetometer data (from GPS headings) and move/no move data (from GPS speed)--
speed_thresh <- 0.2 #threshol to consider a meerkat 'moving'

xs <- HM2017_allX
ys <- HM2017_allY

n_inds <- nrow(xs)
n_times <- ncol(xs)

heads <- moving <- matrix(NA, n_inds, n_times)
for(i in 1:n_inds){
  out <- cocomo::get_heading_and_speed_temporal(xs[i,],ys[i,],t_window = 10)
  heads[i,] <- out$heads
  moving[i,] <- out$speeds > speed_thresh
}

#---Run comove analysis to get events
agg_events_long <- get_comoves(heads, moving, max_head_diff = 0.2, min_comove_duration = 5, min_event_duration = 10)

#---Plot a single event
row <- 6
inds <- agg_events_long$individuals[row][[1]]
t0 <- agg_events_long$start_time[row]
tf <- agg_events_long$stop_time[row]

xmin <- min(xs[,t0:tf],na.rm=T)
ymin <- min(ys[,t0:tf],na.rm=T)
xmax <- max(xs[,t0:tf],na.rm=T)
ymax <- max(ys[,t0:tf],na.rm=T)
cols <- rainbow(n_inds)
cols[which(!(1:n_inds %in% inds))] <- 'gray'
plot(NULL, xlim = c(xmin,xmax), ylim = c(ymin,ymax), xlab = 'Easting', ylab = 'Northing', asp = 1)
for(i in 1:n_inds){
  lines(xs[i,t0:tf], ys[i, t0:tf], col = cols[i], lwd = 2)
}

