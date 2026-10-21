#Make an animation of GPS data and calls
library(cocomo)

setwd('~/Desktop/')

#LOAD
xs <- readRDS('Matrixes_Updated/xs.rds')
ys <- readRDS('Matrixes_Updated/ys.rds')
timestamps <- readRDS('Matrixes_Updated/timestamps.rds')
calls <- readRDS('Matrixes_Updated/calls_df.rds')
ids <- readRDS('Matrixes_Updated/idx.rds')

t <- '2026-04-14 09:22:00'
idx <- which(timestamps==t)

cocomo::generate_movement_and_calls_visualization(xs = xs, ys = ys, timestamps = timestamps, calls = calls, start_time = idx, end_time = idx + 100, output_dir = '~/Desktop/test_anim')


