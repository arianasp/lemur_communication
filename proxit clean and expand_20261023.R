
library(dplyr)
library(readr)
library(tidyr)
library(fuzzyjoin)
library(lubridate)


##COMBINE DATA----------------------------


    #import
      folder_path <- "C:/Users/kwierucka/Dropbox/lemur_communication/rawdata_csv/proxtit"
      xlog <- read.csv("C:/Users/kwierucka/Dropbox/lemur_communication/rawdata_csv/tagging_log_red.csv")
      tagID <- read.csv("C:/Users/kwierucka/Dropbox/lemur_communication/rawdata_csv/Tag_Proxit_Names_20260411.csv")
      
    
      csv_files <- list.files(path = folder_path, 
                            pattern = "\\.csv$", 
                            full.names = TRUE, 
                            recursive = TRUE,
                            include.dirs = T)
    
    #combine
      combined_df <- csv_files %>%
        lapply(function(file) {
          df <- read_csv(file)
          df$source_file <- paste0(basename(dirname(file[1])), "/"  , basename(file[1]))
          return(df)
        }) %>% bind_rows()
    
    

#FILTER for times the animals were wearing the collar----------------------
  
      xlog <- xlog %>% mutate(
          on_datetime  = as.POSIXct(paste(on_year, on_month, on_day, on_time),
                                    format = "%Y %m %d %H:%M", tz = "UTC"),
          off_datetime = as.POSIXct(paste(off_year, off_month, off_day, off_time),
                                    format = "%Y %m %d %H:%M", tz = "UTC")
        )
      
  #build lookup df
      proxit_vesper <- tagID %>% filter(!is.na(Proxit)) %>% select(Proxit, Vesper) %>% distinct()
      
  #vesper tag to Proxit data 
      combined_df <- combined_df %>%
        left_join(proxit_vesper, by = c("Unit1" = "Proxit")) %>% rename(Vesper1 = Vesper) %>%
        left_join(proxit_vesper, by = c("Unit2" = "Proxit")) %>% rename(Vesper2 = Vesper)
      
  #vesper to animal_ID1
      xlog_small <- xlog %>% select(Vesper, on_datetime, off_datetime, animal_ID)
      
      combined_df <- combined_df %>%
        fuzzy_left_join(
          xlog_small,
          by = c("Vesper1" = "Vesper",
                 "Start Time" = "on_datetime",
                 "Start Time" = "off_datetime"),
          match_fun = list(`==`, `>=`, `<=`)
        ) %>% rename(animal_ID1 = animal_ID) %>%
        select(-Vesper, -on_datetime, -off_datetime)
      
  #vesper to animal_ID2
      combined_df <- combined_df %>%
        fuzzy_left_join(
          xlog_small,
          by = c("Vesper2" = "Vesper",
                 "Start Time" = "on_datetime",
                 "Start Time" = "off_datetime"),
          match_fun = list(`==`, `>=`, `<=`)
        ) %>%
        rename(animal_ID2 = animal_ID) %>%
        select(-Vesper, -on_datetime, -off_datetime)
      
  #clean
      combined_df <- combined_df %>%
        filter(!is.na(animal_ID1), !is.na(animal_ID2))
  
  #save
      saveRDS(combined_df, "output/clean_proxit.rds")
      

#EXPAND----------------------------------------------------------------
#to have each second in separate row
      
    #temp_df <- combined_df[1:10000,] #for testing
    
      expanded_df <- combined_df %>%
        mutate(row_id = row_number()) %>%                     
        uncount(`Seconds Duration` + 1, .id = "offset") %>%    #duplicate rows
        mutate(offset = offset - 1,                           #make offset start at 0
               `Start Time` = `Start Time` + offset) %>%
        select(-row_id, -offset)    
    
    #save
     saveRDS(expanded_df, "output/expanded.rds")
    
    
    #clean up duplicates (keeping the last one (based on time))
      cleaned_df <- expanded_df %>%
        group_by(across(-RSSI)) %>%
        slice_tail(n = 1) %>%
        ungroup()
    
    #save
      saveRDS(cleaned_df, "output/cleaned.rds")


# #dyadic RSSI-------------------------------
# #for testing two way interactions between tags. 
# #code works, but isn't useful, sience tags don't log data both ways in the same second.
#       
# 
# dyad_df <- expanded_df %>%
#   mutate(
#     Animal1 = pmin(animal_ID1, animal_ID2),
#     Animal2 = pmax(animal_ID1, animal_ID2),
#     direction = if_else(animal_ID1 == Animal1, "RSSI1", "RSSI2")
#   ) %>%
#   group_by(Animal1, Animal2, `Start Time`, direction) %>%
#   summarise(RSSI = mean(RSSI, na.rm = TRUE), .groups = "drop") %>%   
#   pivot_wider(
#     id_cols = c(Animal1, Animal2, `Start Time`),
#     names_from = direction,
#     values_from = RSSI
#   )
# 
# 
# #check diference in measurements
# dyad_df$RSSI_dif <- dyad_df$RSSI2 - dyad_df$RSSI1 
# mean(dyad_df$RSSI_dif) 