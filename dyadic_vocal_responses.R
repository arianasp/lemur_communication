# Dyadic vocal analysis (23 Oct 2026)
# 1. Data: Calls from 15 Apr 2026, including only 'w' and 'm' calls
# 2. Response window: Start = half mean call duration; End = 90th percentile latency (~4.1 sec)
# 3. Dyadic events: For each initiator call, found first response from another individual within window
# 4. Matrices: Standard, caller-proportional, and responder-proportional call-response matrices
# 5. Analyses: 
#    - Call-type proportions (initial vs. response, per animal & per sex)
#    - Sex-specific responses (counts and likelihood of response by caller/responder sex)
# 6. Key output: dyadic_events dataframe with latency, sex, age, call types


#Last update: 22.10.2026 Lily Broos

library(readr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)



Sys.setenv(TZ='UTC')

#calls_df <- read_csv("C:/Users/Lily/Dropbox/lemur_communication/rawdata_csv/calls_df_new.csv")
calls_df_end <- read_csv("D:/Satyam_Data/Combined_Processed_Data/calls_df_20261021.csv")
ids_matrix <- read_csv("C:/Users/Lily/Dropbox/lemur_communication/rawdata_csv/ids_matrix.csv")

View(calls_df_end)



# Choose the day to visualize
selected_day <- as.Date("2026-04-15")

# Filter data
calls_day <- calls_df_end %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S"),
         animal_ID = as.character(animal_ID)) %>%
  filter(
    call_type != "h",
    as.Date(time) == selected_day
  )%>%
  mutate(
    animal_ID = factor(animal_ID, levels = sort(unique(animal_ID)))
  )


calls_15.04 <- ggplot(
  calls_day,
  aes(x = time, y = as.numeric(animal_ID), color = call_type)
) +
  geom_segment(
    aes(
      xend = time,
      y = as.numeric(factor(animal_ID)) - 0.35,
      yend = as.numeric(factor(animal_ID)) + 0.35
    ),
    linewidth = 0.8
  ) +  
  scale_y_continuous(
    breaks = seq_along(levels(calls_day$animal_ID)),
    labels = levels(calls_day$animal_ID)
  ) +
  labs(
    title = paste("Animal calls on", selected_day),
    x = "Time",
    y = "Animal ID",
    color = "Call type"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(color = "black")
  )


plot(calls_15.04)

# Determine start response window
## Determine average length of calls
calls_df_end <- calls_df_end %>%
  mutate(
    duration_sec = as.numeric(
      difftime(end_time, time, units = "secs")
    )
  )

# Calculate average duration per call type
duration_by_type <- calls_df_end %>%
  group_by(call_type) %>%
  summarise(
    n_calls = sum(!is.na(duration_sec)),
    mean_duration_sec = mean(duration_sec, na.rm = TRUE),
    median_duration_sec = median(duration_sec, na.rm = TRUE),
    sd_duration_sec = sd(duration_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(call_type)




# Define dyadic events------

# Settings
selected_day <- as.Date("2026-04-15")
response_window_end <- 10

# Response window start = half the mean call duration, by call type
response_windows <- duration_by_type %>%
  mutate(
    response_window_start = mean_duration_sec / 2
  ) %>%
  select(call_type, response_window_start)

# Add response window start to calls_df
calls_df_end <- calls_df_end %>% left_join(response_windows, by = "call_type")


# Filter calls for 15 April
calls_15.04 <- calls_df_end %>%
  mutate(
    time = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S"),
    animal_ID = as.character(animal_ID)
  ) %>%
  filter(
    as.Date(time) == selected_day,
    call_type != "h"
  ) %>%
  arrange(time) %>%
  mutate(
    event = row_number(),
    caller_id = animal_ID,
    time_start = time + seconds(response_window_start),
    time_end = time + seconds(response_window_end)
  )
# Find calls by different individuals within the response window
dyadic_events <- calls_15.04 %>%
  transmute(
    event,
    time_init = time,
    time_end,
    caller_id
  ) %>%
  left_join(
    calls_15.04 %>%
      transmute(
        time_resp = time,
        responder_id = animal_ID
      ),
    by = join_by(
      time_init < time_resp,
      time_end >= time_resp
    ),
    relationship = "many-to-many"
  ) %>%
  group_by(event, time_init, caller_id) %>%
  slice_min(time_resp, n = 1, with_ties = FALSE, na_rm = FALSE) %>%
  ungroup() %>%
  mutate(
    latency = as.numeric(
      difftime(time_resp, time_init, units = "secs")
    )
  ) %>%
  select(
    event,
    time_init,
    caller_id,
    responder_id,
    latency
  )

# Calculate optimal response window ----
hist(dyadic_events$latency)

# Keep only events with a detected response
latencies <- dyadic_events %>%
  filter(!is.na(latency)) %>%
  pull(latency)

# Calculate the 95th percentile
quantile(latencies, probs = c(0.50, 0.75, 0.90, 0.95, 0.99))

# function to get 90th percentile for response_window_seconds
response_window_end90 <- quantile(
  dyadic_events$latency,
  probs = 0.90,
  na.rm = TRUE
) %>%
  as.numeric()

print(response_window_end90)

## Then rerun code to get the dyadic events
calls_15.04 <- calls_df_end %>%
  mutate(
    time = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S"),
    animal_ID = as.character(animal_ID)
  ) %>%
  filter(
    as.Date(time) == selected_day,
    call_type != "h"
  ) %>%
  arrange(time) %>%
  mutate(
    event = row_number(),
    caller_id = animal_ID,
    time_start = time + seconds(response_window_start),
    time_end = time + seconds(response_window_end90)
  )

dyadic_events <- calls_15.04 %>%
  transmute(
    event,
    time_init = time,
    time_end,
    caller_id,
    call_type_c = call_type  # Initial caller's call type
  ) %>%
  left_join(
    calls_15.04 %>%
      transmute(
        time_resp = time,
        responder_id = animal_ID,
        call_type_r = call_type  # Responder's call type
      ),
    by = join_by(
      time_init < time_resp,
      time_end >= time_resp
    ),
    relationship = "many-to-many"
  ) %>%
  group_by(event, time_init, caller_id) %>%
  slice_min(time_resp, n = 1, with_ties = FALSE, na_rm = FALSE) %>%
  ungroup() %>%
  mutate(
    latency = as.numeric(
      difftime(time_resp, time_init, units = "secs")
    )
  ) %>%
  select(
    event,
    time_init,
    caller_id,
    responder_id,
    latency,
    call_type_c,
    call_type_r
  )



# Add individual information
ids_df <- ids_matrix %>%
  rename(animal_id = id_code)%>%
  mutate(animal_id = as.factor(animal_id))

dyadic_events <- dyadic_events %>%
  left_join(
    ids_df %>% 
      select(animal_id, sex, age),
    by = c("caller_id" = "animal_id")
  ) %>%
  rename(sex_c = sex, age_c = age)

dyadic_events <- dyadic_events %>%
  left_join(
    ids_df %>% 
      select(animal_id, sex, age),
    by = c("responder_id" = "animal_id")
  ) %>%
  rename(sex_r = sex, age_r = age)



# 3. Inspect results
print(dyadic_events)
hist(dyadic_events$age_c)
hist(dyadic_events$age_r)
hist(dyadic_events$latency)




# Call-response matrices-----

## Standard call-response matrix -----

# Count responses for each caller-responder pair
response_matrix <- dyadic_events %>%
  filter(!is.na(responder_id))%>%
  count(caller_id, responder_id, name = "n_responses") %>%
  complete(
    caller_id,
    responder_id,
    fill = list(n_responses = 0)
  )

dyadic_events %>%
  count(call_type_c, call_type_r)

# Plot the matrix as a heatmap
ggplot(response_matrix,
       aes(x = responder_id, y = caller_id, fill = n_responses)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_responses), size = 3) +
  scale_fill_gradient(low = "white", high = "turquoise4") +
  labs(
    title = "Call-response matrix",
    subtitle = "15 April 2026 | Responses within 10 seconds",
    x = "Responder ID",
    y = "Caller ID",
    fill = "Number of\nresponses"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  ) +
  coord_fixed()
## I.E. how many times each animal's call was followed by a call from another animal within the response window


## Proportional call-response matrix based on caller_ID-----
# Total calls produced by each animal
call_totals <- dyadic_events %>%
  count(caller_id, name = "total_calls") %>%
  rename(animal_id = caller_id)

call_totals

# Count responses for each caller–responder pair and divide by the total number of calls produced by each responder.
response_matrix <- dyadic_events %>%
  filter(!is.na(responder_id)) %>%
  count(caller_id, responder_id, name = "n_responses") %>%
  complete(
    caller_id,
    responder_id,
    fill = list(n_responses = 0)
  ) %>%
  left_join(
    call_totals,
    by = c("caller_id" = "animal_id")
  ) %>%
  mutate(
    proportional_response = n_responses / total_calls
  )

# Plot proportional response matrix
ggplot(
  response_matrix,
  aes(
    x = responder_id,
    y = caller_id,
    fill = proportional_response
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = ifelse(
      is.na(proportional_response),
      "NA",
      sprintf("%.2f", proportional_response)
    )),
    size = 3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "turquoise4",
    na.value = "grey80"
  ) +
  labs(
    title = "Proportional call-response matrix",
    subtitle = "15 April 2026 | Responses within 4.1 seconds",
    x = "Responder ID",
    y = "Caller ID",
    fill = "Responses / \ntotal calls"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  ) +
  coord_fixed()

## e.g. 0.2 would mean that out of all calls produced by an individual, a responder responded to 20% of them


## Proportional call_response matrix based on responder ----
response_matrix <- dyadic_events %>%
  filter(!is.na(responder_id)) %>%
  count(caller_id, responder_id, name = "n_responses") %>%
  group_by(responder_id) %>%
  mutate(
    total_responses = sum(n_responses),
    proportional_response = n_responses / total_responses
  ) %>%
  ungroup() %>%
  complete(
    caller_id,
    responder_id,
    fill = list(
      n_responses = 0,
      total_responses = 0,
      proportional_response = 0
    )
  )


# Plot matrix
ggplot(
  response_matrix,
  aes(
    x = responder_id,
    y = caller_id,
    fill = proportional_response
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.2f", proportional_response)),
    size = 3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "turquoise4"
  ) +
  labs(
    title = "Proportional distribution of responses",
    subtitle = "15 April 2026 | Responses within 4.1 seconds",
    x = "Responder ID",
    y = "Caller ID",
    fill = "Proportion of\nall responses"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  ) +
  coord_fixed()
## E.g. 0.2 would mean that out of all responses given by a responder, 20% were given to this caller_ID



# Response type per animal-----
library(scales)

response_types <- dyadic_events %>%
  filter(
    !is.na(responder_id),
    responder_id != caller_id,
    call_type_r %in% c("w", "m")
  ) %>%
  count(responder_id, call_type_r, name = "n") %>%
  group_by(responder_id) %>%
  mutate(
    total_responses = sum(n),
    proportion = n / total_responses
  ) %>%
  ungroup()

call_types <- dyadic_events %>%
  filter(
    !is.na(caller_id),
    call_type_c %in% c("w", "m")
  ) %>%
  count(caller_id, call_type_c, name = "n") %>%
  group_by(caller_id) %>%
  mutate(
    total_responses = sum(n),
    proportion = n / total_responses
  ) %>%
  ungroup()


# Original (initiating) calls
call_types_plot <- call_types %>%
  transmute(
    animal_id = caller_id,
    call_type = call_type_c,
    n,
    total_calls = total_responses,
    proportion,
    call_category = "Initial calls"
  )

# Response calls
response_types_plot <- response_types %>%
  transmute(
    animal_id = responder_id,
    call_type = call_type_r,
    n,
    total_calls = total_responses,
    proportion,
    call_category = "Response calls"
  )

# Combine the data
combined_types <- bind_rows(
  call_types_plot,
  response_types_plot
)


ggplot(
  combined_types,
  aes(
    x = animal_id,
    y = proportion,
    fill = call_type,
    group = call_type
  )
) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.7
  ) +
  facet_wrap(~ call_category) +
  scale_fill_manual(
    values = c("m" = "turquoise4", "w" = "goldenrod2"),
    labels = c("m" = "Moan", "w" = "Wail")
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = "Initial and response call types by animal",
    x = "Animal ID",
    y = "Proportion of calls",
    fill = "Call type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, size = 7
    ),
    panel.grid.major.x = element_blank()
  )



# Response types per sex------

# Initial calls by sex
dyadic_events$sex_c <- tolower(dyadic_events$sex_c)
dyadic_events$sex_r <- tolower(dyadic_events$sex_r)

call_types_sex <- dyadic_events %>%
  filter(
    !is.na(caller_id),
    !is.na(sex_c),
    call_type_c %in% c("w", "m")
  ) %>%
  count(sex = sex_c, call_type = call_type_c, name = "n") %>%
  group_by(sex) %>%
  mutate(
    total_calls = sum(n),
    proportion = n / total_calls
  ) %>%
  ungroup() %>%
  mutate(call_category = "Initial calls")

# Response calls by sex
response_types_sex <- dyadic_events %>%
  filter(
    !is.na(responder_id),
    responder_id != caller_id,
    !is.na(sex_r),
    call_type_r %in% c("w", "m")
  ) %>%
  count(sex = sex_r, call_type = call_type_r, name = "n") %>%
  group_by(sex) %>%
  mutate(
    total_calls = sum(n),
    proportion = n / total_calls
  ) %>%
  ungroup() %>%
  mutate(call_category = "Response calls")

# Combine the data
combined_sex <- bind_rows(
  call_types_sex,
  response_types_sex
)

# Plot
ggplot(
  combined_sex,
  aes(
    x = sex,
    y = proportion,
    fill = call_type
  )
) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.7
  ) +
  facet_wrap(~ call_category) +
  scale_fill_manual(
    values = c("m" = "turquoise4", "w" = "goldenrod2"),
    labels = c("m" = "Moan", "w" = "Wail")
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = "Initial and response call types by sex",
    x = "Sex",
    y = "Proportion of calls",
    fill = "Call type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank()
  )



plot(x=time_init, y=sex_r, sort=sex_c, data=dyadic_events)




# Sex specific responses -----
library(dplyr)
library(ggplot2)

# Count responses by caller sex and responder sex
responses_by_sex <- dyadic_events %>%
  filter(
    !is.na(responder_id),
    responder_id != caller_id,
    !is.na(sex_c),
    !is.na(sex_r)
  ) %>%
  count(sex_c, sex_r, name = "n")

# Plot
ggplot(
  responses_by_sex,
  aes(
    x = sex_c,
    y = n,
    fill = sex_r
  )
) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "female" = "turquoise4",
      "male" = "goldenrod2"
    ),
    name = "Responder sex"
  ) +
  scale_x_discrete(
    labels = c(
      "female" = "Female caller",
      "male" = "Male caller"
    )
  ) +
  labs(
    title = "Responses by caller and responder sex",
    x = "Sex of calling animal",
    y = "Number of responses"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank()
  )


# Likelihood of getting a response from a specific sex -----
# Summarise each initiated call
call_responses <- dyadic_events %>%
  group_by(event) %>%
  summarise(
    sex_c = first(na.omit(sex_c)),
    n_female = n_distinct(responder_id[sex_r == "female" & !is.na(responder_id)]), #Make sure we exclude self reactions
    n_male = n_distinct(responder_id[sex_r == "male" & !is.na(responder_id)]),
    .groups = "drop"
  ) #%>%
  #filter(!is.na(sex_c))

# Calculate response rates by caller sex
response_rates <- call_responses %>%
  group_by(sex_c) %>%
  summarise(
    total_calls = n(),
    female_response_rate = mean(n_female > 0),
    male_response_rate = mean(n_male > 0),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(female_response_rate, male_response_rate),
    names_to = "responder_sex",
    values_to = "response_rate"
  ) %>%
  mutate(
    responder_sex = recode(
      responder_sex,
      female_response_rate = "female",
      male_response_rate = "male"
    )
  )

# Plot
ggplot(
  response_rates,
  aes(x = sex_c, y = response_rate, fill = responder_sex)
) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 1)
  ) +
  scale_fill_manual(
    values = c("female" = "turquoise4", "male" = "goldenrod2"),
    name = "Responder sex"
  ) +
  labs(
    title = "Probability of receiving a response by sex",
    x = "Caller sex",
    y = "Proportion of calls receiving a response"
  ) +
  theme_minimal()
