#look at season in baleen presence days

# Load necessary libraries

options(scipen = 999)

pacman::p_load(terra,spatstat,sf,dplyr,patchwork, ggplot2)

filepath = "input/2025/baleen_presence_laura_2025.csv"

#projections------
UTM20 <- "EPSG:32620" # CODE FOR UTM Zone 20

# Read baleen whale PAM season data------
baleen_DOY <- read.csv(filepath)
baleen_DOY%>%group_by(species)%>%summarise(count = n())

whale_data <-baleen_DOY%>%mutate(UTC = as.POSIXct(rec_date, format = "%Y-%m-%d"))%>%
  mutate(Date = format(as_date(rec_date), "%Y-%m-%d"))%>%
           mutate(month = month(UTC),  
                  Season = case_when(
               month %in%  c(12, 1, 2) ~ "Winter",
               month %in%  c(3, 4, 5) ~ "Spring",
               month %in%  c(6, 7, 8) ~ "Summer",
               month %in%  c(9,10,11) ~ "Fall",
               TRUE ~ NA_character_
             )
           )%>%group_by(site, Season) %>% 
  mutate(
    season_days = n_distinct(Date)
  ) %>%
  ungroup()

baleen_PA_season = whale_data%>%
  group_by(site, Season, species) %>%
  summarise(
    efort_days = first(season_days),
    detection_days = sum(presence > 0),
    proportion_det = detection_days / efort_days,
    latitude  = first(latitude),
    longitude = first(longitude),
    .groups = "drop"
  )


#check data
season = whale_data%>%group_by(Season, site, )%>%summarise(days = n())

# Date range
whale_data%>%summarise(min(Date), max(Date))
