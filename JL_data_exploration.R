
# Script for data exploration of tick data, eco forecast challenge



# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)


# download data -----------------------------------------------------------
ticks <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ticks/master/Ticks_NEON_Field_Site_Metadata_20210928.csv")



# make site data map ------------------------------------------------------

library(rnaturalearth)
usa <- map_data('state')

usmap <- ggplot() +
  geom_polygon(data = usa,
               aes(x=long, y=lat, group = group),
               fill = "grey95",
               color = "grey30",
               size=.15) 

usmap + 
  geom_point(data = site_data,
             aes(x = field_longitude,
                 y = field_latitude),
             color = "darkmagenta",
             size = 2) +
  ggrepel::geom_text_repel(data = site_data,
             aes(x = field_longitude,
                 y = field_latitude,
                 label = stringr::str_wrap(field_site_name,20)),
             size = 3,
             lineheight=.8) +
  
  coord_map() +
  ggthemes::theme_map()




# color map by range ------------------------------------------------------
tickrange <- ticks %>%
  group_by(siteID) %>%
  summarize(min = min(amblyomma_americanum),
            max = max(amblyomma_americanum),
            median = median(amblyomma_americanum)) %>%
  mutate(range = max - min) %>%
  rename(field_site_id = siteID)


site_data_with_range <- 
  site_data %>% 
  select(field_site_id, field_site_name, field_latitude, field_longitude ) %>%
  left_join(tickrange) 



usmap + 
  geom_point(data = site_data_with_range,
             aes(x = field_longitude,
                 y = field_latitude, 
                 color = range),
             size = 3) +
  ggrepel::geom_text_repel(data = site_data,
                           aes(x = field_longitude,
                               y = field_latitude,
                               label = stringr::str_wrap(field_site_name,20)),
                           size = 2.5,
                           lineheight=.8) +
  
  coord_map() +
  ggthemes::theme_map() +
  theme(legend.position = "right") +
  labs(color = "Tick Count Range") +
  scale_color_viridis_c("plasma")





hist(ticks$amblyomma_americanum)
hist(log1p(ticks$amblyomma_americanum))



# make some plots of tick counts ------------------------------------------
ticks %>%
  mutate(year = lubridate::year(time),
         jday = lubridate::yday(time)) %>%
  ggplot(aes(x=jday, 
             y=amblyomma_americanum, 
             color = siteID)) +
  geom_point() +
  geom_line(size=.15) +
  facet_grid(~year) +
  ggthemes::theme_few()



# see yearly anomalies ----------------------------------------------------
ticks %>%
  mutate(year = lubridate::year(time),
         jday = lubridate::yday(time))  %>%
  group_by(year,siteID) %>%
  ggplot(aes(x = siteID,
             y = amblyomma_americanum,
             color = siteID)) +
  ggdist::stat_halfeye() +
  facet_grid(~year) +
  ggthemes::theme_few()


ticks %>%
  mutate(year = lubridate::year(time),
         jday = lubridate::yday(time))  %>%
  group_by(year) %>%
  ggplot(aes(x=year, y = amblyomma_americanum, group = year)) +
  geom_boxplot()



  






