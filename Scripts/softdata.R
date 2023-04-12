# set working directory----
setwd("your_filepath")

# load libraries----
library(tidyverse)
library(dplyr)

# import data----
hoverflies <- read.csv("Data/sex_matrix_March31.csv")

# count total abundance----
# get data in long format
hoverflies_long <- gather(hoverflies, Species, Count, -c(1:5))
# this shows the total amount of hoverflies found across all species - 1070
hov_total_abundance <- sum(hoverflies_long$Count)

# count individual species----
# this shows the total amount of hoverflies found for each species - View table for results
hoverflies_grouped <- hoverflies_long %>% group_by(Species) %>%
  dplyr::summarise(Abundance = sum(Count))

hoverflies_farm <- hoverflies_long %>% group_by(Farm) %>%
  dplyr::summarise(Abundance = sum(Count))

# play with the group_by variable and you can get counts for all sorts of variables/classifications
# i.e. (Month, Farm, Transect, Pantraps) to see hoverflies in each trap, (Colour)
# to see hoverflies captured for each colour pantrap, etc.

# save new dataframe with species count----
# remove . as separator as it is annoying and looks better 
# (this is just for word table, . is good in R as avoids using spaces)
hoverflies_grouped$Species <- str_replace_all(hoverflies_grouped$Species, "[.]", " ")

# save
write.csv(hoverflies_grouped, "Data/species_counts.csv")





