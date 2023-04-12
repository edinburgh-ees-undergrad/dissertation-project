# set working directory----
setwd("your_filepath")

# load libraries (install packages if you haven't yet)----
library(tidyverse)
library(plyr)
library(permute) # needed to run vegan package
library(lattice) # needed to run vegan package
library(vegan)

# import data----
# in this case we will import it in matrix format and then convert it to long format
hoverflies_sex_matrix <- read.csv("Data/sex_matrix_March31.csv")

# data manipulation ----
# manipulating data so we can have a matrix in wide format that doesn't separate them by sex

# get data in long format
hoverflies_long <- gather(hoverflies_sex_matrix, Species, Count, c(6:58))

# remove rows that are not zero
hoverflies_long_no_zeros<- filter(hoverflies_long, Count > 0)

# save dataframe
write.csv(hoverflies_long_no_zeros,"your_filepath/Data/hoverfly_abundance_sex.csv")

# import dataset edited on Excel (combined male and female counts for each species) for biodiversity calculations
hoverfly_abundance <- read.csv("Data/hoverfly_abundance.csv")

# make this dataset into wide format
hoverflies_abundance_wide <- spread(hoverfly_abundance, Species, Count)

# replace NAs with 0s
hoverflies_abundance_wide[is.na(hoverflies_abundance_wide)] <- 0

# save dataset for easy access later- will manually add pantraps with zero hoverfly abundance in Excel
# note you do not have to do this, you can find the completed matrix in the data folder
write.csv(hoverflies_abundance_wide,"your_filepath/Data/species_matrix_incomplete.csv")

# hoverfly diversity calculations----
# import species matrix NOT divided by sex
hoverflies_matrix <- read.csv("Data/species_matrix.csv")

# hoverfly species richness----
# calculate species richness-- using wide matrix (NOT divided by sex)
hoverfly_diversity <- ddply(hoverflies_matrix, .(Month, Farm, Transect, Pantraps),function(x) {   # group data by individual pantraps
data.frame(hoverfly_richness=sum(x[-c(1:5)]>0))   # calculate richness for each pantrap (sum the amount fo species with an abundance greater than zero)
  })

# hoverfly Shannon's diversity----
# using the diversity function from the vegan package
shannon <- ddply(hoverflies_matrix,.(Month, Farm, Transect, Pantraps),function(x) {
data.frame(hoverfly_shannon=diversity(x[-c(1:5)], index="shannon"))
  })
# if shannon = 0 then there is only one species present or no species present at all

# hoverfly abundance per pantrap----
hoverfly_abundance_pantrap <- ddply(hoverflies_matrix, .(Month, Farm, Transect, Pantraps),function(x) {   
  data.frame(hoverfly_abundance=sum(x[-c(1:5)]))  
})

# Global hoverfly diversity dataset----
# create one common dataset with hoverfly richness, shannon's diversity, and hoverfly abundance
hoverfly_diversity$hoverfly_shannon <- shannon$hoverfly_shannon # extract column from shannon dataframe and paste it into the richness dataframe, making sure it matches!!
hoverfly_diversity$hoverfly_abundance <- hoverfly_abundance_pantrap$hoverfly_abundance

# export dataset for sorting in Excel----
write.csv(hoverfly_diversity,"your_filepath/Data/hoverfly_diversity_incomplete.csv")
# data is exported as "incomplete" because data needs sorting in Excel and functional diversity/community weighted mean
# values calculated in Excel need to be added to the "hoverfly_diversity" data set for statistical analysis

# functional diversity----
# functional diversity and community weighted mean calculations will take place in Microsoft 
# Excel using the macro designed by Leps et al. (2006)
# here we are just checking the distribution of the trait values for species to see whether 
# log transformation of the data is necessary (too skewed/outliers)
# import functional trait matrix
functional_matrix <- read.csv("Data/functional_matrix.csv")

# changing column names
functional_matrix <- dplyr::rename(functional_matrix, body_length = Body.length..mm., 
                            wing_body_ratio = Wing.Body.Ratio, sex = Sex)

# plot histogram of each trait to assess its distribution
# body length
hist(functional_matrix$body_length)
# LOOKS NORMALLY DISTRIBUTED

# wing:body ratio
hist(functional_matrix$wing_body_ratio)
# LOOKS NORMALLY DISTRIBUTED

# sex doesn't matter as it is a binary variable

# floral diversity calculations----
# import dataset
flower_matrix <- read.csv("Data/flower_data_combined.csv")

# flower Shannon's diversity----
# for 2 m pantraps
flower_diversity <- ddply(flower_matrix,.(Month, Out_or_in, Farm, Hedgerow, pan_trap_distance),function(x) {
  data.frame(flower_shannon_2m=diversity(x[c(15:38)], index="shannon"))
})

# flower species richness----
flower_richness <- ddply(flower_matrix, .(Month, Out_or_in, Farm, Hedgerow, pan_trap_distance),function(x) {   
  data.frame(flower_richness_2m=sum(x[c(15:38)]>0))  
})

# compile the floral richness values with the diversity values in one common dataset 
flower_diversity$flower_richness_2m <- flower_richness$flower_richness_2m

# flower abundance----
flower_abundance <- ddply(flower_matrix, .(Month, Out_or_in, Farm, Hedgerow, pan_trap_distance),function(x) {   
  data.frame(flower_abundance_2m=sum(x[c(15:38)]))  
})

# compile the floral abundance values with the diversity values in one common dataset 
flower_diversity$flower_abundance_2m <- flower_abundance$flower_abundance_2m


# export dataset to average in and out diversity measurements for each pan trap 
# to create one global diversity value for each pan trap
# This will be done in Microsoft Excel and compiled with hoverfly diversity data for stats
write.csv(flower_diversity,"your_filepath/Data/flower_diversity_incomplete.csv")


# prepare dataset for analysis of hedgerow quality----
# import data
hoverflies_transect <- read.csv("Data/species_matrix_transect.csv")

# calculate hoverfly diversity per transect in the same way we calculated it for pan traps

# hoverfly species richness
# calculate species richness-- using wide matrix (NOT divided by sex)
hoverfly_diversity_transect <- ddply(hoverflies_transect, .(Month, Farm, Transect),function(x) {   # group data by individual pantraps
  data.frame(hoverfly_richness=sum(x[-c(1:3)]>0))   # calculate richness for each pantrap (sum the amount fo species with an abundance greater than zero)
})

# hoverfly Shannon's diversity
# using the diversity function from the vegan package
shannon_transect <- ddply(hoverflies_transect,.(Month, Farm, Transect),function(x) {
  data.frame(hoverfly_shannon=diversity(x[-c(1:3)], index="shannon"))
})
# if shannon = 0 then there is only one species present or no species present at all

# hoverfly abundance per pantrap
hoverfly_abundance_transect <- ddply(hoverflies_transect, .(Month, Farm, Transect),function(x) {   
  data.frame(hoverfly_abundance=sum(x[-c(1:3)]))  
})

# Global hoverfly diversity dataset
# create one common dataset with hoverfly richness, shannon's diversity, and hoverfly abundance
hoverfly_diversity_transect$hoverfly_shannon <- shannon_transect$hoverfly_shannon # extract column from shannon dataframe and paste it into the richness dataframe, making sure it matches!!
hoverfly_diversity_transect$hoverfly_abundance <- hoverfly_abundance_transect$hoverfly_abundance

# export dataset for sorting in Excel----
write.csv(hoverfly_diversity_transect,"your_filepath/Data/hoverfly_hedgerow_incomplete.csv")
# data is exported as "incomplete" because data needs sorting in Excel and functional diversity/community weighted mean
# values calculated in Excel need to be added to the "hoverfly_hedgerow" data set for statistical analysis


      