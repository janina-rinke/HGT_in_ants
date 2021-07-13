# R Script for plotting LGT candidates from automatic prediction
# Janina Rinke, 09.04.2021


#set working directory
setwd("/Users/Janina/Documents/MASTER/Masterarbeit/0_data")

#load necessary packages
library(ggplot2)
library(devtools)
library(tidyr)
library(export)
library(dplyr)

# load data
candidates <- read.csv("Good_bad_candidates_plot.csv")

######################## Candidate plot ########################
# transform data
total <- candidates$sum_good_candidates + candidates$sum_bad_candidates
cand_df <- data.frame(candidates$species, candidates$sum_good_candidates, candidates$sum_bad_candidates, total)
print(cand_df)

cand_df2 <- rbind(
  data.frame(cand_df$candidates.species, "count" = cand_df$candidates.sum_good_candidates, "type" = "good"),
  data.frame(cand_df$candidates.species, "count" = cand_df$candidates.sum_bad_candidates, "type" = "bad")
  )

# Make a bar plot
p <- ggplot(data=cand_df2, aes(x=cand_df.candidates.species, y=count, fill = forcats::fct_rev(type))) +
     geom_bar(stat="identity", width = 0.4, position = "stack") +
     coord_flip() +
     geom_text(aes(label = count), vjust = 0.5, hjust = 1, color="white", size= 4) +
     labs(x="species", y="number of automatically predicted LGTs", fill="Candidates") +
     scale_fill_manual(values = c("darkred", "steelblue")) +
     theme_minimal() +
     theme(
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           axis.line = element_line(colour = "black")
           )
p


p + theme(
 axis.title.x = element_text(size=12, face="bold"),
 axis.title.y = element_text(size=12, face="bold")
)

# Export to powerpoint
graph2ppt(file="Candidate_plot_classic")

######################## Average length plot ########################
#sort data according to average length
str(candidates)
avg_length_good_candidates <- as.factor(candidates$avg_length_good_candidates)

p1 <- candidates %>%
  arrange(desc(avg_length_good_candidates)) %>%   # First sort by length. This sort the dataframe but NOT the factor levels
  mutate(species=factor(species, levels=species)) %>%
  ggplot(aes(x=species, y=avg_length_good_candidates)) +
  geom_bar(stat="identity", width = 0.4, fill = "steelblue", color = "steelblue") +
  coord_flip() +
  geom_text(aes(label = avg_length_good_candidates), vjust = 0.5, hjust = 1.1, color="black", size= 4) +
  labs(x="species", y="Average length of good candidates") +
  scale_fill_manual(values = c("steelblue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
p1


p1 + theme(
  axis.title.x = element_text(size=12, face="bold"),
  axis.title.y = element_text(size=12, face="bold")
)

# Export to powerpoint
graph2ppt(file="length_of_candidates")

Last changed 2021/05/11
Added to GitHub 
