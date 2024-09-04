# Load required libraries
library(ggplot2)
library(dplyr)  

# Set paths 
path <- file.path("/Users/sara/Desktop/Projects/ExternalCollaborations/1_Fingerprinting_AD/repo/outputs/") # change this accordingly to where your outputs path is located 
setwd(path)


data <- read.csv("chi_square_within_between.csv", header = T,na.strings=c("","NA"))
#replace chi-square highest then 300 with 300 (for better visualization)
data$CU_minus[data$CU_minus > 300] <- 300
data$MCI_plus[data$MCI_plus > 300] <- 300

names(data)[names(data) == "CU_minus"] <- "CU Aβ-"
names(data)[names(data) == "MCI_plus"]<- "MCI Aβ+"
names(data)[names(data) == "Dementia_plus"] <- "AD Dementia"

# Transform the data into long format
data_long <- reshape2::melt(data, id.vars = 'Network', variable.name = 'Group', value.name = 'Value')

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data_long$Network), ncol(data)) )
colnames(to_add) <- colnames(data_long)
to_add$Network <- rep(levels(data_long$Network), each=empty_bar)
data_long <- rbind(data_long, to_add)
data_long <- data_long %>% arrange(factor(Network, levels = c('VIS', 'SMT', 'DA', 'SA', 'L', 'FPN', 'DMN', 'SBC' )))
data_long$id <- seq(1, nrow(data_long ))
data_long$id2 <- rep(c(1,2,3), nrow(data_long)/3)

# Get the name and the y position of each label
label_data <- data_long
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data_long %>% 
  group_by(Network) %>% 
  summarize(start=min(id), end=max(id) ) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
base_data <- base_data%>% arrange(factor(Network, levels = c('VIS', 'SMT', 'DA', 'SA', 'L', 'FPN', 'DMN', 'SBC' )))

# prepare a data frame for grid (scales)
grid_data <- base_data

grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- data_long %>% arrange(factor(Network, levels = c('VIS', 'SMT', 'DA', 'SA', 'L', 'FPN', 'DMN', 'SBC' ))) %>%
  ggplot(aes(x=as.factor(id), y=Value, fill=Network)) +

#p<- ggplot(data_long, aes(x=as.factor(id), y=Value, fill=Network)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=Value, fill=Network), stat="identity", alpha=0.5)+

  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data_long$id),3), y = c(100, 200, 300), label = c("100", "200", "300") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=Value, fill=Network), stat="identity", alpha=0.5) +
  ylim(-100,340) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")) +
  coord_polar() +
  geom_text(data = label_data, aes(x=id, y=Value+10, label=as.character(Group), hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3.7, angle= label_data$angle, inherit.aes = FALSE ) +

# Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) + 
  geom_text(data=base_data, aes(x = title, y = -35, label=Network), colour = "black", alpha=1, size=3.5,  inherit.aes = FALSE)

p
library(svglite)

