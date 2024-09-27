### script for Principal coordinate analysis (PCoA) from sequences (nexus), attributing them different colours and plotting on a map
####################################################################################################
##### we need nexus alignment and lat-lon coordinates of samples (in tab delimited text file)

#####PART 1 calculate PCoA and plot it############
##################################################

library(ape)
library(vegan)
library (ggplot2)
library(dplyr)
library(svglite) # for exporting graphics as svg

alignment <- read.nexus.data ("alciphron_coi_final.nex" )
alignment_dnabin <- as.DNAbin(alignment)

genetic_dist_matrix <- dist.dna(alignment_dnabin, model = "K80", pairwise.deletion=TRUE) #choose model
#pairwise deletion is important in case of missing data. However, too many missing data can spoil the analysis.

##pcoa in ape and easy visualization
#pcoa <- pcoa(genetic_dist_matrix, correction="none", rn=NULL)
#biplot(pcoa, Y=NULL, plot.axes = c(1,2), dir.axis1=1,
#       dir.axis2=1, rn=NULL, main=NULL)


pcoa_results <- cmdscale(genetic_dist_matrix, eig = TRUE, k = 2) #number of computed axes


# Extract PCoA coordinates
pcoa_df <- data.frame(PC1 = pcoa_results$points[, 1], 
                      PC2 = pcoa_results$points[, 2])

pcoa_df$ID <- rownames(pcoa_df)  # Create an ID column using row names

# Normalize PC1 and PC2 to range [0, 1] for color mapping
pcoa_df$PC1_scaled <- (pcoa_df$PC1 - min(pcoa_df$PC1)) / (max(pcoa_df$PC1) - min(pcoa_df$PC1))
pcoa_df$PC2_scaled <- (pcoa_df$PC2 - min(pcoa_df$PC2)) / (max(pcoa_df$PC2) - min(pcoa_df$PC2))

# Colors assignment (red, green, blue)
pcoa_df$color <- rgb(pmax(1 - pcoa_df$PC1_scaled - pcoa_df$PC2_scaled, 0),  # Red
                     pcoa_df$PC1_scaled,  # Green
                     pcoa_df$PC2_scaled)  # Blue

pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = color), size = 3) +
  scale_color_identity() +  
  theme_minimal() +
  labs(title = "PCoA of mitochondrial DNA", 
       x = "Principal Coordinate 1", 
       y = "Principal Coordinate 2") +
  theme(legend.position = "none")  # Remove legend for cleaner plot

print(pcoa_plot)

# Save
ggsave("pcoa_plot.pdf", plot = pcoa_plot, width = 8, height = 6)
ggsave("pcoa_plot.svg", plot = pcoa_plot, width = 12, height = 10)

#####PART 2 display results on a map##############
##################################################

coordinates <- read.table("alciphron_coordinates_baps.txt", header = TRUE, sep = "\t") #it has 3 columns: ID, lat, lon

# 'ID' is a common column to merge
merged_data <- left_join(coordinates, pcoa_df, by = "ID")

# Create the map
world_map <- map_data("world")

# Define limits of the target region (here Palearctic)
xlim_min <- -25  
xlim_max <- 170
ylim_min <- 25   
ylim_max <- 80

# Clip the world_map to the area of interest
clipped_map <- world_map %>%
  filter(long >= xlim_min & long <= xlim_max & lat >= ylim_min & lat <= ylim_max)

pcoa_map <-ggplot() +
  geom_polygon(data = clipped_map, aes(x = long, y = lat, group = group), fill = "lightgray") +  
  # Country boundaries
  geom_path(data = clipped_map, aes(x = long, y = lat, group = group), color = "black", linewidth = 0.3) + 
  # Overlay PCoA points
  geom_point(data = merged_data, aes(x = lon, y = lat, color = color), size = 2, alpha = 0.7) +  # Adjust size and transparency
  scale_color_identity() +  # Use the colors as they are (don't transform them)
  theme_minimal() +
  labs(title = "PCoA Lycaena alciphron", 
       x = "long", 
       y = "lat") +
  theme(legend.position = "none") +  # Remove legend for cleaner plot
  xlim(xlim_min, xlim_max) +  # Set x limits (longitude)
  ylim(ylim_min, ylim_max) +  # Set y limits (latitude)
  coord_fixed()  # Maintain aspect ratio

print(pcoa_map)

# Save
ggsave("pcoa_map_alciphron.pdf", plot = pcoa_map, width = 8, height = 6)
ggsave("pcoa_map_alciphron.svg", width = 10, height = 8, dpi = 300, device = "svg")