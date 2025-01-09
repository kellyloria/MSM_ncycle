##  streams map

site_colors <- c(
  "BWL" = "#054fb9",
  "BWU" = "#A6D7bF",
  "GBL" = "#DD6E42",
  "GBU" = "#EDDA42", 
  "ST848" ="black",
  "ST615" = "black",
  "USGS_BW" ="black",
  "USGS_GB" ="black")

library(tidyverse)
library(ggrepel)
library(ggmap)
library(ggsci)
library(scales)

library(dplyr)
library(ggplot2)
library(ggrepel)

# devtools::install_github('oswaldosantos/ggsn')
library(ggsn)

#devtools::install_github("briatte/tidykml")
library(tidykml)
library(sf)

# get ploygons for watersheds from .kml map 
kml_data <- st_read("/Users/kellyloria/Desktop/TempExtra/Tahoe_Soils_and_Hydro_Data.kml")

glimpse(kml_data)
#filter for blackwood
blackwood_creek <- kml_data %>%
  filter(Name == "BLACKWOOD CREEK")

GB_creek <- kml_data %>%
  filter(Name == "GLENBROOK CREEK")

# Blank watershed Plot BLACKWOOD CREEK
ggplot() +
  geom_sf(data = blackwood_creek)
# Blank watershed
ggplot() +
  geom_sf(data = GB_creek)

# read in site cordinates
dat <- read.csv("/Users/kellyloria/Downloads/TahoeCords_stream.csv", header = T, sep = ",")
unique(dat$site)
range(dat$long)
range(dat$lat)
# 
# Define the bounding box coordinates
bbox <- c(left = -120.2730, bottom = 38.900, right = -119.870, top = 39.278)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = bbox, zoom = 14, maptype = 'stamen_terrain')
# need api run for the first time

#library(ggmap)
# Plot the terrain map
stamen_map <- ggmap(map)
stamen_map

# Whole lake map
lake_map <- ggmap(map) +
  geom_point(data = dat, aes(x = long, y = lat, color = site, shape = station_type), size = 2.5) +
  # geom_label_repel(data = dat, aes(x = long, y = lat, label = site),
  #                  color = "black", fontface = "bold", size = 3) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() +
  scale_color_manual(values = site_colors,
                    guide = guide_legend(override.aes = list(label = ""))) +
  scale_shape_manual(values = c(19, 0, 2)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )
# + scalebar(dat, dist = 4, dist_unit = "km", st.size = 3,
#            transform = TRUE, model = "WGS84",
#            box.fill = c("black", "white"),
#            location = "topright") +
#   north(dat, location = "bottomright", symbol = 15) +
#   coord_cartesian()

              
 

# Extracting coordinates from blackwood_creek
coords <- st_coordinates(blackwood_creek)

# Creating a data frame from the extracted coordinates
blackwood_creek_df <- as.data.frame(coords)[,c(1,2)]
names(blackwood_creek_df) <- c("long", "lat")  # Rename columns
head(blackwood_creek_df)

# Extracting coordinates from blackwood_creek
coords <- st_coordinates(GB_creek)

# Creating a data frame from the extracted coordinates
GB_df <- as.data.frame(coords)[,c(1,2)]
names(GB_df) <- c("long", "lat")  # Rename columns
head(GB_df)

# Plotting the ggmap object with the Watershed polygons
lake_map_WS <- lake_map +
  geom_polygon(data = blackwood_creek_df, aes(x = long, y = lat), color = "#054fb9", fill = NA, linewidth = 0.25, alpha=0.5) +
  geom_polygon(data = GB_df, aes(x = long, y = lat), fill = NA, color = "#DD6E42", linewidth = 0.25, alpha=0.5)

# ggsave(plot = lake_map_WS, filename = paste("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/24_draft_lake_map_z13.png",sep=""),width=10,height=8,dpi=300)

scalebar(data = dat_e, dist = 1, location = "bottomright", dist_unit = "km", transform = TRUE)



#######################################
### Zoom in on the west shore sites ###
#######################################

# Define the bounding box coordinates
wbbox <- c(left = -120.179, bottom = 39.102, right = -120.145, top = 39.1441)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = wbbox,  zoom = 18, maptype = 'stamen_terrain') 

# Plot the terrain map for BW without GPS points 
stamen_map <- ggmap(map)
stamen_map

dat_w <-dat%>%
  filter(site=="BWL"|site=="BWNS1"| site=="BWNS2"|site=="BWNS3"|site=="BWO"|
           site=="SSNS1"| site=="SSNS2"| site=="SSNS3")

library(ggmap)
library(ggplot2)
library(ggsn)
library(dplyr)
library(ggrepel)



Westshore_zoom_map <- ggmap(map) + 
  geom_point(data = dat_w, aes(x = long, y = lat, color = site), size = 2) +
   geom_label_repel(data = dat, aes(x = long, y = lat, label = site), 
                    color = "black", fontface = "bold", size = 3) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() + 
  scale_color_manual(values = c("#4697bd", "#3283a8", "#3283a8", "#3283a8","#4697bd",
                                "#136F63", "#136F63", "#136F63"),
                     guide = guide_legend(override.aes = list(label = ""))) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank() 
  )
Westshore_zoom_map

# 
# 
# ## ### 
# ## old map 
# population_map <- ggmap(map) + 
#   geom_point(data = dat_w, aes(x = long, y = lat, color = site), size = 2) +
#   geom_label_repel(data = dat, aes(x = long, y = lat, label = site), 
#                    color = "black", fontface = "bold", size = 3) +
#   xlab("Longitude") + ylab("Latitude") + theme_bw() + 
#   scale_color_manual(values = c("#4697bd", "#3283a8", "#3283a8", "#3283a8","#4697bd",
#                                 "#136F63", "#136F63", "#136F63"),
#                      guide = guide_legend(override.aes = list(label = ""))) +
#   theme(
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 16, colour = "black", face = "bold"),
#     panel.border = element_rect(size = 1.5, colour = "black"),
#     legend.text = element_text(size = 16),
#     legend.title = element_text(size = 16, face = "bold"),
#     panel.grid = element_blank() 
#   )
# 
# ## ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/figures/23_BWzoom.png",sep=""),width=8,height=6,dpi=300)
# # ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/23_west_big_z13.png",sep=""),width=8,height=6,dpi=300)
# 
# 




############# worflow? 
############

library(ggmap)
library(ggplot2)
library(ggsn)
library(dplyr)
library(ggrepel)

# Make sure dat_w has no missing values in long and lat
dat_w <- dat %>%
  filter(site %in% c("BWL", "BWNS1", "BWNS2", "BWNS3", "BWO",
                     "SSNS1", "SSNS2", "SSNS3")) %>%
  drop_na(long, lat)

# Plot the terrain map with points, scale bar, and north arrow
Westshore_zoom_map <- ggmap(map) +
  geom_point(data = dat_w, aes(x = long, y = lat, color = site), size = 2) +
  geom_label_repel(data = dat_w, aes(x = long, y = lat, label = site),
                   color = "black", fontface = "bold", size = 3) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  scale_color_manual(values = c("#4697bd", "#3283a8", "#3283a8", "#3283a8", "#4697bd",
                                "#136F63", "#136F63", "#136F63"),
                     guide = guide_legend(override.aes = list(label = ""))) +
  # Add the scale bar
  scalebar(data = dat_w, dist = 0.1, dist_unit = "km",
           transform = TRUE, model = "WGS84",
           location = "bottomright", st.size = 3, height = 0.02) +
  # Add the north arrow
  north(data = dat_w, location = "topright", scale = 0.1, symbol = 16) +
  coord_cartesian() +  # Wrap in coord_cartesian to enable annotation
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )

# Display the map
Westshore_zoom_map


############
############


# Define the bounding box coordinates
ebbox <- c(left = -119.952, bottom = 39.08, right = -119.93, top = 39.105)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = ebbox,  zoom = 18, maptype = 'stamen_terrain') 


# Plot the terrain map
stamen_map <- ggmap(map)
stamen_map

dat_e <-dat%>%
  filter(site=="GBL"|site=="GBNS1"| site=="GBNS2"|site=="GBNS3"| site=="GBO"|
         site=="SHNS1"|site=="SHNS2"| site=="SHNS3")

population_map <- ggmap(map) + 
  geom_point(data = dat_e, aes(x = long, y = lat, color = site), size = 2) +
  geom_label_repel(data = dat, aes(x = long, y = lat, label = site), 
                   color = "black", fontface = "bold", size = 3.5) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() + 
  scale_color_manual(values = c("#b8902c", "#a67d17", "#a67d17", "#a67d17", "#b8902c",
                                "#c76640", "#c76640", "#c76640"),
                     guide = guide_legend(override.aes = list(label = ""))) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  ) +
  # Add scale bar
  scalebar(data = dat_e, dist = 1, location = "bottomleft", dist_unit = "km", transform = TRUE)


population_map


# ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/figures/23_GBzoom.png",sep=""),width=8,height=6,dpi=300)

# ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/23_GBzoom.png",sep=""),width=8,height=6,dpi=300)
