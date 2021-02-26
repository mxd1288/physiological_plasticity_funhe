install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

install.packages("rgeos")

library("rgeos")
library("ggplot2")
theme_set(theme_bw())
library("sf")

library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library("maps")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
head(states)
states <- cbind(states, st_coordinates(st_centroid(states)))

library("tools")
states$ID <- as.character(states$ID)
states$ID <- toTitleCase(states$ID)
head(states)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# a flat world map
ggplot(data = world) +
  geom_sf()

# adding labels
ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$NAME)), " countries)"))

# adding color
ggplot(data = world) + 
  geom_sf(color = "black", fill = "lightgreen")

# filling with data
ggplot(data = world) +
  geom_sf(aes(fill = pop_est)) +
  scale_fill_viridis_c(option = "plasma", trans = "sqrt")

# My Mp - OCNGS
NJ_map <- 
  ggplot(data = world) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_sf(data = states, fill = NA) + 
  geom_text(data = states, aes(X, Y, label = ID), size = 3, fontface="bold") +
  coord_sf(xlim = c(-76, -73), ylim = c(38.5, 41.5), expand = FALSE) +
  annotate(geom="rect", xmin=-74.5, xmax=-74, ymin=39.5, ymax=40, alpha=0.2, fill="red") +
  xlab("Longitude") + ylab("Latitude") +
  theme(plot.background = element_rect(fill = "white",colour = "black",size = 1))

NJ_map

library("ggsn")

OC_map <- ggplot(data = world) +
  geom_sf(data = states) +
  coord_sf(xlim = c(-74.5, -74), ylim = c(39.5, 40), expand = FALSE) +
  geom_point(size = 3, aes(x=-74.14, y=39.874, color="North Reference")) +
  geom_point(size = 3, aes(x=-74.180, y=39.809, color="Thermal Effluent")) +
  geom_point(size = 3, aes(x=-74.185, y=39.784, color="South Reference")) +
  scale_color_manual(name="Population",values=c("blue", "purple", "red"), labels=c("North Reference", "South Reference", "Thermal Effluent")) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  theme(legend.position=c(0.2, 0.85)) +
  xlab("Longitude") + ylab("Latitude") +
  theme(plot.background = element_rect(fill = "white",colour = "black",size = 1))
  
OC_map

library(cowplot)
library(gridExtra)

gg_inset_map1 = ggdraw() +
  draw_plot(OC_map) +
  draw_plot(NJ_map, x = 0.45, y = 0.05, width = 0.35, height = 0.4)

ggsave("OC_NJ_map.eps", gg_inset_map1, width=6, height=9)
