### study site map ###
library(ggmap)
library(maptools)
library(mapdata)
library(sp)
library(ggrepel)
library(raster)

# create dataframe of approximate study sites
coords <- data.frame(lat = c(64.3, 69.55, 79.38, 65, 55.89),
	long = c(-21.1, -53.52, 13.43, -150.63, 159.67),
	site = c("Iceland", "Greenland", "Svalbard", "Alaska", "Kamchatka"))

# get international borders
arctic <- borders("world", colour = "grey", fill = "grey")

# create plot
siteMap <- ggplot() +
	arctic +
	coord_map("ortho", orientation= c(90, 0, 0)) +
	geom_point(aes(x = 0, y = 90), shape = 17, size = 2) +
	geom_label_repel(data = coords, aes(x = long, y = lat, label = site),
		size = 6, alpha = 0.5, fontface = "bold") +
	geom_point(data = coords, aes(x = long, y = lat), size = 2.5) +
	ylab("Latitude (decimal degrees)") +
	theme_bw() +
	scale_y_continuous(minor_breaks = c(90, 75, 60, 45, 30, 15),
		breaks = c(90, 75, 60, 45, 30, 15)) +
	theme(axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.y = element_text(size = 16),
		axis.title.y = element_text(size = 18),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank())

pdf("graphics/siteMap.pdf")
siteMap
dev.off()
