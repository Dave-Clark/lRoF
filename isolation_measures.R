# area of sampling site
library(data.table)
library(sf)
library(ggplot2)

dat <- fread("envdata_lrof_2020.csv")

# function to calculate area of minimal spanning convex hull
hullArea <- function(x, xCol, yCol){
  coords <- st_as_sf(x, coords = c(xCol, yCol), crs = 4326)
  area <- st_area(st_convex_hull(st_union(coords)))/1000000
  return(as.numeric(area))
}

maxDist <- function(x, xCol, yCol){
  coords <- st_as_sf(x, coords = c(xCol, yCol), crs = 4326)
  return(as.numeric(max(st_distance(coords)))/1000)
}

isolationDat <- dat[, .(
  max_dist = maxDist(.SD, xCol = "Longitude", yCol = "Latitude"),
  area = hullArea(.SD, xCol = "Longitude", yCol = "Latitude")),
  by = Location]

isolationPlot <- ggplot(isolationDat,
    aes(x = max_dist, y = area, col = Location)) +
  geom_point(size = 4) +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Max. dist between streams (km)",
    y = expression(Area~of~site~(km^2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14))

ggsave("site_isolation.pdf", isolationPlot, height = 4, width = 6,
  device = "pdf")

### Get areas of Island sites first e.g. Iceland, Svalbard, Greenland
islands <- c("Disko Island", "Hengill_13", "Svalbard")

# calculate centroid points for each island site
islandCentroids <- dat[Location %in% islands,
  .(longitude = mean(Longitude), latitude = mean(Latitude)),
  by = Location]

# correct longitude of svalbard data
islandCentroids[Location == "Svalbard", longitude := -1 * longitude]
islandCentroids <- st_as_sf(islandCentroids,
  coords = c("longitude", "latitude"), crs = 4326)

# get polygons for each country
iceland <- ne_countries(scale = "large", returnclass = "sf",
  country = "Iceland")
greenland <- ne_countries(scale = "large", returnclass = "sf",
  country = "Greenland")
svalbard <- ne_countries(scale = "large", returnclass = "sf",
  country = "Norway")

countries <- list(greenland, iceland, svalbard)

countries <- lapply(countries, st_cast, "POLYGON")

getIslandArea <- function(x, poly){
  as.numeric(st_area(poly[as.numeric(st_intersects(x, poly)), ])/1000000)
}

# disko, iceland, svalbard (spitsbergen)
isolationDat[Location %in% islands,
  region_area := unlist(lapply(1:3, function(x)
    getIslandArea(islandCentroids[x, ], countries[[x]]))
)]

# get area of continental sites
bioregions <- st_read("bio_regions/CMEC regions & realms/Regions.shp")

# get centroids of continental sites
# adjust longitudes by 12 deg to match bioregions shapefile
conts <- c("Alaska", "Kamchatka")

# calculate centroid points for each island site
contCentroids <- dat[Location %in% conts,
  .(longitude = mean(Longitude) + 12, latitude = mean(Latitude)),
  by = Location]

contCentroids <- st_as_sf(contCentroids,
  coords = c("longitude", "latitude"), crs = 4326)

bioPolys <- st_cast(bioregions, "POLYGON")
bioPolys <- st_transform(bioPolys, crs = st_crs(contCentroids))

isolationDat[Location %in% conts,
  region_area := getIslandArea(contCentroids, bioPolys)]

fwrite(isolationDat, "lRoF_site_isolation.csv")
