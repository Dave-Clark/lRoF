library(data.table)
library(ecolFudge)
library(hillR)
library(vegan)
library(plotrix)
library(lme4)
library(MASS)
library(modEvA)
library(ggplot2)
library(ggparl) #  boxplot/jitter plot hybrids
library(sf)
library(viridis)
library(ggrepel)
library(ggridges)
library(betapart)
library(patchwork)
library(mvabund)
library(parallel)
library(rnaturalearth)
library(rnaturalearthhires)

arctic <- fread("arcticData.csv")
# data cleaning and preparation

# load OTU tables
otuTabs <- list.files(pattern = "OtuTab")
allOtus <- lapply(otuTabs, fread, quote = "")
names(allOtus) <- c("Archaea", "Bacteria", "Eukarya")

# filter non specific OTUs
allOtus[[1]] <- allOtus[[1]][domain == "Archaea"]
allOtus[[2]] <- allOtus[[2]][domain == "Bacteria"]
# for euks, remove metazoa (animals) and Embryophyceae (land plants)
allOtus[[3]] <- allOtus[[3]][!grepl(
	"Metazoa", V5, ignore.case = T) & !grepl(
			"Embryophyceae", V6, ignore.case = T)
	]

taxSamples <- c("arch", "bac", "euk")

# remove taxon cols
filtOtuTabs <- lapply(1:3, function(x) allOtus[[x]][, .SD,
	.SDcols = c(grep(taxSamples[x], colnames(allOtus[[x]]), value = T), "V1")])

# transpose OTU tables
transOtus <- lapply(filtOtuTabs, transDT, transCol = "V1", rowID = "sample")

# create vector for each df with OTU containing cols
otuCols <- lapply(transOtus, function(x) names(x)[grepl("OTU", names(x))])

# name each df
names(transOtus) <- names(allOtus)

# reformat sample labels to match metadata
transOtus <- lapply(transOtus, function(x)
	x[, sample := sapply(strsplit(sample, "_"), "[[", 1)])

# merge OTU tables with main dataframe
sampNames <- c("archSampName", "bacSampName", "eukSampName")
transOtus <- lapply(1:3, function(x)
	merge(transOtus[[x]], arctic, by.x = "sample", by.y = sampNames[x]))

# remove soil samples
transOtus <- lapply(1:3, function(x) transOtus[[x]][sampleType == "sed", ])

# get library sizes
transOtus <- lapply(1:3, function(x)
	transOtus[[x]][, libSize := rowSums(.SD), .SDcols = otuCols[[x]]])

# inspect rarefaction curves
# par(mfrow = c(1, 3))
# lapply(1:3, function(x) quickRareCurve(transOtus[[x]][, .SD,
#	.SDcols = otuCols[[x]]], label = F))

# inspect sample sizes to determine samples to be discarded
lapply(transOtus, function(x) sort(x$libSize))
smallSamples <- c(600, 20000, 2000)

# remove small samples
transOtus <- lapply(1:3, function(x) transOtus[[x]][libSize > smallSamples[x]])

# randomly rarefy OTU tables to even depth
set.seed(101)
transOtus <- lapply(1:3, function(x) transOtus[[x]][,
	(otuCols[[x]]) := as.data.table(rrarefy(.SD, sample = min(libSize))),
	.SDcols = otuCols[[x]]])

# function to remove OTUs with no occs across entire dataset
remove_empty_otus <- function(table, otuId){
	otuAbunds <- table[, colSums(.SD), .SDcols = otuId]
	emptyOtus <- names(otuAbunds[otuAbunds == 0])
	table <- table[, (emptyOtus) := NULL]
	return(table)
}

# remove empty OTUs
transOtus <- lapply(1:3, function(x)
	remove_empty_otus(table = transOtus[[x]], otuId = otuCols[[x]]))

# SI fig
# get library sizes of remaining samples and plot alongside small sample cutoff
libSizes <- rbindlist(lapply(transOtus, function(x) x[, .(taxon, libSize)]))

cutoffs <- data.table(taxon = as.factor(c("Archaea", "Bacteria", "Eukarya")),
	cutoff = c(600, 20000, 2000))

libPlot <- ggplot(libSizes, aes(x = taxon, y = libSize)) +
	geom_boxjitter(width = 0.7, jitter.alpha = 0.6, outlier.shape = NA) +
	scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000)) +
	geom_segment(data = cutoffs,
		aes(x = as.numeric(taxon) - 0.3,
		 		xend = as.numeric(taxon) + 0.3,
				y = cutoff,
				yend = cutoff),
		linetype = 2, col = "red", size = 1) +
	labs(x = "Taxon", y = "Sequences per sample") +
	theme_bw() +
	theme(axis.text = element_text(size = 16),
		axis.title = element_text(size = 18),
		panel.grid = element_blank())

ggsave("../figures/libSizes.pdf", libPlot, height = 4, width = 4.5,
	device = "pdf")

# update otuCols
otuCols <- lapply(transOtus, function(x) names(x)[grepl("OTU", names(x))])
#############################################
##### DATA READY FOR DIVERSITY ANALYSES #####
#############################################
countries <- unique(transOtus[[1]]$country)

isolationDat <- fread("lRoF_site_isolation.csv")
isolationDat <- isolationDat[Location != "Hengill_12"]
# add my site labels to match
isolationDat[, country := countries]

globalDists <- lapply(1:3, function(x)
	beta.pair(transOtus[[x]][, ifelse(.SD == 0, 0, 1), .SDcols = otuCols[[x]]]))

globalSim <- lapply(globalDists, function(x) x$beta.sim)
globalNes <- lapply(globalDists, function(x) x$beta.sne)
globalBeta  <- lapply(globalDists, function(x) x$beta.sor)

# remove values of 1 by subtracting 0.00001
globalSim <- lapply(globalSim, function(x)
	ifelse(x == 1, x - 0.0001, x))
globalNes <- lapply(globalNes, function(x)
	ifelse(x == 1, x - 0.0001, x))
globalBeta <- lapply(globalBeta, function(x)
	ifelse(x == 1, x - 0.0001, x))

tempDists <- lapply(transOtus, function(x)
	vegdist(x$temp_above_ambient, "euclid"))

simMods <- lapply(1:3, function(x)
	decay.model(globalSim[[x]],
							tempDists[[x]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

nesMods <- lapply(1:3, function(x)
	decay.model(globalNes[[x]],
							tempDists[[x]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

betaMods <- lapply(1:3, function(x)
	decay.model(globalBeta[[x]],
							tempDists[[x]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

# calculate country specific matrices
# calculates matrices across taxa for specific country
countryMatrs <- lapply(countries, function(y)
	lapply(1:3, function(x)
		beta.pair(
			transOtus[[x]][country == y, ifelse(.SD == 0, 0, 1),
				.SDcols = otuCols[[x]]])))

names(countryMatrs) <- countries

# organise beta diversity matrices into taxa and metric type, and remove 1s
archSim <- lapply(lapply(lapply(countryMatrs, "[[", 1), "[[", 1), function(x)
	ifelse(x == 1, x - 0.0001, x))
archNes <- lapply(lapply(lapply(countryMatrs, "[[", 1), "[[", 2), function(x)
	ifelse(x == 1, x - 0.0001, x))
archBeta <- lapply(lapply(lapply(countryMatrs, "[[", 1), "[[", 3), function(x)
	ifelse(x == 1, x - 0.0001, x))
bacSim <- lapply(lapply(lapply(countryMatrs, "[[", 2), "[[", 1), function(x)
	ifelse(x == 1, x - 0.0001, x))
bacNes <- lapply(lapply(lapply(countryMatrs, "[[", 2), "[[", 2), function(x)
	ifelse(x == 1, x - 0.0001, x))
bacBeta <- lapply(lapply(lapply(countryMatrs, "[[", 2), "[[", 3), function(x)
	ifelse(x == 1, x - 0.0001, x))
eukSim <- lapply(lapply(lapply(countryMatrs, "[[", 3), "[[", 1), function(x)
	ifelse(x == 1, x - 0.0001, x))
eukNes <- lapply(lapply(lapply(countryMatrs, "[[", 3), "[[", 2), function(x)
	ifelse(x == 1, x - 0.0001, x))
eukBeta <-lapply( lapply(lapply(countryMatrs, "[[", 3), "[[", 3), function(x)
	ifelse(x == 1, x - 0.0001, x))

# calculate euclidean temperature matrices
tempMatrs <- lapply(countries, function(y)
	lapply(1:3, function(x)
		vegdist(
			transOtus[[x]][country == y, temp_above_ambient], "euclid")))

# rearrange temperature matrices into taxa
names(tempMatrs) <- countries
archTempMats <- lapply(tempMatrs, "[[", 1)
bacTempMats <- lapply(tempMatrs, "[[", 2)
eukTempMats <- lapply(tempMatrs, "[[", 3)

# iterate through model combinations
archSimMods <- lapply(1:5, function(n)
	decay.model(x = archTempMats[[n]],
							y = archSim[[n]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

archNesMods <- lapply(1:5, function(n)
	decay.model(x = archTempMats[[n]],
							y = archNes[[n]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

archBetaMods <- lapply(1:5, function(n)
	decay.model(x = archTempMats[[n]],
							y = archBeta[[n]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

archBetaResults <- data.table(taxon = "Archaea",
															country = rep(countries, times = 3),
															modelType = rep(
																c("Turnover", "Nestedness", "Dissimilarity"),
																each = 5),
															slope = c(
																sapply(archSimMods, function(x) x$b.slope),
																sapply(archNesMods, function(x) x$b.slope),
																sapply(archBetaMods, function(x) x$b.slope)),
															rSqd = c(
																sapply(archSimMods, function(x) x$pse),
																sapply(archNesMods, function(x) x$pse),
																sapply(archBetaMods, function(x) x$pse)),
															pVal = c(sapply(archSimMods, function(x) x$p.val),
																sapply(archNesMods, function(x) x$p.val),
																sapply(archBetaMods, function(x) x$p.val))
															)

bacSimMods <- lapply(1:5, function(n)
	decay.model(x = bacTempMats[[n]],
							y = bacSim[[n]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

bacNesMods <- lapply(1:5, function(n)
	decay.model(x = bacTempMats[[n]],
							y = bacNes[[n]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

bacBetaMods <- lapply(1:5, function(n)
	decay.model(x = bacTempMats[[n]],
							y = bacBeta[[n]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

bacBetaResults <- data.table(taxon = "Bacteria",
															country = rep(countries, times = 3),
															modelType = rep(
																c("Turnover", "Nestedness", "Dissimilarity"),
																each = 5),
															slope = c(
																sapply(bacSimMods, function(x) x$b.slope),
																sapply(bacNesMods, function(x) x$b.slope),
																sapply(bacBetaMods, function(x) x$b.slope)),
															rSqd = c(
																sapply(bacSimMods, function(x) x$pse),
																sapply(bacNesMods, function(x) x$pse),
																sapply(bacBetaMods, function(x) x$pse)),
															pVal = c(sapply(bacSimMods, function(x) x$p.val),
																sapply(bacNesMods, function(x) x$p.val),
																sapply(bacBetaMods, function(x) x$p.val))
															)

eukSimMods <- lapply(1:5, function(n)
	decay.model(x = eukTempMats[[n]],
							y = eukSim[[n]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

eukNesMods <- lapply(1:5, function(n)
	decay.model(x = eukTempMats[[n]],
							y = eukNes[[n]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

eukBetaMods <- lapply(1:5, function(n)
	decay.model(x = eukTempMats[[n]],
							y = eukBeta[[n]],
							model.type = "exp",
							y.type = "diss",
							perm = 1000))

eukBetaResults <- data.table(taxon = "Eukarya",
															country = rep(countries, times = 3),
															modelType = rep(
																c("Turnover", "Nestedness", "Dissimilarity"),
																each = 5),
															slope = c(
																sapply(eukSimMods, function(x) x$b.slope),
																sapply(eukNesMods, function(x) x$b.slope),
																sapply(eukBetaMods, function(x) x$b.slope)),
															rSqd = c(
																sapply(eukSimMods, function(x) x$pse),
																sapply(eukNesMods, function(x) x$pse),
																sapply(eukBetaMods, function(x) x$pse)),
															pVal = c(sapply(eukSimMods, function(x) x$p.val),
																sapply(eukNesMods, function(x) x$p.val),
																sapply(eukBetaMods, function(x) x$p.val))
															)

globalBetaResults <- data.table(taxon = rep(c("Archaea", "Bacteria", "Eukarya"),
																	times = 3),
																country = "Global",
																modelType = rep(
																	c("Turnover", "Nestedness", "Dissimilarity"),
																	each = 3),
																slope = c(
																	sapply(simMods, function(x) x$b.slope),
																	sapply(nesMods, function(x) x$b.slope),
																	sapply(betaMods, function(x) x$b.slope)),
																rSqd = c(
																	sapply(simMods, function(x) x$pse),
																	sapply(nesMods, function(x) x$pse),
																	sapply(betaMods, function(x) x$pse)),
																pVal = c(
																	sapply(simMods, function(x) x$p.val),
																	sapply(nesMods, function(x) x$p.val),
																	sapply(betaMods, function(x) x$p.val))
																)

allBetaResults <- rbindlist(
	list(archBetaResults, bacBetaResults, eukBetaResults, globalBetaResults))

fwrite(allBetaResults, "allBetaResults.csv")

# compare slopes to global values using bootstrapping approach
simBoot <- lapply(simMods, boot.coefs.decay, R = 1000)
nesBoot <- lapply(nesMods, boot.coefs.decay, R = 1000)
betaBoot <- lapply(betaMods, boot.coefs.decay, R = 1000)

globalBoot <- data.table(
	taxon = rep(c("Archaea", "Bacteria", "Eukarya"), each = 1000, times = 3),
	country = "Global",
	modelType = rep(c("Turnover", "Nestedness", "Dissimilarity"), each = 3000),
	betaSlope = c(
		sapply(simBoot, function(x) x$boot.coefs[, 2]),
		sapply(nesBoot, function(x) x$boot.coefs[, 2]),
		sapply(betaBoot, function(x) x$boot.coefs[, 2])
	))

# bootstrap country/taxon specific models
archBootSim <- lapply(archSimMods, boot.coefs.decay, R = 1000)
archBootNes <- lapply(archNesMods, boot.coefs.decay, R = 1000)
archBootBeta <- lapply(archBetaMods, boot.coefs.decay, R = 1000)
bacBootSim <- lapply(bacSimMods, boot.coefs.decay, R = 1000)
bacBootNes <- lapply(bacNesMods, boot.coefs.decay, R = 1000)
bacBootBeta <- lapply(bacBetaMods, boot.coefs.decay, R = 1000)
eukBootSim <- lapply(eukSimMods, boot.coefs.decay, R = 1000)
eukBootNes <- lapply(eukNesMods, boot.coefs.decay, R = 1000)
eukBootBeta <- lapply(eukBetaMods, boot.coefs.decay, R = 1000)

bootResults <- data.table(
	taxon = rep(c("Archaea", "Bacteria", "Eukarya"), each = 15000),
	country = rep(countries, each = 1000, times = 3),
	modelType = rep(
		c("Turnover", "Nestedness", "Dissimilarity"), each = 5000, times = 3),
	betaSlope = c(
		unlist(lapply(list(archBootSim, archBootNes, archBootBeta), function(x)
			sapply(x, function(z) z$boot.coefs[, 2]))),
		unlist(lapply(list(bacBootSim, bacBootNes, bacBootBeta), function(x)
			sapply(x, function(z) z$boot.coefs[, 2]))),
		unlist(lapply(list(eukBootSim, eukBootNes, eukBootBeta), function(x)
			sapply(x, function(z) z$boot.coefs[, 2])))))

allBootstraps <- rbindlist(list(bootResults, globalBoot))
rm(list = c("bootResults", "globalBoot"))
gc()

# compute all possible combinations
combs <- combn(c(countries, "Global"), 2)

calcP <- function(x, y){
	overlap <- sample(x) - sample(y)
	p1 <- mean(overlap < 0)
	p2 <- mean(overlap > 0)
	return(c(p1, p2))
}

meanPvals <- function(x, y, i){
	reps <- replicate(i, calcP(x, y))
	return(rowMeans(reps))
}

archaeaCoefPvals <- lapply(
	c("Turnover", "Nestedness", "Dissimilarity"), function(sim)
		sapply(1:15, function(z)
		meanPvals(allBootstraps[taxon == "Archaea" & modelType == sim &
								country == combs[1, z], betaSlope],
							allBootstraps[taxon == "Archaea" & modelType == sim &
								country == combs[2, z], betaSlope], i = 1000))
	)

bacteriaCoefPvals <- lapply(
	c("Turnover", "Nestedness", "Dissimilarity"), function(sim)
		sapply(1:15, function(z)
		meanPvals(allBootstraps[taxon == "Bacteria" & modelType == sim &
								country == combs[1, z], betaSlope],
							allBootstraps[taxon == "Bacteria" & modelType == sim &
								country == combs[2, z], betaSlope], i = 1000))
	)

eukaryaCoefPvals <- lapply(
	c("Turnover", "Nestedness", "Dissimilarity"), function(sim)
		sapply(1:15, function(z)
		meanPvals(allBootstraps[taxon == "Eukarya" & modelType == sim &
								country == combs[1, z], betaSlope],
							allBootstraps[taxon == "Eukarya" & modelType == sim &
								country == combs[2, z], betaSlope], i = 1000))
	)

meanBetaSlopes <- allBootstraps[country != "Global",
	.(meanSlope = mean(betaSlope), slopeErr = std.error(betaSlope)),
	by = .(taxon, country, modelType)]

betaIsolation <- merge(isolationDat, meanBetaSlopes, by = "country")

isolationMods <- lapply(unique(betaIsolation$modelType), function(x)
	lmer(meanSlope ~ (1 + log10(region_area) | taxon) + log10(region_area),
		data = betaIsolation[modelType == x]))

lapply(isolationMods, MuMIn::r.squaredGLMM)

predictIsolation <- data.table(
	region_area = rep(seq(min(betaIsolation$region_area),
		max(betaIsolation$region_area), 10000), times = 3),
	taxon = rep(unique(betaIsolation$taxon),
		each = length(seq(min(betaIsolation$region_area),
			max(betaIsolation$region_area), 10000))))

predictIsolation[, ":="(
	Turnover = predict(isolationMods[[1]], newdata = predictIsolation),
	Nestedness = predict(isolationMods[[2]], newdata = predictIsolation),
	Dissimilarity = predict(isolationMods[[3]], newdata = predictIsolation))
	]

betaIsolation[, taxon := factor(taxon, levels = c("Eukarya", "Archaea", "Bacteria"))]

regionAreaPlot <- ggplot() +
	geom_point(data = betaIsolation[modelType == "Turnover"],
		aes(x = region_area, y = meanSlope, col = taxon, fill = taxon,
			shape = country),
		size = 4, alpha = 0.6) +
	scale_x_log10() +
	scale_shape_manual(values = 21:25) +
	geom_line(data = predictIsolation,
		aes(x = region_area, y  = Turnover, col = taxon), size = 1.2) +
	labs(x = expression(Region~area~(km^2)), y = "Slope", col = "Taxon",
		shape = "Location", fill = "Taxon") +
	theme_bw() +
	theme(axis.text = element_text(size = 16),
		axis.title = element_text(size = 18),
		legend.text = element_text(size = 14),
		legend.title = element_text(size = 14),
		strip.text.x = element_text(size = 14),
		panel.grid = element_blank())

ggsave("../figures/Region_slopes.pdf", regionAreaPlot, height = 4, width = 6,
	device = "pdf")

# function to calculate prediction curve for range of data
calcDecayFit <- function(x){
  xCoord <- seq(min(x$data[, 1]), max(x$data[, 1]),  0.01)
  pred <- 1 - (1 - x$a.intercept) * exp(-x$b.slope * xCoord)
  decayFit <- data.table(x = xCoord, y = pred)
  return(decayFit)
}

# calc predicted for turnover models
archSimFit <- lapply(archSimMods, calcDecayFit)
bacSimFit <- lapply(bacSimMods, calcDecayFit)
eukSimFit <- lapply(eukSimMods, calcDecayFit)

# get into one df
names(archSimFit) <- names(bacSimFit) <- names(eukSimFit) <- countries
simFits <- lapply(list(archSimFit, bacSimFit, eukSimFit), rbindlist,
	idcol = "country")
names(simFits) <- c("Archaea", "Bacteria", "Eukarya")
simFits <- rbindlist(simFits, idcol = "taxon")
setnames(simFits, old = c("x", "y"), new = c("temp", "turnover"))

# get distance data points into dfs for plotting
names(archSimMods) <- names(bacSimMods) <- names(eukSimMods) <- countries
distData <- lapply(list(archSimMods, bacSimMods, eukSimMods), function(x)
	lapply(x, function(z) z$data))

# bind distance data together and add country col
distData <- lapply(distData, rbindlist, idcol = "country")

# bind all dists together and add taxon col
names(distData) <- c("Archaea", "Bacteria", "Eukarya")
distData <- rbindlist(distData, idcol = "taxon")
setnames(distData, old = c("x", "X1...y"), new = c("temp", "turnover"))

distData[, taxon := factor(taxon, levels = c("Eukarya", "Archaea", "Bacteria"))]
simFits[, taxon := factor(taxon, levels = c("Eukarya", "Archaea", "Bacteria"))]

# plot of turnover vs distance
tvrPlot <- ggplot(distData,
		aes(x = temp, y = turnover, col = country, fill = country,
			shape = country)) +
	geom_point(alpha = 0.2) +
	scale_shape_manual(values = 21:25) +
	facet_wrap(~taxon, nrow = 1) +
	geom_line(data = simFits, aes(x = temp, y = turnover, col = country), size = 1.2, show.legend = F) +
	labs(x = expression(Temperature~difference~(degree*C)), col = "",
		y = "Community turnover", fill = "", shape = "") +
	theme_bw() +
	theme(axis.text = element_text(size = 16),
		axis.title = element_text(size = 18),
		strip.text.x = element_text(size = 14),
		legend.text = element_text(size = 14),
		legend.title = element_text(size = 14),
		panel.grid = element_blank()) +
	guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4)),
		fill = guide_legend(override.aes = list(alpha = 1, size = 4)))

# create inset zoomed map of island sizes
# get polygons that intersect with island coords
# st_cast multipolygon object, plot zoomed map of subset
# get world map
# get world map
world <- ne_countries(returnclass = "sf", scale = "large")

# read in michelle data
dat <- fread("envdata_lrof_2020.csv")

# calculate centroid coordinates for each site
coords <- dat[, .(long = mean(Longitude), lat = mean(Latitude)), by = Location]

# remove previous hengill data
coords <- coords[!Location == "Hengill_12"]

# correct svalbard longitude
coords[Location == "Svalbard", long := long * -1]

coords[, country := countries]

coords <- st_as_sf(coords, coords = c("long", "lat"), crs = 4326)

worldMap <- ggplot() +
	geom_sf(data = world, fill = "lightgrey", col = "lightgrey") +
	geom_sf(data = coords, aes(colour = country, shape = country, fill = country),
		size = 3, alpha = 0.7) +
	ggrepel::geom_label_repel(data = coords,
		aes(label = country, geometry = geometry, col = country),
		stat = "sf_coordinates", nudge_y = -1 * st_coordinates(coords)[, 2],
		size = 5) +
	scale_shape_manual(values = 21:25) +
	theme_void() +
	theme(legend.position = "none")

p1 <- worldMap + tvrPlot + theme(legend.position = "none") +
	plot_layout(ncol = 1, heights = c(1, 0.4))
betaPanel <- p1 - regionAreaPlot +
	plot_layout(widths = c(1, 0.9)) &
	plot_annotation(tag_levels = "A") &
	theme(plot.tag = element_text(size = 24))

### need to reorder taxa, and sort map out
ggsave("../figures/beta_diversity_panel.pdf", betaPanel, height = 6, width = 14,
	device = "pdf")

# reorder taxa
allBootstraps[, taxon := factor(taxon, levels = c("Bacteria", "Archaea", "Eukarya"))]

##### DESIGN SI BETA FIGURES #####
bootPlot <- ggplot(allBootstraps,
		aes(x = betaSlope, fill = country, col = country, y = taxon)) +
	geom_density_ridges(alpha = 0.5, scale = 0.9) +
	labs(x = "Bootstrapped slopes", y = "") +
	facet_wrap(~modelType, scales = "free", ncol = 2) +
	theme_bw() +
	theme(axis.text = element_text(size = 16),
		axis.title = element_text(size = 18),
		legend.text = element_text(size = 14),
		legend.title = element_blank(),
		legend.position = c(0.6, 0.25),
		strip.text.x = element_text(size = 14),
		panel.grid = element_blank())

ggsave("../figures/boostrapped_slopes.pdf", bootPlot, height = 7, width = 8,
	device = "pdf")

# pairwise plots of dissimilarity and nestedness against temp differences

########### alpha diversity analyses ############
transOtus <- lapply(1:3, function(x)
	transOtus[[x]][, ":="(
		richness = hill_taxa(.SD, q = 0),
		entropy = hill_taxa(.SD, q = 1),
		evenness = hill_taxa(.SD, q = 2)), .SDcols = otuCols[[x]]])

transOtus <- lapply(transOtus, function(z)
	merge(z, isolationDat[, !"Location"], by = "country"))

biogeoModels <- function(x){
	model <- glm.nb(
		richness ~ temp_above_ambient * country, data = x)
	return(model)}

richnessMods <- lapply(transOtus, function(z) biogeoModels(z))

richnessAICs <- lapply(richnessMods, AIC)

pR2 <- function(x){
	pseudoR2 <- 1 - (x$deviance / x$null.deviance)
	return(pseudoR2)
}

richnessR2 <- lapply(richnessMods, pR2)

# function to get temp gradients for each country to generate predictions
getPredData <- function(x){
	tempRange <- x[, .(minTemp = min(temp_above_ambient), maxTemp = max(temp_above_ambient)),
	by = country]
	tempSeq <- tempRange[,
		.(temp_above_ambient = seq(minTemp, maxTemp, 0.01)), by = country]
	return(tempSeq)
}

# get prediction data for each taxon
predData <- lapply(transOtus, getPredData)

predData <- lapply(1:3, function(x)
	predData[[x]][, predictedRichness := predict(
			richnessMods[[x]],
			newdata = predData[[x]],
			type = "response")
	])

obsRichness <- lapply(transOtus, function(x)
	x[, .(country, temp_above_ambient, richness)])

names(predData) <- names(obsRichness) <- c("Archaea", "Bacteria", "Eukarya")

allAlphaPreds <- rbindlist(predData, idcol = "taxon")
obsRichness <- rbindlist(obsRichness, idcol = "taxon")

obsRichness[, taxon := factor(taxon, levels = c("Eukarya", "Archaea", "Bacteria"))]
allAlphaPreds[, taxon := factor(taxon, levels = c("Eukarya", "Archaea", "Bacteria"))]
richnessPlot <- ggplot(obsRichness, aes(x = temp_above_ambient, y = richness)) +
	geom_point(aes(shape = country, fill = country, col = country), size = 3, alpha = 0.5) +
	scale_shape_manual(values = 21:25) +
	geom_line(data = allAlphaPreds, aes(x = temp_above_ambient, y = predictedRichness, col = country), size = 1.2) +
	facet_wrap(~taxon, scales = "free_y", ncol = 1) +
	theme_bw() +
	labs(x = expression(Temperature~above~ambient~(degree*C)),
		y = "OTU richness") +
	theme(axis.text = element_text(size = 16),
		axis.title = element_text(size = 18),
		legend.text = element_text(size = 14),
		legend.title = element_blank(),
		panel.grid = element_blank(),
		strip.text.x = element_text(size = 14))

alphaCoefs <- data.table(
	taxon = rep(c("Archaea", "Bacteria", "Eukarya"), each = 5),
	country = rep(countries, times = 3),
	intercept = unlist(lapply(richnessMods, function(x)
		c(coef(x)[1], coef(x)[3:6] + coef(x)[1]))),
	slope = unlist(lapply(richnessMods, function(x)
		c(coef(x)[2], coef(x)[2] + coef(x)[7:10]))))

alphaBiogeo <- merge(alphaCoefs, isolationDat, by = "country")

alphaBiogeo <- melt(alphaBiogeo[, !c("Location", "max_dist")],
	id.vars = c(
		"country", "taxon", "region_area", "sample_area", "mainland_dist"),
	measure.vars = c("slope", "intercept"))

# alpha diversity slope increases with increasing sampled area
# could be a species area relationship?
test <- lmer(value ~ (1 + log10(sample_area) | taxon) + log10(sample_area), data = alphaBiogeo[variable == "slope"])

sampleAreaPreds <- data.table(
	sample_area = rep(
		seq(min(isolationDat$sample_area), max(isolationDat$sample_area), 10), times = 3),
	taxon = rep(c("Archaea", "Bacteria", "Eukarya"), each = length(seq(min(isolationDat$sample_area), max(isolationDat$sample_area), 10)), times = 3)
)

sampleAreaPreds[, predictedSlope := predict(test, newdata = sampleAreaPreds, type = "response")]

richnessSlopes <- ggplot(alphaBiogeo[variable == "slope", ],
		aes(x = sample_area, y = value, col = taxon)) +
	geom_point(aes(shape = country, fill = taxon), size = 4, alpha = 0.5) +
	scale_shape_manual(values = 21:25) +
	geom_line(data = sampleAreaPreds,
		aes(x = sample_area, y = predictedSlope), size = 1.2) +
	scale_x_log10() +
	labs(x = "Site area", y = "Slope") +
	theme_bw() +
	theme(axis.text = element_text(size = 16),
		axis.title = element_text(size = 18),
		legend.text = element_text(size = 14),
		legend.title = element_blank(),
		panel.grid = element_blank())

alphaPanel <- richnessPlot + richnessSlopes + plot_layout(widths = c(0.4, 1)) & plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 20))

ggsave("../figures/alpha_panel.pdf", alphaPanel, height = 6, width = 12, device = "pdf")











#region area plot
regionArea <- ggplot(betaIsolation,
		aes(x = region_area, y = meanSlope, col = country, shape = taxon)) +
	geom_point(size = 3) +
	scale_x_log10() +
	labs(x = expression(Region~area~(km^2)), y = "Slope", col = "", shape = "") +
	facet_wrap(~ modelType) +
	theme_bw() +
	theme(axis.text = element_text(size = 16),
		axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
		axis.title = element_text(size = 18),
		legend.text = element_text(size = 14),
		legend.title = element_text(size = 14),
		strip.text.x = element_text(size = 14),
		panel.grid = element_blank())

sampleArea <- ggplot(betaIsolation,
		aes(x = sample_area, y = meanSlope, col = taxon, shape = country)) +
	geom_point(size = 3) +
	labs(x = expression(Site~area~(km^2)), y = "Slope", col = "", shape = "") +
	scale_x_log10() +
	facet_wrap(~ modelType) +
	theme_bw() +
	theme(axis.text = element_text(size = 16),
		axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
		axis.title = element_text(size = 18),
		legend.text = element_text(size = 14),
		legend.title = element_text(size = 14),
		strip.text.x = element_text(size = 14),
		panel.grid = element_blank())

mainDist <- ggplot(betaIsolation,
		aes(x = mainland_dist, y = meanSlope, col = taxon, shape = country)) +
	geom_point(size = 3) +
	labs(x = expression(Distance~to~mainland~(km)), y = "Slope", col = "",
		shape = "") +
	facet_wrap(~ modelType) +
	theme_bw() +
	theme(axis.text = element_text(size = 16),
		axis.title = element_text(size = 18),
		legend.text = element_text(size = 14),
		legend.title = element_text(size = 14),
		strip.text.x = element_text(size = 14),
		panel.grid = element_blank())

ggsave("../figures/region_area.pdf", regionArea, height = 4, width = 12,
	device = "pdf")
ggsave("../figures/sample_area.pdf", sampleArea, height = 4, width = 12,
	device = "pdf")
ggsave("../figures/mainland_dist.pdf", mainDist, height = 4, width = 12,
	device = "pdf")

################################################################################










# function to calculate prediction curve for range of data
calcDecayFit <- function(x){
  xCoord <- seq(min(x$data[, 1]), max(x$data[, 1]),  0.005)
  pred <- 1 - (1 - x$a.intercept) * exp(-x$b.slope * xCoord)
  decayFit <- data.table(x = xCoord, y = pred)
  return(decayFit)
}



allDists <- lapply(1:3, function(x)
	cbind(betaDists[[x]], tempDists[[x]][, !"country"]))

lapply(allDists, setnames, old = "V1", new = "tempDiff")

# calculate relative contribution of turnover to total beta diversity
lapply(allDists, function(x) x[, rel_tvr := beta.sim/beta.sor])

globalSim <- lapply(allDists, function(x)
	decay.model(as.dist(x$beta.sim),
							as.dist(x$tempDiff),
							model.type = "exponential",
							y.type = "dissimilarities"))
globalNes <-
globalBeta <-































# create cols for sp richness, library size, and taxon, ready for plotting
transOtus <- lapply(1:3, function(x)
	transOtus[[x]][, c("hill_1", "hill_2", "hill_3", "taxon") := list(
		hill_taxa(.SD, q = 0), hill_taxa(.SD, q = 1), hill_taxa(.SD, q = 2),
		, rep(names(transOtus)[x], times = nrow(transOtus[[x]]))),
		.SDcols =
#

# calculate rarefied richness to verify results
transOtus <- lapply(1:3, function(x)
	transOtus[[x]][, rareRich := round(rarefy(.SD, sample = min(libSize))),
	.SDcols = otuCols[[x]]])



















	# calculate contribution of nestedness and turnover
turnContr[, ":="(nesCont = (Nestedness/Similarity) * 100,
		turnCont = (Turnover/Similarity) * 100)]

	# boxplot of % contributions of each component to overall dissimilarity
contrDat <- melt(
		turnContr[, .SD, .SDcols = c("taxon", "country", "nesCont", "turnCont")],
		id.vars = c("taxon", "country"))

contrDat[, variable := factor(variable, levels = levels(variable),
		labels = c("Nestedness", "Turnover"))]

contPlot <- ggplot(contrDat, aes(x = taxon, y = value, col = country)) +
		geom_hline(yintercept = 50, linetype = 2, col = "lightgrey") +
		geom_boxplot() +
		facet_wrap(~ variable) +
		labs(x = "Taxonomic group",
			y = "Percentage contribution to\noverall dissimilarity (%)") +
		theme_bw() +
		scale_color_viridis(discrete = T, begin = 0.8, end = 0) +
		theme(axis.text = element_text(size = 16),
			axis.title = element_text(size = 18),
			panel.grid = element_blank(),
			strip.text.x = element_text(size = 14),
			legend.text = element_text(size = 14),
			legend.title = element_blank())

ggsave("../figures/nestedness_turnover.pdf", contPlot, width = 9, height = 4,
		device = "pdf")

# repeat temperature matrices for mantel tests
tm <- rep(tempMatrices, each = 3)

# apply decay models to all Sorensen matrices
sorMats <- grep("sor", names(bm))

decayModels <- mclapply(sorMats, mc.cores = mc, function(x)
	decay.model(bm[[x]], tm[[x]], model.type = "exp", y.type = "dissim",
		perm = 1000))

	# calculate global pairwise dissimilarities
	# Overall community turnover
	rareTabs <- lapply(1:3, function(x)
		transOtus[[x]][, rrarefy(.SD, sample = min(libSize)), .SDcols = otuCols[[x]]])

	# make rarefied tables binary
	rareTabs <- lapply(rareTabs, makeBinary)

	# compute overall B-diversity matrices for each taxon
	betaTabs <- lapply(rareTabs, beta.pair)

	# extract total dissim matrix only
	betaTabs <-  lapply(betaTabs, function(x) x$beta.sim)

	# run NMDS analyses
	betaNMDS <- lapply(betaTabs, metaMDS, trymax = 250, autotransform = F)

	# get stres vals
	lapply(betaNMDS, function(x) x$stress)

	# add scores back to otu tabs
	transOtus <- lapply(1:3, function(x) transOtus[[x]][, ":="(
		NMDS1 = scores(betaNMDS[[x]])[, 1],
		NMDS2 = scores(betaNMDS[[x]])[, 2])])

	# make plots
	nmdsPlots <- lapply(transOtus, function(a)
		ggplot(a, aes(x = NMDS1, y = NMDS2, col = temp, shape = country)) +
			geom_point(size = 5, alpha = 0.7) +
			labs(x = "NMDS 1", y = "NMDS 2", col = "Temperature", shape = "Site") +
			theme_bw() +
			scale_shape_manual(values = c(15, 17, 18, 19, 8)) +
			scale_color_viridis(option = "magma", begin = 0, end = 0.8,
				breaks = c(5, 15, 25), labels = c("5°C", "15°C", "25°C")) +
			theme(axis.text = element_text(size = 16),
				axis.title = element_text(size = 18),
				panel.grid = element_blank(),
				legend.text = element_text(size = 14),
				legend.title = element_text(size = 14),
				legend.box = "horizontal",
				aspect.ratio = 1)
				)

	nmdsLegend <- get_legend(nmdsPlots[[1]])

	nmdsPlots <- lapply(nmdsPlots, function(x) x + theme(legend.position = "none"))

	# arrange panel plot
	nmdsPanel <- plot_grid(nmdsPlots[[1]], nmdsPlots[[2]], nmdsPlots[[3]],
		nmdsLegend, nrow = 2, labels = c("A", "B", "C", ""), label_size = 18,
		align = "hv", axis = "l")

	ggsave("../figures/Figure_2.pdf", nmdsPanel, height = 7, width = 8,
		device = "pdf")

	#### multivariate abundance modelling
	otuOccs <- mclapply(1:3, mc.cores = mc, function(x) as.data.table(
		transOtus[[x]][, specnumber(.SD, MARGIN = 2), .SDcols = otuCols[[x]]],
		keep.rownames = T))

	# add taxonomy col and collapse into one dt
	otuOccs <- Map(cbind, otuOccs, taxon = paste0(
		c("Archaea\n", "Bacteria\n", "Eukarya\n"), " (max = ",
		sapply(transOtus, nrow, simplify = T), ")"))
	otuOccs <- do.call(rbind, otuOccs)

	# extract OTUs occurring in more than 3 samples
	commonOtus <- lapply(unique(otuOccs$taxon), function(x)
		otuOccs[taxon == x & V2 > 3, V1])

	# create OTU matrix of all otus occuring more than 3 samples
	mvMatrices <- lapply(1:3, function(x)
		as.matrix(transOtus[[x]][, .SD, .SDcols = commonOtus[[x]]]))

	# run temperature model
	tempModels <- mclapply(1:3, mc.cores = mc, function(x)
		manyglm(mvMatrices[[x]] ~ temp + I(temp^2), offset = log(libSize),
			data = transOtus[[x]]))

	bioModels <- mclapply(1:3, mc.cores = mc, function(x)
		manyglm(mvMatrices[[x]] ~ country, offset = log(libSize),
			data = transOtus[[x]]))

	### plot AICs against each other
	aicData <- data.table(
		taxon = c(rep("Archaea", times = ncol(mvMatrices[[1]])),
			rep("Bacteria", times = ncol(mvMatrices[[2]])),
			rep("Eukarya", times = ncol(mvMatrices[[3]]))),
		tempAIC = unlist(lapply(tempModels, AIC)),
		geogAIC = unlist(lapply(bioModels, AIC)))

	# add column to determine which model (if any) AIC supports
	aicData[, support := ifelse(abs(geogAIC - tempAIC) <= 2, "Equal support", ifelse(geogAIC - tempAIC < -2, "Site", "Temperature"))]

	# calculate relative difference
	aicData[, ":="(relDiff = abs(tempAIC - geogAIC)/(tempAIC + geogAIC),
		diff = (tempAIC - geogAIC)/(tempAIC + geogAIC))]

	aicPlot <- ggplot(aicData,
			aes(x = log(geogAIC), y = log(tempAIC), col = diff)) +
		geom_point(size = 3) +
		geom_abline(intercept = 0, slope = 1, linetype = 2, col = "black",
			alpha = 0.9, size = 0.7) +
		facet_wrap(~ taxon, scales = "free") +
		labs(x = expression(ln*(AIC[Site])),
			y = expression(ln*(AIC[Temperature])),
			colour = "AIC support") +
		scale_colour_gradientn(colours = viridis(3, end = 0.85),
			limits = c(-1 * max(aicData$diff), max(aicData$diff)),
			breaks = c(-1 * max(aicData$diff), 0, max(aicData$diff)),
			labels = c("Temperature", "No support", "Site")) +
		coord_cartesian(xlim = range(log(aicData$geogAIC)),
			ylim = range(log(aicData$tempAIC))) +
		theme_bw() +
		theme(axis.text = element_text(size = 16),
			axis.title = element_text(size = 18),
			strip.text = element_text(size = 14),
			legend.text = element_text(size = 14),
			legend.title = element_text(size = 14),
			panel.grid = element_blank(),
			aspect.ratio = 1)

	ggsave("../figures/Figure_3.pdf", aicPlot, height = 3.5, width = 9,
		device = "pdf")

	aicSupport <- aicData[, .N, by = c("taxon", "support")]
	aicSupport[, .((N/sum(N))*100, support), by = taxon]

	aicFig <- ggplot(aicSupport, aes(x = taxon, y = N, fill = support)) +
		geom_col(position = "dodge", width = 0.6) +
		theme_bw() +
		facet_wrap(~taxon, scales = "free") +
		scale_fill_viridis(discrete = T, end = 0.85,
			labels = levels(as.factor(aicSupport$support))) +
		labs(x = "Taxon", y = "Number of OTUs") +
		theme(axis.text = element_text(size = 16),
			axis.title = element_text(size = 18),
			strip.text.x = element_blank(),
			legend.text = element_text(size = 14),
			legend.title = element_blank(),
			panel.grid.minor = element_blank(),
			panel.grid.major = element_blank())

	ggsave("../figures/aic_proportions.pdf", aicFig, height = 3.5, width = 7,
		device = "pdf")

## Richness modelling
# richness as a quadratic function of temperature
# intercept dependent on country
tempOnly <- lapply(transOtus, function(x)
	glm.nb(hill_1 ~ offset(log(libSize)) + country + temp + I(temp^2), data = x))

# richness as a quadratic function of temperature
# temperature coefficients interacting with country
tempInt <- lapply(transOtus, function(x)
	glm.nb(hill_1 ~ offset(log(libSize)) + country * temp + I(temp^2), data = x))

# richness as quadratic function of temperature
# all temp coefficients interacting with country, but not each other
tempFull <- lapply(transOtus, function(x)
	glm.nb(hill_1 ~ offset(log(libSize)) + country * temp + I(temp^2) * country,
	data = x))

# likelihood ratio test to test which hypothesis is more likely
lapply(1:3, function(x) anova(tempOnly[[x]], tempInt[[x]], tempFull[[x]]))

richDev <- data.table(Taxon = names(allOtus),
	none = sapply(tempOnly, Dsquared),
	minimal = sapply(tempInt, Dsquared),
	full = sapply(tempFull, Dsquared))

richDev <- melt(richDev, id.vars = "Taxon", variable.name = "Interaction",
	value.name = "Explained_deviance")

# repeat above analysis for rarefied data
rarTempOnly <- lapply(transOtus, function(x)
	glm.nb(rareRich ~ country + temp + I(temp^2), data = x))

# richness as a quadratic function of temperature
# temperature coefficients interacting with country
rarTempInt <- lapply(transOtus, function(x)
	glm.nb(rareRich ~ country * temp + I(temp^2), data = x))

# richness as quadratic function of temperature
# all temp coefficients interacting with country, but not each other
rarTempFull <- lapply(transOtus, function(x)
	glm.nb(rareRich ~ country * temp + I(temp^2) * country,
	data = x))

# likelihood ratio test to test which hypothesis is more likely
lapply(1:3, function(x)
	anova(rarTempOnly[[x]], rarTempInt[[x]], rarTempFull[[x]]))

# get observed data into similar format
richRates <- lapply(transOtus, function(x)
	x[, .(temp = temp, value = hill_1/libSize, country = country)])

richRates <- rbindlist(
	Map(cbind, richRates, taxon = c("Archaea", "Bacteria", "Eukarya")))

# set up prediction data within temp ranges of original data
predData <- lapply(c("Archaea", "Bacteria", "Eukarya"), function(taxon)
	tempRanges[variable == taxon,
	.(temp = seq(minT, maxT, 0.05),
		taxon = taxon,
		libSize = 1000),
	by = country])

# use models to make predictions
predData <- lapply(1:3, function(x)
	predData[[x]][, ":="(
		value = predict(tempFull[[x]], newdata = predData[[x]], type = "link"),
		se = predict(tempFull[[x]], newdata = predData[[x]], type = "link",
			se.fit = T)$se.fit)])

# calculate prediction 95% conf intervals
predData <- lapply(1:3, function(x)
	predData[[x]][, ":="(uppCI = exp(value + 1.96 * se),
		lowCI = exp(value - 1.96 * se))])

# melt predictions into long format
longPred <- rbindlist(predData)

# exponentiate predictions back to response scale
longPred[, value := exp(value)]

richnessPlots <- ggplot(longPred,
		aes(x = temp, y = value/1000, col = country)) +
	geom_ribbon(
		aes(ymin = lowCI/1000, ymax = uppCI/1000, fill = country,
			col = country), alpha = 0.3, linetype = 2, colour = NA) +
	geom_point(data = richRates,
		aes(x = temp, y = value, col = country), size = 3, alpha = 0.6) +
	scale_color_viridis(discrete = T, begin = 0.8, end = 0) +
	scale_fill_viridis(discrete = T, begin = 0.8, end = 0) +
	coord_cartesian(xlim = c(0, 30)) +
	geom_line() +
	theme_bw() +
	labs(x = "Stream temperature (°C)", y = "Normalised OTU richness") +
	facet_grid(taxon~country, scales = "free_y") +
	theme(legend.position = "none",
		aspect.ratio = 1,
		panel.grid = element_blank(),
		axis.text = element_text(size = 16),
		axis.title = element_text(size = 18),
		strip.text = element_text(size = 14))

ggsave("../figures/Figure_4.pdf", richnessPlots, height = 6, width = 10,
	device = "pdf")

# Do evenness models
evenNone <- lapply(transOtus, function(x)
	glm(hill_3 ~ country + temp + I(temp^2), data = x,
	family = Gamma(link = "log")))

# richness as a quadratic function of temperature
# temperature coefficients interacting with country
evenMin <- lapply(transOtus, function(x)
	glm(hill_3 ~ country * temp + I(temp^2), data = x,
	family = Gamma(link = "log")))

# richness as quadratic function of temperature
# all temp coefficients interacting with country, but not each other
evenFull <- lapply(transOtus, function(x)
	glm(hill_3 ~ country * temp + I(temp^2) * country, data = x,
	family = Gamma(link = "log")))

# lr test
lapply(1:3, function(x)
	anova(evenNone[[x]], evenMin[[x]], evenFull[[x]], test = "F"))

# get R squared and AIC vals
# get AIC and explained deviance for each model for each taxon
evenAic <- data.table(Taxon = names(allOtus),
	none = sapply(evenNone, AIC),
	minimal = sapply(evenMin, AIC),
	full = sapply(evenFull, AIC))

evenAic <- melt(evenAic, id.vars = "Taxon", variable.name = "Interaction",
	value.name = "AIC_score")

evenDev <- data.table(Taxon = names(allOtus),
	none = sapply(evenNone, Dsquared),
	minimal = sapply(evenMin, Dsquared),
	full = sapply(evenFull, Dsquared))

evenDev <- melt(evenDev, id.vars = "Taxon", variable.name = "Interaction",
	value.name = "Explained_deviance")

evenDev[, model := "Evenness"]
richDev[, model := "Richness"]

modDev <- rbindlist(list(evenDev, richDev))

# need to adjust model labels
modDev[, Interaction := factor(Interaction, levels = levels(Interaction),
	labels = c(
		"No site interaction", "Minimal site interaction",
		"Full site interaction"))]

devPlot <- ggplot(modDev,
		aes(x = Interaction, y = Explained_deviance, group = model, col = model)) +
	geom_point() +
	geom_line() +
	facet_wrap(~ Taxon, ncol = 1) +
	labs(x = "Model", y = "Explained\ndeviance", col = NULL) +
	theme_bw() +
	theme(axis.text = element_text(size = 16),
		axis.title = element_text(size = 18),
		axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
		panel.grid = element_blank(),
		legend.text = element_text(size = 14),
		strip.text.x = element_text(size = 14))


evenModelFits <- merge(evenAic, evenDev, by = c("Taxon", "Interaction"))

fwrite(evenModelFits, "Table_S2.csv")

# add model predictions to predData
predData <- lapply(1:3, function(x)
	predData[[x]][, ":="(
		evenness = predict(evenFull[[x]], newdata = predData[[x]], type = "link"),
		evennessStdErr = predict(evenFull[[x]], newdata = predData[[x]],
			se.fit = T, type = "link")$se.fit)])

# calculate prediction 95% conf intervals
predData <- lapply(1:3, function(x)
	predData[[x]][, ":="(evenUppCI = exp(evenness + 1.96 * se),
		evenLowCI = exp(evenness - 1.96 * se))])

# get observed data into similar format
evenData <- lapply(transOtus, function(x)
	x[, .(temp = temp, evenness = hill_3, country = country)])

evenData <- rbindlist(
	Map(cbind, evenData, taxon = c("Archaea", "Bacteria", "Eukarya")))

# melt predictions into long format
evenPred <- rbindlist(predData)

# exponentiate predictions back to response scale
evenPred[, evenness := exp(evenness)]

# design figure + make table
evennessPlots <- ggplot(evenPred,
			aes(x = temp, y = evenness, col = country)) +
		geom_ribbon(
			aes(ymin = evenLowCI, ymax = evenUppCI, fill = country, col = country),
			alpha = 0.3, linetype = 2, colour = NA) +
		geom_point(data = evenData,
			aes(x = temp, y = evenness, col = country), size = 3, alpha = 0.6) +
		scale_color_viridis(discrete = T, begin = 0.8, end = 0) +
		scale_fill_viridis(discrete = T, begin = 0.8, end = 0) +
		coord_cartesian(xlim = c(0, 30)) +
		geom_line() +
		theme_bw() +
		labs(x = "Stream temperature (°C)", y = "Community evenness") +
		facet_grid(taxon ~ country, scales = "free_y") +
		theme(legend.position = "none",
			aspect.ratio = 1,
			panel.grid = element_blank(),
			axis.text = element_text(size = 16),
			axis.title = element_text(size = 18),
			strip.text = element_text(size = 14))

ggsave("../figures/Figure_5.pdf", evennessPlots, height = 6, width = 10,
	device = "pdf")

# create dataframe of approximate study sites
coords <- data.table(lat = c(64.3, 69.55, 79.38, 65, 55.89),
	long = c(-21.1, -53.52, 13.43, -150.63, 159.67),
	site = c("Iceland", "Greenland", "Svalbard", "Alaska", "Kamchatka"))

# reorder factor levels
coords <- coords[order(site), ]

# get international borders


ggsave("../figures/Figure_S1.pdf", siteMap, height = 5, width = 8,
	device = "pdf")
