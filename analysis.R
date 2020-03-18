library(data.table)
library(ecolFudge)
library(hillR)
library(vegan)
library(MASS)
library(modEvA)
library(ggplot2)
library(ggparl) #  boxplot/jitter plot hybrids
library(mapdata)
library(sf)
library(viridis)
library(ggrepel)
library(betapart)
library(cowplot)
library(mvabund)
library(maptools)
library(parallel)

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

# create cols for sp richness, library size, and taxon, ready for plotting
transOtus <- lapply(1:3, function(x)
	transOtus[[x]][, c("hill_1", "hill_2", "hill_3", "libSize", "taxon") := list(
		hill_taxa(.SD, q = 0), hill_taxa(.SD, q = 1), hill_taxa(.SD, q = 2),
		rowSums(.SD), rep(names(transOtus)[x], times = nrow(transOtus[[x]]))),
		.SDcols = otuCols[[x]]])

# reformat sample labels to match metadata
transOtus <- lapply(transOtus, function(x)
	x[, sample := sapply(strsplit(sample, "_"), "[[", 1)])

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

# rarefy dat
par(mfrow = c(1, 3))
lapply(1:3, function(x) quickRareCurve(transOtus[[x]][, .SD,
	.SDcols = otuCols[[x]]], label = F))

# inspect sample sizes to determine samples to be discarded
lapply(transOtus, function(x) sort(x$libSize))
smallSamples <- c(600, 20000, 2000)

# remove small samples
transOtus <- lapply(1:3, function(x) transOtus[[x]][libSize > smallSamples[x]])

# merge OTU tables with main dataframe
sampNames <- c("archSampName", "bacSampName", "eukSampName")
transOtus <- lapply(1:3, function(x)
	merge(transOtus[[x]], arctic, by.x = "sample", by.y = sampNames[x]))

# remove soil samples
transOtus <- lapply(1:3, function(x) transOtus[[x]][sampleType == "sed", ])

# calculate rarefied richness to verify results
transOtus <- lapply(1:3, function(x)
	transOtus[[x]][, rareRich := round(rarefy(.SD, sample = min(libSize))),
	.SDcols = otuCols[[x]]])

	# partition beta diversity into nestedness and turnover components
	# calculate matrices
	# rarefy Otu tables
	# need to coerce to list of matrices.
	mc <- getOption("cores", detectCores())
	rarefiedTables <- lapply(1:3, function(x)
		mclapply(unique(transOtus[[x]]$country), mc.cores = mc, function(z)
			transOtus[[x]][country == z, rrarefy(.SD, sample = min(libSize)),
			.SDcols = otuCols[[x]]]))
	rarefiedTables <- unlist(rarefiedTables, recursive = F)
	# [[1:5]] <- arch
	# [[6:10]] <- bacteria
	# [[11:15]] <- eukarya

	# function to turn abundance matrix into binary matrix
	makeBinary <- function(x){
		x <- ifelse(x > 0, 1, 0)
	}

	# turn all rarefied tables into binary tables
	rarefiedTables <- lapply(rarefiedTables, makeBinary)
	names(rarefiedTables) <- unlist(lapply(c("Archaea", "Bacteria", "Eukarya"),
		paste, unique(transOtus[[1]]$country), sep = "_"))

	# compute beta diversity matrices
	betaMatrices <- lapply(rarefiedTables, beta.pair)
	names(betaMatrices) <- unlist(lapply(c("Archaea", "Bacteria", "Eukarya"),
		paste, unique(transOtus[[1]]$country), sep = "_"))

	# list of 15 elements, each is 3 long e.g. turnover, nestedness, similarity
	newMatrixNames <- unlist(lapply(names(betaMatrices), paste0, c("_Turnover", "_Nestedness", "_Similarity")))

	# make temperature distance matrices
	tempMatrices <- lapply(1:3, function(y)
		lapply(unique(transOtus[[y]]$country), function(x)
			transOtus[[y]][country == x, vegdist(temp, "euclid")]))
	tempMatrices <- unlist(tempMatrices, recursive = F)

	bm <- unlist(betaMatrices, recursive = F)
	matLens <- sapply(bm, length)
	nameCol <- unlist(lapply(1:45, function(x) rep(newMatrixNames[x], times = matLens[x])))
	matrData <- data.table(dist = unlist(bm), temp = unlist(rep(tempMatrices, each = 3)), id = nameCol)

	# add taxon and country cols
	matrData[, c("taxon", "country", "distType") :=
		lapply(1:3, function(x) sapply(strsplit(id, "_"), "[[", x))
	]

	# add matrix element ID in for casting
matrData[, row := 1:.N, by = c("taxon", "country", "distType")]

	# reshape data to calculate % contribution of nestedness + turnover = similarity
turnContr <- dcast(matrData, taxon + country + row ~ distType,
		value.var =  "dist")

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

	# conduct mantel tests
mantelTests <- mclapply(1:45, mc.cores = mc, function(x)
		mantel(bm[[x]], tm[[x]], permutations = 10000))
	names(mantelTests) <- names(bm)

# get slopes and intercepts of beta diversity metrics
betaFits <- matrData[, .(intercept = coef(lm(dist ~ temp))[1],
	slope = coef(lm(dist ~ temp))[2],
	pVal = summary(lm(dist ~ temp))$coefficients[, 4][2],
	minX = min(temp),
	maxX = max(temp)),
	by = c("taxon", "country", "distType")]

	# add mantel P values to coefficient table and calculate maximal y value of
	# best fit line
betaFits[, c("mantelP", "mantelR", "mantel_P_value", "maxY") := list(
		sapply(mantelTests, function(x) ifelse(x$signif < 0.05, T, F)),
		sapply(mantelTests, function(x) x$statistic),
		sapply(mantelTests, function(x) x$signif),
		intercept + (maxX*slope))]

	# reorder factor levels
betaFits[, c("maxY", "distType", "mantelP") := list(
		ifelse(maxY > 1, 1, maxY),
		factor(distType, levels = c("Similarity", "Turnover", "Nestedness")),
		factor(mantelP, levels = c(TRUE, FALSE)))]

	# plot beta diversity indices
betaPlot <- ggplot(matrData[distType != "Similarity"],
			aes(x = temp, y = dist, col = distType)) +
		geom_point(size = 2.5, alpha = 0.2) +
		facet_grid(taxon ~ country) +
		geom_segment(data = betaFits[distType != "Similarity"],
			aes(x = minX, xend = maxX, y = intercept, yend = maxY, col = distType, linetype = mantelP, size = mantelP)) +
		labs(x = expression(Temperature~difference~(degree*C)),
			y = expression(beta-diversity),
			colour = expression(beta-diversity~component)) +
		theme_bw() +
		scale_color_viridis(discrete = T, end = 0.85) +
		scale_size_manual(values = c(1.1, 0.7), labels = c(
			expression(italic(P)<=0.05), expression(italic(P)>=0.05))) +
		scale_linetype_manual(values = c(1, 2), labels = c(
			expression(italic(P)<=0.05), expression(italic(P)>=0.05))) +
		theme(axis.text = element_text(size = 16),
			axis.title = element_text(size = 18),
			panel.grid.minor = element_blank(),
			panel.grid.major = element_blank(),
			strip.text = element_text(size = 14),
			legend.text = element_text(size = 14),
			legend.title = element_blank(),
			aspect.ratio = 1)

	ggsave("../figures/Figure_1.pdf", height = 7, width = 12, device = "pdf")


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
arcticMap <- map("worldHires", fill = T, plot = F)
arcticMap <- st_as_sf(arcticMap)

# create base map
siteMap <- ggplot(arcticMap[arcticMap$ID != "Antarctica", ]) +
	geom_sf(col = "grey", fill = "grey") +
	coord_sf(datum=NA) +
	geom_point(data = coords, aes(x = long, y = lat, col = site), size = 2.5) +
	scale_color_viridis(discrete = T, begin = 0.8, end = 0) +
	geom_label_repel(data = coords,
		aes(x = long, y = lat, label = site, col = site),
		fontface = "bold", segment.size = 0.5, nudge_y = 30 - coords$lat,
		direction = "x", size = 4, alpha = 0.8) +
	theme_void() +
	theme(legend.position = "none",
		panel.grid = element_blank())


ggsave("../figures/Figure_S1.pdf", siteMap, height = 5, width = 8,
	device = "pdf")
