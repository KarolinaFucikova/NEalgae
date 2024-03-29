# desmid community analyses for historical and new floristic data
# all records deposited in iNaturalist with images and geolocations
# read in the data
desmids <- read.csv(file = "data/Species_occurrences.csv", row.names =1)

#leave out the genus column for ordination analysis
species_occ <- desmids[c(1:211),c(1:6)]
# this data set is cleaned up, compared to the xlsx file which contains notes and preliminary analyses
# it also has 0s where the original has blanks for species absences

# first, get species richness for each site and time period
colSums(species_occ)

library(vegan)
# the first part of this code is commented out - it does not make sense to include HB1967 in the analysis
# the commands that compute dissimilarity require the data samples to be rows
# we will transpose 
#transp_occ <- t(species_occ)
# calculate the Jaccard index using
#jac <- 1-vegdist(transp_occ,method="jaccard")
# now we'll make a hierarchical cluster object, which will be the basis of the dendrogram
#occ_hc <-hclust(jac)
# and finally the dendrogram
#plot(occ_hc)

# HC1967 only has 3 species present; the study wasn't really comparable to the others
# excluding HC1967 first
clean_occ <- subset(species_occ, select = -HB1967)

# the commands that compute dissimilarity require the data samples to be rows
# we will transpose 
transp_occ <- t(clean_occ)

# calculate the Jaccard index using
# Jaccard is appropriate for presence/absence data
# the default is to calculate dissimilarity, so similarity is calculated by 1-dissimilarity
jac <- 1-vegdist(transp_occ,method="jaccard")
jac
# the following line will also return the similarity table, if desired:
#1-dist(transp_occ, method="binary")

# bray-curtis (sorensen) is also appropriate but slightly different:
bray <- 1-vegdist(transp_occ, method="bray")
bray

# now we'll make a hierarchical cluster object, which will be the basis of the dendrogram
# here we do want to use the default dissimilarity/distance measures
# calculate them first
jac1 <- vegdist(transp_occ, method="jaccard")
occ_hc <-hclust(jac1)
# and finally the dendrogram
plot(occ_hc, xlab="Sampled Localities", ylab="Jaccard Dissimilarity", main = NULL, cex=1.5, cex.lab=1.3, cex.axis=1.5, lwd=2)

# same for bray-curtis dendrogram
bray1 <- vegdist(transp_occ, method="bray")
occ_hcbray <-hclust(bray1)
plot(occ_hcbray, xlab="Sampled Localities", ylab="Bray-Curtis Dissimilarity", main = NULL, cex=1.5, cex.lab=1.3, cex.axis=1.5, lwd=2)
# the plots are about the same, not surprisingly

# trying detrended correspondence analysis
# probably don't have enough data for it but seems most appropriate according to
# https://www.davidzeleny.net/anadat-r/doku.php/en:ca_dca
# there are a lot of zeroes, so PCA and similar analyses are not ideal
# it's only presence/absence data, not abundance or relative abundance
# again have to exclude HB1967 as the three species that are found there will look like outliers

ord <- decorana(t(clean_occ))
ord
summary(ord)

# plotting the first two DCAs
plot(ord, type = "n")
points(ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(ord, display = "spec", cex=0.7, col="blue")

# attempting to label just particular points, automatically without the locator function
# can get species coordinates from the summary table
# https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/scores
# here the species farthest along axes are:
species_toplot <- c("Desmidium_asymmetricum")
# this list can be changed, added to, etc.)
# it won't automatically plot but makes it easier to find coordinates

table_toplot <- scores(ord, display="species")
# just the first two DCAs"
table_toplot <- table_toplot[,c(1:2)]
# make a table of just the DCA1 and 2 coordinates of just the species on the list
species_toplot_table <- table_toplot[rownames(table_toplot) %in% species_toplot, ]


# apparently can add the clustergram, but it looks awful
plot(ord, type = "n", ylim = c(-2,3), xlim = c(-2,3))
ordicluster(ord,occ_hcbray)
points(ord, display = "sites", pch = 2, col="#56B4E9", cex = 1.5)
points(ord, display = "species", pch = 1, col="#0072B2", cex = 1)
text(ord, display = "sites", cex=1, col="black", adj = c(1.2,0))
#points(-2.2656238, -2.121114, col = "magenta", pch=19)
legend(3, 2, c("sites", "species"),
       col=c("#56B4E9", "#0072B2"), pch=c(2,1), cex=1)
# not adding species because that makes the plot illegible
# xlim and ylim don't work for some reason

# nicer plot with ggvegan:
library("ggvegan")
library("ggplot2")
library("viridis")
library("tidyverse")
autoplot(ord, layers = "sites", arrows = FALSE)
auplot <- autoplot(ord, layers = c("sites","species"), arrows = FALSE)
auplot
aufort <- fortify(ord, display = c("sites","species"))

# extract site coords (points)
fort_sites <- aufort %>% 
  filter(Score == "sites")

# extract species coords
fort_species <- aufort %>% 
  filter(Score == "species")

# selecting just the species whose names I want to show on plot, for example:
#fort_species_filtered <- fort_species[c(3:33),c(1:6)]
# genera: rows 1 Actinotaenium, 2 Bambusina, 3:33 Closterium, 34:76 Cosmarium, 
# 77 Cylindrocystis, 78:84 Desmidium, 85:86 Docidium, 87:107 Euastrum
# 108 Gonatozygon, 109 Haplotaenium, 110 Hyalotheca, 111 Mesotaenium, 
# 112:127 Micrasterias, 128:129 Netrium, 130 Onychonema, 131:133 Penium
# 134:143 Pleurotaenium, 144 Roya, 145 Spinoclosterium, 146:147 Spinocosmarium
# 148:149 Spirotaenia, 150:151 Spondylosium, 152:190 Staurastrum
# 191:200 Staurodesmus, 201:202 Teilingia, 203:204 Tetmemorus, 205:206 Triploceras
# 207:211 Xanthidium

# simple version of plot:
fortify_plot <- ggplot() +
  geom_point(data = fort_species, aes(x = DCA1, y = DCA2)) +
  geom_point(data = fort_sites, aes(x = DCA1, y = DCA2), size = 3, shape = 2) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_text(data = fort_sites, aes(x = DCA1, y = DCA2, label = Label), nudge_x = 0.25, nudge_y = 0.25) +
  #geom_text(data = fort_species_filtered, aes(x = DCA1, y = DCA2, label = Label)) +
  labs(x = "DCA1",
       y = "DCA2",
       title = "DCA - Historical and current desmid communities")
fortify_plot

# can be made fancier version with colorblind friendly palette, if number of genera are reduced
#cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#fortify_plot3 <- fortify_plot2 + scale_colour_manual(values=cbp1)
#fortify_plot3

# trying to plot with labels for notable genera; too many points overlap to be visible
fort_species$genus <- desmids$Genus
# comment out the genus that you WANT to be labeled
fort_species$genus[fort_species$genus == "Actinotaenium"] <- "other" 
fort_species$genus[fort_species$genus == "Bambusina"] <- "other" 
fort_species$genus[fort_species$genus == "Closterium"] <- "other" 
fort_species$genus[fort_species$genus == "Cosmarium"] <- "other" 
fort_species$genus[fort_species$genus == "Cylindrocystis"] <- "other" 
fort_species$genus[fort_species$genus == "Desmidium"] <- "other"
#fort_species$genus[fort_species$genus == "Docidium"] <- "other" 
fort_species$genus[fort_species$genus == "Euastrum"] <- "other" 
fort_species$genus[fort_species$genus == "Gonatozygon"] <- "other" 
fort_species$genus[fort_species$genus == "Haplotaenium"] <- "other" 
#fort_species$genus[fort_species$genus == "Hyalotheca"] <- "other" 
fort_species$genus[fort_species$genus == "Mesotaenium"] <- "other" 
fort_species$genus[fort_species$genus == "Micrasterias"] <- "other"
fort_species$genus[fort_species$genus == "Netrium"] <- "other"
#fort_species$genus[fort_species$genus == "Onychonema"] <- "other" 
fort_species$genus[fort_species$genus == "Penium"] <- "other" 
fort_species$genus[fort_species$genus == "Pleurotaenium"] <- "other" 
fort_species$genus[fort_species$genus == "Roya"] <- "other"
fort_species$genus[fort_species$genus == "Spinoclosterium"] <- "other" 
fort_species$genus[fort_species$genus == "Spinocosmarium"] <- "other" 
fort_species$genus[fort_species$genus == "Spirotaenia"] <- "other" 
fort_species$genus[fort_species$genus == "Spondylosium"] <- "other" 
fort_species$genus[fort_species$genus == "Staurastrum"] <- "other" 
fort_species$genus[fort_species$genus == "Staurodesmus"] <- "other" 
#fort_species$genus[fort_species$genus == "Teilingia"] <- "other" 
fort_species$genus[fort_species$genus == "Tetmemorus"] <- "other" 
fort_species$genus[fort_species$genus == "Triploceras"] <- "other"
fort_species$genus[fort_species$genus == "Xanthidium"] <- "other" 

fortify_plot4 <- ggplot() +
  geom_point(data = fort_species, aes(x = DCA1, y = DCA2)) +
  geom_point(data = fort_sites, aes(x = DCA1, y = DCA2), size = 3, shape = 2) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_text(data = fort_sites, aes(x = DCA1, y = DCA2, label = Label), nudge_x = 0.25, nudge_y = 0.25) +
  #geom_text(data = fort_species, aes(x = DCA1, y = DCA2, label = genus), nudge_x = 0.25, nudge_y = 0.25) +
  geom_point(data = fort_species, aes(x = DCA1, y = DCA2, color = genus)) +
  labs(x = "DCA1",
       y = "DCA2",
       title = "DCA - Historical and current desmid communities")
fortify_plot4
# too many points overlap, can't see the genera
