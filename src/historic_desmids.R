# desmid community analyses for historical and new floristic data
# all records deposited in iNaturalist with images and geolocations
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
jac <- 1-vegdist(transp_occ,method="jaccard")

# now we'll make a hierarchical cluster object, which will be the basis of the dendrogram
occ_hc <-hclust(jac)
# and finally the dendrogram
plot(occ_hc)

# trying detrended correspondence analysis
# probably don't have enough data for it but seems most appropriate according to
# https://www.davidzeleny.net/anadat-r/doku.php/en:ca_dca
# again have to exclude HB1967 as the three species that are found there will look like outliers

ord <- decorana(t(clean_occ))
ord
summary(ord)

# plotting the first two DCAs
plot(ord, type = "n")
points(ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(ord, display = "spec", cex=0.7, col="blue")

# apparently can add the clustergram
plot(ord, type = "n", ylim = c(-2,3), xlim = c(-2,3))
ordicluster(ord,occ_hc)
points(ord, display = "sites", pch = 2, col="#56B4E9", cex = 1.5)
points(ord, display = "species", pch = 1, col="#0072B2", cex = 1)
text(ord, display = "sites", cex=1, col="black", adj = c(1.2,0))
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
fort_species$genus[fort_species$genus == "Actinotaenium"] <- "" 
fort_species$genus[fort_species$genus == "Bambusina"] <- "" 
fort_species$genus[fort_species$genus == "Closterium"] <- "" 
fort_species$genus[fort_species$genus == "Cosmarium"] <- "" 
fort_species$genus[fort_species$genus == "Cylindrocystis"] <- "" 
#fort_species$genus[fort_species$genus == "Desmidium"] <- ""
fort_species$genus[fort_species$genus == "Docidium"] <- "" 
fort_species$genus[fort_species$genus == "Euastrum"] <- "" 
fort_species$genus[fort_species$genus == "Gonatozygon"] <- "" 
fort_species$genus[fort_species$genus == "Haplotaenium"] <- "" 
fort_species$genus[fort_species$genus == "Hyalotheca"] <- "" 
fort_species$genus[fort_species$genus == "Mesotaenium"] <- "" 
fort_species$genus[fort_species$genus == "Micrasterias"] <- ""
fort_species$genus[fort_species$genus == "Netrium"] <- ""
fort_species$genus[fort_species$genus == "Onychonema"] <- "" 
fort_species$genus[fort_species$genus == "Penium"] <- "" 
fort_species$genus[fort_species$genus == "Pleurotaenium"] <- "" 
fort_species$genus[fort_species$genus == "Roya"] <- ""
fort_species$genus[fort_species$genus == "Spinoclosterium"] <- "" 
fort_species$genus[fort_species$genus == "Spinocosmarium"] <- "" 
fort_species$genus[fort_species$genus == "Spirotaenia"] <- "" 
fort_species$genus[fort_species$genus == "Spondylosium"] <- "" 
fort_species$genus[fort_species$genus == "Staurastrum"] <- "" 
fort_species$genus[fort_species$genus == "Staurodesmus"] <- "" 
fort_species$genus[fort_species$genus == "Teilingia"] <- "" 
fort_species$genus[fort_species$genus == "Tetmemorus"] <- "" 
fort_species$genus[fort_species$genus == "Triploceras"] <- ""
fort_species$genus[fort_species$genus == "Xanthidium"] <- "" 

fortify_plot4 <- ggplot() +
  geom_point(data = fort_species, aes(x = DCA1, y = DCA2)) +
  geom_point(data = fort_sites, aes(x = DCA1, y = DCA2), size = 3, shape = 2) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_text(data = fort_sites, aes(x = DCA1, y = DCA2, label = Label), nudge_x = 0.25, nudge_y = 0.25) +
  geom_text(data = fort_species, aes(x = DCA1, y = DCA2, label = genus), nudge_x = 0.25, nudge_y = 0.25) +
  labs(x = "DCA1",
       y = "DCA2",
       title = "DCA - Historical and current desmid communities")
fortify_plot4
