# desmid community analyses for historical and new floristic data
# all records deposited in iNaturalist with images and geolocations
species_occ <- read.csv(file = "data/Species_occurrences.csv", row.names =1)

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
points(ord, display = "sites", cex = 0.8, pch=21, col="yellow", bg="yellow")
text(ord, display = "sites", cex=0.7, col="black")
# not adding species because that makes the plot illegible
# xlim and ylim don't work for some reason

# will continue with PCA test


