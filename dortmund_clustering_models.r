{\rtf1\ansi\ansicpg1252\cocoartf1561\cocoasubrtf100
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Part 1                                                                                   14007080\
\
library(stats, mclust, fpc, cluster)\
\
dortmund <- original <- read.table("Dortmund_G3019ica.dat")\
\
str(dortmund)\
\
# Let's start with some quick exploratory analysis. \
\
# It seems sensible to transform many of the variables into proportions; \
# doing so will allow better comparison but remain easily interpretable.\
\
# In order to construct a meaningful, readily-interpretable cluster analysis\
# of the dataset we might benefit from some dimensional reduction. Let's \
# investigate correlations between variables to see if there are any easy\
# aggregations we might be able to do. \
\
\
dortmund$pop <- rowSums(dortmund[, c("male", "female")])\
\
to_scale <- !(names(dortmund) %in% c("area_buildings", "buildings_until_1900", \
                                     "buildings_1900.1918", "buildings_1919.1948", \
                                     "buildings_1949.1957", "buildings_1958.1962", \
                                     "buildings_1963.1972", "buildings_1973.1982", \
                                     "buildings_1983.1992", "buildings_1993.2001", \
                                     "households", "pop"))\
\
dortmund[, to_scale]  <- dortmund[, to_scale] / dortmund$pop\
\
dortmund <- dortmund[, !(names(dortmund) %in% c("males"))]\
\
cor(dortmund[, c("age_under_26", "age_26.35", "age_36.45", "age_46.55", \
                 "age_56.65", "age_above_65")])\
\
# Interestingly there are more children than age_under_26; presumably the \
# latter of these refers to 18 to 26 year olds (thus our pop refers to \
# adult population not total population, this seems fine).\
\
# There is no indication here of which groups to combine; we will instead \
# remove the age groups 26 through to 65 and leave ourselves with data \
# reflecting the percentage of young and old inhabitants. Later we will\
# check our clusterings against the dropped variables to check they do not\
# disagree with our clusterings. TODO\
\
dortmund <- dortmund[, !(names(dortmund) %in% c("age_26.35", "age_36.45", \
                                                "age_46.55", "age_56.65"))]\
\
dortmund$buildings_total <- rowSums(dortmund[, c("buildings_until_1900", \
                                                 "buildings_1900.1918", \
                                                 "buildings_1919.1948", \
                                                 "buildings_1949.1957", \
                                                 "buildings_1958.1962", \
                                                 "buildings_1963.1972", \
                                                 "buildings_1973.1982", \
                                                 "buildings_1983.1992", \
                                                 "buildings_1993.2001")])\
\
to_scale <- c("buildings_until_1900", "buildings_1900.1918", \
              "buildings_1919.1948", "buildings_1949.1957", \
              "buildings_1958.1962", "buildings_1963.1972", \
              "buildings_1973.1982", "buildings_1983.1992", \
              "buildings_1993.2001")\
\
dortmund[, to_scale] <- dortmund[, to_scale] / dortmund$buildings_total\
\
cor(dortmund[, 10:18])\
\
# Again, none of the variables are particularly linearly correlated so it\
# is not obvious what to do with them. Let's look at correlations across\
# the dataset; this isn't easy given the size of the dataset. \
\
abs(cor(dortmund))>0.85 & cor(dortmund)!=1\
\
# Strongly correlated variables include unemployed and benefits, pop and \
# area_buildings, male and female, moves_in and moves_out and most \
# interestingly births and unemployed. Households is also strongly \
# correlated with a few variables, thus we will drop households, \
# aggregated moves_out and moves_in into net_moves, drop benefits and \
# males (which is perfectly negatively correlated with females after\
# scaling by population).\
\
dortmund$net_moves <- dortmund$moves_in - dortmund$moves_out\
\
dortmund <- dortmund[, !(names(dortmund) %in% c("benefits", "moves_in",\
                                                "moves_out", "male", \
                                                "households"))]\
\
# So we are left with quite a large variable set. We could trim it down for \
# the second clustering but for now let's move on.\
\
#############################################################\
#                                                                     				                  #\
#                                   Clustering 1: Gaussian Mixtures                                #\
#                                                                                                                      #\
#############################################################\
\
# For our first clustering we will use a Gaussian Mixture model. We do this\
# rather than the simpler K-means as there is no reason to assume we have \
# spherical covariances.\
\
clust1.1 <- Mclust(dortmund)\
summary(clust1.1)\
clust1.1$BIC\
\
# So the best Gaussian Mixture models with our current variable set are all \
# single component models... sadly these are useless for our purposes. \
# Let's try again but enforce a minimum of 2 clusters.\
\
clust1.2 <- Mclust(dortmund, G=2:20)\
summary(clust1.2)\
clust1.2$BIC\
\
# This is a big reduction in the BIC, however 6 different components is\
# far more useful in identifying groups of districts. Our new model is \
# diagonal with varying volume and shape. Let's do some visualisation. \
\
pairs(dortmund[, 1:5], col=clust1.2$classification, \
      pch=as.character(clust1.2$classification))\
\
# There is one district has an extraordinarily high death rate. \
\
row.names(dortmund)[which.max(dortmund$deaths)]\
original[which.max(dortmund$deaths), ]\
\
# This seems like it could be an error in the data somehow (with a \
# deathrate close to 2), looking at the original untransformed dataset it\
# has 135 deaths in a population of 69. This district, Rombergpark, does\
# have an elderly population but even this is no explanation. Google offers\
# no quick explanations either... Let's investigate the effect of this \
# district on the model fitting. \
\
clust1.3 <- Mclust(dortmund[!row.names(dortmund)=="Rombergpark", ], G=2:20)\
\
adjustedRandIndex(clust1.3$classification, \
                  clust1.2$classification[!row.names(dortmund)=="Rombergpark"])\
\
# So the Adjust Rand Index between the models with and without Rombergpark\
# is as low as 0.43... Seems we must therefore do something about this. As\
# we can't simply exclude Rombergpark from the clustering, let's instead \
# replace its deathrate with the mean deathrate. \
\
new_deaths <- mean(dortmund[!row.names(dortmund)=="Rombergpark", "deaths"])\
dortmund[row.names(dortmund)=="Rombergpark", "deaths"] <- new_deaths\
\
clust1.4 <- Mclust(dortmund, G=2:20)\
summary(clust1.4)\
clust1.4$BIC\
\
adjustedRandIndex(clust1.2$classification, clust1.4$classification)\
\
# This new clustering, with Rombergpark death rate set to the column mean, \
# is quite different to our previous clustering with the original death\
# rate according to ARI. This is a bit worrying, however it seems much more\
# sensible than using the likely erroneous death rate. The new model has \
# only 3 components; let's see how it looks. \
\
pairs(dortmund[, c(1, 19, 4, 5, 24)], \
      col=clust1.4$classification, \
      #pch=clust1.4$classification)\
      pch=as.character(clust1.4$classification)) # This is really slow and\
                                                  # I can't find clusym...\
\
# This clustering looks pretty reasonable; let's stick with this. \
\
#############################################################\
#                                                                     				                  #\
#                                        Clustering 2: Hierachical                                       #\
#                                                                                                                      #\
#############################################################\
\
# For this clustering let's focus on the demographic variables in order to \
# give a more socially orientated clustering. \
\
dortmund2 <- dortmund[, (names(dortmund) %in% c("unemployed", "births", \
                                                "deaths", "net_moves", \
                                                "children", "female", \
                                                "social_insurance", \
                                                "age_under_26", \
                                                "age_above_65", "pop"))]\
\
# Let's try using a hierachical clustering method. We have to decide upon\
# a method; we will use complete linkage in order to give us more \
# homogeneous clusters. This should be more useful for our purposes: for \
# example if we are able to identify groups of similar districts then we \
# might be able to group them together for deciding policy.\
\
# We also need to choose our dissimilarity; euclidean seems sensible as we\
# have hopefully removed some of the most obvious correlations. We will \
# also try euclidean-squared as this may result in more homogeneity. \
# Hierachical clustering with euclidean distances is sensitve to scaling; \
# currently pop is orders of magnitude greater than the other variables so\
# let's scale the whole matrix. \
\
dortmund2.scaled <- scale(dortmund2)\
\
dortmund2.euc <- dist(dortmund2.scaled)\
dortmund2.euc2 <- dortmund2.euc ** 2\
\
clust2.1k <- hclust(dortmund2.euc, method="complete")\
clust2.2k <- hclust(dortmund2.euc2, method="complete")\
\
plot(clust2.1k, main="Cluster Dendrogram for clust2.1", xlab="")\
\
# From the dendrogram, K=5 seems sensible from the dendrogram. Let's look\
# at the gap statistics too for comparison. To do this we must define a \
# hclust-based function that results in a k cluster.\
\
hclustk <- function(x, k) \{\
  list(cluster=cutree(hclust(dist(x), method="complete"), k=k))\
\}\
\
gap_stat <- clusGap(dortmund2, hclustk, K.max=50, B=100)\
\
print(gap_stat, method="Tibs2001SEmax")\
plot(gap_stat)\
\
# Based on the gap statistic it seems K=5 is the most suitable choice.\
# Interesting that the gap statistic is generally increasing with K up\
# past K=20... For our purposes a smaller number of clusters is \
# probably more useful. Let's look at our clustering.\
\
clust2.1 <- cutree(clust2.1k, 5)\
clust2.2 <- cutree(clust2.2k, 5)\
\
adjustedRandIndex(clust2.1,  clust2.2)\
\
# The two distances give the same clustering at level K=5.\
\
pairs(dortmund[, c(1, 19, 3, 4, 24)], \
      col=clust2.1, \
      #pch=clust2.1)\
      pch=as.character(clust2.1))\
\
pairs(dortmund[, c(17, 18, 20, 21)], \
      col=clust2.1, \
      #pch=clust2.1)\
      pch=as.character(clust2.1))\
\
# This clustering looks quite good and readily interpretible. Let's \
# stick with this as our second clustering. \
\
#############################################################\
#                                                                     				                  #\
#                                                  Comparison                                                 #\
#                                                                                                                      #\
#############################################################\
\
# Let's start by looking at ARI and the clusterings for a few variables. \
\
adjustedRandIndex(clust1.4$classification, clust2.1)\
\
par(mfrow=c(3, 2), oma=c(0, 0, 2, 0))\
\
with(data=dortmund, \{\
  plot(births ~ deaths, col=clust1.4$classification, \
       pch=as.character(clust1.4$classification))\
  plot(births ~ deaths, col=clust2.1, pch=as.character(clust2.1))\
\
  plot(unemployed ~ pop, col=clust1.4$classification,\
       pch=as.character(clust1.4$classification))\
  plot(unemployed ~ pop, col=clust2.1, pch=as.character(clust2.1))\
  \
  plot(cars ~ buildings_total, col=clust1.4$classification,\
       pch=as.character(clust1.4$classification))\
  plot(cars ~ buildings_total, col=clust2.1, pch=as.character(clust2.1))\
\})\
\
mtext("clust1.4", side=3, line=-2, at=grconvertX(-0.75,"npc","nic"), \
      outer=TRUE, cex=2)\
mtext("clust2.1", side=3, line=-2, at=grconvertX(0.5,"npc","nic"), \
      outer=TRUE, cex=2)\
\
table(clust2.1); table(clust1.4$classification)\
\
#############################################################\
#                                                                     				                  #\
#                              Cluster Selection and Interpretation                                #\
#                                                                                                                      #\
#############################################################\
\
# Print the districts and the cluster they're assigned to in our clust2.1\
for(i in 1:dim(dortmund)[1]) \{\
  # seems weird we have to print(sprintf(...)) but that's just how R is\
  print(sprintf("district: %25s %s", row.names(dortmund)[i], clust2.1[i]))\
\}\
\
\
\
\
\
\
\
\
\
\
\
\
\
}