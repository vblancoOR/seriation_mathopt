library(seriation)
library(cluster)

data("Psych24")

## create a dist object and also get rid of the one negative entry in the
## correlation matrix
cor_matrix <- cor(Psych24, method = "pearson")
d <- format(cor_matrix, digits = 4, nsmall = 0)
write.table(d, file = "psych24.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

data("Irish")
# Replace NA values with column means
filled_irish <- apply(Irish, 2, function(x) ifelse(is.na(x), mean(x, na.rm=TRUE), x))
scaled_irish <- scale(filled_irish)
# Compute Euclidean distance on filled data
euclidean_distances <- dist(scaled_irish)
d <- format(euclidean_distances, digits = 4, nsmall = 0)
write.table(d, file = "irish.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Install the package if not installed
#if (!require(vegan)) install.packages("vegan")

# Load the package
library(vegan)

# Compute the Jaccard index (similarity matrix)

data("Munsingen")

# x <- Munsingen[sample(nrow(Munsingen)),]
jaccard_matrix <- vegdist(Munsingen, method = "jaccard")
d <- format(jaccard_matrix, digits = 4, nsmall = 0)
write.table(d, file = "munsingen.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

data(Votes, package = "cba")

Votes <- cba::as.dummy(Votes[-17])

x <- Votes[sample(100),]
#Votes_f <- apply(Votes, 2, function(x) ifelse(is.na(x), mean(x, na.rm=TRUE), x))
jaccard_matrix <- vegdist(x, method = "jaccard")
d <- format(jaccard_matrix, digits = 4, nsmall = 0)
write.table(d, file = "votes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# data("Zoo")
# x <- scale(Zoo[, -17])
# d <- dist(x)
# write.table(d, file = "zoo.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# data("iris")
# x <- as.matrix(iris[-5])
# x <- x[sample(seq_len(nrow(x))), ]
# d <- dist(x)

# write.table(x, file = "iris.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# data("Wood")
# Wood <- Wood[sample(nrow(Wood)), sample(ncol(Wood))]

# write.table(Wood, file = "wood.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
