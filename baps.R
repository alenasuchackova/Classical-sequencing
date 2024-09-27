## Bayesian analysis of population structure
############################################
##This is a script to compute Bayesian analysis of population structure (BAPS) (Corander et al. 2003) through R package rhierbaps (Cheng et al. 2013)
#the method is sorting sequences into clusters without a priori knowledge of their geographic origin
#Corander J., Waldmann P., Sillanpää M.J. Bayesian analysis of genetic differentiation between populations. Genetics. 2003; 163:367–374.
#Cheng L., Connor T.R., Sirén J., Aanensen D.M., Corander J. Hierarchical and spatially explicit clustering of DNA sequences with BAPS software. Mol. Biol. Evol. 2013; 30:1224–1228.

library(rhierbaps)
library(phytools)
library(ape)
#??rhierbaps

######################################################################
#PART 1 BAPS calculation##############################################
#we need a nexus alignment
baps.seq<-read.nexus.data(file="hippothoe_coi_final.nex")
baps.seq<-as.DNAbin(baps.seq)
snp.matice<-load_fasta(baps.seq)
snp.results<-hierBAPS(snp.matice,max.depth = 10, n.pops = 20, assignment.probs=TRUE) #max.depth=levels on which we calculate the clusters, n.pops=maximum number of clusters, n.extra.rounds=how many times to calculate
snp.results
names (snp.results) #partition.df=samples assigned into clusters on different levels, cluster.assignment.prob=probability of assignment of samples into individual populations/clusters, lml.list=log marginal likelihoods

write.table(snp.results$partition.df,file="hippothoe_baps.txt",sep="\t")
write.table(snp.results$cluster.assignment.prob, file="hippothoe_baps_assign_prob.txt",sep="\t")
snp.results$lml.list
#depth0 assumes no population structure

######################################################################
#PART2 calculation of best level (optional, it is wise to use level1)#
##Using log marginal likelihoods (lml.list) to visualize the best number of clusters (highest LML, lowest Delta L)
lml_list <- snp.results$lml.list
lml_means <- numeric(length(lml_list))

# calculate the mean lml (without NA values, these are omitted)
for (i in seq_along(lml_list)) {
  lml_values <- unlist(lml_list[[i]])
  lml_means[i] <- mean(lml_values, na.rm = TRUE)
}

# depth with the highest mean log marginal likelihood
best_depth <- which.max(lml_means)
cat("Mean LML values for each depth:\n")
print(lml_means)
cat("\nBest depth based on highest mean LML:", best_depth - 1, "\n")
## it is wise to visualize the results on a map and decide if to use level/depth1, or higher (if it still makes biological sense). Best lml is usually for depth 2-3

depths <- seq_along(lml_means) - 1
plot(depths, lml_means, type = "b", pch = 19, col = "blue",
     xlab = "Depth", ylab = "Mean Log Marginal Likelihood (LML)",
     main = "Mean LML across Depths",
     xaxt = "n")  
axis(1, at = depths, labels = depths)

delta_L <- diff(lml_means)
depths_delta <- seq_along(delta_L) - 1
df_delta <- data.frame(Depth = depths_delta, DeltaL = delta_L)
plot(depths_delta, delta_L, type = "b", pch = 19, col = "red",
     xlab = "Depth", ylab = "ΔL (Change in Mean LML)",
     main = "Change in Mean LML (ΔL) across Depths")



