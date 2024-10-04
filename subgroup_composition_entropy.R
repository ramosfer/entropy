# Script for data analysis accompanying the manuscript "Quantifying uncertainty 
# due to fission-fusion dynamics as a component of social complexity" by 
# Ramos-Fernandez et al. (resubmitted to Proceedings B, february 2018). 

# The script calculates entropy of subgroup composition from the group and
# individual perspectives, including the observed and bootstrap entropies. It
# also compares both using the KL divergence and JS distance. Finally, the
# script calculates the probability of observing the group in partitions of
# different size.

#### 0. Data load and cleanup ####
library(entropy) # load package for entropy calculations
library(copula) # load package for Stirling number calculations
library(stringr) # load package for doing string operations

# Read database: each row corresponds to a subgroup observation and each column
# corresponds to an individual that can be present (1) or absent (0)
m <- read.csv("subgroup_composition_chimpanzees_2008.csv") # substitute this filename for your file; should be in the same directory as this script
# m <- m[,1:29] # remove females from gelada 2014 dataset
View(m) # check content
str(m) # check types of variables and structure


# Minimum cleanup: make sure there are no NA or rows with all 0s
anyNA(m) # should be FALSE
m[is.na(m)] <- 0 # if not, change NA to 0
sgs <- apply(m, 1, sum)
obs <- apply(m, 2, sum)
m <- m[which(sgs>0), which(obs>0)] # take only the observations with at least one individual and the individuals with at least one observation
str(m) # check that each column consists of integer and the data frame dimensions correspond to what you had in your file

# Turn m into a matrix as it will be necessary for further manipulations
m <- as.matrix(m)

#### 1. Entropy from the whole group perspective ####
# Take each observation as a particular combination of 0s, and 1s, calculating
# empirical entropy for the whole set

# Convert matrix into a string of pasted combinations of 0s and 1s
m_mat <- array(NA, nrow(m))
for (i in 1:nrow(m)){
  m_mat[i] <- paste(m[i,], collapse = "")
}
head(m_mat) # check that this vector contains character values composed of combinations of 0s and 1s

# Count the number of unique compositions
length(unique(m_mat))

# Count the number of occurrences for each unique composition
t <- table(m_mat)
str(t)

# Convert to frequencies and check that they sum to 1 (requires loaded "entropy" library)
freqs.empirical(t)
sum(freqs.empirical(t)) # should be 1

# Calculate the information content of the observed data
estimated_h <- entropy.empirical(t, unit="log2")

# Which should be equal to the original formula for the Shannon entropy (equation 1 in the manuscript)
-sum(freqs.empirical(t) * log2(freqs.empirical(t)))

#### 2. Null bootstrap entropy from the subgroup perspective #### 
# Calculate the maximum entropy that would be expected for the dataset, given
# the observed distribution of subgroup size and the size of the data set

# Retrieve fk, the frequency distribution of subgroups of size k
n <- ncol(m)
hh <- hist(sgs, breaks=c(0:n))
fk <- freqs.empirical(hh$counts) # vector containing the frequencies of observations of subgroup size
Nk <- hh$counts # vector containing the number of observations of each subgroup size
No <- sum(Nk) # should be the same as the number of observations
fia <- vector("list", n)
nk_fi <- array(NA, n)
for (k in 1:n){
  fia[[k]] <- array(NA, Nk[k])
  for (i in 1:Nk[k]){
    fia[[k]][i] <- exp((i*log(1/choose(n,k)))+((Nk[k]-i)*log(1-(1/choose(n,k))))+lchoose(Nk[k],i))*(i/No)*log2(i/No)
  }
  nk_fi[k] <- choose(n,k)*sum(fia[[k]])
}
# Calculate the bootstrap entropy:
boot_h <- -sum(nk_fi, na.rm=T)

# Estimation of the same entropy through an actual bootstrap (for illustration purposes)
# Construct nreps bootstrapped versions of original database, maintaining the size of each 
# observed subgroup constant and only switching 0 and 1 between columns
nreps = 1000
m_boot <- array(NA, dim=c(nreps, nrow(m), ncol(m)))
for (i in 1:nreps){
  for (j in 1:nrow(m)){
    m_boot[i,j,] <- sample(m[j,])
  }
}

# Construct pasted version of all bootstrapped versions ###### THERE WAS AN ERROR HERE, THE LOOP WAS WRONGLY CONSTRUCTED #######
m_boot_paste <- matrix(NA, nreps, nrow(m))
for (i in 1:nreps){
  for (j in 1:nrow(m)){
  m_boot_paste[i,j] <- paste(m_boot[i,j,], collapse = "")
  }
}

# Calculate entropy for all
e_boot <- array(NA, nreps)
for (i in 1:nreps){
  e_boot[i] <- entropy.empirical(table(m_boot_paste[i,]), unit="log2")
}
mean(e_boot) # this should be very close to the analytically calculated bootstrap entropy (boot_h) above
sqrt(var(e_boot)) # and this should give an idea of the variation between bootstrapped versions

#### 2.5 Kullback-Leibler divergence between the group's maximum null (not
#bootstrap) and observed entropy #### index i runs over frequencies with which
#each composition was observed 
###### THERE WAS AN ERROR HERE, added na.rm=T to the sum for kl_div ######
f_comp <- as.array(freqs.empirical(t)) # this will contain the frequencies with which each composition was observed
kl <- array(NA, length(f_comp)) 
for (i in 1:length(kl)){
  k_i <- str_count(names(f_comp[i]), "1") # finds the number of occurrences of "1" in the pasted subgroup compositions
  kl[i] <- f_comp[[i]] * log2((f_comp[[i]]*(choose(n, k_i)))/fk[k_i]) # for each composition, do a "p*log(p/q)" calculation 
}
kl_div <- sum(kl, na.rm=T) # the sum of those "p*log(p/q)" calculations is equal to the KL divergence

#### 2.6 Jensen-Shannon distance between bootstrap and observed entropy #### 
###### THERE WAS AN ERROR HERE, start loop in 2 ######
# (this computation may take longer than most others)
# index j runs over compositions, i runs over number of subgroups observed of the size of the jth composition
js <- array(NA, length(f_comp))
js_l <- vector("list", length(f_comp))
for (j in 2:length(f_comp)){
  k_a <- str_count(names(f_comp[j]), "1") # check that no zero value is possible here
  Nk_a <- Nk[k_a]
  js_l[[j]] <- array(NA, Nk_a)
  for (i in 1:Nk_a){
    js_l[[j]][i] <- (1/(choose(n, k_a))^i)*((1-(1/choose(n, k_a)))^(Nk_a-i))*(choose(Nk_a, i))*
      ((f_comp[j]*log2(1+(i/(No*f_comp[j]))))+((i/No)*log2(1+(No*(f_comp[j]/i)))))
  }
  js[j] <- sum(js_l[[j]], na.rm=T) 
}
js_dist <- 1-(0.5*sum(js, na.rm=T))

#### 3. Entropy from each individual perspective #### 
# Take each observation as a particular combination of 0s, and 1s, calculating 
# the random and empirical information content for the set of associations of
# each individual.

# Create a list for the data set of each individual's associations
m_inds <- vector("list", ncol(m))
for (j in 1:ncol(m)){
  m_inds[[j]] <- m[which(m[,j]==1), -j]
}

head(m_inds[[1]]) # check that each of the ncol(m) individuals has its own data set in the list, changing the number in brackets from 1 to n
nrow(m_inds[[17]])

# Remove observations where each individual is alone (because those will not be
# relevant from the individual perspective)
for (j in 1:length(m_inds)){
  m_inds[[j]] <- m_inds[[j]][which(rowSums(m_inds[[j]])>0),]
}

# Convert each observed subgroup to a single character variable composed of 0 and 1
m_inds_p <- vector("list", length(m_inds))
for (j in 1:length(m_inds)){
  m_inds_p[[j]] <- array(NA, nrow(m_inds[[j]]))
  for (i in 1:nrow(m_inds[[j]])){
    if (nrow(m_inds[[j]])==0)
      m_inds_p[[j]][i] <- 0
    else
      m_inds_p[[j]][i] <- paste(m_inds[[j]][i,], collapse = "")
  }
}

# Retrieve fk_inds, the frequency distributions of subgroups of size k for each individual
sgs_inds <- vector("list", length(m_inds))
n_inds <- array(NA, ncol(m))
Nk_inds <- vector("list", length(m_inds))
fk_inds <- vector("list", length(m_inds))
for (j in 1:length(sgs_inds)){
  sgs_inds[[j]] <- apply(m_inds[[j]], 1, sum)
  n_inds[j] <- ncol(as.data.frame(m_inds[[j]]))
  Nk_inds[[j]] <- hist(sgs_inds[[j]], breaks=c(0:n_inds[j]), plot=FALSE)$counts
  fk_inds[[j]] <- freqs.empirical(Nk_inds[[j]])
}

# Count the number of unique compositions
m_unique <- array(NA, length(m_inds_p))
for (i in 1:length(m_unique)){
  m_unique[i] <- length(unique(m_inds_p[[i]]))
}

# Count the number of occurrences for each unique composition, calculate
# frequencies, check the sum to 1 of each set and calculate empirical entropies
# for each set (according to equation 1 in the article and to the entropy
# package)

t_inds <- vector("list", length(m_inds_p))
f_inds <- vector("list", length(m_inds_p))
check_sum_f_inds <- array(NA, length(m_inds_p))
estimated_h_inds <- array(NA, length(m_inds_p))
check_ent <- array(NA, length(m_inds_p))
set_size <- array(NA, length(m_inds_p))

for (j in 1:length(t_inds)){
  t_inds[[j]] <- table(m_inds_p[[j]])
  f_inds[[j]] <- freqs.empirical(t_inds[[j]])
  check_sum_f_inds[j] <- sum(f_inds[[j]])
  estimated_h_inds[j] <- entropy.empirical(t_inds[[j]], unit="log2") # observed individual entropy values
  check_ent[j] <- -sum(freqs.empirical(t_inds[[j]]) * log2(freqs.empirical(t_inds[[j]]))) # should be equal to the above
  set_size[j] <- nrow(m_inds[[j]])
}  

#### 4. Null bootstrap entropy from the individual perspective #### 
# --> Analytical calculation only, avoiding a bootstrapped set for each individual <---
# Construct nreps bootstrapped versions of original database, maintaining the size of each 
# observed subgroup constant and only switching 0 and 1 between columns

# Add 0 for the largest, possibly unobserved subgroup sizes (depends on an exploration of the histogram of fk above)
Nk[(length(Nk)+1):ncol(m)] <- 0 

No_inds <- array(NA, length(Nk_inds))
for (i in 1:length(Nk_inds)){
  No_inds[i] <- sum(Nk_inds[[i]])
}

ni <- ncol(m)-1
fia <- vector("list", ni)
nk_fi <- matrix(NA, ncol(m), ni) # has to be a matrix
for (j in 1:ncol(m)){
  for (k in 1:ni){
    fia[[k]] <- array(NA, Nk_inds[[j]][k])
    for (i in 1:Nk_inds[[j]][k]){
      fia[[k]][i] <- exp((i*log(1/choose(ni,k)))+((Nk_inds[[j]][k]-i)*log(1-(1/choose(ni,k))))+lchoose(Nk_inds[[j]][k],i))*(i/No_inds[j])*log2(i/No_inds[j])
    }
    nk_fi[j,k] <- choose(ni,k)*sum(fia[[k]])
  }
}

# Calculate the bootstrap entropy for each individual:
boot_h_inds <- array(NA, ncol(m))
for (i in 1:ncol(m)){
  boot_h_inds[i] <- -sum(nk_fi[i,], na.rm=T)
}

#### 4.5 Kullback-Leibler divergence between maximum null (not bootstrap) and observed individual entropy ####

# j runs over individuals, i runs over compositions observed for each individual
f_comp_inds <- f_inds # the frequencies with which each composition was observed, for each individual
kl_div_inds <- array(NA, n)
kl_inds <- vector("list", n)
for (j in 1:ncol(m)){
  for (i in 1:length(f_comp_inds[[j]])){
    k_i_inds <- str_count(names(f_comp_inds[[j]][i]), "1") # finds the number of occurrences of "1" in the pasted subgroup compositions
    kl_inds[[j]][i] <- f_comp_inds[[j]][i] * log2((f_comp_inds[[j]][i]*(choose(ni, k_i_inds)))/fk_inds[[j]][k_i_inds]) # for each composition, do a "p*log(p/q)" calculation 
  }
  kl_div_inds[j] <- sum(kl_inds[[j]]) # the sum of those "p*log(p/q)" calculations is equal to the KL divergence
}

#write.table(kl_div_inds, file="kl_divergence_individuals_geladas2015_onlymales.csv", row.names=FALSE, sep=",", append=TRUE, col.names="geladas_2015")

#### 4.6 Jensen-Shannon distance between bootstrap and observed individual entropy ####
# (this computation may take longer than most others): index h runs
# over individuals, [...] j runs over compositions, i runs over number of
# subgroups observed of the size of the jth composition START HERE
js_inds <- vector("list", n)
js_l_inds <- vector("list", n)
js_dist_inds <- array(NA, n)
for (h in 1:n){
  js_inds[[h]] <- array(NA, length(f_comp_inds[[h]]))
  js_l_inds[[h]] <- vector("list", length(f_comp_inds[[h]]))
  for (j in 1:length(f_comp_inds[[h]])){
    k_a_inds <- str_count(names(f_comp_inds[[h]][j]), "1") # check that no zero value is possible here
    Nk_a_inds <- Nk[k_a_inds]
    for (i in 1:Nk_a_inds){
      js_l_inds[[h]][[j]][i] <- (1/(choose(n, k_a_inds))^i)*((1-(1/choose(n, k_a_inds)))^(Nk_a_inds-i))*(choose(Nk_a_inds, i))*
        ((f_comp_inds[[h]][j]*log2(1+(i/(No*f_comp_inds[[h]][j]))))+((i/No)*log2(1+(No*(f_comp_inds[[h]][j]/i)))))
    }
    js_inds[[h]][j] <- sum(js_l_inds[[h]][[j]], na.rm=T) 
  }
  js_dist_inds[h] <- 1-(0.5*sum(js_inds[[h]], na.rm=T))
}

#### 5. Null model of partition size #### 
# This section estimates the probability of different partition sizes given the
# subgroup size and composition observations.

# First variant of fit (power law)
n <- ncol(m) # group size

# the parameter of the power law in the estimation of partition size
# probability. Here it is necessary to select different ranges for this
# parameter, depending on an exploration of the plots for g vs. dist.ks and
# dist.kl below
g <- seq(-13, -9) # these values should be changed depending on whether the plot below 
# plot(g, dist.kl) shows a minimum or not. 

Kp <- array(NA, length(g)) # constant to normalize the subgroup size distribution
k <- 1:(n-1) # k=n will be estimated separately

fp <- matrix(NA, length(g), length(k))
fdp <- matrix(0, length(g), length(k))

for (i in seq_along(g)){
  for (j in seq_along(k)){
    fdp[i,j] <- choose(n,k[j])*sum(Stirling2.all(n-k[j])*(seq(2,n-k[j]+1)^(g[i]-1)))
  }
  Kp[i] <- 1/(1+sum(fdp[i,])) # here is the constant Kp which is specific to each value of gamma
  for (j in seq_along(k)){
    fp[i,j] <- Kp[i]*choose(n,k[j])*sum(Stirling2.all(n-k[j])*(seq(2,n-k[j]+1)^(g[i]-1)))
  }
}

# Check that the Kp values coincide with an alternative expression (see Denis' notes)
Kp_alt <- array(NA, length(g))
for (i in seq_along(g)){
  Kp_alt[i] <- 1/sum((seq(1,n)^g[i])*Stirling2.all(n)) # in this case we can use a direct method to calculate Stirling numbers of the second order
}

# Add the value of fp for n=k
fp <- cbind(fp, Kp)

rowSums(fp) # check that the sum of rows of f is always 1

# Calculate KL distance and KS test
dist.ks <- array(NA, length(g))
dist.kl <- array(NA, length(g))
for (i in 1:nrow(fp)){
  dist.ks[i] <- ks.test(fk, fp[i,])$statistic # if there is no clear minimum, you will get a warning: "cannot compute exact p-value with ties"
  dist.kl[i] <- KL.plugin(fk, fp[i,]) 
} 
plot(g, dist.ks)
plot(g, dist.kl) # if this plot shows a minimum, continue. Otherwise, modify the values of g above

# Identify minima for g and Kp using the KL distance criterion
gmin <- g[which(dist.kl == min(dist.kl))] 
Kp_min <- Kp[which(g==gmin)]

# Plot the approximation and the empirical data, with the best fit according to
# the KL distance criterion (may need graphical adjustments)
plot(fp[1,], type="n", ylim=c(0,0.4), ylab="f(k)", xlab="k", xlim=c(1,n))
for (i in 1:nrow(fp)){
  lines(fp[i,], col="dark grey")
}
lines(fk, col="red", lwd=1.5)
legend(10, 0.4, c("fk empirical", "fk approximated (power law)", "best fit (KL distance)"), cex=0.9, col=c("red", "dark gray", "black"),
       pch=c("-","-"), bty="n", lwd=c(1.5, 1, 1.5))
lines(fp[which(g == gmin),], lwd=1.5)

# Depending on how the best fit curve fits the observed distribution of subgroup
# size, you may choose to try the second variant or modify the values of the g parameter above

# Second variant of fit: exponential
a <- seq(2, 4, .1) # the parameter of the exponential function in the estimation of partition size probability
# these values should be changed according to the plot below plot(a, dist.kl.e)

Ke <- array(NA, length(a)) # constant to normalize the subgroup size distribution

fe <- matrix(NA, length(a), length(k))
fde <- matrix(0, length(a), length(k))

for (i in seq_along(a)){
  for (j in seq_along(k)){
    fde[i,j] <- choose(n,k[j])*sum(Stirling2.all(n-k[j])*(exp(-a[i]*seq(2,n-k[j]+1))/seq(2,n-k[j]+1)))
  }
  Ke[i] <- 1/(exp(-a[i])+sum(fde[i,])) # here is the constant Ke which is specific to each value of gamma
  for (j in seq_along(k)){
    fe[i,j] <- Ke[i]*choose(n,k[j])*sum(Stirling2.all(n-k[j])*exp(-a[i]*seq(2,n-k[j]+1))/seq(2,n-k[j]+1))
  }
}

# Alternative calculation of K
Ke_alt <- array(NA, length(a))
for (i in seq_along(a)){
  Ke_alt[i] <- 1/sum(exp(-a[i]*seq(1,n))*Stirling2.all(n)) # in this case we can use a direct method to calculate Stirling numbers of the second order
}

# Add last column, corresponding to k=n
fe <- cbind(fe, Ke*exp(-a))

# check that the sum across rows of f is always 1
rowSums(fe) 

# Calculate KL distance and KS test
dist.ks.e <- array(NA, length(a))
dist.kl.e <- array(NA, length(a))
for (i in 1:nrow(fe)){
  dist.ks.e[i] <- ks.test(fk, fe[i,])$statistic # you could get a warning: "cannot compute exact p-value with ties"
  dist.kl.e[i] <- KL.plugin(fk, fe[i,]) # this test might not be working well for this particular case
} 
plot(a, dist.ks.e)
plot(a, dist.kl.e) # if this plot shows a minimum, continue. Otherwise, modify the value of a above.

# Identify minima for a and Ke using the KL distance criterion
amin <- a[which(dist.kl.e == min(dist.kl.e))] # this minimum value of a should show in the plots above as a clear minima
Ke_min <- Ke[which(a==amin)] # as well

### Plot the approximation with the empirical and the best fit according to the
### KL distance criterion (may need graphical adjustments)
plot(fe[1,], type="n", ylim=c(0,0.4), ylab="f(k)", xlab="k", xlim=c(0,n))
for (i in 1:nrow(fe)){
  lines(fe[i,], col="dark grey")
}
lines(fk, col="red", lwd=1.5)
legend(10, 0.4, c("fk empirical", "fk approximated (exponential)", "best fit (KL distance)"), cex=0.9, col=c("red", "dark gray", "black"),
       pch=c("-","-"), bty="n", lwd=c(1.5, 1, 1.5))
lines(fe[which(a == amin),], lwd=1.5)


# Finally, use the estimated values for K to estimate a probability distribution of partition size
# For the power law fit
Pm_p <- Kp_min*(seq(1:n)^gmin) # this is the probability of finding a given partition, simply based on the frequency of observations of subgroups of different size
Nm_p <- Pm_p*Stirling2.all(n) # But this considers the possible numbers of partitions of each size m (Stirling numbers of the second kind)
sum(Nm_p) # check that it adds to 1

# For the exponential fit
Pm_e <- Ke_min*exp(-amin*seq(1:n))
Nm_e <- Pm_e * Stirling2.all(n)
sum(Nm_e) # should add to 1

# Plot the probability N(m) of finding a partition of size m (may need graphical adjustments)
plot(seq(1:n), Nm_p, ylab="Frequency of partitions of size m", xlab="m", ylim=c(0,1), type="n")
lines(seq(1:n), Nm_e, type="b", col="blue", pch=1)
lines(seq(1:n), Nm_p, type="b", col="red", pch=2)
legend(n-5, 0.9, c("power-law fit", "exponential fit"), col=c("red", "blue"), pch=c("-","-"), bty="n", cex=0.9)
