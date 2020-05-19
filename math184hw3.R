#Initial set up
setwd("Downloads")
data <- read.table("hcmv.txt", header=TRUE)

#Generating 5 random samples with 296 locations
set.seed(22343)
n <- 296
for ( i in 1:5) {
sample[[i]] <- runif(n, min = 0, max = 229354)
}
#Sorting the sample (because the palindrome locations in the DNA are sorted)
for (i in 1:5) {
  sample[[i]] = sort(sample[[i]])
}

#Comparing DNA vs Random Sample distribution using histograms
par(mfrow = c(2,3))
for ( i in 1:5) {
hist(sample[[i]], breaks = 50, probability = TRUE, col = i+1, main=paste("Random sample", i, sep=" "))
}
hist(data$location, breaks = 50, probability = TRUE, col = 4, main = "CMV DNA palindrome locations")


#Calculate the means of the 5 samples to compare spacing
sampleMean <- rowMeans(cbind(sample[[1]],sample[[2]],sample[[3]],sample[[4]],sample[[5]]))
samplegap = c()
samplegap2 = c()
samplegap3 = c()
#Consecutive
gap = c()
for(i in 1:nrow(data))
{
  gap[i] = data$location[i+1]-data$location[i]
}
for(i in 1:length(sampleMean)-1)
{
samplegap[i] = sampleMean[i+1]-sampleMean[i]
}

#Consecutive w/ 1 in b/w
gap2 = c()
for(i in 1:nrow(data))
{
  gap2[i] = data$location[i+2]-data$location[i]
}
for(i in 1:length(sampleMean)-2)
{
samplegap2[i] = sampleMean[i+2]-sampleMean[i]
}

#Consecutive w/ 2 in b/w
gap3 = c()
for(i in 1:nrow(data))
{
  gap3[i] = data$location[i+3]-data$location[i]
}
for(i in 1:length(sampleMean)-3)
{
samplegap3[i] = sampleMean[i+3]-sampleMean[i]
}

#Visualizing differences in spacing
par(mfrow=c(2,3))
hist(gap, breaks = 30, col = 2, main = "DNA - Gap b/w consecutive palindromes")
hist(gap2, breaks = 30, col = 2, main = "DNA - Gap b/w consecutive palindromes w/ 1 in between")
hist(gap3, breaks = 30, col = 2, main = "DNA - Gap b/w consecutive palindromes w/ 2 in between")
hist(samplegap, breaks = 30, col = 4, main = "Simulated data - Gap b/w consecutive palindromes")
hist(samplegap2, breaks = 30, col = 4, main = "Simulated data - Gap b/w consecutive palindromes w/ 1 in between")
hist(samplegap3, breaks = 30, col = 4, main = "Simulated data - Gap b/w consecutive palindromes w/ 2 in between")
par(mfrow=c(3,1))
plot(gap, main="DNA : Spacing between consecutive palindromes",xlab="Palindrome #", ylab="Spacing", pch=19, col=1)
points(samplegap, pch=19, col=2)
legend( "topright", c("DNA","Random sample"), text.col=c(1,2) )
plot(gap2, main="DNA : Spacing between 2 palindromes w/ 1 in between",xlab="Palindrome #", ylab="Spacing", pch=19, col=1)
points(samplegap2,pch=19, col=2)
legend( "topright", c("DNA","Random sample"), text.col=c(1,2) )
plot(gap3, main="DNA : Spacing between 2 palindromes w/ 2 in between",xlab="Palindrome #", ylab="Spacing", pch=19, col=1)
points(samplegap3, pch=19, col=2)
legend( "topright", c("DNA","Random sample"), text.col=c(1,2) )
tab <- table(cut(data$location, breaks = seq(0,230000, length.out = 57), include.lowest = TRUE))
int = table(count)
count = as.vector(tab)
chisq.test(int,p=Pois)
Pois[8] = 1-sum(Pois)
Pois <- dpois(0:6, lambda = mean(count))

tab <- table(cut(data$location, breaks = seq(92400,92850, length.out = 10), include.lowest = TRUE))
tab <- table(cut(data$location, breaks = seq(92500,92850, length.out = 10), include.lowest = TRUE))
tab <- table(cut(data$location, breaks = seq(92500,92800, length.out = 10), include.lowest = TRUE))
tab <- table(cut(data$location, breaks = seq(92600,93000, length.out = 10), include.lowest = TRUE))
tab <- table(cut(data$location, breaks = seq(92600,93000, length.out = 57), include.lowest = TRUE))
tab <- table(cut(data$location, breaks = seq(92600,93000, length.out = 57), include.lowest = TRUE))
tab <- table(cut(data$location, breaks = seq(0,230000, length.out = 57), include.lowest = TRUE))
tab2 <- table(cut(data$location, breaks = seq(193000,196000, length.out = 10), include.lowest = TRUE))
tab2 <- table(cut(data$location, breaks = seq(193000,196000, length.out = 6), include.lowest = TRUE))
tab2 <- table(cut(data$location, breaks = seq(195000,195300, length.out = 6), include.lowest = TRUE))
tab2 <- table(cut(data$location, breaks = seq(19500,195300, length.out = 6), include.lowest = TRUE))
tab2 <- table(cut(data$location, breaks = seq(194800,195400, length.out = 6), include.lowest = TRUE))
tab2 <- table(cut(data$location, breaks = seq(193000,196000, length.out = 6), include.lowest = TRUE))
pois230 = rpois(230,lambda = mean(count230))
pois150 = rpois(150,lambda = mean(count150))
pois100 = rpois(100,lambda = mean(count100))
poisValues = rpois(57,lambda = mean(counts))
count230 = as.vector(tab230)
tab230 <- table(cut(data$location, breaks = seq(0, 230000, length.out = k+1), include.lowest = TRUE))
count150 = as.vector(tab150)
tab150 <- table(cut(data$location, breaks = seq(0, 230000, length.out = k+1), include.lowest = TRUE))
count100 = as.vector(tab100)
tab100 <- table(cut(data$location, breaks = seq(0, 230000, length.out = k+1), include.lowest = TRUE))
hist(counts100, breaks = 10, col = rgb(1,0,0,0.5), probability = TRUE, xlab = "number of points inside an interval", ylim = c(0,0.2), main="Palindrome counts in intervals of size 2300")
hist(count150, breaks = 10, col = rgb(1,0,0,0.5), probability = TRUE, xlab = "number of points inside an interval", ylim = c(0,0.5), main="Palindrome counts in intervals of size 1533")
hist(count230, breaks = 10, col = rgb(1,0,0,0.5), probability = TRUE, xlab = "number of points inside an interval", ylim = c(0,0.8), main="Palindrome counts in intervals of size 1000")
lines(density(count230, adjust=2), col = rgb(1,0,0,1))
lines(density(pois230, adjust=2), col = rgb(0,0,0,1))
lines(density(count150, adjust=2), col = rgb(1,0,0,1))
lines(density(pois150, adjust=2), col = rgb(0,0,0,1))
lines(density(counts100, adjust=2), col = rgb(1,0,0,1))
lines(density(pois100, adjust=2), col = rgb(0,0,0,1))
barplot(tab2, col=(rgb(0,0,1,0.5)), main="Palindrome counts in the interval (193,000, 196,000)")
barplot(tab3, col=(rgb(0,1,0,0.5)), main="Palindrome counts in the interval (184,000, 199,000)")
barplot(tab1, col=(rgb(1,0,0,0.5)), main="Palindrome counts in the interval (0, 230,000)")

