# first observe our data set
# must have hcmv.txt in current directory
data <- read.table("hcmv.txt", header=TRUE) 
summary(data)
hist(data$location, breaks = 50, xlab="Location", main="Palindromes Found on CMV DNA Strand")

#Comparing DNA to uniform distribution using lattice
N <- 229354
n <- 296
set.seed(21)
gene <- seq(1, N)
site.random <- sample.int(N, size=n)
library(lattice)
stripplot(site.random, pch=16, cex=0.25, xlab='Location', main = 'Uniform Distribution')
stripplot(data$location, pch=16, cex=0.25, xlab='Location', main ='CMV DNA')

#generating table of 100 uniform distributions 
normal <- function(x,y){sample.int(x, size=y)}
table1 <- replicate(100, normal(N,n))

#Method to split the regions
regionsplit <- function(n.region, gene, site){
  count.int <- table(cut(site, breaks = seq(1, length(gene), length.out=n.region+1), include.lowest=TRUE))
  count.vector <- as.vector(count.int)
  count.tab <- table(count.vector)
  return (count.tab)
}

#Observe counts and histograms of some of the simulated uniform distributions w/
#50 split regions
n.region <- 50
regionsplit(n.region, gene, table1[1,])
hist(table1[1,], breaks=n.region)
regionsplit(n.region, gene, table1[2,])
hist(table1[2,], breaks=n.region)

#Comparing max # of palindromes in intervals for 100 uniform distributions and our data
#break size is 50
n.region <- 50
avheight <- 0
for(i in 1:ncol(table1)-1){avheight<- avheight +length(regionsplit(n.region, gene, table1[i,])) }
avheight <- avheight/100 #average height of largest region for uniform distributions
avheight
length(regionsplit(n.region, gene, data$location)) #height of largest region for data

#Aggregating the Table
numBins <- 50
binSize <- N/numBins
counter <- 0

aggArray <- c(1:numBins)
aggArray[] <- 0

for (c in 1:100) {
  for (r in 1: 296) {
    val <- table1[r,c]
    counter <- 1
    
    while (val - binSize >= 0) {
      val <- val - binSize
      counter <- counter + 1
    }
    
    aggArray[counter] <- aggArray[counter] + 1
  }
}
aggArray[] <- aggArray[]/100
plot(aggArray, xlab="Bin[i], with each bin sized 
     9174.16 pairs", ylab="Number hits", main="100 
     Uniform Scatters aggregated by bins")


#Chi-Squared Table Function
#50P (k palindromes in an interval of ~length 4587) = 50e^−λ[1+λ+λ^2/2..+λ^k/k]
chisqtable <- function(n.region, site, N){
  n <- length(site)
  lambda.est <- n/n.region
  count.int <- table(cut(site, breaks = seq(1, length(gene), length.out=n.region+1), include.lowest=TRUE))
  count.vector <- as.vector(count.int)
  count.range <- max(count.vector) - min(count.vector) + 1
  
  table <- matrix(rep(NA, count.range*3), count.range, 3)
  for (i in 1:count.range){
    offset <- min(count.vector) - 1
    table[i, 1] <- i + offset
    table[i, 2] <- sum(count.vector == i + offset)
    if ((i + offset == min(count.vector)) && (min(count.vector) != 0))
      table[i, 3] <- ppois(i+offset, lambda.est)*n.region
    else if (i + offset == max(count.vector))
      table[i, 3] <- 1 - ppois(i + offset - 1, lambda.est)
    else
      table[i, 3] <- (ppois(i+offset, lambda.est) - ppois(i + offset - 1, lambda.est))*n.region
  }
  
  return (table)
}
#Chi-Squared Test, DNA vs Uniform
counts <- table(cut(data$location, breaks = seq(1, length(gene), length.out=50+1), include.lowest=TRUE))
counts <- as.vector(counts)
E_i <- 296/50
chi_2 <- sum((counts - E_i)^2/E_i)
chi_2 # = 66.5
chi2_compare <- qchisq(p = 0.95, df = 49)
chi2_compare
p_value <- 1 - pchisq(chi_2, df = 49)
p_value


#Chi-Squared Test, DNA vs Poisson, determining counts
hist(counts, breaks = 15, col = rgb(1,0,0,0.5), probability = TRUE, xlab = "number of points inside an interval", ylim = c(0,0.2))
Pois <- rpois(294, lambda = mean(counts))
hist(Pois, breaks = 50, col = rgb(0,0,1,0.5), probability = TRUE, add = TRUE)
lines(density(counts, adjust = 2), col = rgb(1,0,0,0.5))
lines(density(Pois, adjust = 2), col = rgb(0,0,1,0.5))
legend(x = 10, y = 0.15, legend = c("sample", "Poisson"), lty = c(1,1), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
site.data.tabtemp <- chisqtable(50, data$location, N)
site.data.tabtemp
site.data.tab <- matrix(rep(NA, 7*2), 7, 2)
site.data.tab[1,] <- colSums(site.data.tabtemp[1:4, 2:3])
site.data.tab[2:6,] <- site.data.tabtemp[5:9, 2:3]
site.data.tab[7,] <- colSums(site.data.tabtemp[10:16, 2:3])
site.data.tab
site.data.stats <- sum((site.data.tab[,2] - site.data.tab[,1])^2/site.data.tab[,2])
pchisq(site.data.stats, 7 - 2, lower.tail=FALSE) # == 0.7981743
#Keep in mind that the null hypothesis of the testing procedures above assumes the site locations are formed by a random scatter. 
#Thus, the conclusions show that with a significance level of >0.05, we accept the null hypothesis.

#Standard Residual Plot
Residuals <- (counts - E_i) / sqrt(E_i)
plot(Residuals, type = 'h', ylab = "standardized residuals", xlab = "interval index")

#Initialize Data for Next Part
mydata=read.table("hcmv.txt",header=TRUE)
head(mydata)
#K-means clustering

# Prepare Data
mydata <- na.omit(mydata) # listwise deletion of missing
mydata <- scale(mydata) # standardize variables

# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 5) # 5 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)

# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")

# Model Based Clustering
install.packages("mclust")
library(mclust)
fit <- Mclust(mydata)
plot(fit) # plot results 
summary(fit) # display the best model

# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 5)

# Cluster Plot against 1st 2 principal components

# Vary parameters for most readable graph
installed.packages("cluster")
library(cluster) 
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
install.packages("fpc")
library(fpc)
library(cluster)
plotcluster(mydata, fit$cluster)

#Explore Function
explore <- function(x,interval,threshold,print=T,samples=1,par=1){
  #bins <- seq(1,(interval)*ceiling(229354/(interval))+1,interval)
  bins <- c(((0: (floor(229354/interval)-1))*interval)+1,229354)
  histinfo <- hist(x,breaks=bins,plot=F)
  if(par==2){
    info <- hist(histinfo$counts,breaks=c(0:15))
    barplot(height = (info$counts)/samples,width = 1,names.arg = c(0:14),col = "lightgreen",las=2,ylim=c(0, max(info$counts/samples)),ylab = "# of intervals",main = paste("# of intervals vs palindrome count (interval length =", interval,")"))
    return (info)
  }
  if(par==1){
    counts_d <- sort(histinfo$counts,decreasing = T)
    counts_i <- unique(counts_d[which(counts_d>threshold)])
    barplot(height = (histinfo$counts)/samples,width = diff(bins),names.arg = histinfo$mids,col = "lightgreen",las=2,ylim=c(0, max(histinfo$counts/samples)),ylab = "# of palindromes encountered",main = paste("# of palindromes vs location of occurence (interval length =", interval,")"))
    if(print){
      print("Midpoints of intervals with maximum number of palindromes clustered are:")
      for (j in counts_i/samples){
        print(paste(histinfo$mids[(which(histinfo$counts/samples==j))],":",j,"palindromes"))
      }
    }
    l <- list("hi"=histinfo, "ci"= counts_i/samples, "c"=histinfo$counts/samples)
    return(l)
  }
  
  
}



interval1 = 4000
threshold1 = 6
bins <- c(((0: (floor(229354/interval1)-1))*interval1)+1,229354)
ret1 <- explore(dat$location,interval1,threshold1,print=T)
#Intervals with 0-14 palindromes
ret2 <- explore(dat$location,200,par=2)



####################################      PART A2      #################################
# generating a simulated sample of 296 palindrom locations from 229354 possible sites
dat_simulated1 <- sample(1:229354,296*15,replace = T)
ret21<- explore(dat_simulated1,interval1,threshold1,print=T,samples = 15)
#Intervals with 0-14 palindromes
ret22<- explore(dat_simulated1,200,par=2)


########################    Uniform Distribution   ####################
dat<-read.table("hcmv.txt",header=TRUE)
uniform_bins <- round(0:10*229354/10)
uniform_bins2 <- round(0:20*229354/20)
hist_uni <- hist(dat$location,breaks=uniform_bins,col="gray", freq = T, xlab="No. of Palindromes",ylab="frequency",main="Histogram for Uniform distribution analysis ")
hist_uni2 <- hist(dat$location,breaks=uniform_bins2,col="gray", freq = F, xlab="No. of Palindromes",ylab="frequency",main="Histogram for Uniform distribution analysis")
p_uni <- rep(1/length(hist_uni$counts),length(hist_uni$counts))
p_uni2 <- rep(1/length(hist_uni2$counts),length(hist_uni2$counts))
chi_uni <- chisq.test(hist_uni$counts,p=p_uni)
chi_uni2 <- chisq.test(hist_uni2$counts,p=p_uni2)
library(ggplot2)
ggplot(data.frame(cbind(hist_uni$stdres,hist_uni$counts)),aes(x=1:length(hist_uni$counts),y=0))+geom_errorbar(aes(ymin = chi_uni$stdres,ymax=0))+labs(x="Levels of No. of Palindromes",y="Standardized residual",title="Standardized residuals for the 10 different intervals")

print("Uniform Distribution for 10 intervals Goodness-of-Fit test:")
print(paste("The value of chi-squared statistic is",chi_uni$statistic))
print(paste("The p-value is",chi_uni$p.value))
cat("The standard residuals are as follows: \n", chi_uni$residuals)

ggplot(data.frame(cbind(hist_uni2$stdres,hist_uni2$counts)),aes(x=1:length(hist_uni2$counts),y=0))+geom_errorbar(aes(ymin = chi_uni2$stdres,ymax=0))+labs(x="Levels of No. of Palindromes",y="Standardized residual",title="Standardized residuals for the 20 different intervals")
print("Uniform Distribution for 20 intervals Goodness-of-Fit test:")
print(paste("The value of chi-squared statistic is",chi_uni2$statistic))
print(paste("The p-value is",chi_uni2$p.value))
cat("The standard residuals are as follows: \n", chi_uni2$residuals)

######################## Poisson Distribution for Actual data  #################

poisson_bins <- c(0,2,3,4,5,6,7,8,max(ret1$c))
lambda = 296/(length(poisson_bins)-1)
hist_p <- hist(ret1$c,breaks=poisson_bins,col="gray", freq = T, xlab="No. of Palindromes",ylab="frequency",main="Histogram for Poisson distribution analysis")
p_poisson <- diff(ppois(poisson_bins,lambda))
chi_ppoisson <- (p_poisson+((1-sum(p_poisson))/length(p_poisson)))
chi_poisson <- chisq.test(hist_p$counts,p=chi_ppoisson)

print("Poisson Goodness-of-Fit test:")
print(paste("The value of chi-squared statistic for is",chi_poisson$statistic))
print(paste("The p-value is",chi_poisson$p.value))

###############  Poisson Distribution for Simulated Random data  ##################

hist_pr <- hist(ret21$c,breaks=poisson_bins,col="gray", freq = T, xlab="No. of Palindromes",ylab="frequency",main="Histogram for Poisson distribution analysis of simulated data")
chi_poissonr <- chisq.test(hist_pr$counts,p=chi_ppoisson)
print("Poisson Goodness-of-Fit test on Simulated Random Data:")
print(paste("The value of chi-squared statistic for is",chi_poissonr$statistic))
print(paste("The p-value is",chi_poissonr$p.value))

########################  Exponential Distribution   #################

plot(sort(diff(dat$location),decreasing = TRUE),xlab = "index of palindromes",ylab = "Difference between palindromes",main ="Sorted distances among adjacent Palindromes")
exp_bins <- c(0,1000,1500,2000,2500,3000,3500,4000,5500)
hist_exp <- hist(diff(dat$location),breaks=exp_bins,col="gray", freq = T, xlab="No. of Palindromes",ylab="frequency",main="Histogram for Exponential distribution analysis")
p_exp <- diff(pexp(exp_bins,rate=295/229354))
chi_pexp <- ((1-sum(p_exp))/length(p_exp))+p_exp
chi_exp <- chisq.test(hist_exp$counts,p=chi_pexp)
print("Exponential distribution Goodness-of-Fit test:")
print(paste("The value of chi-squared statistic for is",chi_exp$statistic))
print(paste("The p-value is",chi_exp$p.value))

########################  Gamma Distribution   #################

gamma_bins <- c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,6000)
hist_gamma <- hist(diff(dat$location,lag=2),breaks=gamma_bins,col="gray", freq = T, xlab="No. of Palindromes",ylab="frequency",main="Histogram for Gamma distribution analysis")
p_gamma <- diff(pgamma(gamma_bins,2,rate=229354/295))
chi_pgamma <- ((1-sum(p_gamma))/length(p_gamma))+p_gamma
chi_gamma <- chisq.test(hist_gamma$counts,p=chi_pgamma)
print("Gamma distribution Goodness-of-Fit test:")
print(paste("The value of chi-squared statistic for is",chi_gamma$statistic))
print(paste("The p-value is",chi_gamma$p.value))


