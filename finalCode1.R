# Filename: finalCode.R
# Name: Himanshu Makhija - A09845605, Math/CS, B.S.
#       Hanna Goldman - A11436767, Math/CS, B.S.
#       Amey Paranjape - A53218045, ECE, M.S.
#       Shagun Gupta - A91068956, BENG-Bioinformatics, B.S.
# Date: 1/24/2016
# Description: This file contains all the relevant code needed to plot any charts, graphs, or tables generated in
#       our text file. Before running the script, PLEASE HAVE BABIES.TXT in your working directory. 

data <- read.table("babies.txt", header=TRUE)

# Cleaning the data in "babies.txt"
babies <- babies[babies$gestation != 999,]
babies <- babies[babies$age != 99,]
babies <- babies[babies$height != 99,]
babies <- babies[babies$height != 999,]
babies <- babies[babies$weight != 999,]
babies <- babies[babies$smoke != 9,]

# Sort the data
smoker.ind <- which(babies['smoke'] == 1) 
nonsmoker.ind <- which(babies['smoke'] == 0)
babies.smoker <- babies[smoker.ind,]
babies.nonsmoker <- babies[nonsmoker.ind,]

# BOX PLOTS OF UNCLEANED/CLEANED
# Unclean
uncleanBabies <- read.table("babies.txt", header=TRUE)
head(uncleanBabies)
uncleanBwt = uncleanBabies$bwt
uncleanSmoke = uncleanBabies$smoke
boxplot(uncleanBwt~uncleanSmoke, uncleanBabies, main = "Boxplot of unfiltered data", xlab="smoking status", 
        ylab="birthweight in ounces")

#Clean
cleanBwt = babies$bwt
cleanSmoke = babies$smoke
boxplot(cleanBwt~cleanSmoke, babies, main = "Boxplot of cleaned data", xlab="smoking status", 
        ylab="birthweight in ounces")

# --------------

# Presenting summary of dataset in smokers and nonsmokers
summary(babies.nonsmoker)
summary(babies.smoker)

# --------------

# Scatter plot: Gestational age and Birth Weight
plot(babies$gestation, babies$bwt, main ="Relationship Between Gestational Age and Babies Birth Weight", 
     xlab="Gestational Days", ylab="Birth Weight(oz)")
reg1 <- lm(babies$bwt~babies$gestation) 
abline(reg1, col='red')

# --------------

# Scatter plot: BMI and Birth Weight
babies$bmi <-((babies$weight)*.45)/((babies$height)*.025)^2

plot(babies$bmi, babies$bwt, xlab= "Mother's BMI", ylab="Baby's Birth Weight", 
     main="Effect of Mothers BMI on Babies Birth Weight")
abline(h=88, col="blue")
abline(h=141, col="blue")
abline(v=19, col="red")
abline(v=25, col="red")
legend(x = "topright", c("Range of Healthy BMI's","Range of Healthy Birth Weights"),col = c("red","blue"),
       lwd = c(2,2,2),cex=0.5)

# --------------

# Plots: Parity and Birth weight
# Scatter for nonsmokers
al <- babies[babies$parity == 0,]
bw = al$bwt
plot(density(bw), type="l", main="Plot of bwt frequency for parity = 0", xlab="Weight in ounces ", ylab="Frequency ")
abline(a = NULL, b = NULL, h = NULL, v = 88, col="blue") #Limit line
abline(a = NULL, b = NULL, h = NULL, v = 141, col="red") #Limit line

#Scatter for smokers
an <- babies[babies$parity == 1,]
bt = an$bwt
plot(density(bt), type="l", main="Plot of bwt frequency for parity = 1", xlab="Weight in ounces ", ylab="Frequency ")
abline(a = NULL, b = NULL, h = NULL, v = 88, col="blue") #Limit line
abline(a = NULL, b = NULL, h = NULL, v = 141, col="red") #Limit line

#Q-Q Plot: smokers vs. nonsmokers
qqplot(bt, bw, xlab="parity=1", ylab="parity=0", main = "q-q plot of bwt for parity=1 vs parity=0") 
abline(c(0,1))

# --------------

# Q-Q Plots: Observing the Distribution of the Birthweights
qqnorm(babies$bwt)
qqline(babies$bwt, col='red')
qqnorm(babies.smoker$bwt)
qqline(babies.smoker$bwt, col='red')
qqnorm(babies.nonsmoker$bwt)
qqline(babies.nonsmoker$bwt, col='red')

# --------------

# Histograms: Birth Weights in Smokers and Nonsmokers
# Nonsmokers
h<-hist(babies.nonsmoker$bwt,breaks=20, col ="grey", xlab="Birth Weight [oz]", 
        main="Frequency of Birth Weights of Babies of Nonsmokers", xlim = c(50,180))
x<-babies.nonsmoker$bwt
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)
abline(v=median(babies.nonsmoker$bwt), col="red")
abline(v=88, col="darkgreen")
legend(x = "topright", c("Normal Curve","Median","Small"),col = c("blue","red","darkgreen"),lwd = c(2,2,2),cex=0.5)

# Smokers
h<-hist(babies.smoker$bwt,breaks=15, col ="grey", xlab="Birth Weight [oz]", 
        main="Frequency of Birth Weights of Babies of Smokers")
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
abline(v=median(babies.smoker$bwt), col="red")
abline(v=88, col="darkgreen")
legend(x = "topright", c("Normal Curve","Median","Small"),col = c("blue","red","darkgreen"),lwd = c(2,2,2),cex=0.5)

# ---------------

# Histograms: Effects of Smoking on Gestation Period
# Smokers
h<-hist(babies.smoker$bwt,breaks=15, col ="grey", xlab="Birth Weight [oz]", 
        main="Frequency of Birth Weights of Babies of Smokers")
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
abline(v=median(babies.smoker$bwt), col="red")
abline(v=88, col="darkgreen")
legend(x = "topright", c("Normal Curve","Median","Small"),col = c("blue","red","darkgreen"),lwd = c(2,2,2),cex=0.5)

# Nonsmokers
h<-hist(babies.nonsmoker$bwt,breaks=20, col ="grey", xlab="Birth Weight [oz]", 
        main="Frequency of Birth Weights of Babies of Nonsmokers", xlim = c(50,180))
x<-babies.nonsmoker$bwt
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)
abline(v=median(babies.nonsmoker$bwt), col="red")
abline(v=88, col="darkgreen")
legend(x = "topright", c("Normal Curve","Median","Small"),col = c("blue","red","darkgreen"),lwd = c(2,2,2),cex=0.5)

# ---------------

# MCS: Kurtosis Theoretical vs. Data Sample Kurtosis
# Script: monteCarlo_cs1
# Name: Himanshu Makhija
# Date: 1/22/2017
# Description: This script will generate a number of n sized samples and plot their kurtosesis with a box plot
# or, if you uncomment, with a histogram as well.

### NonSmokers
# Include the moments calculating library
library(moments)

# Initialize index, and the matrix to store our many samples
sampleKurtosis <- matrix(0, 1000, 1)
index <- 1
rnormVec <- 0

# Generate n sized samples, m times, each time storing the resulting kurtosis in the matrix
while(index < 1001) {
  rnormVec <- rnorm(715)
  
  randomSampleKurtosis <- kurtosis(rnormVec)
  
  sampleKurtosis[index, 1] <- randomSampleKurtosis

  index <- index + 1
}

# Plot to show distribution
hist(sampleKurtosis, 30, , , , , 10, , , , "MCS of Nonsmoker Kurtosis Distribution", , , 
     "Kurtosis")
abline(v=3, col="red")
abline(v=4.013421, col="blue")

## Smokers
library(moments)

# Initialize index, and the matrix to store our many samples
sampleKurtosis <- matrix(0, 1000, 1)
index <- 1
rnormVec <- 0

# Generate n sized samples, m times, each time storing the resulting kurtosis in the matrix
while(index < 1001) {
  rnormVec <- rnorm(459)
  
  randomSampleKurtosis <- kurtosis(rnormVec)
  
  sampleKurtosis[index, 1] <- randomSampleKurtosis
  
  index <- index + 1
}

# Plot to show distribution
hist(sampleKurtosis, 30, , , , , 10, , , , "MCS of Smoker Kurtosis Distribution", , , 
     "Kurtosis")
abline(v=3, col="red")
abline(v=2.947538, col="blue")

