# Filename: assignment2
# Name: Himanshu Makhija - A09845605, Math/CS, B.S.
#       Hanna Goldman - A11436767, Math/CS, B.S.
#       Amey Paranjape - A53218045, ECE, M.S.
#       Shagun Gupta - A91068956, BENG-Bioinformatics, B.S.
# Date: 2/07/2016
# Description: This file contains all the relevant code needed to plot any charts, graphs, or tables generated in
#       our text file. Before running the script, PLEASE HAVE videodata.txt in your working directory. 
install.packages("moments")
library(moments)

#Initialize the data
data <- read.table("videodata.txt", header=TRUE)
# Clean the data
data[data == 99] <- NA
sum(is.na(data))

#Scenario 1
boot.population <- rep(data$time, length.out = 314)
set.seed(12345)
B = 500
boot.sample <- array(dim = c(B, 91))
for (i in 1:B) {boot.sample[i, ] <- sample(boot.population, size = 91, replace = FALSE)}
boot.playnumber <- apply(X = boot.sample, MARGIN = 1, FUN = function(x) sum(x !=0))
boot.playpercentage <- boot.playnumber/91
hist(boot.playpercentage, breaks = 20, probability = TRUE, density = 20, col = 3, border = 3, xlab = "Percentage of People", main="")
lines(density(boot.playpercentage, adjust = 2), col = 2)
par(pty = 's')
qqnorm(boot.playpercentage)
qqline(boot.playpercentage)
skewness(boot.playpercentage)
kurtosis(boot.playpercentage)
n.sample <- rnorm(n=500)
skewness(n.sample)
kurtosis(n.sample)
boot.sd <- sd(boot.playpercentage)
mean.boot <- mean(boot.playpercentage)
mean.boot
int.boot <- c(quantile(boot.playpercentage, 0.025), quantile(boot.playpercentage, 0.975))
int.boot

#Scenario 2
p.daily <- data[which(data$freq==1),]
p.weekly <- data[which(data$freq==2),]
p.monthly <- data[which(data$freq==3),]
p.semesterly <- data[which(data$freq==4),]

summary(p.daily$time)
summary(p.weekly$time)
summary(p.monthly$time)
summary(p.semesterly$time)

sd(p.daily$time)
sd(p.weekly$time)
sd(p.monthly$time)
sd(p.semesterly$time)

#Scenario 3
boot.population <- rep(data$time, length.out = 314)
set.seed(12345)
B = 500
boot.sample <- array(dim = c(B, 91))
for (i in 1:B) { boot.sample[i, ] <- sample(boot.population, size = 91, replace = FALSE) }
boot.mean <- apply(X = boot.sample, MARGIN = 1, FUN = mean) 
hist(boot.mean, breaks = 20, probability = TRUE, density = 20, col = 3, border = 3, xlab="Average Play Time [hrs]", main="Average Play Time for Bootstrap Samples")
lines(density(boot.mean, adjust = 2), col = 2)
par(pty = 's')
qqnorm(boot.mean)
qqline(boot.mean)
shapiro.test(boot.mean)
skewness(boot.mean)
kurtosis(boot.mean)
n.sample <- rnorm(n=500)
skewness(n.sample)
kurtosis(n.sample)
mean.boot <- mean(boot.sample)
mean.boot
int.boot <- c(quantile(boot.mean, 0.025), quantile(boot.mean, 0.975))
int.boot

#Scenario 4
data <- read.table("videodata.txt", header=TRUE)

#Cleaning Data
data <- data[data$like != 99,]
data <- data[data$time != 99,]
data <- data[data$where != 99,]
data <- data[data$freq != 99,]
data <- data[data$busy != 99,]
data <- data[data$educ != 99,]
data <- data[data$sex != 99,]
data <- data[data$age != 99,]
data <- data[data$home != 99,]
data <- data[data$math != 99,]
data <- data[data$work != 99,]
data <- data[data$cdrom != 99,]
data <- data[data$email != 99,]
data <- data[data$grade != 99,]

data['dis_like'] <- rep(NA, dim(data)[1])
for(i in 1:dim(data)[1]){
  like <- data[i, 'like']
  if(like==0 || like==4 || like==5){
    data[i, 'dis_like'] = 0
  }else{
    data[i, 'dis_like'] = 1
  }
}


#Comparison in tree package
install.packages("tree")
library(tree)
data.tree <- tree(dis_like~educ+sex+age+home+math+work+own+cdrom+grade, data=data)
plot(data.tree, type="uniform")
text(data.tree)


# Regression Tree Example
install.packages("rpart")
library(rpart)

# grow tree for like
fit <- rpart(dis_like~educ+sex+age+home+math+work+own+cdrom+grade, data=data)

printcp(fit) # display the results 
plotcp(fit) # visualize cross-validation results 
summary(fit) # detailed summary of splits

# create additional plots 
par(mfrow=c(1,2)) # two plots on one page 
rsq.rpart(fit) # visualize cross-validation results  	

# plot tree 
plot(fit, uniform=TRUE, 
     main="Regression Tree for Like ")
text(fit, use.n=TRUE, all=TRUE, cex=.8)
     
#Scenario 5
data <- read.table("videodata.txt", header=TRUE)
dat1p = na.omit(data)
dat1p[which(dat1p$like==1 | dat1p$like==4 | dat1p$like==5),]$like =1
dat1p[which(dat1p$like==2 | dat1p$like==3),]$like =0

m_nl<- dat1p$sex[which(dat1p$sex==1 & dat1p$like==1)]
m_l <- dat1p$sex[which(dat1p$sex==1 & dat1p$like==0)]
f_nl <- dat1p$sex[which(dat1p$sex==0 & dat1p$like==1)]
f_l <- dat1p$sex[which(dat1p$sex==0 & dat1p$like==0)]
f<- cbind(v1=c(length(m_l),length(f_l)),v2=c(length(m_nl),length(f_nl)))
barplot(f, main="Group Barplots: Gender comparison", ylab= "Count",names.arg=c("Like","Don't Like"),beside=TRUE, col=c("green","darkblue"))
legend("topright", fill= c("green","darkblue"), legend=c("Male", "Female")  )
dev.off()

c_nl<- dat1p$home[which(dat1p$home==1 & dat1p$like==1)]
c_l <- dat1p$home[which(dat1p$home==1 & dat1p$like==0)]
nc_nl <- dat1p$home[which(dat1p$home==0 & dat1p$like==1)]
nc_l <- dat1p$home[which(dat1p$home==0 & dat1p$like==0)]
f<- cbind(v1=c(length(c_l),length(nc_l)),v2=c(length(c_nl),length(nc_nl)))
barplot(f, main="Group Barplots: Computer at home comparison", ylab= "Count",names.arg=c("Like","Don't Like"),beside=TRUE, col=c("green","darkblue"))
legend("topright", fill= c("green","darkblue"), legend=c("Have computer at home", "Don't have computer at home")  )

o_nl<- dat1p$own[which(dat1p$own==1 & dat1p$like==1)]
o_l <- dat1p$own[which(dat1p$own==1 & dat1p$like==0)]
no_nl <- dat1p$own[which(dat1p$own==0 & dat1p$like==1)]
no_l <- dat1p$own[which(dat1p$own==0 & dat1p$like==0)]
f<- cbind(v1=c(length(o_l),length(no_l)),v2=c(length(o_nl),length(no_nl)))
barplot(f, main="Group Barplots: Own PC comparison", ylab= "Count",names.arg=c("Like","Don't Like"),beside=TRUE, col=c("green","darkblue"))
legend("topright", fill= c("green","darkblue"), legend=c("Own PC", "Dont Own PC")  )

e_nl<- dat1p$educ[which(dat1p$educ==1 & dat1p$like==1)]
e_l <- dat1p$educ[which(dat1p$educ==1 & dat1p$like==0)]
ne_nl <- dat1p$educ[which(dat1p$educ==0 & dat1p$like==1)]
ne_l <- dat1p$educ[which(dat1p$educ==0 & dat1p$like==0)]
f<- cbind(v1=c(length(e_l),length(ne_l)),v2=c(length(e_nl),length(ne_nl)))
barplot(f, main="Group Barplots: Find Playing educational? comparison", ylab= "Count",names.arg=c("Like","Don't Like"),beside=TRUE, col=c("green","darkblue"))
legend("topright", fill= c("green","darkblue"), legend=c("Find Educational", "Don't find Educational")  )

d_nl<-dat1p$cdrom[which(dat1p$cdrom==1 & dat1p$like==1)]
d_l <- dat1p$cdrom[which(dat1p$cdrom==1 & dat1p$like==0)]
nd_nl <- dat1p$cdrom[which(dat1p$cdrom==0 & dat1p$like==1)]
nd_l <- dat1p$cdrom[which(dat1p$cdrom==0 & dat1p$like==0)]
f<- cbind(v1=c(length(d_l),length(nd_l)),v2=c(length(d_nl),length(nd_nl)))
barplot(f, main="Group Barplots: PS has CD-ROM?  comparison", ylab= "Count",names.arg=c("Like","Don't Like"),beside=TRUE, col=c("green","darkblue"))
legend("topright", fill= c("green","darkblue"), legend=c("Has CD-ROM", "No CD-ROM")  )

w_nl<-dat1p$work[which(dat1p$work>0 & dat1p$like==1)]
w_l <- dat1p$work[which(dat1p$work>0 & dat1p$like==0)]
nw_nl <- dat1p$work[which(dat1p$work==0 & dat1p$like==1)]
nw_l <- dat1p$work[which(dat1p$work==0 & dat1p$like==0)]
f<- cbind(v1=c(length(w_l),length(nw_l)),v2=c(length(w_nl),length(nw_nl)))
barplot(f, main="Group Barplots: Working or Non-Working?  comparison", ylab= "Count",names.arg=c("Like","Don't Like"),beside=TRUE, col=c("green","darkblue"))
legend("topright", fill= c("green","darkblue"), legend=c("Working", "Non-working")  )

#Scenario 6

# Get distribution of expected grades for sample students
summary(data$grade)

# ---- -----

# BOOSTRAP ON EXPECTED GRADES FROM SAMPLE #

# Bootstrap both theoretical and sample distributions to N=314
boot.Spopulation <- rep(data$grade, length.out = 314)

# Now we have Bootstrap pop, and draw 400 samples from each
runs = 400
boot.sample = array(dim = c(runs, 91))

for (i in 1:runs) {
  boot.sample[i, ] <- sample(boot.Spopulation, size=91, replace=FALSE)
}

# Calculate sample mean for each row of boot.Spop..
boot.Smean <- apply(X=boot.sample, MARGIN=1, FUN=mean)

# Plot histogram of sample Bootstrap means
hist(boot.Smean, breaks=20, probability=TRUE, density=20, col=3, border=3, main="Bootstrap Distribution of Sample (n= 91) Population's Expected Grade")
lines(density(boot.Smean, adjust=2), col=2)

# REPEAT PROCEDURE FOR THEORETICAL DISTRIBUTION

# Build Theoretical Grade distribution
TheoreticalGrade <- matrix(c(4,4,3,3,3,2,2,2,2,1), nrow=10, ncol=1)
boot.Tpopulation <- rep(TheoreticalGrade, length.out=314)
boot.theoretical = array(dim=c(runs, 91))

for (i in 1:runs) {
  boot.theoretical[i, ] <- sample(boot.Tpopulation, size=91, replace=FALSE)
}

# Calculate sample mean for each row of boot.theoretical
boot.TMean <- apply(X=boot.theoretical, MARGIN=1, FUN=mean)

# Plot histogram of theoretical Bootstrap means
hist(boot.TMean, breaks=20, probability=TRUE, density=20, col=3, border=3, main="Bootstrap Distribution of Theoretical (n = 91) Grade Distribution for Expected Grade")
lines(density(boot.TMean, adjust=2), col=2)

# Finally, do Q-Q plots of both distributions
par(pty='s')
qqnorm(boot.Smean, main="Normal Q-Q Plot of Sample")
qqline(boot.Smean)

qqnorm(boot.TMean, main="Normal Q-Q Plot of Theoretical")
qqline(boot.TMean)

# Perform S-W Normality Test
shapiro.test(boot.Smean)
shapiro.test(boot.TMean)

# ---- ----

# ADD 4 MORE DATAPOINTS TO ORIGINAL DATA AND REPEAT
newrow1 <- 0
appendedVideodata <- rbind(data, newrow1)
appendedVideodata <- rbind(appendedVideodata, newrow1)
appendedVideodata <- rbind(appendedVideodata, newrow1)
appendedVideodata <- rbind(appendedVideodata, newrow1)

# Set new rows' grades to 1
appendedVideodata[92, 'grade'] <- 1
appendedVideodata[93, 'grade'] <- 1
appendedVideodata[94, 'grade'] <- 1
appendedVideodata[95, 'grade'] <- 1

# Bootstrap to N=314
boot.S95population <- rep(appendedVideodata$grade, length.out = 314)

# Now we have Bootstrap pop, and draw 400 samples from each
runs = 400
boot.95sample = array(dim = c(runs, 95))

for (i in 1:runs) {
  boot.95sample[i, ] <- sample(boot.S95population, size=95, replace=FALSE)
}

# Calculate sample mean for each row of boot.Spop..
boot.S95mean <- apply(X=boot.95sample, MARGIN=1, FUN=mean)

# Plot histogram of sample Bootstrap means
hist(boot.S95mean, breaks=20, probability=TRUE, density=20, col=3, border=3, main="Bootstrap Distribution of Sample (n = 95) Population's Expected Grade")
lines(density(boot.S95mean, adjust=2), col=2)

mean(boot.S95mean)

# QQ Plots for n=95 bootstrap
par(pty='s')
qqnorm(boot.S95mean, main="Normal Q-Q Plot of Sample (n = 95)")
qqline(boot.S95mean)

# Perform S-W Normality Test
shapiro.test(boot.Smean)

# ---- ----

# CHECK GAMERS VS NON-GAMERS
gamers.ind <- which(data['freq'] <= 2)
nongamers.ind <- which(data['freq'] > 2)
video.gamers <- data[nongamers.ind,]
video.nongamers <- data[gamers.ind,]


# Get Bootstrap pop for both groups
boot.Gpopulation <- rep(video.gamers$grade, length.out = 314)
boot.Npopulation <- rep(video.nongamers$grade, length.out = 314)

# Now we have Bootstrap pop, and draw 400 samples from each
runs = 400
boot.gamer = array(dim = c(runs, 41))
boot.nongamer = array(dim = c(runs, 37))

for (i in 1:runs) {
  boot.gamer[i, ] <- sample(boot.Gpopulation, size=41, replace=FALSE)
  boot.nongamer[i, ] <- sample(boot.Npopulation, size = 37, replace = FALSE)
}

# Calculate sample mean for each row of boot.Spop..
boot.Gmean <- apply(X=boot.gamer, MARGIN=1, FUN=mean)
boot.Nmean <- apply(X=boot.nongamer, MARGIN=1, FUN=mean)

# Plot histogram of sample Bootstrap means
hist(boot.Gmean, breaks=20, probability=TRUE, density=20, col=3, border=3)
lines(density(boot.Gmean, adjust=2), col=2)

hist(boot.Nmean, breaks=20, probability = TRUE, density=20, col=3, border=3)
lines(density(boot.Nmean, adjust=2), col=2)

# Perform S-W test and QQplots
par(pty='s')
qqnorm(boot.Gmean, main="Normal Q-Q Plot of Gamers (n = 41)")
qqline(boot.Gmean)

qqnorm(boot.Nmean, main="Normal Q-Q Plot of Non-gamers (n = 37)")
qqline(boot.Nmean)

# Perform S-W Normality Test
shapiro.test(boot.Gmean)
shapiro.test(boot.Nmean)

