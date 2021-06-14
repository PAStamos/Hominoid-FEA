library("ggplot2")
library("rstan")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(rethinking)

setwd("C:/Users/peter/Google Drive/R Files/FEA_Project")


#read in data
csv.data <- read.csv(file = "FEA_results.csv", stringsAsFactors = TRUE)


##For Scaled_VM_stress_98 or Scaled_VM_stress_Mean
#Subset data to just the rows without NA in Scaled_VM_stress_98_to_UKEN_SC026 column
d.sub <- csv.data[is.na(csv.data$Scaled_VM_stress_98_to_UKEN_SC026)==FALSE , ]
#Subset data to just the rows without NA in Scaled_VM_stress_Mean
d.sub <- csv.data[is.na(csv.data$Scaled_VM_stress_Mean)==FALSE , ]

##4-way model
d.sub$GenusLoco_index <-as.integer(d.sub$GenusLoco)



#Rename column for ease of reference
d.sub$ScaledVM <- d.sub$Scaled_VM_stress_98_to_UKEN_SC026
#or
d.sub$ScaledVM <- d.sub$Scaled_VM_stress_Mean

#Subset the data to only the parameters used in the model. Stan likes to throw errors otherwise
d.sub2 <- list(ScaledVM=d.sub$ScaledVM, GenusLoco_index=d.sub$GenusLoco_index)


mstan.ScaledVM <- ulam(
  alist(
    ScaledVM ~ dnorm(mu , sigma),
    mu <- a[GenusLoco_index],
    a[GenusLoco_index] ~ dnorm (0 , sigma_genusloco),
    sigma ~ dunif(0 , 1),
    sigma_genusloco ~ dunif(0 , 1)
  ) , data=d.sub2, iter=10000, chains=4,
  contro=list(adapt_delta = 0.99))

labels <- paste("a[", 1:4, "]:" , levels(d.sub$GenusLoco) , sep="" )
#labels without numbered prefix
labels <- paste("" , levels(d.sub$GenusLoco) , sep="" )

plot(precis(mstan.ScaledVM , depth=2 , digits = 3, prob=0.95, pars = "a") , labels=labels , pch = 16, col = col.alpha("blue",1.0),
     xlab="Scaled von Mises stress")

#postcheck(mstan.ScaledVM, prob = 0.95)

#Contrasts between Homo and Pan
poststan.ScaledVM <- extract.samples(mstan.ScaledVM , n=1e7)
poststan.ScaledVM$diff_HFE <- poststan.ScaledVM$a[,2] - poststan.ScaledVM$a[,1]
poststan.ScaledVM$ratio_HFE <- poststan.ScaledVM$a[,2] / poststan.ScaledVM$a[,1]
poststan.ScaledVM$diff_PFE <- poststan.ScaledVM$a[,4] - poststan.ScaledVM$a[,3]
poststan.ScaledVM$ratio_PFE <- poststan.ScaledVM$a[,4] / poststan.ScaledVM$a[,3]
poststan.ScaledVM$diff_HPE <- poststan.ScaledVM$a[,1] - poststan.ScaledVM$a[,3]
poststan.ScaledVM$ratio_HPE <- poststan.ScaledVM$a[,1] / poststan.ScaledVM$a[,3]
poststan.ScaledVM$diff_HPF <- poststan.ScaledVM$a[,2] - poststan.ScaledVM$a[,4]
poststan.ScaledVM$ratio_HPF <- poststan.ScaledVM$a[,2] / poststan.ScaledVM$a[,4]
precis(poststan.ScaledVM , depth=2 , digits = 3, prob=0.95)

dens(poststan.ScaledVM$diff_HFE , show.HPDI = 0.95 , show.zero = TRUE) #Difference between Homo Climb and Homo Walk
dens(poststan.ScaledVM$diff_PFE , show.HPDI = 0.95 , show.zero = TRUE) #Difference between Pan Climb and Pan Walk
dens(poststan.ScaledVM$diff_HPE , show.HPDI = 0.95 , show.zero = TRUE) #Difference between Homo Walk and Pan Walk
dens(poststan.ScaledVM$diff_HPF , show.HPDI = 0.95 , show.zero = TRUE) #Difference between Homo Climb and Pan Climb

dens(poststan.ScaledVM$ratio_HFE , show.HPDI = 0.95 , show.zero = TRUE)
dens(poststan.ScaledVM$ratio_PFE , show.HPDI = 0.95 , show.zero = TRUE)
dens(poststan.ScaledVM$ratio_HPE , show.HPDI = 0.95 , show.zero = TRUE)
dens(poststan.ScaledVM$ratio_HPF , show.HPDI = 0.95 , show.zero = TRUE)



#modeling relationship between ariaDNE and von Mises stress
##2-way models
d.sub3 <- d.sub
d.sub3$Genus_index <- as.integer(d.sub3$Species)
d.sub3$Loco_index <- as.integer(d.sub3$Motion)


#Rename column for ease of reference
d.sub3$ScaledVM_Max <- d.sub3$Scaled_VM_stress_98_to_UKEN_SC026
d.sub3$ScaledVM_Mean <- d.sub3$Scaled_VM_stress_Mean



#Climbing models
d.sub3.climb <- d.sub3[d.sub3$Motion=='climb' , ]

##ScaledVM_Max
#Subset the data to only the parameters used in the model. Stan likes to throw errors otherwise
d.sub3.climb <- list(ScaledVM_Max=d.sub3.climb$ScaledVM_Max, ariaDNE_bw_10=d.sub3.climb$ariaDNE_bw_10)


mstan.ariaDNE.climb <- ulam(
  alist(
    ScaledVM_Max ~ dnorm(mu , sigma),
    mu <- a + b * ariaDNE_bw_10  ,
    a ~ dnorm( 0 , 1),
    b ~ dnorm( 0 , 1),
    sigma ~ dunif(0 , 1)
  ) , data=d.sub3.climb, iter=10000, chains=4,
  contro=list(adapt_delta = 0.99))

precis(mstan.ariaDNE.climb , depth=2 , digits = 3, prob=0.95)

##Make figures with shaded CI regions. Based on Edition 2 R code 4.57 and 4.58
#ClimbMax
poststan.mstan.ariaDNE.climb <- extract.samples(mstan.ariaDNE.climb, n=1e6)
mu.link <- function(ariaDNE_bw_10) poststan.mstan.ariaDNE.climb$a +poststan.mstan.ariaDNE.climb$b*(ariaDNE_bw_10)
ariaDNE.seq <- seq(from=0.01 , to=0.08, by=0.002)
mu <- sapply (ariaDNE.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI , prob=0.95)

plot(ScaledVM_Max ~ ariaDNE_bw_10, data=d.sub3.climb, col = 'black', lwd= 2, cex = 1.1, xlab = "ariaDNE", ylab = "scaled stress", main = "Flexed Maximum von Mises stress") 
lines(ariaDNE.seq, mu.mean,col = 'black', lwd= 1, cex = 1,)
shade(mu.PI, ariaDNE.seq, col = col.alpha("blue",0.20))


##ScaledVM_Mean
#Subset the data to only the parameters used in the model. Stan likes to throw errors otherwise
d.sub3.climb <- list(ScaledVM_Mean=d.sub3.climb$ScaledVM_Mean, ariaDNE_bw_10=d.sub3.climb$ariaDNE_bw_10)

mstan.ariaDNE.climb <- ulam(
  alist(
    ScaledVM_Mean ~ dnorm(mu , sigma),
    mu <- a + b * ariaDNE_bw_10  ,
    a ~ dnorm( 0 , 1),
    b ~ dnorm( 0 , 1),
    sigma ~ dunif(0 , 1)
  ) , data=d.sub3.climb, iter=10000, chains=4,
  contro=list(adapt_delta = 0.99))

precis(mstan.ariaDNE.climb , depth=2 , digits = 3, prob=0.95)

#make figure
poststan.mstan.ariaDNE.climb <- extract.samples(mstan.ariaDNE.climb, n=1e6)
mu.link <- function(ariaDNE_bw_10) poststan.mstan.ariaDNE.climb$a +poststan.mstan.ariaDNE.climb$b*(ariaDNE_bw_10)
ariaDNE.seq <- seq(from=0.01 , to=0.08, by=0.002)
mu <- sapply (ariaDNE.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI , prob=0.95)

plot(ScaledVM_Mean ~ ariaDNE_bw_10, data=d.sub3.climb, col = 'black', lwd= 2, cex = 1.1, xlab = "ariaDNE", ylab = "scaled stress", main = "Flexed Mean von Mises stress") 
lines(ariaDNE.seq, mu.mean)
shade(mu.PI, ariaDNE.seq, col = col.alpha("blue",0.20))



##Walking models
#
d.sub3.walk <- d.sub3[d.sub3$Motion=='walk' , ]

#ScaledVM_Max
#Subset the data to only the parameters used in the model. Stan likes to throw errors otherwise
d.sub3.walk <- list(ScaledVM_Max=d.sub3.walk$ScaledVM_Max, ariaDNE_bw_10=d.sub3.walk$ariaDNE_bw_10)

mstan.ariaDNE.walk <- ulam(
  alist(
    ScaledVM_Max ~ dnorm(mu , sigma),
    mu <- a + b * ariaDNE_bw_10  ,
    a ~ dnorm( 0 , 1),
    b ~ dnorm( 0 , 1),
    sigma ~ dunif(0 , 1)
  ) , data=d.sub3.walk, iter=10000, chains=4,
  contro=list(adapt_delta = 0.99))

precis(mstan.ariaDNE.walk , depth=2 , digits = 3, prob=0.95)

#make figure
poststan.mstan.ariaDNE.walk <- extract.samples(mstan.ariaDNE.walk, n=1e6)
mu.link <- function(ariaDNE_bw_10) poststan.mstan.ariaDNE.walk$a +poststan.mstan.ariaDNE.walk$b*(ariaDNE_bw_10)
ariaDNE.seq <- seq(from=0.01 , to=0.08, by=0.002)
mu <- sapply (ariaDNE.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI , prob=0.95)

plot(ScaledVM_Max ~ ariaDNE_bw_10, data=d.sub3.walk, col = 'black', lwd= 2, cex = 1.1, xlab = "ariaDNE", ylab = "scaled stress", main = "Extended Maximum von Mises stress") 
lines(ariaDNE.seq, mu.mean)
shade(mu.PI, ariaDNE.seq, col = col.alpha("blue",0.20))

#ScaledVM_Mean
#Subset the data to only the parameters used in the model. Stan likes to throw errors otherwise
d.sub3.walk <- list(ScaledVM_Mean=d.sub3.walk$ScaledVM_Mean, ariaDNE_bw_10=d.sub3.walk$ariaDNE_bw_10)

mstan.ariaDNE.walk <- ulam(
  alist(
    ScaledVM_Mean ~ dnorm(mu , sigma),
    mu <- a + b * ariaDNE_bw_10  ,
    a ~ dnorm( 0 , 1),
    b ~ dnorm( 0 , 1),
    sigma ~ dunif(0 , 1)
  ) , data=d.sub3.walk, iter=10000, chains=4,
  contro=list(adapt_delta = 0.99))

precis(mstan.ariaDNE.walk , depth=2 , digits = 3, prob=0.95)

#make figure
poststan.mstan.ariaDNE.walk <- extract.samples(mstan.ariaDNE.walk, n=1e6)
mu.link <- function(ariaDNE_bw_10) poststan.mstan.ariaDNE.walk$a +poststan.mstan.ariaDNE.walk$b*(ariaDNE_bw_10)
ariaDNE.seq <- seq(from=0.01 , to=0.08, by=0.002)
mu <- sapply (ariaDNE.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI , prob=0.95)

plot(ScaledVM_Mean ~ ariaDNE_bw_10, data=d.sub3.walk, col = 'black', lwd= 2, cex = 1.1, xlab = "ariaDNE", ylab = "scaled stress", main = "Extended Mean von Mises stress") 
lines(ariaDNE.seq, mu.mean)
shade(mu.PI, ariaDNE.seq, col = col.alpha("blue",0.20))







