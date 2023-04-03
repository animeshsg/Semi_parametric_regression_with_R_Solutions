############ R script: BostMortGAMfit ##########

# For plotting the estimated function components of
# a generalized additive model fit to the Boston mortgage
# data.

# Last changed: 06 APR 2021

# Load required packages:

library(mgcv) ; library(HRW)

# Load in data:

data(BostonMortgages)
BostonMortgages$denyBinary <- as.numeric(BostonMortgages$deny == "yes")

# Obtain GAM fit:

fit1GAMBostMort <- gam(denyBinary ~ black + s(dir) + s(lvr)
                                    + pbcr + self + single + as.factor(ccs),
                                    family = binomial,data = BostonMortgages)

# Plot the estimated GAM components:

par(mfrow = c(2,2),bty = "l",cex.lab = 1.5,
     cex.main = 1.5,col.main = "navy",lwd = 2,
     mai = c(0.9,0.7,0.35,0.05))

plot(fit1GAMBostMort,shade = TRUE,shade.col = "palegreen",
     select = 1,xlim = c(0,1),ylab = "effect on logit(probab. of denial)",
     xlab = "debt payment to income ratio",
     main = "link scale",rug = FALSE)
rug(BostonMortgages$dir,col = "dodgerblue",quiet = TRUE)

plot(fit1GAMBostMort,shade = TRUE,shade.col = "palegreen",
     xlim = c(0,1), select = 2, ylab = "effect on logit(probab. of denial)",
     xlab = "loan size to property value ratio",main = "link scale",rug = FALSE)
rug(BostonMortgages$lvr,col = "dodgerblue",quiet = TRUE)

# Set function for obtaining the mode of a vector:

modalValue <- function(x)
   return(unique(x)[which.max(tabulate(match(x,unique(x))))])

# Set grids for `dir' and `lvg':

ng <- 401 ; dirg <- seq(0,1,length = ng) ; lvrg <- seq(0,1,length = ng)

# Obtain and plot slice of the probability surface
# in the `dir' direction corresponding to the modal
# values of categorical predictors and the mean of 
# other continuous predictors:

newdataDF <- data.frame(dir = dirg, 
                        black = modalValue(BostonMortgages$black),
                        lvr = mean(BostonMortgages$lvr), 
                        pbcr = modalValue(BostonMortgages$pbcr), 
                        self = modalValue(BostonMortgages$self),
                        single = modalValue(BostonMortgages$single),
                        ccs = modalValue(BostonMortgages$ccs))

predObjdir <- predict(fit1GAMBostMort,newdata = newdataDF,
                      type = "response",se.fit = TRUE)
etahatdirg <- predObjdir$fit
lowdirg <- etahatdirg - qnorm(0.975)*predObjdir$se.fit
uppdirg <- etahatdirg + qnorm(0.975)*predObjdir$se.fit
lowdirg[lowdirg<0] <- 0 ; lowdirg[lowdirg>0.5] <- 0.5

plot(0,type = "n",xlim = c(0,1),ylim = c(0,0.5),xlab  =  "debt payments to income ratio",
     ylab  =  "probability of denial",main  =  "response scale")
polygon(c(dirg,rev(dirg)),c(lowdirg,rev(uppdirg)),col = "palegreen",border = FALSE)
lines(dirg,etahatdirg,col = "darkgreen",lwd = 2)
rug(BostonMortgages$dir,col = "dodgerblue",quiet = TRUE)

# Obtain and plot slice of the probability surface
# in the `lvr' direction corresponding to the modal
# values of categorical predictors and the mean of 
# other continuous predictors:

newdataDF <- data.frame(dir = mean(BostonMortgages$dir), 
                        black = modalValue(BostonMortgages$black),
                        lvr = lvrg,
                        pbcr = modalValue(BostonMortgages$pbcr), 
                        self = modalValue(BostonMortgages$self),
                        single = modalValue(BostonMortgages$single),
                        ccs = modalValue(BostonMortgages$ccs))

predObjlvr <- predict(fit1GAMBostMort,newdata = newdataDF,
                      type = "response",se.fit = TRUE)
etahatlvrg <- predObjlvr$fit
lowlvrg <- etahatlvrg - qnorm(0.975)*predObjlvr$se.fit
upplvrg <- etahatlvrg + qnorm(0.975)*predObjlvr$se.fit
lowlvrg[lowlvrg<0] <- 0

plot(0,type = "n",xlim = c(0,1),ylim = c(0,0.5),xlab  =  "loan size to property value ratio",
     ylab  =  "probability of denial",main  =  "response scale")
polygon(c(lvrg,rev(lvrg)),c(lowlvrg,rev(upplvrg)),col = "palegreen",border = FALSE)
lines(lvrg,etahatlvrg,col = "darkgreen",lwd = 2)
rug(BostonMortgages$lvr,col = "dodgerblue",quiet = TRUE)

############ End of BostMortGAMfit ############


