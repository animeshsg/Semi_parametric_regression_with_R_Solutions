######### R script: WarsawAptsSimpFacByCurv ##########

# For fitting a district-specific curves model to the
# Warsaw apartment data. The district division is
# Srodmiescie and non-Srodmiescie.

# Last changed: 29 AUG 2018

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Load required packages:

library(HRW);  library(rstan)

# Set MCMC parameters:

nWarm <- 1000           # Length of burn-in.
nKept <- 1000           # Size of the kept sample.
nThin <- 1              # Thinning factor. 

# Load Warsaw apartment data:

data(WarsawApts)

# Set graphics parameters:

cex.axisVal <- 1.8 ; cex.labVal <- 1.8
cex.legendVal <- 1.5 ; lwdVal <- 2

# Standardize both predictor and response variable
# and set hyperparameter values:

xOrig <- WarsawApts$construction.date
yOrig <- WarsawApts$areaPerMzloty
mean.x <- mean(xOrig)  ; sd.x <- sd(xOrig)
mean.y <- mean(yOrig)  ; sd.y <- sd(yOrig)
x <- (xOrig - mean.x)/sd.x
y <- (yOrig - mean.y)/sd.y
treatIndic <- as.numeric(as.character(WarsawApts$district) != "Mokotow")

sigmaBeta <- 1e5 ; Au <- 1e5 ; Aeps <- 1e5

# Obtain linear and spline basis design matrices (X and Z):

X <- cbind(rep(1,length(y)),x,treatIndic,treatIndic*x)
aOrig <- min(xOrig) ; bOrig <- max(xOrig)
a <- (aOrig - mean.x)/sd.x  ;  b <- (bOrig - mean.x)/sd.x
numIntKnots <- 25
intKnots <-  quantile(unique(x),seq(0,1,length = numIntKnots+2)
                      [-c(1,numIntKnots+2)])
Z <- ZOSull(x,intKnots = intKnots,range.x = c(a,b))

# Set dimension variables:

n <- length(y)
ncZ <- ncol(Z)

# Specify model in Stan:

WarsawAptSimpFacByCurvModel <- 
   'data
   {
      int<lower=1> n;           int<lower=1> ncZ;
      vector[n] y;              matrix[n,4] X;        
      matrix[n,ncZ] Z ;         vector[n] treatIndic;              
      real<lower=0> sigmaBeta;  
      real<lower=0> Aeps;       real<lower=0> Au;           
   }
   parameters 
   {
      vector[4] beta;       
      vector[ncZ] uControl;           vector[ncZ] uTreatmt;            
      real<lower=0> sigmaUcontrol;    real<lower=0> sigmaUtreatmt;
      real<lower=0> sigmaEps;   
   }
   model 
   {
      for (i in 1:n)
         y[i] ~ normal((dot_product(beta,X[i])
                       + dot_product(uControl,((1-treatIndic[i])*Z[i]))
                       + dot_product(uTreatmt,(treatIndic[i]*Z[i]))),sigmaEps);
      uControl ~ normal(0,sigmaUcontrol) ; uTreatmt ~ normal(0,sigmaUtreatmt); 
      beta ~ normal(0,sigmaBeta);
      sigmaEps ~ cauchy(0,Aeps); 
      sigmaUcontrol ~ cauchy(0,Au); sigmaUtreatmt ~ cauchy(0,Au); 
   }'

# Set up input data:

allData <- list(n = n,ncZ = ncZ,y = y,X = X,Z = Z,treatIndic = treatIndic,
                sigmaBeta = 1e5,Aeps = 1e5,Au = 1e5)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = WarsawAptSimpFacByCurvModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = WarsawAptSimpFacByCurvModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 100,
                 fit = stanCompilObj)

# Extract relevant MCMC sampes:

betaMCMC <- NULL
for (j in 1:4)
{
   charVar <- paste("beta[",as.character(j),"]",sep = "") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar,permuted = FALSE))
}
uControlMCMC <- NULL  ;  uTreatmtMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("uControl[",as.character(k),"]",sep = "") 
   uControlMCMC <- rbind(uControlMCMC,extract(stanObj,charVar,permuted = FALSE))
                   
   charVar <- paste("uTreatmt[",as.character(k),"]",sep = "") 
   uTreatmtMCMC <- rbind(uTreatmtMCMC,extract(stanObj,charVar,permuted = FALSE))
}
sigmaEpsMCMC <- as.vector(extract(stanObj,"sigmaEps",permuted = FALSE))
sigmaUcontrolMCMC <- as.vector(extract(stanObj,"sigmaUcontrol",permuted = FALSE))
sigmaUtreatmtMCMC <- as.vector(extract(stanObj,"sigmaUtreatmt",permuted = FALSE))

# Plot data and fitted curves:

ng <- 201
cex.labVal <- 1.8
xLow <- min(x) ; xUpp <- max(x)
xg <- seq(xLow,xUpp,length = ng)

Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,intKnots = intKnots,range.x = c(a,b))

par(mai = c(1.02,0.9,0.82,0.42))
plot(xOrig,yOrig,type = "n",bty = "l",xlab = "construction date (year)",
     ylab = "area (square meters) per million zloty",
     cex.lab = cex.labVal,cex.axis = cex.axisVal)

points(xOrig[treatIndic == 0],yOrig[treatIndic == 0],col = "dodgerblue")

points(xOrig[treatIndic == 1],yOrig[treatIndic == 1],col = "deeppink")

fhatTreatMCMC <- Xg%*%(betaMCMC[c(1,2),]+ betaMCMC[c(3,4),]) + Zg%*%uTreatmtMCMC

fCntrlg <- apply(fhatCntrlMCMC,1,mean)
credLowCntrl <- apply(fhatCntrlMCMC,1,quantile,0.025)
credUppCntrl <- apply(fhatCntrlMCMC,1,quantile,0.975)

fTreatg <- apply(fhatTreatMCMC,1,mean)
credLowTreat <- apply(fhatTreatMCMC,1,quantile,0.025)
credUppTreat <- apply(fhatTreatMCMC,1,quantile,0.975)

x <- (xOrig - mean.x)/sd.x

# Convert grids to original units:

xgOrig <- sd.x*xg + mean.x

fCntrlgOrig <- sd.y*fCntrlg + mean.y
credLowCntrlOrig <- sd.y*credLowCntrl + mean.y
credUppCntrlOrig <- sd.y*credUppCntrl + mean.y

fTreatgOrig <- sd.y*fTreatg + mean.y
credLowTreatOrig <- sd.y*credLowTreat + mean.y
credUppTreatOrig <- sd.y*credUppTreat + mean.y

lines(xgOrig,fCntrlgOrig,col = "blue",lwd = 2)
lines(xgOrig,credLowCntrlOrig,col = "blue",lty = 2,lwd = 2)
lines(xgOrig,credUppCntrlOrig,col = "blue",lty = 2,lwd = 2)

lines(xgOrig,fTreatgOrig,col = "darkmagenta",lwd = 2)
lines(xgOrig,credLowTreatOrig,col = "darkmagenta",lty = 2,lwd = 2)
lines(xgOrig,credUppTreatOrig,col = "darkmagenta",lty = 2,lwd = 2)

legend(1971,65,legend = c("Mokotow","other districts"),
       col = c("deeppink","dodgerblue"),pch = rep(1,2),cex = cex.legendVal,
       pt.cex = rep(1,2))

# Plot contrast curve fit:

estFunCol <- "darkgreen"; 
varBandCol <- "palegreen"

ContrastMCMC <- Xg%*%betaMCMC[c(3,4),] + Zg%*%(uTreatmtMCMC-uControlMCMC)

credLower <- apply(ContrastMCMC,1,quantile,0.025)
credUpper <- apply(ContrastMCMC,1,quantile,0.975)
Contrastg <- apply(ContrastMCMC,1,mean)

# Convert to original units:

ContrastgOrig <- sd.y*Contrastg
credLowerOrig <- sd.y*credLower 
credUpperOrig <- sd.y*credUpper 

par(mfrow = c(1,1),mai = c(1.02,0.9,0.82,0.42))
plot(0,0,type = "n",bty = "l",xlim = range(xgOrig),
     range(c(credLowerOrig,credUpperOrig)),
     xlab = "construction date (year)",
     ylab = " difference in area (square meters) per million zloty",
     cex.lab = cex.labVal,cex.axis = cex.axisVal)

polygon(c(xgOrig,rev(xgOrig)),c(credLowerOrig,rev(credUpperOrig)),
        col = varBandCol,border = FALSE)

lines(xgOrig,ContrastgOrig,lwd = 2,col = estFunCol)
abline(0,0,col = "slateblue",lwd = 2)

############ End of WarsawAptsSimpFacByCurv ############

