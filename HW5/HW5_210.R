########## R script: WarsawAptsBayes ##########

# For obtaining a Bayesian penalized spline
# fit to the Warsaw apartment data using Stan
# via the package rstan. Approximate Bayesian
# inference is achieved via Markov chain
# Monte Carlo (MCMC).

# Last changed: 05 JAN 2019

# Load required packages:

library(HRW) ; data(TreasuryRate)
date <- (as.Date(TreasuryRate$date,"%m/%d/%Y")[!is.na(TreasuryRate$rate)])
r <- TreasuryRate$rate[!is.na(TreasuryRate$rate)]
plot(date,r,type = "l",bty = "l",col = "darkgreen",xlab = "date",ylab="U.S. Treasury rate")

# Standardize both predictor and response variable
# and set hyperparameter values:

yOrig=r[-1]-r[-length(r)]
xOrig=r[-length(r)]
mean.x <- mean(xOrig)  ; sd.x <- sd(xOrig)
mean.y <- mean(yOrig)  ; sd.y <- sd(yOrig)
x <- (xOrig - mean.x)/sd.x
y <- (yOrig - mean.y)/sd.y
sigmaBeta <- 1e5 ; Au <- 1e5 ; Aeps <- 1e5
sigmaGamma<-1e-5; Bv<-1e-5;
# Obtain linear and spline basis design matrices (X and Z):

X <- cbind(rep(1,length(y)),x)
aOrig <- min(xOrig) ; bOrig <- max(xOrig)
a <- (aOrig - mean.x)/sd.x  ;  b <- (bOrig - mean.x)/sd.x
numIntKnots <- 10
intKnots <-  quantile(unique(x),seq(0,1,length=numIntKnots+2)
                      [-c(1,numIntKnots+2)])
Z <- ZOSull(x,intKnots=intKnots,range.x=c(a,b))
ncZ <- ncol(Z)

numIntKnotsB <- 100
intKnotsB <-  quantile(unique(x),seq(0,1,length=numIntKnots+2)
                      [-c(1,numIntKnots+2)])
W <- ZOSull(x,intKnots=intKnots,range.x=c(a,b))
ncW <- ncol(W)

# Specify model in Stan:

npRegModel <- 
   'data
   {
      int<lower=1> n;         int<lower=1> ncZ;
      vector[n] y;            matrix[n,2] X;
      matrix[n,ncZ] Z;        real<lower=0> sigmaBeta;
      real<lower=0> Au;       real<lower=0> Aeps;
      int<lower=1> ncW;       matrix[n,ncW] W;
      real<lower=0> sigmaGamma;
      real<lower=0> Bv;
   }
   parameters 
   {
      vector[2] beta;          vector[ncZ] u;
      real<lower=0> sigmaeps;  real<lower=0> sigmau;
      vector[2] gamma;          vector[ncW] v;
      real<lower=0> sigmav;
   }
   transformed parameters 
   {
      vector[n] f; // f function
      vector[n] g; // g function
      f = X*beta + Z*u;
      g = exp(X*gamma + W*v);
   }
   model 
   {
      y ~ normal(f,sqrt(g));
      u ~ normal(0,sigmau); beta ~ normal(0,sigmaBeta); 
      sigmaeps ~ cauchy(0,Aeps); sigmau ~ cauchy(0,Au);
      sigmav ~ cauchy(0,Bv); gamma ~ normal(0,sigmaGamma);
      v ~ normal(0,sigmav);
   }'

# Store data in a list in format required by Stan: 

allData <- list(n=length(x),ncZ=ncZ,y=y,X=X,Z=Z,
                sigmaBeta=sigmaBeta,Au=Au,Aeps=Aeps,
                sigmaGamma=sigmaGamma,Bv=Bv,
                W=W,ncW=ncW)

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code=npRegModel,data=allData,
                         iter=1,chains=1)

# Set MCMC sample size parameters:

nWarm <- 1000        # Length of warm-up.
nKept <- 2000        # Size of the kept sample.
nThin <- 2           # Thinning factor. 

# Obtain MCMC samples for each parameter using Stan:

initFun <- function()
   return(list(sigmau=1,sigmaeps=0.7,beta=rep(0,2),u=rep(0,ncZ),sigmav=1,gamma=rep(0,2),v=rep(0,ncW)))

stanObj <-  stan(model_code=npRegModel,data=allData,warmup=nWarm,
                 iter=(nWarm+nKept),chains=1,thin=nThin,refresh=100,
                 fit=stanCompilObj,init=initFun,seed=13)

# Extract relevant MCMC samples:

betaMCMC <- NULL
for (j in 1:2)
{
   charVar <- paste("beta[",as.character(j),"]",sep="") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar,permuted=FALSE))
}
gammaMCMC <- NULL
for (l in 1:2)
{
  charVar <- paste("gamma[",as.character(l),"]",sep="") 
  gammaMCMC <- rbind(gammaMCMC,extract(stanObj,charVar,permuted=FALSE))
}
uMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("u[",as.character(k),"]",sep="") 
   uMCMC <- rbind(uMCMC,extract(stanObj,charVar,permuted=FALSE))
}
vMCMC <- NULL
for (m in 1:ncW)
{
  charVar <- paste("v[",as.character(m),"]",sep="") 
  vMCMC <- rbind(vMCMC,extract(stanObj,charVar,permuted=FALSE))
}
sigmaepsMCMC <- as.vector(extract(stanObj,"sigmaeps",permuted=FALSE))
sigmauMCMC <- as.vector(extract(stanObj,"sigmau",permuted=FALSE))
sigmavMCMC <- as.vector(extract(stanObj,"sigmav",permuted=FALSE))

# Obtain MCMC samples of regression curves over a fine grid:

ng <- 101
xgOrig <- seq(aOrig,bOrig,length=ng)
xg <- (xgOrig - mean.x)/sd.x
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,intKnots=intKnots,range.x=c(a,b))
Wg <- ZOSull(xg,intKnots=intKnotsB,range.x=c(a,b))
fhatMCMC <- Xg%*%betaMCMC + Zg%*%uMCMC
ghatMCMC <- exp(Xg%*%gammaMCMC+Wg%*%vMCMC)

# Convert fhatMCMC matrix to original scale:

fhatMCMCOrig <- fhatMCMC*sd.y + mean.y
fhatgOrig <- apply(fhatMCMCOrig,1,mean)
credLower <- apply(fhatMCMCOrig,1,quantile,0.025)
credUpper <- apply(fhatMCMCOrig,1,quantile,0.975)

ghatMCMCOrig <- sqrt(ghatMCMC)*sd.y + mean.y
ghatgOrig <- apply(ghatMCMCOrig,1,mean)
credLowerg <- apply(ghatMCMCOrig,1,quantile,0.025)
credUpperg <- apply(ghatMCMCOrig,1,quantile,0.975)

# Display the fit:

par(mai=c(1,1.1,0.1,0.1))
cex.labVal <- 2   ;   cex.axisVal <- 1.5
plot(xOrig,yOrig,type="n",xlab="construction date (year)",
     ylab="area (square meters) per million zloty",
     bty="l",xlim=range(xgOrig),ylim=range(c(credLower,credUpper,yOrig)),
     cex.lab=cex.labVal,cex.axis=cex.axisVal)
polygon(c(xgOrig,rev(xgOrig)),c(credLowerg,rev(credUpperg)),
        col="pink",border=FALSE)
polygon(c(xgOrig,rev(xgOrig)),c(credLower,rev(credUpper)),
        col="palegreen",border=FALSE)
points(xOrig,yOrig,col="dodgerblue")
lines(xgOrig,fhatgOrig,col="darkgreen",lwd=2)
lines(xgOrig,ghatgOrig,col="red",lwd=2)

abline(v=quantile(xOrig,0.25),lty=2,col="darkorange")
abline(v=quantile(xOrig,0.50),lty=2,col="darkorange")
abline(v=quantile(xOrig,0.75),lty=2,col="darkorange")

## Standardized residuals test

ng <- 2987
xgOrig <- seq(aOrig,bOrig,length=ng)
xg <- (xgOrig - mean.x)/sd.x
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,intKnots=intKnots,range.x=c(a,b))
Wg <- ZOSull(xg,intKnots=intKnotsB,range.x=c(a,b))
fhatMCMC <- Xg%*%betaMCMC + Zg%*%uMCMC
ghatMCMC <- exp(Xg%*%gammaMCMC+Wg%*%vMCMC)
fhatMCMCOrig<-apply(fhatMCMC,1,mean)
ghatMCMCOrig<-apply(ghatMCMC,1,mean)
EhatSq<-((yOrig-fhatMCMCOrig)/sqrt(ghatMCMCOrig))**2
xg<-seq(min(date),max(date),length=2987)
plot(xg,EhatSq,col="white")
lines(xg,EhatSq,col="blue")