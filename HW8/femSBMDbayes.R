########## R script: femSBMDbayes ##########

# For conducting a Bayesian additive mixed
# model analysis on the spinal bone mineral
# density data using Markov chain Monte Carlo 
# and Stan.

# Last changed: 29 AUG 2018

# Set flag for code compilation (needed if 
# running script first time in current session):

compileCode <- TRUE

# Set MCMC sample size paramaters:

nWarm <- 1500
nKept <- 1500
nThin <- 1

# Load required packages:

library(HRW)   ;   library(rstan)   ;   library(lattice)

# Load in the female spinal bone data and extract the 
# component variables:

data(femSBMD)
femSBMD$asian<-ifelse(femSBMD$ethnicity=='Asian',1,0)
idnum <- femSBMD$idnum 
x1 <- femSBMD$black 
x2 <- femSBMD$hispanic 
x3 <- femSBMD$white
x0 <- femSBMD$asian
x4Orig <- femSBMD$age 
yOrig <- femSBMD$spnbmd

# Standardise data for Bayesian analysis:

mean.x4 <- mean(x4Orig) ; sd.x4 <- sd(x4Orig)
mean.y <- mean(yOrig)   ; sd.y <- sd(yOrig)
x4 <- (x4Orig - mean.x4)/sd.x4
y <- (yOrig - mean.y)/sd.y

# Set up matrices for additive model:

numObs <- length(y)
X <- cbind(rep(1,numObs),x0,x1,x2,x3,x4)
ncX <- ncol(X)

numIntKnots <- 15
intKnots <- quantile(unique(x4),seq(0,1,length = 
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])
range.x4 <- c(1.01*min(x4)-0.01*max(x4),1.01*max(x4)-0.01*min(x4))
Zspl <- ZOSull(x4,intKnots = intKnots,range.x = range.x4)
ncZ <- ncol(Zspl)
Zspl_asian<-X[,2]*Zspl;Zspl_black<-X[,3]*Zspl;
Zspl_hispanic<-X[,4]*Zspl;Zspl_white<-X[,5]*Zspl;

numGrp <- length(unique(idnum))
numObs <- length(y)

# Set hyperparameters:

sigmaBeta <- 1e5 ; AU <- 1e5 ; Au_asian <- 1e5 ; Au_black <- 1e5 ; Au_hispanic <- 1e5 ; Au_white <- 1e5 ; Aeps <- 1e5

# Specify model in Stan:

addMixModModel <- 
'data
{
   int<lower=1> numObs;         int<lower=1> numGrp;
   int<lower=1> ncX;            int<lower=1> ncZ;        
   real<lower=0> sigmaBeta;     real<lower=0> AU;
   real<lower=0> Aeps;          real<lower=0> Au_asian;
   vector[numObs] y;            int<lower=1> idnum[numObs];
   matrix[numObs,ncX] X;        matrix[numObs,ncZ] Zspl_asian;
   real<lower=0> Au_black;           real<lower=0> Au_hispanic;
   real<lower=0> Au_white;      matrix[numObs,ncZ] Zspl_black;
   matrix[numObs,ncZ] Zspl_hispanic;  matrix[numObs,ncZ] Zspl_white;
}
parameters
{
   vector[ncX] beta;            vector[numGrp] U;            
   vector[ncZ] u_asian;         real<lower=0> sigmaU;        
   real<lower=0> sigmau_asian;        real<lower=0> sigmaEps;    
   vector[ncZ] u_black;         vector[ncZ] u_hispanic;
   vector[ncZ] u_white;         real<lower=0> sigmau_black;
   real<lower=0> sigmau_hispanic;        real<lower=0> sigmau_white; 
}
model
{
   y ~ normal(X*beta + U[idnum] + Zspl_asian*u_asian+Zspl_black*u_black+Zspl_hispanic*u_hispanic+Zspl_white*u_white,sigmaEps);
   U ~ normal(0,sigmaU);           u_asian ~ normal(0,sigmau_asian);
   beta  ~ normal(0,sigmaBeta) ;   sigmaEps ~ cauchy(0,Aeps);
   sigmaU ~ cauchy(0,AU) ;   sigmau_asian ~ cauchy(0,Au_asian);
   u_black ~ normal(0,sigmau_black);  sigmau_black ~ cauchy(0,Au_black);
   u_hispanic ~ normal(0,sigmau_hispanic);  sigmau_hispanic ~ cauchy(0,Au_hispanic);
   u_white ~ normal(0,sigmau_white);  sigmau_white ~ cauchy(0,Au_white);
}'
# Zspl*u_black*X[,3]+Zspl*u_hispanic*X[,3]+Zspl*u_white*X[,4]
# Fit model using MCMC via Stan:

allData <- list(numObs = numObs,numGrp = numGrp,ncX = ncX,ncZ = ncZ,
                idnum = idnum,X = X,y = y,Zspl_asian = Zspl_asian,sigmaBeta = sigmaBeta,
                AU = AU,Au_asian = Au_asian,Au_black = Au_black,Au_hispanic = Au_hispanic,
                Au_white = Au_white,Aeps = Aeps,Zspl_black = Zspl_black,
                Zspl_hispanic = Zspl_hispanic,Zspl_white = Zspl_white)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = addMixModModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = addMixModModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 100,
                 fit = stanCompilObj)

# Save and extract relevant MCMC samples:

beta0MCMC <- as.vector(extract(stanObj,"beta[1]",permuted = FALSE))
beta1MCMC <- as.vector(extract(stanObj,"beta[2]",permuted = FALSE))
beta2MCMC <- as.vector(extract(stanObj,"beta[3]",permuted = FALSE))
beta3MCMC <- as.vector(extract(stanObj,"beta[4]",permuted = FALSE))
beta4MCMC <- as.vector(extract(stanObj,"beta[5]",permuted = FALSE))

sigmaUMCMC <- as.vector(extract(stanObj,"sigmaU",permuted = FALSE))
sigmaEpsMCMC <- as.vector(extract(stanObj,"sigmaEps",permuted = FALSE))


UMCMC <- NULL
for (k in 1:numGrp)
{
   charVar <- paste("U[",as.character(k),"]",sep = "") 
   UMCMC <- rbind(UMCMC,extract(stanObj,charVar,permuted = FALSE))
}

uMCMC_asian <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("u_asian[",as.character(k),"]",sep = "") 
   uMCMC_asian <- rbind(uMCMC_asian,extract(stanObj,charVar,permuted = FALSE))
}

uMCMC_black <- NULL
for (k in 1:ncZ)
{
  charVar <- paste("u_black[",as.character(k),"]",sep = "") 
  uMCMC_black <- rbind(uMCMC_black,extract(stanObj,charVar,permuted = FALSE))
}

uMCMC_hispanic <- NULL
for (k in 1:ncZ)
{
  charVar <- paste("u_hispanic[",as.character(k),"]",sep = "") 
  uMCMC_hispanic <- rbind(uMCMC_hispanic,extract(stanObj,charVar,permuted = FALSE))
}

uMCMC_white <- NULL
for (k in 1:ncZ)
{
  charVar <- paste("u_white[",as.character(k),"]",sep = "") 
  uMCMC_white <- rbind(uMCMC_white,extract(stanObj,charVar,permuted = FALSE))
}

# Convert to parameters original scale:

beta1MCMCorig <- beta1MCMC*sd.y
beta2MCMCorig <- beta2MCMC*sd.y
beta3MCMCorig <- beta3MCMC*sd.y
beta4MCMCorig <- beta4MCMC*(sd.y/sd.x4)

beta0MCMCorig <-  mean.y + sd.y*beta0MCMC -  mean.x4*beta4MCMCorig

sigmaUMCMCorig <- sigmaUMCMC*sd.y
sigmaEpsMCMCorig <- sigmaEpsMCMC*sd.y

# Do parameters plot:

parms <- list(cbind(beta1MCMCorig,beta2MCMCorig,beta3MCMCorig,
                    sigmaUMCMCorig,sigmaEpsMCMCorig))
parNamesVal <- list(c("Black","vs. Asian"),c("Hispanic","vs. Asian"),
                      c("White","vs. Asian"),c(expression(sigma[U])),
                      c(expression(sigma[epsilon])))
summMCMC(parms,parNames = parNamesVal)

# Do fitted curves using lattice graphics:

ng <- 101
ylim.val <- c(min(yOrig),max(yOrig))
x4g <- seq(min(x4),max(x4),length = ng)
Xg <- cbind(rep(1,ng),x4g)
Zg <- ZOSull(x4g,intKnots = intKnots,range.x = range.x4)

midCurvsg <- list()
lowCurvsg <- list()
uppCurvsg <- list()

for (ip in 1:4)
{
   if (ip == 1){
     betaSplMCMC <- rbind(beta0MCMC,beta1MCMC)

     fMCMC <-   Xg%*%betaSplMCMC + Zg%*%uMCMC_asian
   } 
   if (ip == 2){
     betaSplMCMC <- rbind(beta0MCMC+beta2MCMC,beta1MCMC)
     fMCMC <-   Xg%*%betaSplMCMC + Zg%*%uMCMC_black
     
   } 
   if (ip == 3){
     betaSplMCMC <- rbind(beta0MCMC+beta3MCMC,beta1MCMC)

     fMCMC <-   Xg%*%betaSplMCMC + Zg%*%uMCMC_hispanic
     
   }
   if (ip == 4){
     betaSplMCMC <- rbind(beta0MCMC+beta4MCMC,beta1MCMC)
     fMCMC <-   Xg%*%betaSplMCMC + Zg%*%uMCMC_white
     
   }

   
   fMCMCorig <- fMCMC*sd.y + mean.y
   credLowerOrig <- apply(fMCMCorig,1,quantile,0.025)
   credUpperOrig <- apply(fMCMCorig,1,quantile,0.975)
   fhatgOrig <- apply(fMCMCorig,1,mean)

   x4gOrig <- x4g*sd.x4 + mean.x4

   midCurvsg[[ip]] <- apply(fMCMCorig,1,mean)
   lowCurvsg[[ip]] <- apply(fMCMCorig,1,quantile,0.025)
   uppCurvsg[[ip]] <- apply(fMCMCorig,1,quantile,0.975)
}

fitFig <- xyplot(yOrig~x4Orig|ethnicity,groups = idnum,
                  as.table = TRUE,data = femSBMD,
                  strip = strip.custom(par.strip.text = list(cex = 1.5)),
                  par.settings = list(layout.heights = list(strip = 1.6)),
                  scales = list(cex = 1.25),
                  xlab = list("age (years)",cex = 1.5),
                  ylab = list(expression(paste(
                  "spinal bone mineral density (g/c",m^2,")")),cex = 1.5),
                  subscripts = TRUE,
                  panel = function(x,y,subscripts,groups)
                  {
                     panel.grid() 
                     if (any(femSBMD$ethnicity[subscripts] == "Asian"))    panNum <- 1
                     if (any(femSBMD$ethnicity[subscripts] == "Black"))    panNum <- 2
                     if (any(femSBMD$ethnicity[subscripts] == "Hispanic")) panNum <- 3
                     if (any(femSBMD$ethnicity[subscripts] == "White"))    panNum <- 4

                     panel.superpose(x,y,subscripts,groups,type = "b",pch = 16)

                     panel.polygon(c(x4gOrig,rev(x4gOrig)),c(lowCurvsg[[panNum]],rev(uppCurvsg[[panNum]])),
                                   col = "palegreen",border = FALSE)

                     panel.xyplot(x4gOrig,midCurvsg[[panNum]],lwd = 2,type = "l",col = "darkgreen")     
                })
   
print(fitFig)

########## End of femSBMDbayes ##########


estFunCol <- "darkgreen"; 
varBandCol <- "palegreen"

betaSplMCMC <- beta0MCMC+beta1MCMC
ContrastMCMC1 <- Xg%*%betaSplMCMC + Zg%*%(uMCMC_black-uMCMC_asian)
credLower1 <- apply(ContrastMCMC1,1,quantile,0.025)
credUpper1 <- apply(ContrastMCMC1,1,quantile,0.975)
Contrastg1 <- apply(ContrastMCMC1,1,mean)

# Convert to original units:

ContrastgOrig1 <- sd.y*Contrastg1
credLowerOrig1 <- sd.y*credLower1 
credUpperOrig1<- sd.y*credUpper1 

par(mfrow = c(1,1),mai = c(1.02,0.9,0.82,0.42))
plot(0,0,type = "n",bty = "l",xlim = range(xgOrig),
     range(c(credLowerOrig1,credUpperOrig1)),
     main="Black-Asian Contrast Curve Contrast Curve",
     xlab = "Age",
     ylab = " Spinal Bone Mineral Density",
     cex.lab = cex.labVal,cex.axis = cex.axisVal)

polygon(c(xgOrig,rev(xgOrig)),c(credLowerOrig1,rev(credUpperOrig1)),
        col = varBandCol,border = FALSE)

lines(xgOrig,ContrastgOrig1,lwd = 2,col = estFunCol)
abline(0,0,col = "slateblue",lwd = 2)