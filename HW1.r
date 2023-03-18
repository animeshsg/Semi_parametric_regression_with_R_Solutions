library(matlib)
#Q2.1z
ng <- 101 ; xg <- seq(0,1,length = ng)
T1g <- rep(1,ng) ; T2g <- xg ;
T3g <- (xg - 0.5)*(xg - 0.5>0)
B1g <- (1 - 2*xg)*(1 - 2*xg>0)
B2g <- 1 - abs(2*xg - 1) ; B3g <- 2*T3g
par(mfrow = c(2,1))
plot(0,type = "n",xlim = c(0,1),ylim = c(0,1),xlab = "x",ylab = "",bty = "l")
lines(xg,T1g,col = 1) ; lines(xg,T2g,col = 2)
lines(xg,T3g,col = 3)
text(0.1,0.8,expression(T[1]),col = 1)
text(0.4,0.5,expression(T[2]),col = 2)
text(0.8,0.2,expression(T[3]),col = 3)
plot(0,type = "n",xlim = c(0,1),ylim = c(0,1),xlab = "x",ylab = "",bty = "l")
lines(xg,B1g,col = 4) ; lines(xg,B2g,col = 5)
lines(xg,B3g,col = 6)
text(0.1,0.9,expression(B[1]),col = 4)
text(0.4,0.9,expression(B[2]),col = 5)
text(0.9,0.6,expression(B[3]),col = 6)

#Q2.1b
B1gh<-T1g-2*T2g+2*T3g
B2gh<-2*T2g-4*T3g
B3gh<-2*T3g

all.equal(B1gh,B1g)
all.equal(B2gh,B2g)
all.equal(B3gh,B3g)
B4g<-B1gh+B2gh+B3gh
all.equal(B4g,T1g)

#Q2.1c
L<-rbind(c(1,-2,2),c(0,2,-4),c(0,0,2))
#Q2.1d
iL<-solve(L)
dL<-det(L)
print(iL)
print(dL)

#Q2.1e
par(mfrow=c(1,1))
set.seed(1) ; n <- 100 ; x <- sort(runif(100))
y <- cos(2*pi*x) + 0.2*rnorm(n)
plot(x,y,col = "dodgerblue",bty = "l")
XT <- cbind(rep(1,n),x,(2*x - 1)*(2*x - 1>0))
XB <- cbind((1 - 2*x)*(1 - 2*x>0),1 - abs(2*x - 1),
              + (2*x - 1)*(2*x - 1>0))
fitT <- lm(y~-1+XT) ; fitB <- lm(y~-1+XB)
lines(x,fitted(fitT),col = "orange",lwd = 6)
lines(x,fitted(fitB),col = "darkgreen",lwd = 2)

#2.1f
#The alternate spline bases will give the same mathematical fit is when the column functions of XB are an alternative bases of column functions of XT

#2.1g
vB<-vcov(fitB)
vT<-vcov(fitT)
summary(vB)
summary(vT)
#The linear B-spline basis has lower covariance. This is desirable because having highly correlated spline basis would not able to explain the variability in dataset and would be redundant basis'. 
#Hence having spline basis with lower covariance will be able to model and capture the variability more efficiently by individually contributing more to model complexity. 

#2.1h
B4g<-B1gh+B2gh+B3gh
all.equal(B4g,T1g)
#the plot in 2.1a captures this effect efficiently with each Basis function capturing a certain interval of X. b1 from [0,0.5], b3 from [0.5,1]

#2.2a

library(HRW); data(WarsawApts)
x <- WarsawApts$construction.date
y <- WarsawApts$areaPerMzloty
fitCubic <- lm(y ~ poly(x,3,raw = TRUE))
ng <- 101 ; xg <- seq(1.01*min(x) - 0.01*max(x),1.01*max(x) - 0.01*min(x),length = ng)
fHatCubicg <- as.vector(cbind(rep(1,ng),xg,xg^2,xg^3)%*%fitCubic$coef)
plot(x,y,col = "dodgerblue")
lines(xg,fHatCubicg,col = "darkgreen",lwd = 2)
plot(fitted(fitCubic),residuals(fitCubic),col = "dodgerblue")
abline(0,0,col = "slateblue",lwd = 2)

#2.2b
trLin <- function(x,kappa) return((x-kappa)*(x>kappa))

#2.2c
knots <- seq(min(x),max(x),length = 5)[-c(1,5)]
X <- cbind(1,x)
for (k in 1:3) X <- cbind(X,trLin(x,knots[k]))
fitTLQ <- lm(y ~ -1 + X)
Xg <- cbind(1,xg)
for (k in 1:3) Xg <- cbind(Xg,trLin(xg,knots[k]))
fHatTLQg <- as.vector(Xg%*%fitTLQ$coef)
plot(x,y,col = "dodgerblue")
lines(xg,fHatTLQg,col = "darkgreen",lwd = 2)
plot(fitted(fitTLQ),residuals(fitTLQ),col = "dodgerblue")
abline(0,0,col = "slateblue",lwd = 2)

#2.2d
knots <- seq(min(x),max(x),length = 22)[-c(1,22)]
X <- cbind(1,x)
for (k in 1:20) X <- cbind(X,trLin(x,knots[k]))
fitTLQ <- lm(y ~ -1 + X)
Xg <- cbind(1,xg)
for (k in 1:20) Xg <- cbind(Xg,trLin(xg,knots[k]))
fHatTLQg <- as.vector(Xg%*%fitTLQ$coef)
plot(x,y,col = "dodgerblue") 
lines(xg,fHatTLQg,col = "darkgreen",lwd = 2)
plot(fitted(fitTLQ),residuals(fitTLQ),col = "dodgerblue")
abline(0,0,col = "slateblue",lwd = 2)

#2.2e
fitSS <- smooth.spline(x,y,lambda=100)
fitSS

plot(x,y,col="blue")
lines(fitSS,lwd=2,col="red")
plot(fitted(fitSS),residuals(fitSS),col = "dodgerblue")
abline(0,0,col = "slateblue",lwd = 2)




