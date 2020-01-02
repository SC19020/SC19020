## ----eval=FALSE---------------------------------------------------------------
#  bootstrapR <- function(X, n) {
#    m <- length(X)
#    out <- matrix(nrow = n, ncol = m)
#    for(i in 1:n){
#      out[i,] <- sample(X, m, replace = TRUE)
#    }
#    return(out)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(SC19020)
X <- c(96,87,88,90,85)
outer <- bootstrapR(X, 5)
outer

## ----eval=FALSE---------------------------------------------------------------
#  KDER <- function(X, h, x) {
#    n <- length(X)
#    out <- (1/h)*mean(dnorm((x-X)/h))
#    return(out)
#  }

## ----eval=TRUE----------------------------------------------------------------
X <- rnorm(100,1,2)
x <- seq(-1, 1, by = 0.01)
h <- 0.1
y <- numeric(length(x))
for(i in 1:length(x)){
  y[i] <- KDER(X, h, x[i])
}
plot(x, y, type = "l")

## -----------------------------------------------------------------------------
x<-c(20,30,40,45,60)
y<-c(16,20,27,40,60)
plot(x,y,type = "l",col="blue")
title(main="An example")

## -----------------------------------------------------------------------------
math<-c(95,80,74,66,84)
engl<-c(80,77,68,92,85)
phys<-c(90,82,65,74,90)
chem<-c(88,96,70,84,80)
mydata<-data.frame(math,engl,phys,chem)
mydata

## -----------------------------------------------------------------------------
sigma <- 1;n <- 100000
u <- runif(n)
x <- sqrt(-2*(sigma^2)*log(1-u)) # F(x) = 1-exp(-x^2/(2*sigma^2)), x>=0,sigma>0
hist(x,prob =TRUE, breaks=100,xlim=c(0,6),main=expression(f(x)==x*exp(-x^2/2)))
y <- seq(0, 6, .01)
lines(y, y/(sigma^2)*exp(-y^2/(2*(sigma^2))))

## -----------------------------------------------------------------------------
sigma <- 2; n <- 100000
u <- runif(n)
x <- sqrt(-2*(sigma^2)*log(1-u)) # F(x) = 1-exp(-x^2/(2*sigma^2)), x>=0,sigma>0
hist(x,prob =TRUE,breaks=100,xlim=c(0,8), main=expression(f(x)==x/4*exp(-x^2/8)))
y <- seq(0, 8, 0.01)
lines(y, y/(sigma^2)*exp(-y^2/(2*(sigma^2))))

## -----------------------------------------------------------------------------
n <- 1e4
X1 <- rnorm(n,0,1)
X2 <- rnorm(n,3,1) 
r <- sample(c(0,1),n,replace=TRUE,prob=c(0.25,0.75))
Z <- r*X1+(1-r)*X2
hist(Z,prob = TRUE,breaks = 100,xlim = c(-4,6))
x<-seq(-4,6,0.01)
y<-0.75*dnorm(x,0,1)+0.25*dnorm(x,3,1)
lines(x,y)

## -----------------------------------------------------------------------------
n <- 1e4
p <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

for (i in 1:9){
  X1 <- rnorm(n,0,1)
  X2 <- rnorm(n,3,1) 
  r <- sample(c(0,1),size=n,replace=TRUE,prob=c(p[i],1-p[i]))
  Z<- r*X1+(1-r)*X2
  hist(Z,breaks = 100,xlim = c(-4,6))
}
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
n<-10; d<-5
X<-numeric(d)
for(i in 1:d){
  X[i]<-sqrt(rchisq(1,df=n-i+1))
  }
A <- diag(X)#生成符合题意的对角阵
B <- matrix(0,nrow=d,ncol=d)
for(i in 1:d)#生成下三角元服从正态分布的矩阵
  for(j in 1:i-1){
    B[i,j]<-rnorm(1)
    }
T <- A+B#符合题意的下三角阵
Sigma <- matrix(0,nrow=d,ncol=d)
for(i in 1:d)#输入多元正态分布的协方差矩阵
  for(j in 1:d){
    if (i==j)
   Sigma[i,j]<-4
    else
   Sigma[i,j]<-0.8}
L<-chol(Sigma)#求出Sigma矩阵的Choleski factorization
W<-t(L)%*%T%*%t(T)%*%L#利用Bartlett’s decomposition生成符合题意的样本
W

## -----------------------------------------------------------------------------
m <- 1e4
x <- runif(m, min=0, max=pi/3) #generate random numbers
theta.hat <- mean(sin(x)) * pi/3 #use sample mean to approximate the integration
print(c(theta.hat,cos(0) - cos(pi/3))) #compare the approximation and the true value

## -----------------------------------------------------------------------------
m <- 1e4
x1 <- runif(m)
x2 <- runif(m/2) #generate m/2 samples to ensure the total are the same with the 1st method
MC_standard<-numeric(m); MC_antithetic<-numeric(m/2)
MC_standard<-exp(-x1)/(1+x1^2)
MC_antithetic<-1/2*(exp(-x2)/(1+x2^2)+exp(-1+x2)/(1+(1-x2)^2))
theta.stan <- mean(MC_standard) #the standard Mento Carlo method 
se.stan<- sd(MC_standard)
theta.anti <- mean(MC_antithetic) #the Monte Carlo method with antithetic variables
se.anti <- sd(MC_antithetic)
reduc_sd <- 1-sd(MC_antithetic)/sd(MC_standard) #the reduction of sd
cbind(theta.stan,theta.anti,se.stan,se.anti,reduc_sd)

## -----------------------------------------------------------------------------
m <- 10000
theta.hat <- se <- numeric(5)
g <- function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1) }

x <- runif(m) #using f0
fg <- g(x)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)

x <- rexp(m, 1) #using f1
fg <- g(x) / exp(-x)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)

x <- rcauchy(m) #using f2
i <- c(which(x > 1), which(x < 0))
x[i] <- 2 #to catch overflow errors in g(x)
fg <- g(x) / dcauchy(x)
theta.hat[3] <- mean(fg)
se[3] <- sd(fg)

u <- runif(m) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat[4] <- mean(fg)
se[4] <- sd(fg)

u <- runif(m) #f4, inverse transform method
x <- tan(pi * u / 4)
fg <- g(x) / (4 / ((1 + x^2) * pi))
theta.hat[5] <- mean(fg)
se[5] <- sd(fg)

rbind(theta.hat, se)

## -----------------------------------------------------------------------------
m <- 10000 #number of replicates
k <- 5 #number of strata
r <- m/k #replicates per stratum
#compute the 0%, 20%, 40%, 60%, 80% and 100%percentiles
alpha <- numeric(6)
for(j in 1:6){
  alpha[j] <- -log(1-(1/5)*(j-1)*(1-exp(-1)))
}
x <- matrix(nrow=2000,ncol=5)
g <- function(x)(1-exp(-1))/(5*(1+x^2))*(x>0)*(x<1)
est <- numeric(5)
#the inverse transform method to generate samples 
for (j in 1:k){
  u <- runif(r)
  x[,j] <- - log(1 - (1/5)*(u+j-1) * (1 - exp(-1)))
  est[j] <- mean(g(x[,j]))
  }
theta.hat<-sum(est); se<-sd(g(x))
cbind(theta.hat,se)

## -----------------------------------------------------------------------------
set.seed(123)
covp <- function(m){
  n <- 20; alpha <- 0.05
  LCL_mean <- numeric(m)
  UCL_mean <- numeric(m)
  UCL_var_chisq <- UCL_var_normal <- numeric(m)
  for(i in 1:m){
    x <- rchisq(n,2) #samples from chisq with mean=2, var=4
    #the lower and upper confidence limit of mean
    LCL_mean[i] <- mean(x)-sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
    UCL_mean[i] <- mean(x)+sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
    #the one-side upper confidence limit of variance
    UCL_var_chisq[i] <- (n-1)*var(x)/qchisq(alpha, df=n-1)
    y <- rnorm(n) #samples from norm with mean=0, var=1
    #the one-side upper confidence limit of variance
    UCL_var_normal[i] <- (n-1)*var(y)/qchisq(alpha, df=n-1)
  }
  #coverage probability of mean from chisq distribution
  covp_mean_chisq <- mean(LCL_mean<2&2<UCL_mean)
  #coverage probability of variance from chisq distribution
  covp_var_chisq <- mean(4<UCL_var_chisq)
  #coverage probability of variance from normal distribution
  covp_var_normal <- mean(1<UCL_var_normal)
  cbind(covp_mean_chisq,covp_var_chisq,covp_var_normal)
}
m<-1e3; covp(m)
m<-1e4; covp(m)
m<-1e5; covp(m)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 50 #sample sizes
m <- 1000 #replicates to compute a sample quantile once
k <- 100 #replicates of Mento Carlo experiment
alpha <- c(0.025,0.05,0.95,0.975) #different probabilities of quantiles
skew <- numeric(m)
se_quantile.hat <- quantile.appro <- quantile.hat <- numeric(4)
quantile <- matrix(nrow=k,ncol=4)
#we compute k sample quantiles of every probability
for(i in 1:k){
  for(j in 1:m){
    x <- rnorm(n)
    skew[j] <- mean((x-mean(x))^3)/(mean((x-mean(x))^2))^(3/2)  #sample skewness
    }
  quantile[i,] <- quantile(skew, probs=alpha) #a sample quantiles of skewness
}
#use the mean of sample quantiles to estimate the sample quantiles
quantile.hat <- colMeans(quantile)
for(i in 1:4){
#quantiles of the large sample approximation
  quantile.appro[i] <- qnorm(alpha[i],0,sqrt(6/n))
#the standard error of the estimates
  se_quantile.hat[i] <- sqrt((alpha[i]*(1-alpha[i])/(n*(dnorm(quantile.appro[i],0,sqrt(6/n))^2))))
}
cbind(alpha,quantile.hat,quantile.appro,se_quantile.hat)   

## -----------------------------------------------------------------------------
#skewness test of normality against symmetric Beta distribution
alpha<-0.05; n<-100
a<-seq(1,50,1)
# critical value for the skewness test
cv<-qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
# the function to compute skewness of samples
sk <- function(x) mean((x-mean(x))^3)/(mean((x-mean(x))^2))^(3/2)
pwr <- numeric(length(a))
m <- 10000 # num. repl. each sim.
for (i in 1:length(a)) {
  sktests<-numeric(m)
  for (j in 1:m) {
    x <-  rbeta(n,a[i],a[i]) # test decision is 1 (reject) or 0
    sktests[j] <- as.integer(abs(sk(x)) >= cv) 
    }
  pwr[i] <- mean(sktests) # proportion rejected, that is the power of test
}
plot(a,pwr,xlab = "parameter of Beta distribution",ylab = "power of test")

## -----------------------------------------------------------------------------
# skewness test of normality against heavy-tailed symmetric alternatives such as t distribution
alpha<-0.05; n<-100
a<-seq(1,30,1)
cv<-qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
sk <- function(x) mean((x-mean(x))^3)/(mean((x-mean(x))^2))^(3/2)
pwr <- numeric(length(a))
m <- 10000 #num. repl. each sim.
for (i in 1:length(a)) {
  sktests<-numeric(m)
  for (j in 1:m) {
    x <- rt(n,a[i]) #test decision is 1 (reject) or 0
    sktests[j] <- as.integer(abs(sk(x)) >= cv)
    }
  pwr[i] <- mean(sktests) #proportion rejected
}
plot(a,pwr,xlab = "df of t distribution",ylab = "power of test")

## -----------------------------------------------------------------------------
#H0:mu=1
m <- 1e4; n <- 1000
set.seed(123)
p.valx <- p.valy <- p.valz <- numeric(m)
for(i in 1:m){
  # generate samples from $\chi^{2}(1)$ and compute p-value
  x <- rchisq(n,1)
  p.valx[i] <- 2*(1-pt(abs(sqrt(n)*(mean(x)-1)/sd(x)),n-1))
  # generate samples from uniform(0,2) and compute p-value
  y<-runif(n,0,2)
  p.valy[i] <- 2*(1-pt(abs(sqrt(n)*(mean(y)-1)/sd(y)),n-1))
  # generate samples from Exp(1) and compute p-value
  z<-rexp(n,1)
  p.valz[i] <- 2*(1-pt(abs(sqrt(n)*(mean(z)-1)/sd(z)),n-1))
}
#print the t1es
print(c(mean(p.valx<=0.05),mean(p.valy<=0.05),mean(p.valz<=0.05)))

## -----------------------------------------------------------------------------
set.seed(12345)
library(boot)
data(scor, package = "bootstrap")
scor<-as.matrix(scor)
# display the scatter plots for each pair of test scores
pairs(scor)
# compute the correlations
cor(scor)
# bootstrap estimates of the correlations and their standard errors
b.cor <- function(x,i) cor(x[i,1],x[i,2])
obj <- boot(data=scor,statistic=b.cor,R=10000)
round(c(rho12=obj$t0,se=sd(obj$t)),3)

b.cor <- function(x,i) cor(x[i,3],x[i,4])
obj <- boot(data=scor,statistic=b.cor,R=10000)
round(c(rho34=obj$t0,se=sd(obj$t)),3)

b.cor <- function(x,i) cor(x[i,3],x[i,5])
obj <- boot(data=scor,statistic=b.cor,R=10000)
round(c(rho35=obj$t0,se=sd(obj$t)),3)

b.cor <- function(x,i) cor(x[i,4],x[i,5])
obj <- boot(data=scor,statistic=b.cor,R=10000)
round(c(rho45=obj$t0,se=sd(obj$t)),3)

## -----------------------------------------------------------------------------
set.seed(12345)
# the package to compute sample skewness and bootstrap
library(moments); library(boot)
n<-1e2; m<-1e3
boot.skew <- function(x,i) skewness(x[i])
ci_norm<-ci_chisq<-matrix(NA,m,2)

# mento carlo methods to compute coverage probability
for(i in 1:m){
  # generate original sample
  U<-rnorm(n)
  V<-rchisq(n,5)
  # bootstrap to compute skewness
  out1 <- boot(U,statistic=boot.skew, R = 999)
  out2 <- boot(V,statistic=boot.skew, R = 999)
  # compute the confidence interval
  ci_norm[i,] <- boot.ci(out1,type="norm")$norm[2:3]
  ci_chisq[i,] <- boot.ci(out2,type="norm")$norm[2:3]
}
# print the results
cbind(cp_norm =mean(ci_norm[,1]<=0&ci_norm[,2]>=0),prop_left=mean(0<ci_norm[,1]),prop_right=mean(0>ci_norm[,2]))

cbind(cp_chisq =mean(ci_chisq[,1]<=sqrt(8/5)&ci_chisq[,2]>=sqrt(8/5)),prop_left=mean(sqrt(8/5)<ci_chisq[,1]),prop_right=mean(sqrt(8/5)>ci_chisq[,2]))

## -----------------------------------------------------------------------------
data(scor, package = "bootstrap")
lambda.hat <- eigen(cov(scor))$values
theta.hat <- max(lambda.hat)/sum(lambda.hat)
n <- nrow(scor)
theta.j <- numeric(n)
#leave one out
for (i in 1:n) {
  x <- scor [-i,]
  lambda <- eigen(cov(x))$values
  theta.j[i] <- max(lambda)/sum(lambda)
}
bias.jack <- (n-1)*(mean(theta.j)-theta.hat)
se.jack <- sqrt((n-1)^2/n*var(theta.j))
#print bias and se
round(c(bias.jack=bias.jack, se.jack=se.jack),5)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)
#for n-fold cross validation
for (i in 1:n) {
  y <- magnetic[-i]
  x <- chemical[-i]
#Fit the model(s) using only the n−1 observations in the training set
#Compute the predicted response for the test point
#Compute the prediction error
  
  M1 <- lm(y ~ x) #linear model
  yhat1 <- M1$coef[1]+M1$coef[2]*chemical[i]
  e1[i] <- magnetic[i] - yhat1
  
  M2 <- lm(y ~ x + I(x^2)) #quadratic polynomial model
  yhat2 <- M2$coef[1]+M2$coef[2]*chemical[i]+M2$coef[3]*chemical[i]^2
  e2[i] <- magnetic[i] - yhat2
  
  M3 <- lm(log(y) ~ x) #exponential model
  logyhat3 <- M3$coef[1]+M3$coef[2]*chemical[i]
  yhat3 <- exp(logyhat3)
  e3[i] <- magnetic[i] - yhat3

  M4 <- lm(y ~ x + I(x^2)+I(x^3)) #cubic polynomial model
  yhat4 <- M4$coef[1]+M4$coef[2]*chemical[i]+M4$coef[3]*chemical[i]^2+M4$coef[4]*chemical[i]^3
  e4[i] <- magnetic[i] - yhat4
}
#calculate the average squared prediction error
 c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
#print the adjusted R squares
 c(summary(M1)$adj.r.squared, summary(M2)$adj.r.squared, summary(M3)$adj.r.squared, summary(M4)$adj.r.squared)

## -----------------------------------------------------------------------------
set.seed(12345)
## count5 function
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X)) 
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
## generate samples with unequal sizes and same variance
n1 <- 20; n2 <- 30; N <- n1 + n2
mu1 <- mu2 <- 0; sigma1 <- sigma2 <- 1 
x <- rnorm(n1,mu1,sigma1); y <- rnorm(n2,mu2,sigma2)

#use permutation to estimate t1e with sample size of 20 and 30
z <- c(x,y); R <- 9999; t <- numeric(R)
for (i in 1:R) {
  k <- sample(1:N, size = N/2, replace = FALSE)
  x1 <- z[k]; y1 <- z[-k] # a new group of X and Y with the same size
  t[i] <- count5test(x1,y1)
}
p1 <- mean(t)# t1e

#generate another group of samples with sample size of 20 and 50
n1 <- 20; n2 <- 50; N <- n1+n2
mu1 <- mu2 <- 0; sigma1 <- sigma2 <- 1
x <- rnorm(n1,mu1,sigma1); y <- rnorm(n2,mu2,sigma2)

#use permutation to estimate t1e 
z <- c(x,y); R <- 9999; t <- numeric(R)
for (i in 1:R) {
  k <- sample(1:N, size = N/2, replace = FALSE)
  x1 <- z[k]; y1 <- z[-k] # a new group of X and Y with the same size
  t[i] <- count5test(x1,y1)
}
p2 <- mean(t)# t1e

# use Mento carlo method to estimate the empirical t1e with sample size of 20 and 50
p <- numeric(100)
for(j in 1:100){
  x <- rnorm(n1,mu1,sigma1); y <- rnorm(n2,mu2,sigma2)
  
  z <- c(x,y); R <- 9999; t <- numeric(R)
  for (i in 1:R) {
    k <- sample(1:N, size = N/2, replace = FALSE)
    x1 <- z[k]; y1 <- z[-k]
    t[i] <- count5test(x1,y1)
  }
  p[j] <- mean(t)
}
p3 <- mean(p)# empirical t1e

# print the result
round(c(p1=p1,p2=p2,p3=p3),3)

## -----------------------------------------------------------------------------
set.seed(12345)
library(mvtnorm); library(Ball); library(boot)
seed <- 12345; m <- 1e2
mean <- c(0,0); sigma <- matrix(c(1,0,0,1),nrow =2)
## the statistic of distance correlation test
DCOR <- function(z,ix,dims) {
  p <- dims[1]
  q <- dims[2]
  d <- p + q
  x <- as.matrix(z[ , 1:p])
  y <- as.matrix(z[ix, -(1:p)])
  n <- nrow(x)
  m <- nrow(y)
  if (n != m || n < 2) stop("Sample sizes must agree")
  if (! (all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d); M <- mean(d)
    a <- sweep(d, 1, m); b <- sweep(a, 2, m)
    b + M
  }
  A <- Akl(x)
  B <- Akl(y)
  dCov <- mean(A * B)
  dVarX <- sqrt(mean(A * A))
  dVarY <- sqrt(mean(B * B))
  dCor <- dCov / sqrt(dVarX * dVarY)
  return(dCor)
}

N <- seq(10, 150, by=10)# different sample sizes
power11 <- power12 <- power21 <- power22 <- numeric(length(N))
alpha <- 0.05

## model 1
## Mento Carlo Method
for(i in 1:length(N)){
  X1 <- rmvnorm(N[i],mean,sigma); e1 <- rmvnorm(N[i],mean,sigma)
  Y1 <- X1/4+e1 #generate samples
  z <- cbind(X1,Y1); z1 <- as.matrix(z)
  
  p.values <- matrix(NA,m,2)
  
  for(j in 1:m){
    boot.obj <- boot(data = z1, statistic = DCOR, R = 99, sim ="permutation", dims = c(2, 2))# distance correlation test with permutation
    tb <- c(boot.obj$t0, boot.obj$t)
    p.values[j,1] <- mean(tb >= tb[1])
    
    p.values[j,2] <- bcov.test(X1,Y1,R=99,seed=j*seed)$p.value# ball covariance test
    }
  power11[i] <- colMeans(p.values < alpha)[1]
  power12[i] <- colMeans(p.values < alpha)[2]
}
## draw the graph about powers of two testing methods with model 1
par <- par(no.readonly = TRUE)
plot(N,power11,type = "b",xlab="n",ylab="power",main="Model 1",lty=2,col=2)
lines(N,power12,type="b",lty=2,col=3)
legend("topleft", legend = c("dcor", "ball"), col = 2:3,lty = c(2,2))


## model 2 (the same process with model 1)
for(i in 1:length(N)){
  X2 <- rmvnorm(N[i],mean,sigma); e2 <- rmvnorm(N[i],mean,sigma)
  Y2 <- (X2/4)*e2
  z <- cbind(X2,Y2); z2 <- as.matrix(z)
  
  p.values <- matrix(NA,m,2)
  for(j in 1:m){
    boot.obj <- boot(data = z2, statistic = DCOR, R = 99, sim = "permutation", dims = c(2, 2))
    tb <- c(boot.obj$t0, boot.obj$t)
    p.values[j,1] <- mean(tb>=tb[1])
    
    p.values[j,2] <- bcov.test(X2,Y2,R=99,seed=j*seed)$p.value
    }
  power21[i] <- colMeans(p.values < alpha)[1]
  power22[i] <- colMeans(p.values < alpha)[2]
}

plot(N,power21,type = "b",xlab="n",ylab="power",main="Model 2",lty=2,col=2)
lines(N,power22,type="b",lty=2,col=3)
legend("topleft", legend = c("dcor", "ball"), col = 2:3,lty = c(2,2))

## -----------------------------------------------------------------------------
set.seed(12345)
## random walk Metropolis sampler
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= ((1/2)*exp(-abs(y))) / ((1/2)*exp(-abs(x[i-1]))))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

## samples generated by random walk Metropolis sampler with proposal distributions of different variances
N <- 2000; sigma <- c(0.05, 0.5, 2, 16); x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

## plots of four chains
index <- 1:2000
refline <- c(log(2*0.025), -log(2*(1-0.975)))
for(i in 1:4){
  y <- rw.Metropolis(sigma[i], x0, N)$x
  plot(index,y,type = "l",ylab="x")
  abline(h=refline,lwd=1,col=2)
  title(paste("sigma=",sigma[i]))
}

## compare the quantiles
alpha <- seq(0.05,0.95,by=0.05)
Q <- c(log(2*alpha[1:10]), -log(2*(1-alpha[11:19])))# the true quantiles of Lapalace disribution
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
mc <- rw[401:N, ]
Qrw <- apply(mc,2,function(x) quantile(x, alpha))# the sample quantiles of the four chains
qq <- data.frame(round(cbind(Q, Qrw), 3))
names(qq) <- c('True','sigma=0.05','sigma=0.5','sigma=2','sigma=16')
print(qq)

## acceptance rates
print(c(acceptRate1=1-rw1$k/N, acceptRate2=1-rw2$k/N, acceptRate3=1-rw3$k/N, acceptRate4=1-rw4$k/N))

## -----------------------------------------------------------------------------
isTRUE(exp(log(0.1))==log(exp(0.1))); isTRUE(exp(log(0.1))==0.1); isTRUE(log(exp(0.1))==0.1)

isTRUE(all.equal(exp(log(0.1)),log(exp(0.1)))); isTRUE(all.equal(exp(log(0.1)),0.1)); isTRUE(all.equal(log(exp(0.1)),0.1))

## -----------------------------------------------------------------------------
## answer of 11.4
k <- c(4:25,100,500,1000)# different values of parameter k (df of t distribution)
solution1 <- numeric(length(k))
for(i in 1:length(k)){
# the one-dimensional nonlinear function 
  f <- function(a) {pt(sqrt(a^2*(k[i]-1)/(k[i]-a^2)),k[i]-1)-pt(sqrt(a^2*k[i]/(k[i]+1-a^2)),k[i])}
  solution1[i] <- uniroot(f,lower=0+1e-6,upper=sqrt(k[i])-1e-6)$root
}
round(solution1,4)

## answer of 11.5
k <- c(4:25,100,500,1000)
solution2 <- numeric(length(k))
for(i in 1:length(k)){
# the upper value of integrate function
  ck1 <- function(a) sqrt(a^2*(k[i]-1)/(k[i]-a^2))
# the left function of the equation 
  g1 <- function(a) 2*exp(lgamma(k[i]/2)-log(sqrt(pi*(k[i]-1)))-lgamma((k[i]-1)/2))*(integrate(function(u) {(1+u^2/(k[i]-1))^(-k[i]/2)},lower = 1e-6,upper = ck1(a)-1e-6)$value)
  ck2 <- function(a) sqrt(a^2*(k[i])/(k[i]+1-a^2))
  g2 <- function(a) 2*exp(lgamma((k[i]+1)/2)-log(sqrt(pi*k[i]))-lgamma(k[i]/2))*(integrate(function(u) {(1+u^2/(k[i]))^(-(k[i]+1)/2)},lower = 1e-6,upper = ck2(a)-1e-6)$value)

  # the one-dimensional nonlinear function
  f <- function(a) {g1(a)-g2(a)}
  solution2[i] <- uniroot(f,lower=0+1e-2,upper=sqrt(k[i])-1e-2)$root
}
round(solution2,4)


## -----------------------------------------------------------------------------
library("stats4")
## the log likelihood function
ll <- function(p,nAA,nBB) -{nAA*log(p[1]^2)+nBB*log(p[2]^2)+41*log((1-p[1]-p[2])^2)+(28-nAA)*log(2*p[1]*(1-p[1]-p[2]))+(24-nBB)*log(2*p[2]*(1-p[1]-p[2]))+70*log(2*p[1]*p[2])}
tol <- 1e-10# convergence condition
N <- 100# max numbers of iterations
lmlv <- numeric(N)# log-maximum likelihood values
p <- c(0.1,0.1)# initial estimate of the target parameter
par <- matrix(nrow= N ,ncol = 2); par[1,] <- p
for(i in 1:N){
  # the E step
  nAA <- 28*p[1]^2/(p[1]^2+2*p[1]*(1-p[1]-p[2]))
  nBB <- 24*p[2]^2/(p[2]^2+2*p[2]*(1-p[1]-p[2]))
  lmlv[i]<-ll(p,nAA,nBB)
  # the M step
  par[i+1,] <- optim(par=p,fn=ll,nAA=nAA,nBB=nBB)$par
  p <- par[i+1,]
  if (sum(abs(par[i+1,]-par[i,])) < tol) break
}
print(round(c(p_mle=p[1],q_mle=p[2]),5))# MLE of p and q
index <- c(1:i)
plot(index,-lmlv[index],type="b",ylab="log-maximum likelihood values",xlab="iteration numbers")

## -----------------------------------------------------------------------------
attach(mtcars)
formulas <- list( mpg ~ disp, mpg ~ I(1 / disp), mpg ~ disp + wt, mpg ~ I(1 / disp) + wt)

## for loops
for(i in seq_along(formulas)){
  print(lm(formulas[[i]]))
}

## lapply
lapply(formulas, lm)

detach(mtcars)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
  })
## with an anonymous function
# usefor loops
for(i in seq_along(bootstraps)){
  print(lm(bootstraps[[i]]$mpg ~ bootstraps[[i]]$disp))
}
# use lapply
lapply(bootstraps, lm, formula = mpg ~ disp)

## use lapply without an anonymous function
fit <- function(x){
  lm(mpg ~ disp,data=x)
}
lapply(bootstraps, fit)

## -----------------------------------------------------------------------------
attach(mtcars)
rsq <- function(mod) summary(mod)$r.squared
## R2 in ex3
formulas <- list(mpg ~ disp, mpg ~ I(1 / disp), mpg ~ disp + wt, mpg ~ I(1 / disp) + wt)

# for loops
for(i in seq_along(formulas)){
  print(rsq(lm(formulas[[i]])))
}
# lapply
lapply(lapply(formulas,lm),rsq)

## R2 in ex4
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
  })

# for loops
for(i in seq_along(bootstraps)){
  print(rsq(lm(bootstraps[[i]]$mpg ~ bootstraps[[i]]$disp)))
}
# lapply
lapply(lapply(bootstraps, lm, formula = mpg ~ disp),rsq)
detach(mtcars)

## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
## use sapply() and an anonymous function
pvalue <- function(test) test$p.value
sapply(trials, pvalue)

## use [[ directly
sapply(trials, '[[', 3)

## -----------------------------------------------------------------------------
library(parallel)
cl <- makeCluster(getOption("cl.cores", 2))

bootstraps <- lapply(1:100000, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, c(1,3)]
})
## for sapply()
system.time(sapply(bootstraps, lm))
## for parSapply
parSapply(cl, 1:10, sqrt)
system.time(parSapply(cl,bootstraps, lm))

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
set.seed(135)

## random walk Metropolis sampler in R
rwMetropolisR <- function(sigma, x0, N) {
  x <- numeric(N)
  x[0] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= ((1/2)*exp(-abs(y))) / ((1/2)*exp(-abs(x[i-1]))))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

## random walk Metropolis sampler in C++
cppFunction('List rwMetropolisC(double sigma, double x0, int N) {
            NumericVector x(N);
            x[0] = x0;
            int k = 0;
            for(int i = 1; i < N; i++){
              double u = as<double>(runif(1));
              double y = as<double>(rnorm(1, x[i-1], sigma));
              if(u <= exp(-abs(y)+abs(x[i-1]))) 
              x[i] = y;
              else {
                x[i] = x[i-1];
                k += 1;
              }
            }
            return List::create(x,k);
          }')

N <- 2000; sigma <- c(0.05, 0.5, 2, 16); x0 <- 25

# plots of four chains in C++
index <- 1:2000; refline <- c(log(2*0.025), -log(2*(1-0.975)))
for(i in 1:4){
  y <- unlist(rwMetropolisC(sigma[i], x0, N)[1])
  plot(index, y, type = "l", ylab="x")
  abline(h = refline, lwd = 1, col = 2)
  title(paste("sigma=", sigma[i]))
}

# acceptance rates in C++
acceptRateC <- numeric(4)
for(i in 1:4){
  acceptRateC[i] <- 1-(as.numeric(rwMetropolisR(sigma[i], x0, N)[2]))/N
}
print(c(acceptRateC1 = acceptRateC[1], acceptRateC2 = acceptRateC[2], acceptRateC3 = acceptRateC[3], acceptRateC4 = acceptRateC[4]))

# qqplot
xR <- rwMetropolisR(2, 25, 10000)$x[1000:10000]
xC <- unlist(rwMetropolisC(2, 25, 10000)[1])[1000:10000]
qqplot(xR, xC, main="qqplot", xlab = "rwMetropolisR", ylab = "rwMetropolisC")
abline(a = 0, b = 1)

# computation time
ts <- microbenchmark(xR = rwMetropolisR(2, 25, 10000), xC = rwMetropolisC(2, 25, 10000))
summary(ts)[, c(1,3,5,6)]

