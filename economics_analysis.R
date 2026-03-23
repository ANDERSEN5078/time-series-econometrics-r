
Sys.setenv(GITHUB_PAT = "ton token GitHub")
remotes::install_github("ccolonescu/PoEdata")
data("usa")



library(tseries) # for ADF unit root tests
library(dynlm)
library(nlWaldTest) # for the `nlWaldtest()` function
library(lmtest) #for `coeftest()` and `bptest()`.
library(broom) #for `glance(`) and `tidy()`
library(PoEdata) #for PoE4 datasets
library(car) #for `hccm()` robust standard errors
library(sandwich)
library(knitr) #for kable()
library(forecast)


data("usa", package="PoEdata")
usa.ts <- ts(usa, start=c(1984,1), end=c(2009,4),
             frequency=4)
Dgdp <- diff(usa.ts[,1])
Dinf <- diff(usa.ts[,"inf"])
Df <- diff(usa.ts[,"f"])
Db <- diff(usa.ts[,"b"])
usa.ts.df <- ts.union(gdp=usa.ts[,1], # package tseries
                      inf=usa.ts[,2], 
                      f=usa.ts[,3],
                      b=usa.ts[,4],
                      Dgdp,Dinf,Df,Db,
                      dframe=TRUE)

plot(usa.ts.df$gdp)
plot(usa.ts.df$Dgdp)
plot(usa.ts.df$inf)
plot(usa.ts.df$Dinf)
plot(usa.ts.df$f)
plot(usa.ts.df$Df)
plot(usa.ts.df$b)
plot(usa.ts.df$Db)


kable(head(usa.ts.df), 
      caption="Time series data frame constructed with 'ts.union'")



N <- 500
a <- 1
l <- 0.01
rho <- 0.7

set.seed(246810)
v <- ts(rnorm(N,0,1))

y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- rho*y[t-1]+v[t]
}
plot(y,type='l', ylab="rho*y[t-1]+v[t]")
abline(h=0)

y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- a+rho*y[t-1]+v[t]
}
plot(y,type='l', ylab="a+rho*y[t-1]+v[t]")
abline(h=0)

y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- a+l*time(y)[t]+rho*y[t-1]+v[t]
}
plot(y,type='l', ylab="a+l*time(y)[t]+rho*y[t-1]+v[t]")
abline(h=0)

y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- y[t-1]+v[t]
}
plot(y,type='l', ylab="y[t-1]+v[t]")
abline(h=0)

a <- 0.1
y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- a+y[t-1]+v[t]
}
plot(y,type='l', ylab="a+y[t-1]+v[t]")
abline(h=0)

y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- a+l*time(y)[t]+y[t-1]+v[t]
}
plot(y,type='l', ylab="a+l*time(y)[t]+y[t-1]+v[t]")
abline(h=0)



###Regression falacieuse

T <- 1000
set.seed(1357)
y <- ts(rep(0,T))
vy <- ts(rnorm(T))
for (t in 2:T){
  y[t] <- y[t-1]+vy[t]
}

set.seed(4365)
x <- ts(rep(0,T))
vx <- ts(rnorm(T))
for (t in 2:T){
  x[t] <- x[t-1]+vx[t]
}
y <- ts(y[300:1000])
x <- ts(x[300:1000])
ts.plot(y,x, ylab="y and x")


spurious.ols <- lm(y~x)
summary(spurious.ols)


plot(x, y, type="p", col="grey")


plot(usa.ts.df$f)
Acf(usa.ts.df$f)

adf.test(usa.ts.df$f, k=10)

plot(usa.ts.df$b)

Acf(usa.ts.df$b)


adf.test(usa.ts.df$b, k=10)

f <- usa.ts.df$f
f.dyn <- dynlm(d(f)~L(f)+L(d(f)))
tidy(f.dyn)


b <- usa.ts.df$b
b.dyn <- dynlm(d(b)~L(b)+L(d(b)))
tidy(b.dyn)


df <- diff(usa.ts.df$f)
plot(df)
Acf(df)
adf.test(df, k=2)

db <- diff(usa.ts.df$b)
plot(db)
Acf(db)
adf.test(db, k=1)

df.dyn <- dynlm(d(df)~L(df)-1)
db.dyn <- dynlm(d(db)~L(db)-1)
tidy(df.dyn)

tidy(db.dyn)
ndiffs(f)

ndiffs(b)

###Cointegration


fb.dyn <- dynlm(b~f)
ehat.fb <- resid(fb.dyn)
ndiffs(ehat.fb) #result: 1

output <- dynlm(d(ehat.fb)~L(ehat.fb)+L(d(ehat.fb))-1) #no constant
foo <- tidy(output)
foo

bfx <- as.matrix(cbind(b,f), demean=FALSE)
po.test(bfx)


###Le modele de correction d'erreur

b.ols <- dynlm(L(b)~L(f))
b1ini <- coef(b.ols)[[1]]
b2ini <- coef(b.ols)[[2]]
d.ols <- dynlm(b~L(b)+f+L(f))
aini <- 1-coef(d.ols)[[2]]
d0ini <- coef(d.ols)[[3]]
d1ini <- coef(d.ols)[[4]]
Db <- diff(b)
Df <- diff(f)
Lb <- lag(b,-1)
Lf <- lag(f,-1)
LDf <- lag(diff(f),-1)
bfset <- data.frame(ts.union(cbind(b,f,Lb,Lf,Db,Df,LDf)))
formula <- Db ~ -a*(Lb-b1-b2*Lf)+d0*Df+d1*LDf
bf.nls <- nls(formula, na.action=na.omit, data=bfset,
              start=list(a=aini, b1=b1ini, b2=b2ini, 
                         d0=d0ini, d1=d1ini))
kable(tidy(bf.nls), 
      caption="Parameter estimates in the error correction model")

ehat <- bfset$Lb-coef(bf.nls)[[2]]-coef(bf.nls)[[3]]*bfset$Lf
ehat <- ts(ehat)
ehat.adf <- dynlm(d(ehat)~L(ehat)+L(d(ehat))-1)
kable(tidy(ehat.adf), 
      caption="Stationarity test within the error correction model")



foo <- tidy(ehat.adf)
