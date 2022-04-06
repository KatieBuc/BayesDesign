#devtools::install_github("eliocamp/metR")

library(ggplot2)
library(mgcv)
library(proxy)
library(metR)

set.seed(1)

# same as search_windows_utility_grid.R
n.monte.draws=1500
utility.function.name='KL.divergence'
output.path = paste('results/windows', utility.function.name, sep='/')

# define ordered vectors of neighbourhoods
x.n2 = c("28848", "28563", "28564", "28565", "28566", "28567")
y.n1 = c('29073', '29074', '29075','29049', '28603', '28604', '28605', '29039')

get.res <- c()
for(pt.i in y.n1){
  file.name=paste('sampling_network_clearwarer', pt.i, 'ndraws', n.monte.draws,'wilcoxmedian_n2.Rdata', sep='')
  load(paste(output.path, file.name, sep='/'))
  tmp <- data.frame(row.names=n2,ut1)[x.n2,]
  get.res <- rbind(get.res, cbind(x.n2, i, tmp))
}

# not proud of this code, but it arranges utility corresponding to expand.grid
grid <- expand.grid(x.n2, y.n1)
colnames(grid) <- c('x.n2', 'y.n1')
dat <- as.data.frame(get.res)
colnames(dat) <- c( 'x.n2', 'y.n1', 'u')
dat$u <- as.numeric(dat$u)
dat$n1 <- as.factor(dat$y.n1)
dat$n2 <- as.factor(dat$x.n2)
get.ut.n <- function(n1.pt, n2.pt){
  dat[(dat$n1==n1.pt)&(dat$n2==n2.pt),]$u
}
y <- apply(grid, 1, function(x){ get.ut.n(x[2], x[1])})

# neighbourhoods can have differing number of points
n <- length(y)
n.x1 = length(x.n2)
n.x2 = length(y.n1)

# but we use a proxy for dist
x <- expand.grid(x1=seq(0,1,length=n.x1), x2=seq(0,1,length=n.x2))
x1<- x$x1
x2<- x$x2

d1 <- dist(x1,method = "manhattan") # proxy for stream dist, in 1D same as Euc anyway
d2 <- dist(x2,method = "manhattan") 

distdd <- array(0,c(2,n,n))
distdd[1,1:n,1:n] <- as.matrix(d1)
distdd[2,1:n,1:n] <- as.matrix(d2)

# cross-validation function to minimise
CV<-function(paras,yp){
  delta<-paras[1:2]
  lambda<-paras[3]
  K<-matrix(0,nrow=n,ncol=n)
  for(uu in 1:2){
    K<-K-delta[uu]*distdd[uu,,]}
  K<-exp(K)
  PHI<-K+lambda*diag(n)
  iPHI<-solve(PHI)
  S<-K%*%iPHI
  yhat<-as.vector(S%*%yp)
  out <- sum(((yp-yhat)/(1-diag(S)))^2)
  if(!is.finite(out)|is.nan(out)){out <- 200}
  out}

# start params and minimise
paras <- c(1, 1, 0.1)
opt <- optim(paras,CV,method="L-BFGS-B",yp=y,lower=c(0.0001,0.0001,0.0001),upper=c(5,5,5))


############################## CHECK ##############################

# delta<-opt$par[1:2]
# lambda<-opt$par[3]
# K<-matrix(0,nrow=n,ncol=n)
# for(uu in 1:2){
#   K<-K-delta[uu]*distdd[uu,,]}
# K<-exp(K)
# PHI<-K+lambda*diag(n)
# iPHI<-solve(PHI)
# S<-K%*%iPHI
# yhat<-as.vector(S%*%y)
# plot(y,yhat)

