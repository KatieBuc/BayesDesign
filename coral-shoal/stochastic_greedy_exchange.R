if(getwd() == "/home/virtualreef/katie/BayesDesign/coral-shoal"){
  library(doFuture, lib.loc="~/katie/packages")
  library(parallel, lib.loc="~/katie/packages")
  library(mvnfast, lib.loc="~/katie/packages")
  library(geodist, lib.loc="~/katie/packages")
  library(sp, lib.loc="~/katie/packages")
  library(raster, lib.loc="~/katie/packages")
  library(psych, lib.loc="~/katie/packages")
  library(MASS, lib.loc="~/katie/packages")
} else{
  library(doFuture)
  library(parallel)
  library(tidyverse)
  library(mvnfast)
  library(geodist)
  library(raster)
  library(psych)
  library(MASS)
}

cl = parallel::makeCluster(16) # number of threads, can slightly overload number of cores
doFuture::registerDoFuture()
plan(cluster, workers = cl)

#j=1
utility.function.name='KL.divergence'
n.monte.draws = 350
year = c(2013)
shoal='Barracouta East'
n=20
dim.grid = 10
meters.apart = 5
length.transect = 500
number.transects = 3
refit.prior=FALSE
N=50
B=600
K=15

source('coral_data.R')
set.seed(2*j + j)

################################ OUTPUTS #######################################

output.path = paste('results', number.transects, utility.function.name, sep='/')
if(!dir.exists(output.path)) dir.create(output.path, recursive = TRUE)
file.name=paste('optimal_design_greedy_exchange_stochastic', j,'ndraws', n.monte.draws, '_median.Rdata', sep='')

######################## EXPECTED UTILITY MONTE CARLO ##########################

get.expected.utility <- function(utility.function, prior.draws, c.point, n.monte.draws){
  
  # get transect points and interpolate
  c.transects <- get.transects(c.point)
  
  # unpack to vector form
  pts <- do.call("rbind", lapply(c.transects, `[[`, 1))
  depth <- do.call("c", lapply(c.transects, `[[`, 2))
  fishnets <- do.call("c", lapply(c.transects, `[[`, 3))
  x <- cbind(1, depth, depth**2)
  
  utility.c <- c()
  postr.l <- list()
  y.data <- list()
  s <- Sys.time()
  
  print('simulating data...')
  
  # do the data stuff first, in sequence
  for(iii in 1:n.monte.draws){
    
    # draw params from priors
    p <- prior.draws[iii, ]
    
    # generate data (y depends on params, regenerate x)
    y.data[[iii]] <- generate.y(p, h, x, fishnets)
    
  }
  
  s.N.all <- c()
  for(ii in 1:(N*n.monte.draws)){
    s.i <- as.vector(mvnfast::rmvn(1, rep(0, dim(h)[1]), Sigma.model(p, h)))
    s.N.all <- rbind(s.N.all, s.i)
  }
  print('finding utility in parallel...')
  
  out <- foreach(i = 1:n.monte.draws, .export=ls(.GlobalEnv), .combine='c',.packages=c('mvnfast', 'psych'))%dopar%{
    
    # read in priors and data draws
    p <- prior.draws[i, ]
    y.sim <- y.data[[i]]
    s.N <- s.N.all[((i-1)*N+1):(i*N),]
    
    tryCatch({
      
      # approx posterior
      postr <- laplace(p, x, y.sim, fishnets, prior.mu=prior.mu, prior.cov=prior.cov, s.N=s.N)
      
      # calculate utility
      utility <- utility.function(postr[[2]], prior.mu=prior.mu, prior.cov=prior.cov, postr.mu=postr[[1]], n=length(p))
      
      # return
      list(list(utility=utility, postr=postr))
      
    }, error=function(e){cat(i, "ERROR :",conditionMessage(e), "\n");list(list(utility=0,postr=NA))},finally = {})
    #, warning=function(w){cat(i, "WARN :",conditionMessage(w), "\n");list(list(utility=0,postr=NA))},
  }
  
  utility.c <-  unlist(lapply(out, function(x){x$utility}))
  postr.l <- lapply(out, function(x){x$postr})
  
  f <- Sys.time()
  time <- as.numeric(f-s)
  unts <- units(f-s)
  
  return(list(utility.c, postr.l, time, unts))
}

############################## GREEDY EXCHANGE #################################

utility.function <- match.fun(utility.function.name)

# random start deign
angles = c(0, 45, 90, 135)
x.design.space = x.center[c(3,6,9)]
y.design.space = y.center[c(3,6,9)]

start.alpha = sample(angles, 3, replace=TRUE)
start.x <- sample(x.design.space, 3, replace=FALSE)
start.y <- sample(y.design.space, 3, replace=FALSE)

# initialise design
d <- d.init <- list(c(start.x[1], start.y[1], start.alpha[1]), c(start.x[2], start.y[2], start.alpha[2]), c(start.x[3], start.y[3], start.alpha[3]))
all.res <- list()

# initial utility (saved as a distribution)
prior.draws <- get.prior.draws(B, prior.mu, prior.cov) 
res <- get.expected.utility(utility.function, prior.draws, d.init, B)
all.res <- c(all.res, res)
U.best.dist <- res[[1]] 

# utility results
U.utility = c(median(U.best.dist, na.rm=TRUE))
U.design = list(d)

# the utility that's been accepted/not accepted
trace = c() 
trace.proposals <- c()
acceptance <- c()

for(kk in 1:K){
  
  # this loop takes about 60h I think (now it should take 8h)
  for(transect.num in 1:number.transects){
    other.transect.nums <- (1:number.transects)[-(transect.num)]
    
    for(coord.num in 1:3){ #x, y, alpha
      if(coord.num==1) not.d = setdiff(x.design.space, d[[transect.num]][coord.num])
      if(coord.num==2) not.d = setdiff(y.design.space, d[[transect.num]][coord.num])
      if(coord.num==3) not.d = setdiff(angles, d[[transect.num]][coord.num])
      
      for(jj in 1:length(not.d)){
        
        # take a clean copy of d, and swap the ith element in d for the jth element in not.d
        d.swap <- d
        d.swap[[transect.num]][coord.num] = not.d[jj]

        # if swapped coordinates make for an x,y duplocate with existing design, next
        test_condition = all(as.numeric(d.swap[[transect.num]][1:2]) == as.numeric(d[[other.transect.nums[1]]][1:2]))|
                      all(as.numeric(d.swap[[transect.num]][1:2]) == as.numeric(d[[other.transect.nums[2]]][1:2]))
          
        if (test_condition) {
          next
        }
        
        # calculate utility for new design
        prior.draws <- get.prior.draws(n.monte.draws, prior.mu, prior.cov) 
        res <- get.expected.utility(utility.function, prior.draws, d.swap, n.monte.draws)
        all.res <- c(all.res, res)
        
        U <- median(res[[1]], na.rm=TRUE)
        
        # store results
        U.utility <- c(U.utility, U)
        U.design <- c(U.design, list(d.swap))
        
        # save
        save.image(file=paste(output.path, file.name.cont, sep='/'))
      }
      
      max.idx <- which.max(tail(U.utility, length(not.d)))
      d.max = U.design[[length(U.utility) - length(not.d) + max.idx]]
      
      # recalculate U.i.max.dist with higher B value here...
      prior.draws <- get.prior.draws(B, prior.mu, prior.cov)
      res <- get.expected.utility(utility.function, prior.draws, d.max, B)
      all.res <- c(all.res, res)
      U.i.max.dist <- res[[1]]
      
      out <- wilcox.test(U.i.max.dist, U.best.dist, alternative = c("greater"))
      acceptance.prob = 1 - out$p.value
      
      if (runif(1) < acceptance.prob) 
      {
        d = d.max
        U.best.dist = U.i.max.dist
      }
      
      trace.proposals <- c(trace.proposals, median(U.i.max.dist, na.rm=TRUE))
      trace <- c(trace, median(U.best.dist, na.rm=TRUE))
      acceptance <- c(acceptance, acceptance.prob)
      
      save.image(file=paste(output.path, file.name.cont, sep='/'))
    }
  }
}

