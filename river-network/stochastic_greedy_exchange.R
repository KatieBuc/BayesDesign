############################### INFO ################################
# - Runs the stochastic coordinate exchange algorithm with proposed
#   non-parametric acceptance criteria (Wilcoxon Rank Sum test)
# - Outputs .Rdata to corresponding sub directory, where
#   a separate Notebook is used to analyse results from outputs
# - This code is intended to run on the HPC, see .sh script
#####################################################################

library(SSN)
library(mvtnorm)
library(dplyr)
library(DescTools)
library(psych)
library(MASS)
library(doFuture)
library(parallel)

# uncomment when testing, these variables are defined in the .sh
# j=1
# n.monte.draws=10
# n.choose=5
# utility.function.name='KL.divergence'
# postr.approx.function.name='laplace'
# setwd("Z:/clearwater")

B = 15000
K=15

source('clearwater_data.R')
source('expected_utility.R')

############################## OUTPUTS ##############################

output.path = paste('results', n.choose, utility.function.name, sep='/')
if(!dir.exists(output.path)) dir.create(output.path, recursive = TRUE)

file.name=paste('optimal_design_stochastic_greedy_exchange', j, 'ndraws', n.monte.draws, priors, 'wilcox_median.Rdata', sep='')
save.image(file=paste(output.path, file.name, sep='/'))

##################### STOCHASTIC GREEDY EXCHANGE #####################

utility.function <- match.fun(utility.function.name)
postr.function <- match.fun(postr.approx.function.name)

# initialise design
set.seed(2*j + j)
all.design.pids <- row.names(parseDFobs)[!row.names(parseDFobs)%in%c("163","176","193")]
d <- d.init <- sample(all.design.pids, size=n.choose)

# initial utility
prior.draws <- get.prior.draws(B, prior.mu, prior.sd, params)
res <- get.expected.utility(utility.function, postr.function, prior.draws, n, sort(d.init), B)
U.best.dist <- res[[1]] # as a distribution

# utility results, for each jj
U.utility = c(median(U.best.dist, na.rm=TRUE))
U.design = list(d.init)

# the utility that's been accepted/not accepted
all.res <- list()
trace = c() 
acceptance <- c()
trace.proposals <- c()

for(kk in 1:K)
{
  
  for (i in 1:n.choose)
  {
    # store utility results for coordinate i swap
    U.jj = rep()
    D.jj = list()
    U.dist.jj = list()
    
    # find the points not in the design....
    not.d <- setdiff(all.design.pids, d)
    
    for (jj in 1:(length(all.design.pids)-n.choose))
    {
      
      # take a clean copy of d, and swap the ith element in d for the jth element in not.d
      d.swap <- d
      d.swap[i] <- not.d[jj]
      
      # calculate utility for new design
      prior.draws <- get.prior.draws(n.monte.draws, prior.mu, prior.sd, params)
      res <- get.expected.utility(utility.function, postr.function, prior.draws, n, sort(d.swap), n.monte.draws)
      all.res <- c(all.res, res)
      
      U <- median(res[[1]], na.rm=TRUE)
      
      # store results
      U.utility <- c(U.utility, U)
      U.design <- c(U.design, list(d.swap))
      
      U.jj <- c(U.jj, U)
      D.jj <- c(D.jj, list(d.swap))
      U.dist.jj <- c(U.dist.jj, list(res[[1]]))
    }
    
    # accept the best swapped design with probability...
    d.max <- D.jj[[which.max(U.jj)]]
    #U.i.max.dist <- U.dist.jj[[which.max(U.jj)]]
    
    # recalculate U.i.max.dist with higher B value here...
    prior.draws <- get.prior.draws(B, prior.mu, prior.sd, params)
    res <- get.expected.utility(utility.function, postr.function, prior.draws, n, d.max, B)
    U.i.max.dist <- res[[1]]
    
    # calculate acceptance probability
    
    # the below code is the ACE acceptance criteria
    # v = (sum((median(U.best.dist, na.rm=TRUE) - U.best.dist)**2) + sum((median(U.i.max.dist, na.rm=TRUE) - U.i.max.dist)**2))/(2*B-2)
    # xx = (median(U.i.max.dist, na.rm=TRUE) - median(U.best.dist, na.rm=TRUE))/sqrt(2 * B * v)
    # acceptance.prob = 1 - pt(-xx, 2*B-2)
    
    # the proposed acceptance criteria
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

    save.image(file=paste(output.path, file.name, sep='/'))
  }
}  

