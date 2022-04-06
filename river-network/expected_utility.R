############################### INFO ################################
# - Defines the likelihood model (exponential tail-down covariance)
# - Laplace approximation to the posterior, as well as MCMC 
#   (Metropolis-Hastings) for comparison
# - Defines utility functions
# - Defines expected utility function
#####################################################################

cl = parallel::makeCluster(16) # number of threads, can slightly overload number of cores
doFuture::registerDoFuture()
plan(cluster, workers = cl)

########################## POSTERIOR APPROX #########################

Sigma.model <- function(p, h) (exp(p['theta.1']) * exp(- 3*h/exp(p['theta.2'])) + exp(p['theta.0']) * diag(dim(h)[1])) 

Mean.model <- function(p, x) {
  betas = p[grepl("beta", names(p))]
  beta.x <- as.data.frame(mapply("*", x, betas))
  x.mean <- apply(beta.x, 1, sum)
  x.mean
  }

log.posterior <- function(p, x, y, h) {
  log.likelihood <- dmvnorm(y, Mean.model(p, x), Sigma.model(p, h), log=T)
  log.prior <- sum(dnorm(p, prior.mu, prior.sd, log=T))
  log.likelihood + log.prior
}

# laplace aproximation to posterior
laplace <- function(x0, f, x, y, h, ...){
  
  fit <- optim(x0, f, method="L-BFGS-B", lower = prior.mu - 3*prior.sd, upper = prior.mu + 3*prior.sd, control = list(fnscale = -1, maxit=20), hessian = TRUE, x=x, y=y, h=h)
  posterior.mu <- fit$par
  posterior.var <- solve(-fit$hessian) #hessian is negative fisher
  
  return(list(posterior.mu, posterior.var))
}

# mcmc
get.proposal <- function(p)
  p + rnorm(length(p), 0, rwsd)

# accept new theta with probability 
step <- function(p, f, x, y, h) {
  
  # pick new point
  p.new <- get.proposal(p)
  
  # acceptance probability
  alpha <- min(0, f(p.new, x, y, h) - f(p, x, y, h))
  
  # accept new point with probability alpha
  if (runif(1) < exp(alpha))
    p <- p.new
  
  # returning the point
  p
}

mcmc <- function(x0, f, x, y, h, ...) {
  args <- list(...)
  
  for(i in 1:length(args)) {
    assign(x = names(args)[i], value = args[[i]])
  }
  
  res <- matrix(NA, nsteps, length(x0))
  for (k in seq_len(nsteps))
    res[k,] <- x0 <- step(x0, f, x, y, h)
  drop(res)
  
  burnin <- 5000
  samples.mcmc <- res[-(1:burnin),]
  
  posterior.mu <- apply(samples.mcmc, 2, mean)
  posterior.var <- var(samples.mcmc)
  
  # acceptance rate
  accepted.tmp <- rep(0)
  for (i in 2:dim(samples.mcmc)[1]) accepted.tmp[i] <- all(samples.mcmc[i-1,] != samples.mcmc[i,])
  accepted.ratio <- sum(accepted.tmp)/iter.mcmc
  
  # check autocorrelation
  acf.1 = acf(samples.mcmc,plot=FALSE)$acf[, , 1][2,]
  
  # check linear fit
  linear.fit.coeff <- c()
  for(ii in 1:ncol(samples.mcmc)){
    glm.df <- data.frame(cbind(samples.mcmc[,ii],1:nrow(samples.mcmc)))
    names(glm.df) = c('val', 'iter')
    glm.fit <- glm(glm.df)
    v = as.numeric(glm.fit$coefficients['iter'])
    linear.fit.coeff <- c(linear.fit.coeff, v)
  }
  
  return(list(posterior.mu, posterior.var, accepted.ratio, acf.1, linear.fit.coeff))
}

############################## UTILITY ##############################

D.posterior.precision <- function(postr.cov, ...) -log(det(postr.cov))

KL.divergence <- function(postr.cov, ...) {
  args <- list(...)
  
  for(i in 1:length(args)) {
    assign(x = names(args)[i], value = args[[i]])
  }
  
  kl <- as.numeric((t(prior.mu - postr.mu)%*%solve(prior.cov)%*%(prior.mu - postr.mu) + tr(solve(prior.cov)*postr.cov) + log(det(prior.cov)/det(postr.cov)) - n)/2)
  if(is.na(kl)) return(0)
  return(kl)
}

#################### EXPECTED UTILITY MONTE CARLO ####################

get.expected.utility <- function(utility.function, postr.function, prior.draws, n, c.point, n.monte.draws, generate.y.function.name = NA){
  
  if(is.na(generate.y.function.name)) generate.y.function.name = 'generate.y'
  generate.y.function <- match.fun(generate.y.function.name)
  
  utility.c <- c()
  postr.l <- list()
  y.data <- c()
  s <- Sys.time()
  
  # do the data stuff first, in sequence
  for(i in 1:n.monte.draws){
    
    # draw params from priors
    p <- prior.draws[i, ]
    
    # generate data (y depends on params, regenerate x)
    y.sim <- generate.y.function(n, p) 
    y.data <- rbind(y.data, y.sim)
    
  }
  
  out <- foreach(i = 1:n.monte.draws, .export=ls(.GlobalEnv), .combine='c',.packages=c('mvtnorm','MASS', 'dplyr'))%dopar%{

    # read in priors and data draws
    p <- prior.draws[i, ]
    y.sim <- y.data[i, ]
    
    # filter based on location of design
    y.sim <- y.sim[c.point]; x.subset <- x[c.point,]; h.subset <- h[c.point,c.point]
    
    tryCatch({
      
    # approx posterior
    postr <- postr.function(p, log.posterior, x.subset, y.sim, h.subset, prior.sd=prior.sd, nsteps=iter.mcmc, rwsd=rwsd)
    
    # calculate utility
    utility <- utility.function(postr[[2]], prior.mu=prior.mu, prior.cov=prior.cov, postr.mu=postr[[1]], n=length(p))
    
    # return
    list(list(utility=utility, postr=postr))
    
    }, warning=function(w){cat(i, "WARN :",conditionMessage(w), "\n");list(list(utility=0,postr=NA))}, error=function(e){cat(i, "ERROR :",conditionMessage(e), "\n");list(list(utility=0,postr=NA))},finally = {})
    
  }
  
  utility.c <-  unlist(lapply(out, function(x){x$utility}))
  postr.l <- lapply(out, function(x){x$postr})
  
  f <- Sys.time()
  time <- as.numeric(f-s)
  unts <- units(f-s)
  
  return(list(utility.c, postr.l, y.data, time, unts))
}

