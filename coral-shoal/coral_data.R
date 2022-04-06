################################### INFO #######################################
# - Reads in coral data from csv
# - Defines models, creates prior distributions, creates transect objects
################################################################################

get.coral.data <- function(year,shoal){
  coral.all <- read.csv('coral_data.csv')
  cols <- c('Easting', 'Northing', 'FishNet', 'Depth', 'CoralHab', 'HardCoral')
  data.frame(coral.all[(coral.all$Year%in%year)&(coral.all$Shoal==shoal), cols])
}

get.fishnet.grid <- function(dim.grid){
  mins <- apply(coral[,c('Easting', 'Northing')], 2, min)
  maxs <- apply(coral[,c('Easting', 'Northing')], 2, max)
  x.grid <- seq(mins[1], maxs[1], length.out = dim.grid+1)
  y.grid <- seq(mins[2], maxs[2], length.out = dim.grid+1)
  return(list(x.grid, y.grid))
}

get.fishnet <- function(point.x, point.y){
  paste0(min(which(x.grid[-(1)] >= point.x)), '-', min(which(y.grid[-(1)] >= point.y)))
}

Sigma.model <- function(p, h) (exp(p['theta.1']) * exp(- 3*h/exp(p['theta.2'])) + exp(p['theta.0']) * diag(dim(h)[1])) 

Mean.model <- function(p, x) {
  betas = p[grepl("beta", names(p))]
  beta.x <- as.data.frame(mapply("*", data.frame(x), betas)) ##### t(t(x) * betas) 
  x.mean <- apply(beta.x, 1, sum)
  x.mean
}

get.log.likelihood <- function(p, x, y, fishnets, s.N){
  log.likelihood.N <- c()
  
  for(ii in 1:N){
    s.i <- s.N[ii,]
    names(s.i) = colnames(h)
    l_p <- Mean.model(p, x) + as.numeric(s.i[fishnets]) # linear predictor
    pr <- exp(l_p)/(1+exp(l_p)) # probability of coral
    pr <- replace(pr, pr==1, 0.999999999) # avoid -Inf
    pr <- replace(pr, pr==0, 0.000000001)
    tmp <- sum(y * log(pr) + (n-y) * log(1 - pr)) # log of Bernoulli distribution
    log.likelihood.N <- c(log.likelihood.N, tmp)
  }
  mean(log.likelihood.N, na.rm=TRUE)
}

log.posterior <- function(p, x, y, fishnets, prior.mu, prior.cov, s.N) { 
  log.likelihood <- get.log.likelihood(p, x, y, fishnets, s.N) # take MC, integrate out s
  log.prior <- mvnfast::dmvn(p, prior.mu, prior.cov, log=T) # Bernstein-von Mises theorem
  log.posterior.val = log.likelihood + log.prior
  #if(is.infinite(log.posterior.val)) return(-1e10)
  #if(is.nan(log.posterior.val)) return(-1e10)
  log.posterior.val
}

################################# DATA #########################################

# read the data
coral <- get.coral.data(year,shoal)
rownames(coral) <- 1:nrow(coral)

# make fishnets
xy.grid <- get.fishnet.grid(dim.grid)
x.grid <- xy.grid[[1]]; y.grid <- xy.grid[[2]]

# save fishnet string in dataframe
coral[, 'new.fishnet'] <- apply(coral[,c('Easting', 'Northing')], 1, function(x) get.fishnet(x[1], x[2]))

# find easting/northing coordinate centers of fishnets
x.center = x.grid[-(1)] - 0.5*abs(x.grid[2] - x.grid[1])
y.center = y.grid[-(1)] - 0.5*abs(y.grid[2] - y.grid[1])

# centre of fishet to make h -> mapping to h
h = as.matrix(dist(expand.grid(x.center,y.center)))
colnames(h) <- rownames(h) <- apply(expand.grid(1:dim.grid,1:dim.grid), 1 , paste , collapse = "-" )

################################ PRIORS ########################################

# this takes some time, so only do if necessary i.e. something needs to change
if(refit.prior){
  d2 = coral$Depth**2
  glm.fit = glm(cbind(coral$HardCoral, 20-coral$HardCoral) ~ coral$Depth + d2, data=coral, family = binomial) 
  beta.0.glm = as.numeric(glm.fit$coefficients[1])
  beta.1.glm = as.numeric(glm.fit$coefficients[2])
  beta.2.glm = as.numeric(glm.fit$coefficients[3])
  
  x0 = prior.mu.glm = c('beta.0'=beta.0.glm, 'beta.1' = beta.1.glm, 'beta.2' = beta.2.glm, 'theta.0' = 0, 'theta.1' = 1, 'theta.2'=1)
  prior.sd.glm = c('beta.0'= 1, 'beta.1' = 0.25, 'beta.2' = 0.15, 'theta.0' = 1, 'theta.1' = 1, 'theta.2'= 8)
  prior.cov.glm = prior.sd.glm^2 * diag(length(prior.mu.glm))
  
  # fit the model, posterior becomes prior for design, generate thetas and betas
  y.data = coral$HardCoral # binomial n = 20 
  x.data = cbind(1, coral[,'Depth'],coral[,'Depth']**2)
  fishnets.data = coral$new.fishnet
  fit.prior <- optim(x0, log.posterior, control = list(fnscale = -1), method="L-BFGS-B", lower=c(-100,-100,-100, -Inf, -Inf, -Inf), 
                     upper=c(100, 100, 100, 2, 2, 11), hessian = TRUE, x=x.data, y=y.data, fishnets=fishnets.data, prior.mu=prior.mu.glm, prior.cov=prior.cov.glm)
  
  prior.mu=fit.prior[[1]]
  prior.cov=solve(-fit.prior$hessian)
  names(prior.mu) = names(prior.mu.glm)
  colnames(prior.cov) = rownames(prior.cov) = names(prior.mu.glm)
  save(glm.fit,prior.mu.glm,prior.cov.glm, fit.prior,prior.mu,prior.cov, file = "prior_fit.RData")
  
} else {
  load("prior_fit.RData")
}

# not positive definite
prior.cov = prior.cov+diag(nrow(prior.cov))*1e-7

################################## TRANSECTS ###################################

# given two centre points and an angle, find xy of all points in transect
get.transect.locations.xy <- function(point.x, point.y, angle){
  length.out <- 0.5*length.transect/meters.apart
  if(angle==0){
    # vertical points
    new.x = point.x
    new.y <- c(seq(point.y-meters.apart, by=-meters.apart, length.out=length.out), seq(point.y, by=meters.apart, length.out=length.out))
  } else if(angle==45){
    # points at 45 degree
    new.x <- c(seq(point.x-sqrt(meters.apart**2/2), by=-sqrt(meters.apart**2/2), length.out=length.out), seq(point.x, by=sqrt(meters.apart**2/2), length.out=length.out))
    new.y <- c(seq(point.y-sqrt(meters.apart**2/2), by=-sqrt(meters.apart**2/2), length.out=length.out), seq(point.y, by=sqrt(meters.apart**2/2), length.out=length.out))
  } else if(angle==90){
    # horizontal points
    new.x <- c(seq(point.x-meters.apart, by=-meters.apart, length.out=length.out), seq(point.x, by=meters.apart, length.out=length.out))
    new.y = point.y
  } else if(angle==135){
    # points at 135 degree
    new.x <- c(seq(point.x-sqrt(meters.apart**2/2), by=-sqrt(meters.apart**2/2), length.out=length.out), seq(point.x, by=sqrt(meters.apart**2/2), length.out=length.out))
    new.y <- c(seq(point.y, by=sqrt(meters.apart**2/2), length.out=length.out), seq(point.y-sqrt(meters.apart**2/2), by=-sqrt(meters.apart**2/2), length.out=length.out))
  } else{
    NA
  }
  return(cbind(new.x, new.y))
}

# given distances of all observations to proposed transect points, find interpolated depth for transect points
interpolate.depth <- function(dist.obs.to.transect.pt){
  names(dist.obs.to.transect.pt) <- 1:length(dist.obs.to.transect.pt)
  nearest.neighbours <- sort(dist.obs.to.transect.pt)[1:3]
  depth <- coral[as.numeric(names(nearest.neighbours)),'Depth']
  nearest.weight = as.numeric(1/as.numeric(nearest.neighbours)/sum(1/as.numeric(nearest.neighbours)))
  #c(nearest.neighbours, depth, sum(nearest.weight * depth))
  sum(nearest.weight * depth)
}

# populate current transects list
get.transects <- function(c.point, r=NA){
  
  c.transects <- list()
  for(i in 1:number.transects){
    # get all transect points, given current point
    transect.pts <- get.transect.locations.xy(c.point[[i]][1], c.point[[i]][2], c.point[[i]][3])
    
    # if radius, then add random noise
    if(!is.na(r)){transect.pts <- transect.pts + cbind(runif(dim(transect.pts)[1],-r[i],r[i]),runif(dim(transect.pts)[1],-r[i],r[i]))}
    
    # make sure wthin bounds (and fishnets)
    transect.pts <- transect.pts[(transect.pts[,2] < max(y.grid))&(transect.pts[,2] >min(y.grid))&(transect.pts[,1] < max(x.grid))&(transect.pts[,1] >min(x.grid)),] # keep inside domain defined
    c.transects[[i]] <- list(pts = transect.pts)
    
    # nearest neighbors to get depth value
    dist.transect.obs <- raster::pointDistance(coral[,c('Easting', 'Northing')], transect.pts, longlat=FALSE)
    interpolated.depth.transect.pt <- apply(dist.transect.obs, 2, function(x) interpolate.depth(c(x)))
    c.transects[[i]][['depth']] <- interpolated.depth.transect.pt
    
    # fishnet string
    fishnets.new <- apply(transect.pts, 1, function(x) get.fishnet(x[1],x[2]))
    c.transects[[i]][['fishnets']] <- fishnets.new
  }
  
  c.transects
}

############################### INFERENCE ######################################

laplace <- function(x0, x, y, fishnets, prior.mu, prior.cov, s.N){
  
  fit <- optim(x0, log.posterior, control = list(fnscale = -1), method="L-BFGS-B", lower=c(-100,-100, -100, -Inf, -Inf, -Inf), 
               upper=c(100, 100, 100, 2, 2, 11), hessian = TRUE, x=x, y=y, fishnets=fishnets, prior.mu=prior.mu, prior.cov=prior.cov, s.N=s.N)
  
  posterior.mu <- fit$par
  posterior.var <- solve(-fit$hessian) # hessian is negative fisher
  
  return(list(posterior.mu, posterior.var))
}

D.posterior.precision <- function(postr.cov, ...) -log(det(postr.cov))

KL.divergence <- function(postr.cov, ...) {
  args <- list(...)
  
  for(i in 1:length(args)) {
    assign(x = names(args)[i], value = args[[i]])
  }
  
  kl <- as.numeric((t(prior.mu - postr.mu)%*%solve(prior.cov)%*%(prior.mu - postr.mu) + psych::tr(solve(prior.cov)*postr.cov) + log(det(prior.cov)/det(postr.cov)) - n)/2)
  if(is.na(kl)) return(0)
  return(kl)
}

generate.y <- function(p, h, x, fishnets){
  s.i <- as.vector(mvnfast::rmvn(1, rep(0, dim(h)[1]), Sigma.model(p, h))) 
  names(s.i) = colnames(h)
  l_p <- Mean.model(p, x) + as.numeric(s.i[fishnets]) # linear predictor
  pr <- exp(l_p)/(1+exp(l_p)) # probability of coral
  y.sim <- rbinom(length(pr), 20, pr) # bionomial distribution
  y.sim
}

get.prior.draws<- function(n.monte.draws, prior.mu, prior.cov){
  prior.draws <- MASS::mvrnorm(n.monte.draws, prior.mu, prior.cov, tol=1e-15) 
  prior.draws
}

