############################## INFO ####################################
# - This file is run a the beginning of main files: stochastic_greedy_exchange
#   and search_windows_utility_grid
# - Reads .ssn object, parses covariates, creates prior distributions
##########################################################################

date.str = '2013-7-01'
iter.mcmc = 50000
rwsd = c(beta.0=0.75, beta.1=0.2, beta.2=0.22, beta.3=0.2, beta.4=0.22, theta.0=0.2, theta.1=0.12, theta.2=0.75)

cov.models = c("Exponential.taildown")
params = c(intercept='beta.0', slope='beta.1', elev='beta.2', cumdrainag='beta.3', air_temp='beta.4', nugget='theta.0', psill='theta.1', range='theta.2')

############################## CLEARWATER DATA ##############################

# import SpatialStreamNetwork object
ssn.obj.path = 'clearwater.ssn'
n <- importSSN(ssn.obj.path, predpts = "preds")
dist.junc <- readRDS(paste(ssn.obj.path, '/distance/obs/dist.net2.RData', sep=''))

# import data.frame containing observations and covariates
merge.obj.rds.data <- function(ssn, obs.or.preds.str, date.str){
  
  covariate.cols <- c('slope', 'elev', 'cumdrainag','air_temp')
  raw.df <- getSSNdata.frame(ssn, obs.or.preds.str)
  colnames(raw.df)[tolower(colnames(raw.df)) %in% covariate.cols] <- tolower(colnames(raw.df)[tolower(colnames(raw.df)) %in% covariate.cols])
  
  clear <- readRDS(paste("clear_", tolower(obs.or.preds.str), ".RDS", sep=''))
  clear <- clear[clear$date==date.str,]
  colnames(clear)[tolower(colnames(clear)) %in% covariate.cols] <- tolower(colnames(clear)[tolower(colnames(clear)) %in% covariate.cols])
  
  # select columns needed
  cols = c('NEAR_X', 'NEAR_Y', 'slope', 'elev', 'cumdrainag','air_temp')
  if(obs.or.preds.str=='Obs') cols = c(cols, 'temp')
  clear.slice <- clear[,cols]
  
  # normalise the data
  for(col in covariate.cols){
    clear.slice[,col] = (clear.slice[,col] - mean(clear.slice[,col], na.rm=TRUE))/sqrt(var(clear.slice[,col]))
  }
  
  parse.df <- raw.df[!(names(raw.df) %in% cols[-(1:2)])] %>% left_join(clear.slice, by = c('NEAR_X', 'NEAR_Y'))
  row.names(parse.df) = raw.df$pid
  parse.df
}

parseDFobs = merge.obj.rds.data(n, "Obs", date.str)
parseDFpreds = merge.obj.rds.data(n, "preds", date.str)

# putting it back to the ssn object
n <- putSSNdata.frame(parseDFobs, n, Name = 'Obs')
n <- putSSNdata.frame(parseDFpreds, n, Name = 'preds')

# fit the ssn object, linear model
glmssn.out <- glmssn(temp ~ slope + elev + cumdrainag + air_temp, n,
                     CorModels = cov.models,
                     addfunccol = "afvArea")

# weights matrix
b.mat <- pmin(dist.junc, base::t(dist.junc))
flow.con.mat <- 1 - (b.mat > 0) * 1
addfunccol <- "afvArea"
afv <- parseDFobs[c('locID', addfunccol)]
n.sites=dim(parseDFobs)[1]

w.matrix <- sqrt(pmin(outer(afv[, addfunccol],rep(1, times = n.sites)),
                      base::t(outer(afv[, addfunccol],rep(1, times = n.sites) ))) /
                   pmax(outer(afv[, addfunccol],rep(1, times = n.sites)),
                        base::t(outer(afv[, addfunccol], rep(1, times = n.sites))))) * flow.con.mat

y.true <- parseDFobs$temp
x <- cbind(1, parseDFobs[,c('slope', 'elev', 'cumdrainag','air_temp')])
h <- (dist.junc + t(dist.junc)) #* w.matrix

############################## PRIORS ##############################

# empirical prior fit
empirical.priors <- function(glmssn.out){
  thetas = log(glmssn.out$estimates$theta) # note: on the log scale
  betas = glmssn.out$estimates$betahat
  prior.mu = c(as.vector(betas), thetas[3], thetas[1], thetas[2])
  prior.sd = c(as.vector(sqrt(diag(glmssn.out$estimates$covb))),log(c(2, 2.5, 1000))) # no std dev reported, manual correction
  names(prior.mu) <- names(prior.sd) <- params
  list(prior.mu, prior.sd)
}

generate.y <- function(n, p){
  obs_df <- getSSNdata.frame(n, "Obs")
  sim.out <- SimulateOnSSN(n, ObsSimDF = obs_df,
                           formula = ~ slope + elev + cumdrainag + air_temp,  coefficients = p[grepl("beta", names(p))], 
                           CorModels = cov.models, use.nugget = TRUE,
                           CorParms = c(exp(p['theta.1']), exp(p['theta.2']), exp(p['theta.0'])), addfunccol = "afvArea")
  
  sim.ssn <- sim.out$ssn.object
  simDFobs <- getSSNdata.frame(sim.ssn, "Obs")
  y.sim <- simDFobs$Sim_Values
  names(y.sim) <- simDFobs$pid
  y.sim
}

# read in priors
tmp <- empirical.priors(glmssn.out)
prior.mu = tmp[[1]]; prior.sd <- tmp[[2]]

# override for priors
if(exists("priors")){
  if(priors=='frequentist') {prior.sd = c(prior.sd[1:(length(prior.sd)-1)]*0.9, 0.9)}
} else{
  priors='' # if doesn't exist yet, introduce variable for saving filename
}

prior.cov = as.matrix(prior.sd^2 * diag(length(prior.sd)))

get.prior.draws<- function(n.monte.draws, prior.mu, prior.sd, params){
  prior.draws<-MASS::mvrnorm(n.monte.draws*1.1, prior.mu, prior.sd^2 * diag(length(params))) #grab more draws
  prior.draws <- prior.draws[exp(prior.draws[,length(prior.sd)])<4*max(h),][1:n.monte.draws,] #filter for max length, range parameter only
  if(dim(prior.draws)[1] < n.monte.draws){prior.draws <- rbind(prior.draws,prior.draws)[1:n.monte.draws,]}#repeat if necessary
  prior.draws
}
