################################ INFO ###############################
# - This script defines neighbourhoods on network (from dense set of 
#   predictions assumed generated in .ssn) and computes expected utility 
# - For a given point in neighbourhood 1, loop through and calculate 
#   utility for neighbourhood2
# - Outputs of this script form inputs for Gaussian Process emulator
# - This code is intended to run on the HPC, see .sh script
#####################################################################

library(rgdal)
library(sp)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(dplyr)
library(DescTools)
library(psych)
library(MASS)
library(parallel)
library(doFuture)
library(SSN)

n.monte.draws=1500
utility.function.name='KL.divergence'
postr.approx.function.name='laplace'

source('clearwater_data.R')
source('expected_utility.R')
source('ssn_functions.R')

############################## OUTPUTS ##############################

output.path = paste('results/windows', utility.function.name, sep='/')
if(!dir.exists(output.path)) dir.create(output.path, recursive = TRUE)

file.name=paste('sampling_network_clearwarer', pt.1, 'ndraws', n.monte.draws,'wilcoxmedian_n2.Rdata', sep='')
save.image(file=paste(output.path, file.name, sep='/'))

############################# SSN DATA ##############################

# read in the densely gridded prediction sites needed for windows
# created by generate_prediction.grid.R
ss = SSN::importSSN(ssn.obj.path, predpts = "predsSystematic500")

# get all that we know, locations and covariates
cols = c('NEAR_X', 'NEAR_Y', 'slope', 'elev', 'cumdrainag','air_temp')
obs.preds.covariates = rbind(parseDFobs[,cols], parseDFpreds[,cols])

# get Euclidean matrix from obs.preds.covariates to unknown gridded predictions
preds.grid.pt.data = as.SpatialPointsDataFrame(ss, data = "predsSystematic500")
preds.grid.pt.df <- preds.grid.pt.data@data
pids <- row.names(preds.grid.pt.data@coords)

# nearest neighbours average for elevation and temp covariate values
for(pid in pids){
  pt = preds.grid.pt.data@coords[pid,]
  coords.pt.grid = rbind(pt, obs.preds.covariates[c('NEAR_X', 'NEAR_Y')])
  nearest.neighbours <- sort(as.matrix(dist(coords.pt.grid))[-(1),1])[1:3]
  if(nearest.neighbours[1] < 0.1) {nearest.weight = c(1,0,0)
  } else {nearest.weight = as.numeric(1/nearest.neighbours/sum(1/nearest.neighbours))}
  
  for(col in c('elev', 'air_temp')){
    new.val = sum(obs.preds.covariates[names(nearest.neighbours),col] * nearest.weight)
    preds.grid.pt.df[pid,col] = new.val
  }
}

# add segment based covariate values, cumdrainag
rids <- preds.grid.pt.df$rid
cumdrainag = ss@data[rids,'CUMdrainAG']
preds.grid.pt.df[,'cumdrainag'] = (cumdrainag - mean(cumdrainag))/sd(cumdrainag)

# add segment based covariate values, slope (from file, not in original .ssn)
obs_preds_slope = read.csv('obs_preds_data.csv')
other_slope = read.csv('segments_with_slope.csv')
all_slope = unique(rbind(unique(obs_preds_slope[c('rid', 'slope')]), unique(other_slope[c('rid', 'slope')])))
rownames(all_slope) <- as.character(all_slope$rid)
slope <- all_slope[rids,'slope']
preds.grid.pt.df[,'slope'] = (slope - mean(slope))/sd(slope)
head(preds.grid.pt.df)

# save in ssn object
ss <- putSSNdata.frame(preds.grid.pt.df, ss, Name = "predsSystematic500")
obs_df <- getSSNdata.frame(n, "Obs")
ss <- putSSNdata.frame(obs_df, ss, Name = "Obs")

# get obs and new dense preds, combine on common cols
obs_df <- getSSNdata.frame(ss, "Obs")
preds_df <- getSSNdata.frame(ss, "predsSystematic500")
obs_df <- obs_df[,colnames(obs_df) %in% colnames(preds_df)]
preds_df <- preds_df[,colnames(preds_df) %in% colnames(obs_df)]
all_df <- rbind(obs_df, preds_df)

###################### GENERATE PREDICTIONS #########################

output.path = paste('results/windows', utility.function.name, sep='/')
if(!dir.exists(output.path)) dir.create(output.path, recursive = TRUE)

file.name=paste('sampling_network_clearwarer', pt.1, 'ndraws', n.monte.draws,'wilcoxmedian_n2.Rdata', sep='')
save.image(file=paste(output.path, file.name, sep='/'))


utility.function <- match.fun(utility.function.name)
postr.function <- match.fun(postr.approx.function.name)

generate.y.with.preds <- function(ss, p){
  
  sim.out <- SimulateOnSSN_fixed(ss, ObsSimDF = obs_df, PredSimDF = preds_df,  PredID = "predsSystematic500",
                                 formula = ~ slope + elev + cumdrainag + air_temp,  coefficients = p[grepl("beta", names(p))], 
                                 CorModels = cov.models, use.nugget = TRUE,
                                 CorParms = c(exp(p['theta.1']), exp(p['theta.2']), exp(p['theta.0'])), addfunccol = NULL)
  
  sim.ssn <- sim.out$ssn.object
  simDFobs <- getSSNdata.frame(sim.ssn, "Obs")
  simDFpreds <- getSSNdata.frame(sim.ssn, "predsSystematic500")
  
  y.sim <- c(simDFobs$Sim_Values, simDFpreds$Sim_Values)
  names(y.sim) <- all_df$pid
  y.sim
}

get.dist.matrix <- function(ssn.obj.path, netID = 2, PredID = "predsSystematic500"){
  
  # get distance matrix, assume created
  obs2obs <- readRDS(paste(ssn.obj.path, '/distance/obs/dist.net',netID,'.RData', sep=''))
  pred2pred <- readRDS(paste(ssn.obj.path, '/distance/', PredID, '/dist.net', netID, '.RData', sep=''))
  pred2obs <- readRDS(paste(ssn.obj.path, '/distance/', PredID, '/dist.net',netID,'.a.RData', sep=''))
  obs2pred <- readRDS(paste(ssn.obj.path, '/distance/', PredID, '/dist.net',netID,'.b.RData', sep=''))
  
  dist.junc <- cbind(rbind(obs2obs, obs2pred), rbind(pred2obs, pred2pred))
  
  # weights matrix
  # update h with * w.matrix for a tail down model
  
  h <- dist.junc + t(dist.junc)
  h
}

get.neighbours <- function(i, h, r){
  pid <- neighbourhood.pids[i]
  col.mask <- h[pid, ] < r[i]
  names(h[pid,col.mask])
}

############################## WINDOWS ##############################

#overwrite h, x global
h <- get.dist.matrix(ssn.obj.path)
x <- cbind(1, all_df[,c('slope', 'elev', 'cumdrainag','air_temp')])

# optimal design, utility is conditional on these points as existing
added.pids <- c("167", "169", "172", "174", "183")

# centre points for windows
neighbourhood.pids <- c("163","176","193")

# radius to search for neighbours
r= c(2000, 2000, 2000)
n2 <- get.neighbours(2, h, r)

ut1 <- c()
prior.draws <- get.prior.draws(n.monte.draws, prior.mu, prior.sd, params)
all.res <- list()

# outer loop across n1 done in shell script....
for(pt.2 in n2){
  print(c(pt.1, pt.2))
  
  c.point <- c(pt.1, pt.2)
  
  u1 <- get.expected.utility(utility.function, postr.function, prior.draws, ss, c.point, n.monte.draws,
                                                 generate.y.function.name = 'generate.y.with.preds')
  ut1 <- c(ut1, mean(u1[[1]], na.rm=TRUE))
  all.res <- c(all.res, list(u1))
  
  save.image(file=paste(output.path, file.name, sep='/'))
}

############################## CHECK ##############################

#get.neighbours will return all points within the radius of neighbourhood.pids[i]
#n1 <- get.neighbours(1, h, r)

#make sure it's an actual line segment, that is,
#n1 <- c('29073', '29074', '29075', '29049', '28603', '28604', '28605', '29039')

#plot(as.SpatialLines(n))
#points(as.SpatialPoints(n)[added.pids], col='red')

# plot(as.SpatialLines(n))
# points(as.SpatialPoints(ss, data='predsSystematic500')[n1[-(1)]], col='blue')
# 
# subset <- n2[-(1)]
# text(as.numeric(as.SpatialPoints(ss, data='predsSystematic500')[subset]$coords.x1+500),
#      as.numeric(as.SpatialPoints(ss, data='predsSystematic500')[subset]$coords.x2+200),
#      labels=names(as.SpatialPoints(ss, data='predsSystematic500')[subset]$coords.x1), cex=0.7)
# 
# points(as.SpatialPoints(ss, data='predsSystematic500')[n2[-(1)]], col='blue')
# points(as.SpatialPoints(ss, data='predsSystematic500')[n3[-(1)]], col='blue')
