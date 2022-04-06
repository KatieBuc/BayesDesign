##https://www.kaggle.com/edsans/ssnbayes

library(SSN)
library(mvtnorm)
library(dplyr)
library(DescTools)
library(psych)
library(MASS)
library(doFuture)
library(parallel)
library(viridis)
library(ggplot2)

source('clearwater_data.R')

# greate the densely gridded prediction sites needed for windows
ss = SSN::importSSN(ssn.obj.path, predpts = "predsSystematic500")

# get all that we know, locations and covariates
cols = c('NEAR_X', 'NEAR_Y', 'slope', 'elev', 'cumdrainag','air_temp')
obs.preds.covariates = rbind(parseDFobs[,cols], parseDFpreds[,cols])

# get euclidean matrix from obs.preds.covariates to unknown gridded predictions
preds.grid.pt.data = as.SpatialPointsDataFrame(ss, data = "predsSystematic500")
preds.grid.pt.df <- preds.grid.pt.data@data
pids <- row.names(preds.grid.pt.data@coords)

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

# add NEAR_X, NEAR_Y
all(preds.grid.pt.df[,c('pid')]==rownames(preds.grid.pt.data@coords))
preds.grid.pt.df$NEAR_X = as.numeric(preds.grid.pt.data@coords[,1])
preds.grid.pt.df$NEAR_Y = as.numeric(preds.grid.pt.data@coords[,2])

# save in ssn object
ss <- putSSNdata.frame(preds.grid.pt.df, ss, Name = "predsSystematic500")
obs_df <- parseDFobs[,cols]
ss <- putSSNdata.frame(obs_df, ss, Name = "Obs")

# get obs and new dense preds, combine on common cols
obs_df <- getSSNdata.frame(ss, "Obs")
preds_df <- getSSNdata.frame(ss, "predsSystematic500")
obs_df <- obs_df[,colnames(obs_df) %in% colnames(preds_df)]
preds_df <- preds_df[,colnames(preds_df) %in% colnames(obs_df)]
all_df <- rbind(obs_df, preds_df)

########################### PLOT ############################

y.n1 = c('29073', '29074', '29075', '29049', '28603', '28604', '28605', '29039')
x.n2 = c("28848", "28563", "28564", "28565", "28566", "28567")
optimal.d <- c("167", "169", "171", "174", "183") # median
col <- 'gray'

clear <- all_df[c(y.n1, x.n2, optimal.d),]
clear.df <- SSNbayes::collapse(ss)
clear.df$addfunccol_cat <- cut(clear.df$computed_afv, 
                               breaks = seq(min(clear.df$computed_afv),
                                            max(clear.df$computed_afv),
                                            length.out=5),
                               labels = 1:4,
                               include.lowest = T)

clear$optimal.col <- NA
clear[optimal.d,c('optimal.col')]=1

set.seed(1)
p1 <- ggplot(clear.df, aes(ymin)) + 
  geom_path(aes(X1, X2, group = slot, size = addfunccol_cat), lineend = 'round', linejoin = 'round', col = col)+
  geom_point(data = clear,
             aes(x = NEAR_X, y = NEAR_Y, col = optimal.col), size = 2, show.legend = FALSE)+
  geom_text(data = clear[optimal.d,],
            aes(x = NEAR_X-100, y = NEAR_Y+700, label = pid),size = 2,position=position_jitter(width=700,height=0))+
  scale_size_manual(values = seq(0.2,2,length.out = 5))+
  scale_color_viridis(option = 'C')+
  scale_shape_manual(values = c(16))+
  ylab("Latitude") +
  xlab("Longitude")+
  coord_fixed()+
  theme_bw()+
  guides(size = F)+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=13),
        legend.text=element_text(size=13),
        legend.title=element_text(size=13),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.background =element_rect(fill='white'))

p2 <- ggplot(clear.df) + 
  geom_path(aes(X1, X2, group = slot, size = addfunccol_cat), lineend = 'round', linejoin = 'round', col = col)+
  geom_point(data = clear,
             aes(x = NEAR_X, y = NEAR_Y, col = optimal.col), size = 1.5, show.legend = FALSE)+
  geom_text(data = clear,
            aes(x = NEAR_X-400, y = NEAR_Y+200, label = c("")),size = 2)+
  scale_size_manual(values = seq(0.2,2,length.out = 5))+
  scale_color_viridis(option = 'C')+
  scale_shape_manual(values = c(16))+
  coord_fixed()+
  theme_bw()+
  guides(size = F)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background =element_rect(fill='white'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlim(1377000, 1382000)+
  ylim(1824600, 1828300)

p3 <- ggplot(clear.df) + 
  geom_path(aes(X1, X2, group = slot, size = addfunccol_cat), lineend = 'round', linejoin = 'round', col = col)+
  geom_point(data = clear,
             aes(x = NEAR_X, y = NEAR_Y, col = optimal.col), size = 1.5, show.legend = FALSE)+
  geom_text(data = clear,
            aes(x = NEAR_X-400, y = NEAR_Y+200, label = ""),size = 2)+
  scale_size_manual(values = seq(0.2,2,length.out = 5))+
  scale_color_viridis(option = 'C')+
  scale_shape_manual(values = c(16))+
  coord_fixed()+
  theme_bw()+
  guides(size = F)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background =element_rect(fill='white'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlim(1378300, 1383500)+
  ylim(1835800, 1839500)


p1 + annotation_custom(ggplotGrob(p3), xmin = 1370000, xmax = 1390000, #top left
                       ymin = 1840000, ymax = 1850000) + 
  annotation_custom(ggplotGrob(p2), xmin = 1382000, xmax = 1400000, 
                    ymin = 1820000, ymax = 1830000) #bottom right


ggsave('clearwater_network_windows_inlayed.png', width = 6, height = 4.5)


