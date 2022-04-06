#https://www.kaggle.com/edsans/ssnbayes

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

col <- 'gray'
optimal.d <- c("167", "169", "171", "174", "183") # median

clear <- parseDFobs
clear.df <- SSNbayes::collapse(n)
clear.df$addfunccol_cat <- cut(clear.df$computed_afv, 
                               breaks = seq(min(clear.df$computed_afv),
                                            max(clear.df$computed_afv),
                                            length.out=5),
                               labels = 1:4,
                               include.lowest = T)


# restrict clear to the design point domain
all.design.pids <- row.names(parseDFobs)[!row.names(parseDFobs)%in%c("163","176","193")]
clear <- clear[clear$pid%in%all.design.pids,]


################################ PLOT ###################################

set.seed(2010) # for random jitter on labels
clear$optimal.col <- 16
clear[optimal.d,c('optimal.col')]=17

ggplot(clear.df) + 
  geom_path(aes(X1, X2, group = slot, size = addfunccol_cat), lineend = 'round', linejoin = 'round', col = col)+
  geom_point(data = clear,
             aes(x = Xcoord, y = Ycoord, col = temp, shape=as.factor(optimal.col)), size = 2)+
  ggrepel::geom_text_repel(data = clear,
            aes(x = Xcoord-200, y = Ycoord+500, label = pid),size = 2.3,box.padding=0.1,direction="x",force=0.15)+
  scale_size_manual(values = seq(0.2,2,length.out = 5))+
  scale_color_viridis(option = 'C')+
  scale_shape_identity()+
  scale_shape_manual(name = "Design",
                     labels = c("non-optimal", "optimal"),
                     values = c(1,16)) + 
  ylab("Latitude") +
  xlab("Longitude")+
  coord_fixed()+
  theme_bw()+
  guides(size = F)+
  labs(shape="Optimal",colour = "Temperature(°C)")+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=13),
        legend.text=element_text(size=13),
        legend.title=element_text(size=13),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.background =element_rect(fill='white'))

ggsave('clearwater_network_optimal_colour.png', width = 8, height = 4)

