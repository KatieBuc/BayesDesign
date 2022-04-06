
source('gaussian_process.R')

label.x.n2 = c(expression(n[21]),expression(n[22]),expression(n[23]),expression(n[24]),expression(n[25]),expression(n[26]))
label.y.n1 = c(expression(n[11]),expression(n[12]),expression(n[13]),expression(n[14]),expression(n[15]),expression(n[16]),expression(n[17]),expression(n[18]))

muhat<-function(z,delta,u){
  out<-rep(0,nrow(z))
  for(j in 1:nrow(z)){
    out[j]<-sum(u*exp(-delta[1]*abs(x1-z1[j])-delta[2]*abs(x2-z2[j]))) # Manhattan/stream distance
  }
  out}

############################## PHASE 1 ##############################

# utility grid
df_check <- data.frame(cbind(x, y))
ggp1 <- ggplot(df_check, aes(x1, x2, fill= y)) +
  geom_tile()+
  scale_fill_viridis_c(name = expression(tilde(U)(X)))+  
  theme_bw()+
  scale_x_continuous(breaks = seq(0,1,length=n.x1), labels = label.x.n2) +
  scale_y_continuous(breaks = seq(0,1,length=n.x2), labels = label.y.n1)+
  xlab("Neighbourhood 2 (N2)")+ ylab("Neighbourhood 1 (N1)")+
  ggtitle("Phase I: Utility Grid")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))


############################## PHASE 2 ##############################

# smooth GP
# df_check <- data.frame(cbind(x, yhat))
# ggplot(df_check, aes(x1, x2, fill= yhat)) +
#   geom_tile()+
#   scale_fill_viridis_c()

# interpolate, new grid points
z <- expand.grid(z1=seq(0,1, length=100), z2=seq(0,1, length=100))
z1<- z$z1
z2<- z$z2
z.n <- dim(z)[1]

# apply GP interpolation
u<-as.vector(iPHI%*%y)
y.interp <- muhat(z,delta,u)

df_check <- data.frame(cbind(z1, z2, y.interp, v4=y.interp))
ggp2 <- ggplot(df_check, aes(z1, z2, fill= v4)) +
  geom_tile()+
  scale_fill_viridis_c(name = expression(f(X["*"])))+  
  theme_bw()+
  scale_x_continuous(breaks = seq(0,1,length=n.x1), labels = label.x.n2) +
  scale_y_continuous(breaks = seq(0,1,length=n.x2), labels = label.y.n1)+
  xlab("Neighbourhood 2 (N2)")+ ylab("Neighbourhood 1 (N1)")+
  ggtitle("Phase II: Utility Surface")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

############################## PHASE 3 ##############################

df_check <- data.frame(cbind(z1, z2, y.interp, v4=y.interp/max(y.interp)))
ggp3 <- ggplot(df_check,  aes(z1, z2, z= v4)) +
  geom_contour_fill(breaks = seq(0,1,by=0.02)) +
  geom_contour(color = "white", alpha=0.2, breaks = seq(0,1,by=0.02))+
  metR::geom_text_contour(aes(z = v4),breaks = seq(0,1,by=0.02), size=3, skip=0)+
  scale_fill_viridis_c(name = expression(Eff(X["*"])))+  
  theme_bw()+
  scale_x_continuous(breaks = seq(0,1,length=n.x1), labels = label.x.n2) +
  scale_y_continuous(breaks = seq(0,1,length=n.x2), labels = label.y.n1)+
  xlab("Neighbourhood 2 (N2)")+ ylab("Neighbourhood 1 (N1)")+
  ggtitle("Phase III: Design Efficiency")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))


ggpl4 <- grid.arrange(ggp1, ggp2, ggp3, ncol = 3)
ggsave(file='windows_3.png', ggpl4, width = 10, height = 3.5)

