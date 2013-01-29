library(R.matlab)
library(ggplot2)
library(scales)
library(grid)

theme_set(theme_grey(16))
fig_2a <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/fig_2a_R_short.mat")

PP<-fig_2a$fig.2a.R

Model<-PP[,4];
Num<-PP[,2];
Num[Num==8]<-"8  Neurons"
Num[Num==32]<-"32  Neurons"
Num[Num==64]<-"64  Neurons"
Num[Num==100]<-"100  Neurons"

N_list<-PP[,3]
Prob<-PP[,1]
Prob[Prob<1e-8]<-1e-9

fig_data_orig <- data.frame(Prob,Num,N_list,Model)
fig_data <- subset(fig_data_orig,Num!="100  Neurons")

fig_data$Num <- ordered(fig_data$Num,levels=rev(c("64  Neurons","32  Neurons","8  Neurons")))

p1<-ggplot(fig_data,aes(N_list,Prob,fill=Model))
p1<-p1+geom_bar(stat="identity",position="identity",alpha=0.75)
p1<-p1+scale_fill_continuous("Model",breaks=c(1,2),labels=c("EIF","Ising"))
p1<-p1+facet_grid(Num~.,scales="free")
p1<-p1+xlab("Population spike count")+ylab("Probability")
p1<-p1+opts(legend.position=c(.95, .95),legend.background = theme_rect(fill="gray90", size=0.4))+opts(legend.title = theme_text(size=14),legend.text = theme_text(size = 12))

##########

fig_data2 <- subset(fig_data_orig,Num=="100  Neurons")

p2<-ggplot(fig_data2,aes(N_list,Prob,group=Model,color=Model))
p2<-p2+geom_line(size=1.5)
p2<-p2+scale_color_gradient("Model",breaks=c(1,2),labels=c("EIF","Ising"))
p2<-p2+scale_y_log10(limits=c(1e-05,1),breaks=c(1e-1,1e-2,1e-3,1e-4,1e-5),labels=c(expression(10^-1),expression(10^-2),expression(10^-3),expression(10^-4),expression(10^-5)))
p2<-p2+scale_x_discrete(limits=c(0,100))
#p2<-p2+coord_cartesian(x=c(0,100),y=c(1e-04,1))
p2<-p2+xlab("Population spike count")+ylab("Probability")
p2<-p2+facet_grid(Num~.)+opts(legend.position="none")


p2lin<-ggplot(fig_data2,aes(N_list,Prob,group=Model,color=Model))
p2lin<-p2lin+geom_line()+labs(x=NULL,y=NULL)
p2lin<-p2lin+opts(legend.position="none")
p2lin<-p2lin+opts(plot.background=theme_rect(size=0.5))

########

fig_2b <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/fig_2b_R.mat")

PP_2b<-fig_2b$fig.2b.R

Prob_e<-PP_2b[,4]
Rho   <-PP_2b[,2]
N_list_2b<-PP_2b[,1]
Prob_m<-PP_2b[,3]

fig_data_2b <- data.frame(Prob_m,Prob_e,N_list_2b,Rho)

p2b<-ggplot(fig_data_2b,aes(N_list_2b,Prob_m,group=Rho,color=Rho))
p2b<-p2b+geom_line()
p2b<-p2b+scale_color_gradient("Rho",breaks=c(0.05,0.1,0.25))
p2b<-p2b+xlab("Number of neurons")+ylab("Kullback-Leibler divergence")
p2b<-p2b+opts(legend.position=c(.1, .9),legend.background = theme_rect(fill="gray90", size=0.4))+opts(legend.title = theme_text(size=14),legend.text = theme_text(size = 12),legend.key.width = unit(2.5, "cm"))


########

fig_2c <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/fig_2c_R.mat")

PP_2c<-fig_2c$fig.2c.R

Probc_e<-PP_2c[,4]
Mu   <-PP_2c[,2]
N_list_2c<-PP_2c[,1]
Probc_m<-PP_2c[,3]

fig_data_2c <- data.frame(Probc_m,Probc_e,N_list_2c,Mu)

captionc<-paste(strwrap("Varying the mean firing rate at a constant correlation value",40),collapse="\n")

p2c<-ggplot(fig_data_2c,aes(N_list_2c,Probc_m,group=Mu,color=Mu))
p2c<-p2c+geom_line()
p2c<-p2c+scale_color_gradient("Mu",breaks=c(0.1,0.2,0.3))
p2c<-p2c+xlab("Number of neurons")+ylab("Kullback-Leibler divergence")
p2c<-p2c+scale_y_continuous(breaks=c(0.5,1.0,1.5,2),labels=c("0.5","1.0","1.5","2"))
p2c<-p2c+coord_cartesian(y=c(0,3.1))
p2c<-p2c+opts(legend.position=c(.1, .9),legend.background = theme_rect(fill="gray90", size=0.4))+opts(legend.title = theme_text(size=14),legend.text = theme_text(size = 12),legend.key.width = unit(2.5, "cm"))


#########


#pdf("/Users/dleen/Dropbox/Research/Simple_HOC_draft/R/polishing.pdf",width=10,height=5)
subvp<-viewport(width=0.15,height=0.1,x=0.35,y=0.18)
#print(p2)
#print(p2lin,vp=subvp)
#dev.off()

pdf("/Users/dleen/Dropbox/Research/Simple_HOC_draft/R/fig_2/fig_2a.pdf",width=20,height=15)
grid.newpage()
#pushViewport(viewport(layout=grid.layout(4,1)))
pushViewport(viewport(layout=grid.layout(4,2)))

vplayout<-function(x,y)
	viewport(layout.pos.row=x, layout.pos.col=y)
	
#print(p1,vp=vplayout(1:3,1))
#print(p2,vp=vplayout(4,1))
#print(p2lin,vp=subvp)
print(p1,vp=vplayout(1:3,1))
print(p2,vp=vplayout(4,1))
print(p2lin,vp=subvp)
print(p2b,vp=vplayout(1:2,2))
print(p2c,vp=vplayout(3:4,2))
dev.off()