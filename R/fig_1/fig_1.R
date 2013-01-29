library(R.matlab)
library(ggplot2)
library(scales)
library(grid)

rm(list = ls())
theme_set(theme_grey(7))

fig_1_volt <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/fig_1_volt.mat")
fig_1_bins <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/bin_times.mat")

VV<-fig_1_volt$volt.R.data
Vbins<-fig_1_bins$bin.times

Neuron = VV[,3]
Time = VV[,2]
Voltage = VV[,1]

fig_data<-data.frame(Voltage,Time,Neuron)
fig_data_short<-subset(fig_data,640<Time & Time<750)

Start = Vbins[,1]
End = Vbins[,2]
Bins = Vbins[,3]
Bins[Bins==1]<-"A"
Bins[Bins==2]<-"B"

fig_bin_data<-data.frame(Start,End,Bins)
fig_bin_short<-subset(fig_bin_data,640<End & End<750)

xrng = range(fig_data_short$Time)

#pv<-ggplot(fig_data_short,aes(Time,Voltage,color=Neuron,group=Neuron))
pv<-ggplot(fig_data_short,aes(Time,Voltage))
pv<-pv+geom_line(aes(color=Neuron,group=Neuron),size=0.3)
pv<-pv+geom_rect(aes(NULL,NULL,xmin=Start,xmax=End,ymin=-75,ymax=-10,fill=Bins),data=fig_bin_short,alpha=0.25)
pv<-pv+coord_cartesian(x=c(xrng[1],xrng[2]),y=c(-75,-10))
pv<-pv+scale_fill_manual(values=c(NULL,"#55B1F7"))
pv<-pv+opts(legend.position="none")
pv<-pv+opts(axis.text.y=theme_text(angle=90, vjust=0.3,size=5),axis.title.y=theme_text(angle=90,vjust=0.55,size=10))
pv<-pv+scale_y_discrete(breaks=c(-70,-60,-52,-20,-13),labels=c(-70,-60,-53,0,20))
pv<-pv+scale_x_discrete(breaks=NULL)
pv<-pv + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
pv<-pv + opts(panel.background = theme_rect(colour = NA)) 
pv<-pv+xlab("Time (ms)")+opts(axis.text.x = theme_text(size=5), axis.title.x=theme_text(size=10,vjust=1.4))
pv<-pv+opts(plot.margin = unit(c(0,0,-0.25,-0.15), "cm"))


pdf("/Users/dleen/Dropbox/Research/Simple_HOC_draft/R/fig_1/fig_1_newest.pdf",width=3.375,height=1.2)
print(pv)
dev.off()
