library(R.matlab)
library(ggplot2)
library(scales)
library(grid)

rm(list = ls())
theme_set(theme_grey(7))

fig_1_filter <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Linear_filter/fig_data/filter_R.mat")
fig_1_nl <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Linear_filter/fig_data/nln_R.mat")
fig_1_rates <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Linear_filter/fig_data/firing_R.mat")

filt<-fig_1_filter$A.t

Filter = filt[,1]
Filter = 1000*Re(Filter)
Time = filt[,2]
Time = Re(Time)
Rho_vals = filt[,3]
Rho_vals = Re(Rho_vals)

Rho_vals[Rho_vals == 0.05] <- "A"
Rho_vals[Rho_vals == 0.1] <- "B"
Rho_vals[Rho_vals == 0.25] <- "C"

fig_data<-data.frame(Filter,Time,Rho_vals)
fig_data_short<-subset(fig_data,-10<Time & Time<50)

fig_data_short$Rho_vals <- ordered(fig_data_short$Rho_vals, levels = rev(c("C", "B", "A")))

pv <- ggplot(fig_data_short,aes(Time,Filter))
pv <- pv+geom_line(aes(color=Rho_vals,group=Rho_vals),size=0.3)

pv <- pv + scale_color_brewer("",palette="Blues",labels = c(expression("" * rho * " = 0.05"), expression("" * rho * " = 0.1  "), expression("" * rho * " = 0.25")))
pv <- pv + xlab("Time (ms)") + ylab(expression("Linear Filter (Hz/mV)"))
pv <- pv + opts(legend.key = theme_blank())
pv <- pv + scale_x_continuous(limits = c(-10,40), breaks = c(-10, 0, 10, 20, 30, 40))
pv <- pv + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
pv <- pv + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + opts(panel.grid.major = theme_blank()) + opts(axis.line = theme_segment())
	
pv <- pv + opts(plot.margin = unit(c(0.05,0,0,-0.2), "cm"))	
pv <- pv + opts(legend.position = c(0.7, 0.85)) + opts(legend.text = theme_text(size = 9), legend.key.width = unit(0.5,"cm"), legend.key.height = unit(0.25, "cm"))
pv <- pv + opts(axis.text.x = theme_text(size = 7))
pv <- pv + opts(axis.text.y = theme_text(size = 7))
pv <- pv + opts(axis.title.y = theme_text(angle = 90, vjust = 0.6, size = 8))
pv <- pv + opts(axis.title.x = theme_text(size = 8, vjust=0.6))

#############
#############

nonlin <- fig_1_nl$nln.R

F_nl = 1000*nonlin[,1]
Freq = 1000*nonlin[,2]
Rho_vals = nonlin[,3]

Rho_vals[Rho_vals == 0.05] <- "A"
Rho_vals[Rho_vals == 0.1] <- "B"
Rho_vals[Rho_vals == 0.25] <- "C"

fig_data_nl<-data.frame(F_nl,Freq,Rho_vals)

fig_data_nl$Rho_vals <- ordered(fig_data_nl$Rho_vals, levels = rev(c("C", "B", "A")))

pv_nl <- ggplot(fig_data_nl,aes(Freq,F_nl))
pv_nl <- pv_nl+geom_line(aes(color=Rho_vals,group=Rho_vals),size=0.3)

pv_nl <- pv_nl + scale_color_brewer("",palette="Blues",labels = c(expression("" * rho * " = 0.05"), expression("" * rho * " = 0.1  "), expression("" * rho * " = 0.25")))
pv_nl <- pv_nl + xlab("Frequency (Hz)") + ylab(expression("Nonlinearity (Hz)"))
pv_nl <- pv_nl + opts(legend.key = theme_blank())
pv_nl <- pv_nl + scale_x_continuous(limits = c(-20,60),breaks=c(-20,0,20,40,60))
pv_nl <- pv_nl + scale_y_continuous(limits = c(-5,300))
pv_nl <- pv_nl + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
pv_nl <- pv_nl + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + opts(panel.grid.major = theme_blank()) + opts(axis.line = theme_segment())
	
pv_nl <- pv_nl + opts(plot.margin = unit(c(0.05,0.1,0,-0.1), "cm"))	
#pv_nl <- pv_nl + opts(legend.position = c(0.7, 0.25)) + opts(legend.text = theme_text(size = 9), legend.key.width = unit(0.5,"cm"), legend.key.height = unit(0.25, "cm"))
pv_nl <- pv_nl + opts(legend.position = 'none')
pv_nl <- pv_nl + opts(axis.text.x = theme_text(size = 7))
pv_nl <- pv_nl + opts(axis.text.y = theme_text(size = 7))
pv_nl <- pv_nl + opts(axis.title.y = theme_text(angle = 90, vjust = 0.6, size = 8))
pv_nl <- pv_nl + opts(axis.title.x = theme_text(size = 8, vjust=0.6))

#############
#############

rates <- fig_1_rates$firing.rate

Frate = 1000*rates[,1]
Time = rates[,2]
Type = rates[,3]

Type[Type == 1] <- "A"
Type[Type == 2] <- "B"

Zero_line = rep(0,length(Time))
Z_data = data.frame(Zero_line,Time)

fig_data_rates<-data.frame(Frate,Time,Type)
fig_data_rates_short<-subset(fig_data_rates,200<Time & Time<500)


pv_rate <- ggplot(fig_data_rates_short,aes(Time,Frate))
pv_rate <- pv_rate+geom_line(aes(color=Type,group=Type),size=0.3)
pv_rate <- pv_rate+geom_line(data=Z_data,aes(x=Time,y=Zero_line),linetype=2,size=0.2)

pv_rate <- pv_rate + scale_color_brewer("",palette="Paired",labels = c("Linear","LNL"))
pv_rate <- pv_rate + xlab("Time (ms)") + ylab(expression("Firing rate (Hz)"))
pv_rate <- pv_rate + opts(legend.key = theme_blank())
pv_rate <- pv_rate + scale_x_continuous(limits = c(200,500),breaks=c(200,250,300,350,400,450,493),labels=c(200,250,300,350,400,450,500))
pv_rate <- pv_rate + scale_y_continuous(limits = c(-10,50),breaks=c(0,10,20,30,40),labels=c(0,10,20,30,40))
pv_rate <- pv_rate + coord_cartesian(x=c(200,500),y=c(-10,50))
pv_rate <- pv_rate + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
pv_rate <- pv_rate + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + opts(panel.grid.major = theme_blank()) + opts(axis.line = theme_segment())
	
pv_rate <- pv_rate + opts(plot.margin = unit(c(-0.15,0.1,-0.1,-0.2), "cm"))	
pv_rate <- pv_rate + opts(legend.position = c(0.7, 0.85)) + opts(legend.text = theme_text(size = 9), legend.key.width = unit(0.5,"cm"), legend.key.height = unit(0.25, "cm"))
pv_rate <- pv_rate + opts(axis.text.x = theme_text(size = 7))
pv_rate <- pv_rate + opts(axis.text.y = theme_text(size = 7))
pv_rate <- pv_rate + opts(axis.title.y = theme_text(angle = 90, vjust = 0.6, size = 8))
pv_rate <- pv_rate + opts(axis.title.x = theme_text(size = 8, vjust=0.6))

#############
#############

subvp <- viewport(width = 0.3, height = 0.12, x = 0.83, y = 0.55)
subvp_1 <- viewport(width = 0.3, height = 0.12, x = 0.34, y = 0.55)
subvp_2 <- viewport(width = 0.3, height = 0.12, x = 0.83, y = 0.88)

pdf("/Users/dleen/Dropbox/Research/Simple_HOC_draft/R/fig_1/fig_1_right_half.pdf", width = 3.375, height = 2.755)
grid.newpage()

pushViewport(viewport(layout = grid.layout(2, 2)))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

print(pv_rate, vp = vplayout(2, 1:2))
print(pv, vp = vplayout(1, 1))
print(pv_nl, vp = vplayout(1, 2))
#print(p2,vp=subvp)
#print(p2_1,vp=subvp_1)
#print(p2_2,vp=subvp_2)

dev.off()
