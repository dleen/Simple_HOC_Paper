library(R.matlab) # Reads from matlab .mat files
library(ggplot2) # Great plots
library(scales)
library(grid) # Multiple plots on one figure

rm(list = ls())
theme_set(theme_grey(7))

fig_4 <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/L_func_comp.mat")

L_comp<-fig_4$L.function

L_s = L_comp[,1]
s_var = L_comp[,2]
Type = L_comp[,3]

L_s_short = L_comp[6603:7203,1]
s_var_short = L_comp[6603:7203,2]

Type[Type == 1] <- "A"
Type[Type == 2] <- "B"
Type[Type == 3] <- "C"

fig_data<-data.frame(L_s,s_var,Type)
fig_data_1<-subset(fig_data,Type != "C")
fig_data_2<-subset(fig_data,Type == "C")
#fig_data_2<-data.frame(L_s_short,s_var_short,Type=="C")

pv_b <- ggplot(fig_data_1,aes(s_var,L_s))
pv_b <- pv_b+geom_line(aes(color=Type,group=Type),size=0.3)

pv_b <- pv_b + scale_color_brewer("",palette="Paired",labels = c("DG","LNL"))
pv_b <- pv_b + xlab("s") + ylab(expression("L(s)"))
pv_b <- pv_b + opts(legend.key = theme_blank())
pv_b <- pv_b + scale_x_continuous(limits = c(-2,2.75))
pv_b <- pv_b + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
pv_b <- pv_b + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + opts(panel.grid.major = theme_blank()) + opts(axis.line = theme_segment())
	
pv_b <- pv_b + opts(plot.margin = unit(c(0.05,0,0,-0.2), "cm"))	
pv_b <- pv_b + opts(legend.position = c(0.7, 0.85)) + opts(legend.text = theme_text(size = 9), legend.key.width = unit(0.5,"cm"), legend.key.height = unit(0.25, "cm"))
pv_b <- pv_b + opts(axis.text.x = theme_text(size = 7))
pv_b <- pv_b + opts(axis.text.y = theme_text(size = 7))
pv_b <- pv_b + opts(axis.title.y = theme_text(angle = 90, vjust = 0.6, size = 8))
pv_b <- pv_b + opts(axis.title.x = theme_text(size = 8, vjust=0.6))

pv_b <- pv_b + geom_area(data=fig_data_2,aes(s_var_short,L_s_short,fill=Type,alpha=0.1))
pv_b <- pv_b + scale_fill_manual(values="gray60")

################
################

fig_4a <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/fig_4a_R.mat")

EIF_DG<-fig_4a$EIF.DG.comp

Prob = EIF_DG[,1]
Nums = EIF_DG[,2]
Rho = EIF_DG[,3]
Type_a = EIF_DG[,4]

Rho[Rho==0.05] <- "Small"
Rho[Rho==0.1] <- "Med"
Rho[Rho==0.25] <- "Large"

Type_a[Type_a == 1] <- "EIF"
Type_a[Type_a == 2] <- "DG"

fig_data_a<-data.frame(Prob,Nums,Rho,Type_a)
fig_data_a$Rho <- ordered(fig_data_a$Rho, levels = rev(c("Large", "Med", "Small")))

test <-subset(fig_data_a,Type_a=="EIF")
test1 <-subset(fig_data_a,Type_a=="DG")

pv_a <- ggplot(test1,aes(Nums,Prob))
pv_a <- pv_a+geom_line(data=test1,aes(color=Rho,group=Rho),size=0.3)

pv_a <- pv_a+geom_line(data=test,aes(color=Rho,group=Rho),size=0.3)
#pv_a <- pv_a+scale_linetype_manual(values="dashed")

pv_a <- pv_a + scale_y_log10(limits=c(1e-4,1),breaks = c(0.1, 0.01,0.001, 1e-04), labels = c(expression(10^-1),expression(10^-2), expression(10^-3), expression(10^-4)))

pv_a <- pv_a + scale_color_brewer("",palette="Blues",labels = c(expression("" * rho * " = 0.05"), expression("" * rho * " = 0.1  "), expression("" * rho * " = 0.25")))
pv_a <- pv_a + xlab("Population spike count") + ylab("P(k)")
pv_a <- pv_a + opts(legend.key = theme_blank())
pv_a <- pv_a + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 99.9), labels = c(0, 25, 50, 75, 100))
pv_a <- pv_a + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
pv_a <- pv_a + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + 
	opts(panel.grid.major = theme_blank()) + opts(axis.line = theme_segment())
	
pv_a <- pv_a + opts(plot.margin = unit(c(0,0,0,0), "cm"))	
pv_a <- pv_a + opts(legend.position = c(0.33, 0.86)) + opts(legend.text = theme_text(size = 9), legend.key.width = unit(0.5,"cm"), legend.key.height = unit(0.25, "cm"))
pv_a <- pv_a + opts(axis.text.x = theme_text(size = 7))
pv_a <- pv_a + opts(axis.text.y = theme_text(size = 7))
pv_a <- pv_a + opts(axis.title.y = theme_text(angle = 90, vjust = 0.45, size = 8))
pv_a <- pv_a + opts(axis.title.x = theme_text(size = 8, vjust=0.3))

################
################

fig_4c <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/heat_cap_R.mat")

heat<-fig_4c$heat.cap

heat_cap = heat[,1]
Nums = heat[,2]
Type_c = heat[,3]

Type_c[Type_c == 1] <- "DG"
Type_c[Type_c == 2] <- "EIF"
Type_c[Type_c == 3] <- "LNL"
Type_c[Type_c == 4] <- "Ising"

fig_data_c<-data.frame(heat_cap,Nums,Type_c)
fig_data_c$Type_c <- ordered(fig_data_c$Type_c, levels = rev(c("DG", "EIF", "LNL", "Ising")))

pv_c <- ggplot(fig_data_c,aes(Nums,heat_cap))
pv_c <- pv_c+geom_line(aes(color=Type_c,group=Type_c),size=0.3)

pv_c <- pv_c + scale_color_brewer("",palette="Paired")
pv_c <- pv_c + xlab("Population size N") + ylab("Heat capacity")
pv_c <- pv_c + opts(legend.key = theme_blank())
pv_c <- pv_c + scale_x_continuous(limits = c(4, 100), breaks = c(4, 25, 50, 75, 99.9), labels = c(4, 25, 50, 75, 100))
pv_c <- pv_c + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
pv_c <- pv_c + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + 
	opts(panel.grid.major = theme_blank()) + opts(axis.line = theme_segment())
	
pv_c <- pv_c + opts(plot.margin = unit(c(0,0,0,0), "cm"))	
pv_c <- pv_c + opts(legend.position = c(0.33, 0.86)) + opts(legend.text = theme_text(size = 9), legend.key.width = unit(0.5,"cm"), legend.key.height = unit(0.25, "cm"))
pv_c <- pv_c + opts(axis.text.x = theme_text(size = 7))
pv_c <- pv_c + opts(axis.text.y = theme_text(size = 7))
pv_c <- pv_c + opts(axis.title.y = theme_text(angle = 90, vjust = 0.45, size = 8))
pv_c <- pv_c + opts(axis.title.x = theme_text(size = 8, vjust=0.3))


pdf("/Users/dleen/Dropbox/Research/Simple_HOC_draft/R/fig_4/fig_4_newest.pdf", width = 6.75, height = 1.666)
grid.newpage()

pushViewport(viewport(layout = grid.layout(1,3)))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

print(pv_a,vp = vplayout(1,1))
print(pv_b,vp = vplayout(1,2))
print(pv_c,vp = vplayout(1,3))

dev.off()