library(R.matlab)
library(ggplot2)
library(scales)
library(grid)

theme_set(theme_grey(7))
fig_3a <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/fig_3a_R.mat")

PP <- fig_3a$fig.3a.R

Model <- PP[, 4]
Num <- PP[, 2]
Num[Num == 8] <- "8  Neurons"
Num[Num == 32] <- "32  Neurons"
Num[Num == 64] <- "64  Neurons"
Num[Num == 100] <- "100  Neurons"

#Model[Model==1] <- "EIF"
#Model[Model==2] <- "Ising"

N_list <- PP[, 3]
Prob <- PP[, 1]
Prob[Prob < 1e-08] <- 1e-09

fig_data_orig <- data.frame(Prob, Num, N_list, Model)
#fig_data <- subset(fig_data_orig,Num!="100  Neurons")
fig_data <- fig_data_orig

fig_data$Num <- ordered(fig_data$Num, levels = rev(c("100  Neurons", "64  Neurons", "32  Neurons", "8  Neurons")))

p1 <- ggplot(fig_data, aes(N_list, Prob, fill = Model))
p1 <- p1 + geom_bar(stat = "identity", position = "identity", alpha = 0.75)
p1 <- p1 + scale_fill_continuous("", breaks = c(1, 2), labels = c("EIF", "LNL Th."))
#p1 <- p1 + scale_fill_brewer(palette="Dark2")
p1 <- p1 + facet_wrap(~Num, scales = "free", ncol = 2)
p1 <- p1 + xlab("Population spike count") + ylab("Probability")
p1 <- p1 + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank())
p1 <- p1 + opts(legend.position = c(0.28, 0.9)) + opts(legend.text = theme_text(size = 9), legend.key.width = unit(0.25,"cm"), legend.key.height = unit(0.25, "cm"))
p1 <- p1 + opts(strip.text.x = theme_text(size = 7),strip.background = theme_rect(colour="white",fill="#E8E8E8"))
p1 <- p1 + opts(axis.ticks = theme_segment(size = 0.05), axis.ticks.length = unit(0.05, "cm"))
p1 <- p1 + opts(axis.text.x = theme_text(size = 7))
p1 <- p1 + opts(axis.text.y = theme_text(size = 7))
p1 <- p1 + opts(axis.title.y = theme_text(angle = 90, size = 8,vjust=0.55))
p1 <- p1 + opts(axis.title.x = theme_text(size = 8,vjust=0.35))
p1 <- p1 + opts(plot.margin = unit(c(0,0.05,0,-0.15), "cm"))+opts(panel.margin = unit(0.1,"cm"))


################
################

fig_data2_100 <- subset(fig_data_orig, Num == "100  Neurons")

p2 <- ggplot(fig_data2_100, aes(N_list, Prob, group = Model, color = Model))
p2 <- p2 + geom_line(size = 0.5)
p2 <- p2 + scale_color_gradient("Model", breaks = c(1, 2))
p2 <- p2 + scale_y_log10(limits = c(1e-04, 0.2), breaks = c(0.1, 0.01,0.001, 1e-04), labels = c(expression(10^-1),expression(10^-2), expression(10^-3), expression(10^-4)))
p2 <- p2 + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), labels = c("0", 
	"25", "50", "75", "100"))
p2 <- p2 + xlab(NULL) + ylab(NULL)
p2 <- p2 + opts(legend.position = "none")
p2 <- p2 + scale_x_continuous(limits = c(0, 100), breaks = c(0.1, 25, 50, 75, 99.9), labels = c(0, 25, 50, 75, 100))
p2 <- p2 + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
p2 <- p2 + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + 
	opts(axis.line = theme_segment())
p2 <- p2 + opts(plot.margin = unit(c(0,0,0,0), "cm"))
	

fig_data2_64 <- subset(fig_data_orig, Num == "64  Neurons")

p2_1 <- ggplot(fig_data2_64, aes(N_list, Prob, group = Model, color = Model))
p2_1 <- p2_1 + geom_line(size = 0.5)
p2_1 <- p2_1 + scale_color_gradient("Model", breaks = c(1, 2))
p2_1 <- p2_1 + scale_y_log10(limits = c(1e-04, 0.2), breaks = c(0.1, 0.01,0.001, 1e-04), labels = c(expression(10^-1),expression(10^-2), expression(10^-3), expression(10^-4)))
p2_1 <- p2_1 + xlab(NULL) + ylab(NULL)
p2_1 <- p2_1 + opts(legend.position = "none")
p2_1 <- p2_1 + scale_x_continuous(limits = c(0, 64), breaks = c(0.1, 16, 32, 48, 64), labels = c(0, 16, 32, 48, 64))
p2_1 <- p2_1 + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
p2_1 <- p2_1 + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + 
	opts(axis.line = theme_segment())
p2_1 <- p2_1 + opts(plot.margin = unit(c(0,0,0,0), "cm"))


fig_data2_32 <- subset(fig_data_orig, Num == "32  Neurons")

p2_2 <- ggplot(fig_data2_32, aes(N_list, Prob, group = Model, color = Model))
p2_2 <- p2_2 + geom_line(size = 0.5)
p2_2 <- p2_2 + scale_color_gradient("Model", breaks = c(1, 2))
p2_2 <- p2_2 + scale_y_log10(limits = c(1e-04, 0.2), breaks = c(0.1, 0.01,0.001, 1e-04), labels = c(expression(10^-1),expression(10^-2), expression(10^-3), expression(10^-4)))
p2_2 <- p2_2 + xlab(NULL) + ylab(NULL)
p2_2 <- p2_2 + opts(legend.position = "none")
p2_2 <- p2_2 + scale_x_continuous(limits = c(0, 32), breaks = c(0.1, 8, 16, 24, 32), labels = c(0, 8, 16, 24, 32))
p2_2 <- p2_2 + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
p2_2 <- p2_2 + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + 
	opts(axis.line = theme_segment())
p2_2 <- p2_2 + opts(plot.margin = unit(c(0,0,0,0), "cm"))




###############
###############

fig_3b <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/fig_3b_R_logN_DJS.mat")

PP_2b <- fig_3b$fig.3b.R

Prob_e <- PP_2b[, 4]
Rho <- PP_2b[, 2]
N_list_2b <- PP_2b[, 1]
Prob_m <- PP_2b[, 3]

Rho[Rho==0.05] <- "Small"
Rho[Rho==0.1] <- "Med"
Rho[Rho==0.25] <- "Large"

fig_data_2b <- data.frame(Prob_m, Prob_e, N_list_2b, Rho)
fig_data_2b$Rho <- ordered(fig_data_2b$Rho, levels = rev(c("Large", "Med", "Small")))

p2b <- ggplot(fig_data_2b, aes(N_list_2b, Prob_m, group = Rho, color = Rho))
p2b <- p2b + geom_line(size = 0.75)
p2b <- p2b + scale_color_brewer("",palette="Blues",labels = c(expression("" * rho * " = 0.05"), expression("" * rho * " = 0.1  "), expression("" * rho * " = 0.25")))
p2b <- p2b + xlab("Population size N") + ylab(expression("JS-div / log N   "* (10^{-3})*""))
p2b <- p2b + opts(legend.key = theme_blank())
p2b <- p2b + scale_x_continuous(limits = c(4, 100), breaks = c(4, 25, 50, 75, 99.9), labels = c(4, 25, 50, 75, 100))
#p2b <- p2b + scale_y_continuous(labels=scientific_format())
p2b <- p2b + scale_y_continuous(limits=c(0,10e-03),breaks=c(0,2e-03,4e-03,6e-03,8e-03,10e-3),labels=c(0,2,4,6,8,10))
p2b <- p2b + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
p2b <- p2b + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + 
	opts(panel.grid.major = theme_blank()) + opts(axis.line = theme_segment())
	
p2b <- p2b + opts(plot.margin = unit(c(0.25,0.05,0,-0.06), "cm"))	
p2b <- p2b + opts(legend.position = c(0.31, 0.88)) + opts(legend.text = theme_text(size = 9), legend.key.width = unit(0.35,"cm"), legend.key.height = unit(0.25, "cm"))
p2b <- p2b + opts(axis.text.x = theme_text(size = 7))
p2b <- p2b + opts(axis.text.y = theme_text(size = 7))
p2b <- p2b + opts(axis.title.y = theme_text(angle = 90, vjust = 0.45, size = 8))
p2b <- p2b + opts(axis.title.x = theme_text(size = 8, vjust=0.3))


###############
###############

fig_3c <- readMat("/Users/dleen/Dropbox/Research/Simple_HOC_draft/codes/Figure_code/fig_3c_R_logN_DJS.mat")

PP_2c <- fig_3c$fig.3c.R

Probc_e <- PP_2c[, 4]
Mu <- PP_2c[, 2]
N_list_2c <- PP_2c[, 1]
Probc_m <- PP_2c[, 3]

Mu[Mu==0.1] <- "Small"
Mu[Mu==0.2] <- "Med"
Mu[Mu==0.3] <- "Large"

fig_data_2c <- data.frame(Probc_m, Probc_e, N_list_2c, Mu)
fig_data_2c$Mu <- ordered(fig_data_2c$Mu, levels = rev(c("Large", "Med", "Small")))

p2c <- ggplot(fig_data_2c, aes(N_list_2c, Probc_m, group = Mu, color = Mu))
p2c <- p2c + geom_line(size = 0.75)
p2c <- p2c + scale_color_brewer("", palette="Blues", labels = c(expression("" * mu * " = 0.1"), expression("" * mu * " = 0.2"), expression("" * mu * " = 0.3")))
p2c <- p2c + xlab("Population size N") + ylab(NULL)
p2c <- p2c + opts(legend.key = theme_blank())
p2c <- p2c + scale_x_continuous(limits = c(4, 100), breaks = c(4, 25, 50, 75, 100), labels = c(4, 25, 50, 75, 100))
p2c <- p2c + scale_y_continuous(limits=c(0,5e-03),breaks=c(0,1e-03,2e-03,3e-03,4e-03,5e-03),labels=c(0,1,2,3,4,5))
p2c <- p2c + opts(axis.ticks = theme_segment(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
p2c <- p2c + opts(panel.background = theme_rect(colour = NA)) + opts(panel.grid.minor = theme_blank()) + opts(panel.grid.major = theme_blank()) + 
	opts(axis.line = theme_segment())
	
p2c <- p2c + opts(plot.margin = unit(c(0.25,0.1,0,0.2), "cm"))
p2c <- p2c + opts(legend.position = c(0.33, 0.86)) + opts(legend.text = theme_text(size = 9), legend.key.width = unit(0.5,"cm"), legend.key.height = unit(0.25, "cm"))
p2c <- p2c + opts(axis.text.x = theme_text(size = 7))
p2c <- p2c + opts(axis.text.y = theme_text(size = 7))
p2c <- p2c + opts(axis.title.y = theme_text(angle = 90, vjust = 0.2, size = 8))
p2c <- p2c + opts(axis.title.x = theme_text(size = 8, vjust=0.3))


#############
#############

subvp <- viewport(width = 0.3, height = 0.12, x = 0.83, y = 0.55)
subvp_1 <- viewport(width = 0.3, height = 0.12, x = 0.34, y = 0.55)
subvp_2 <- viewport(width = 0.3, height = 0.12, x = 0.83, y = 0.88)

pdf("/Users/dleen/Dropbox/Research/Simple_HOC_draft/R/fig_3/fig_3a_test.pdf", width = 3.375, height = 5)
grid.newpage()
grid.text(paste("margin"))

pushViewport(viewport(layout = grid.layout(3, 2)))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)



print(p1, vp = vplayout(1:2, 1:2))
#print(p2lin,vp=subvp)
print(p2b, vp = vplayout(3, 1))
print(p2c, vp = vplayout(3, 2))
print(p2,vp=subvp)
print(p2_1,vp=subvp_1)
print(p2_2,vp=subvp_2)

dev.off()
