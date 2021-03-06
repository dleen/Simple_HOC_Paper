library(R.matlab) # Read matlab matrix files
library(ggplot2) # For the plotting
library(scales) 
library(grid) # For plotting multiple plots in a grid

theme_set(theme_bw(7))
fig_4a <- readMat(paste("/Users/dleen/Dropbox/Research/",
                        "Simple_HOC_draft/codes/Figure_code/",
                        "fig_3a_R.mat", sep=""))

PP <- fig_4a$fig.3a.R

Model <- PP[, 4]
Num <- PP[, 2]
Num[Num == 8] <- "8  Neurons"
Num[Num == 32] <- "32  Neurons"
Num[Num == 64] <- "64  Neurons"
Num[Num == 100] <- "100  Neurons"

N_list <- PP[, 3]
Prob <- PP[, 1]
Prob[Prob < 1e-08] <- 1e-09

fig_data_orig <- data.frame(Prob, Num, N_list, Model)
fig_data <- fig_data_orig

fig_data$Num <- ordered(fig_data$Num, 
                        levels = rev(c("100  Neurons", "64  Neurons", "32  Neurons", "8  Neurons")))

p1 <- ggplot(fig_data, aes(N_list, Prob, fill = factor(Model)))
p1 <- p1 + geom_bar(stat = "identity", position = "identity", alpha = 0.75)

p1 <- p1 + scale_fill_manual("Model", values = c("#132B43","#55B1F7"), labels = c("EIF", "LNL Th."))

p1 <- p1 + facet_wrap(~Num, scales = "free", ncol = 2)
p1 <- p1 + xlab("Population spike count") + ylab("Probability")

p1 <- p1 + theme(panel.background = element_rect(colour = NA)) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank())
p1 <- p1 + theme(legend.position = c(0.28, 0.9)) + 
  theme(legend.text = element_text(size = 9), 
        legend.key.width = unit(0.25,"cm"), 
        legend.key.height = unit(0.25, "cm"))
p1 <- p1 + theme(strip.text.x = element_text(size = 7),
                 strip.background = element_rect(colour="white",fill="#E8E8E8"))
p1 <- p1 + theme(axis.ticks = element_line(size = 0.05), 
                 axis.ticks.length = unit(0.05, "cm"))
p1 <- p1 + theme(axis.text.x = element_text(size = 7))
p1 <- p1 + theme(axis.text.y = element_text(size = 7))
p1 <- p1 + theme(axis.title.y = element_text(angle = 90, size = 8,vjust=0.55))
p1 <- p1 + theme(axis.title.x = element_text(size = 8,vjust=0.35))
p1 <- p1 + theme(plot.margin = unit(c(0,0.05,0,-0.15), "cm")) + 
  theme(panel.margin = unit(0.1,"cm"))


################
################

fig_data2_100 <- subset(fig_data_orig, Num == "100  Neurons")

p2 <- ggplot(fig_data2_100, aes(N_list, Prob, group = Model, color = Model))
p2 <- p2 + geom_line(size = 0.5)
p2 <- p2 + scale_color_gradient("Model", breaks = c(1, 2))
p2 <- p2 + scale_y_log10(limits = c(1e-04, 0.2), 
                         breaks = c(0.1, 0.01,0.001, 1e-04), 
                         labels = c(expression(10^-1),expression(10^-2), expression(10^-3), expression(10^-4)))
p2 <- p2 + xlab(NULL) + ylab(NULL)
p2 <- p2 + theme(legend.position = "none")
p2 <- p2 + scale_x_continuous(limits = c(0, 100), 
                              breaks = c(0.1, 25, 50, 75, 99.9), 
                              labels = c(0, 25, 50, 75, 100))
p2 <- p2 + theme(axis.ticks = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
p2 <- p2 + theme(panel.background = element_rect(colour = NA)) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) + 
	theme(axis.line = element_line())
p2 <- p2 + theme(plot.margin = unit(c(0,0,0,0), "cm"))
	

fig_data2_64 <- subset(fig_data_orig, Num == "64  Neurons")

p2_1 <- ggplot(fig_data2_64, aes(N_list, Prob, group = Model, color = Model))
p2_1 <- p2_1 + geom_line(size = 0.5)
p2_1 <- p2_1 + scale_color_gradient("Model", breaks = c(1, 2))
p2_1 <- p2_1 + scale_y_log10(limits = c(1e-04, 0.2), 
                             breaks = c(0.1, 0.01,0.001, 1e-04), 
                             labels = c(expression(10^-1),
                                        expression(10^-2), 
                                        expression(10^-3), 
                                        expression(10^-4)))
p2_1 <- p2_1 + xlab(NULL) + ylab(NULL)
p2_1 <- p2_1 + theme(legend.position = "none")
p2_1 <- p2_1 + scale_x_continuous(limits = c(0, 64), 
                                  breaks = c(0.1, 16, 32, 48, 64), 
                                  labels = c(0, 16, 32, 48, 64))
p2_1 <- p2_1 + theme(axis.ticks = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
p2_1 <- p2_1 + theme(panel.background = element_rect(colour = NA)) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) + 
	theme(axis.line = element_line())
p2_1 <- p2_1 + theme(plot.margin = unit(c(0,0,0,0), "cm"))

#
# Subfigure
#
fig_data2_32 <- subset(fig_data_orig, Num == "32  Neurons")

p2_2 <- ggplot(fig_data2_32, aes(N_list, Prob, group = Model, color = Model))
p2_2 <- p2_2 + geom_line(size = 0.5)
p2_2 <- p2_2 + scale_color_gradient("Model", breaks = c(1, 2))
p2_2 <- p2_2 + scale_y_log10(limits = c(1e-04, 0.2), 
                             breaks = c(0.1, 0.01,0.001, 1e-04), 
                             labels = c(expression(10^-1),
                                        expression(10^-2), 
                                        expression(10^-3), 
                                        expression(10^-4)))
p2_2 <- p2_2 + xlab(NULL) + ylab(NULL)
p2_2 <- p2_2 + theme(legend.position = "none")
p2_2 <- p2_2 + scale_x_continuous(limits = c(0, 32), 
                                  breaks = c(0.1, 8, 16, 24, 32), 
                                  labels = c(0, 8, 16, 24, 32))
p2_2 <- p2_2 + theme(axis.ticks = element_line(size = 0.1), 
                     axis.ticks.length = unit(0.05, "cm"))
p2_2 <- p2_2 + theme(panel.background = element_rect(colour = NA)) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) + 
	theme(axis.line = element_line())
p2_2 <- p2_2 + theme(plot.margin = unit(c(0,0,0,0), "cm"))


###############
###############

fig_4b <- readMat(paste("/Users/dleen/Dropbox/Research/",
                        "Simple_HOC_draft/codes/Figure_code/",
                        "fig_3b_R_logN_DJS.mat", sep=""))

PP_2b <- fig_4b$fig.3b.R

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
p2b <- p2b + scale_color_brewer("",palette="Blues",
                                labels = c(expression("" * rho * " = 0.05"), 
                                           expression("" * rho * " = 0.1  "), 
                                           expression("" * rho * " = 0.25")))
p2b <- p2b + xlab("Population size N") + 
  ylab(expression("JS-div / log N   "* (10^{-3})*""))
p2b <- p2b + theme(legend.key = element_blank())
p2b <- p2b + scale_x_continuous(limits = c(4, 100), 
                                breaks = c(4, 25, 50, 75, 99.9), 
                                labels = c(4, 25, 50, 75, 100))
p2b <- p2b + scale_y_continuous(limits=c(0,10e-03),
                                breaks=c(0,2e-03,4e-03,6e-03,8e-03,10e-3),
                                labels=c(0,2,4,6,8,10))
p2b <- p2b + theme(axis.ticks = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"))
p2b <- p2b + theme(panel.background = element_rect(colour = NA)) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) + 
	theme(panel.grid.major = element_blank()) +
  theme(axis.line = element_line())
	
p2b <- p2b + theme(plot.margin = unit(c(0.25,0.05,0,-0.06), "cm"))	
p2b <- p2b + theme(legend.position = c(0.31, 0.88)) +
  theme(legend.text = element_text(size = 9), 
        legend.key.width = unit(0.35,"cm"), 
        legend.key.height = unit(0.25, "cm"))
p2b <- p2b + theme(axis.text.x = element_text(size = 7))
p2b <- p2b + theme(axis.text.y = element_text(size = 7))
p2b <- p2b + theme(axis.title.y = element_text(angle = 90, vjust = 0.45, size = 8))
p2b <- p2b + theme(axis.title.x = element_text(size = 8, vjust=0.3))


###############
###############

fig_3c <- readMat(paste("/Users/dleen/Dropbox/Research/",
                        "Simple_HOC_draft/codes/Figure_code/",
                        "fig_3c_R_logN_DJS.mat", sep=""))

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
p2c <- p2c + scale_color_brewer("", palette="Blues",
                                labels = c(expression("" * mu * " = 0.1"),
                                           expression("" * mu * " = 0.2"),
                                           expression("" * mu * " = 0.3")))
p2c <- p2c + xlab("Population size N") + ylab(NULL)
p2c <- p2c + theme(legend.key = element_blank())
p2c <- p2c + scale_x_continuous(limits = c(4, 100), 
                                breaks = c(4, 25, 50, 75, 100), 
                                labels = c(4, 25, 50, 75, 100))
p2c <- p2c + scale_y_continuous(limits=c(0,5e-03),
                                breaks=c(0,1e-03,2e-03,3e-03,4e-03,5e-03),
                                labels=c(0,1,2,3,4,5))
p2c <- p2c + theme(axis.ticks = element_line(size = 0.1), 
                   axis.ticks.length = unit(0.05, "cm"))
p2c <- p2c + theme(panel.background = element_rect(colour = NA)) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) + 
	theme(axis.line = element_line())
	
p2c <- p2c + theme(plot.margin = unit(c(0.25,0.1,0,0.2), "cm"))
p2c <- p2c + theme(legend.position = c(0.33, 0.86)) + 
  theme(legend.text = element_text(size = 9), 
        legend.key.width = unit(0.5,"cm"), 
        legend.key.height = unit(0.25, "cm"))
p2c <- p2c + theme(axis.text.x = element_text(size = 7))
p2c <- p2c + theme(axis.text.y = element_text(size = 7))
p2c <- p2c + theme(axis.title.y = element_text(angle = 90, vjust = 0.2, size = 8))
p2c <- p2c + theme(axis.title.x = element_text(size = 8, vjust=0.3))


#############
#############

subvp <- viewport(width = 0.3, height = 0.12, x = 0.83, y = 0.55)
subvp_1 <- viewport(width = 0.3, height = 0.12, x = 0.34, y = 0.55)
subvp_2 <- viewport(width = 0.3, height = 0.12, x = 0.83, y = 0.88)

pdf(paste("/Users/dleen/Dropbox/Research/",
          "Simple_HOC_draft/R/fig_4/fig_4a_new.pdf", sep=""), width = 3.375, height = 5)

grid.newpage()
grid.text(paste("margin"))

pushViewport(viewport(layout = grid.layout(3, 2)))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)



print(p1, vp = vplayout(1:2, 1:2))
print(p2b, vp = vplayout(3, 1))
print(p2c, vp = vplayout(3, 2))
print(p2,vp=subvp)
print(p2_1,vp=subvp_1)
print(p2_2,vp=subvp_2)

dev.off()