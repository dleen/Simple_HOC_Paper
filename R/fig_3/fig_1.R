# This is one section from figure 1 plot.
# This file produces a plot of the voltage from a simulation
# of some neurons (doesn't matter what type) and plots it 
# against time in milliseconds. It also plots an alternating
# light transparent blue and clear, vertical checkboard
# pattern behind the voltage to indicate the bin size 
# relative to the time.

library(R.matlab) # Load matlab .mat files
library(ggplot2) # Great graphs
library(scales)
library(grid)

rm(list = ls())
theme_set(theme_bw(7)) # ggplot blackwhite theme

########
# Data #
########

# Load the voltage values
fig_1_volt <- readMat(paste("/Users/dleen/Dropbox/Research/Simple_HOC_draft/",
                            "codes/Figure_code/fig_1_volt.mat",sep=""))
# Load the bin times
fig_1_bins <- readMat(paste("/Users/dleen/Dropbox/Research/Simple_HOC_draft/",
                            "codes/Figure_code/bin_times.mat", sep=""))

# Just take the voltage readings from the matlab dataframe
VV <- fig_1_volt$volt.R.data
# Just take the bin times from the other matlab frame
Vbins <- fig_1_bins$bin.times

# 3rd column indicates which neuron the voltage belongs to
# neuron 1, 2, 3, 4
Neuron = VV[,3]
# The time in milliseconds
Time = VV[,2]
# The voltage reading in mV 
Voltage = VV[,1]

# Put these variables into a data frame
fig_data<-data.frame(Voltage,Time,Neuron)

# Don't display all 80,000 rows, take a small indicative sample
# The values 640, 750 show 1 spike per bin, and also 
# 2 spikes occuring in the same bin.
fig_data_short<-subset(fig_data,640<Time & Time<750)

# This keeps track of where the 10ms bins start and end
# This is done in a separate file to make the plotting easier.
Start = Vbins[,1]
End = Vbins[,2]
Bins = Vbins[,3]
Bins[Bins==1]<-"A" # Bin to color blue
Bins[Bins==2]<-"B" # Bin to leave transparent

# Dataframe for the bins
fig_bin_data<-data.frame(Start,End,Bins)
# Pick a subset like before
fig_bin_short<-subset(fig_bin_data,640<End & End<750)

# The min and max of the time 
xrng = range(fig_data_short$Time)

############
# Plotting #
############

# Start the plot using the shortened data, plotting voltage vs time
pv <- ggplot(fig_data_short, aes(Time,Voltage))

# Connect the points using a line and group according to neuron number
pv <- pv + geom_line(aes(color=Neuron, group = Neuron), size=0.3)

# Make the boxes corresponding to the time bins
pv <- pv + geom_rect(aes(NULL, NULL, xmin=Start, xmax=End, ymin=-75, ymax=-10, fill=Bins),
                     data=fig_bin_short, alpha=0.25)

# Tighen up the axes and make them fit to the data ranges
pv <- pv + coord_cartesian(x=c(xrng[1], xrng[2]), y=c(-75,-10))

# Color for the boxes, one blue, one empty
pv <- pv + scale_fill_manual(values = c(NA, "#55B1F7"))

# Turn off legend, don't think it's needed here
pv <- pv + theme(legend.position = "none")

# Making the y-axis text and title fit just right
pv <- pv + theme(axis.text.y=element_text(angle=90, vjust=0.3,size=5),
                 axis.title.y=element_text(angle=90,vjust=0.55,size=10))

# Put breaks along the y-axis at values which are important, like the reset voltage
pv <- pv + scale_y_discrete(breaks=c(-70,-60,-52,-20,-13),labels=c(-70,-60,-53,0,20))

# Remove breaks along the x-axis as these are redundant with the boxes
# marking off the bins anyway
pv <- pv + scale_x_discrete(breaks=NULL)

# Adjust the tick size and length as necessary
pv <- pv + theme(axis.ticks = element_line(size = 0.1), 
                 axis.ticks.length = unit(0.05, "cm"))

# Add x label and adjust its text size
pv <- pv+xlab("Time (ms)")+theme(axis.text.x=element_text(size=5), 
                                 axis.title.x=element_text(size=10,vjust=1.4))

# Adjust the margins for snugness
pv <- pv+theme(plot.margin = unit(c(0,0,-0.25,-0.15), "cm"))

############
# Printing #
############

# Open pdf to which to write
pdf(paste("/Users/dleen/Dropbox/Research/Simple_HOC_draft/",
          "R/fig_1/fig_1_newest.pdf", sep=""),width=3.375,height=1.2)
print(pv) # Print output to pdf 
dev.off() # close file after use
