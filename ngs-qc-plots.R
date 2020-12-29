args <- commandArgs(TRUE)
srcPath <- args[1]
samplename=args[2]



#----------------------------------------------------------------------------------------
### Read Length Distribution ###
#Import Library
library("ggplot2")

#Read data from file
dat = read.table(file.path(srcPath, "length-dist.txt"), header=F)

# make V1 an ordered factor
dat$V1 <- factor(dat$V1, levels = dat$V1)

#Plot
qplot(V1, V2, data=dat, group=1, geom="line", colour="blue") + scale_color_manual(values=c("blue")) + labs(x = "Read Length (bp)", y = "Read Count", title = paste0(samplename,': ',"Read Length Distribution")) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size=8)) 

#Save to PNG
ggsave(file=file.path(srcPath, paste0(samplename,'-',"length-dist.png")), width=10, height=4)


#----------------------------------------------------------------------------------------
### Read Quality Distribution ###
#Import Library
#library("ggplot2")

#Read data from file
dat = read.table(file.path(srcPath, "read-quality.txt"), header=F)

# make V1 an ordered factor
dat$V1 <- factor(dat$V1, levels = dat$V1)

#Plot
qplot(V1, V2, data=dat, group=1, geom="line", colour="blue") + scale_color_manual(values=c("blue")) + labs(x = "Phred Quality Score", y = "Number of Reads", title = paste0(samplename,': ',"Read Quality Distribution")) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size=8)) + geom_vline(xintercept = 20, linetype = "longdash") + geom_vline(xintercept = 28, linetype = "longdash")

#Save to PNG
ggsave(file=file.path(srcPath, paste0(samplename,'-',"read-quality.png")), width=10, height=4)

#----------------------------------------------------------------------------------------
### Base Quality Distribution ###
#Import Library

#Read data from file
dat = read.table(file.path(srcPath, "base-quality.txt"), header=F)

# make V1 an ordered factor
dat$V1 <- factor(dat$V1, levels = dat$V1)

#Generate QC group 
group = dat$V3>28
cols=c("TRUE"="#66CC00", "FALSE"="#FFB266")

#Plot
ggplot(dat, aes(x = factor(V1))) + geom_boxplot(aes(lower = V4, upper = V5, middle = V3, ymin = V6, ymax = V7, colour=group), stat = "identity") + geom_hline(yintercept = 28, linetype = "longdash", colour = "#66CC00", size=0.3) + geom_hline(yintercept = 20, linetype = "longdash", colour = "red", size=0.3)  + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size=8)) + labs(y = "Phred Quality Score", x = "Base Position", title = paste0(samplename,': ',"Base Quality Distribution"))+scale_colour_manual(values=cols)

#Save to PNG
ggsave(file=file.path(srcPath, paste0(samplename,'-',"base-quality.png")), width=10, height=4)

#----------------------------------------------------------------------------------------
