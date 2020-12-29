library(ballgown)

args <-commandArgs(trailingOnly=TRUE)
dir <-args[1]
pattern<-args[2]
measure <- args[3]
file <- args[4]
bg = ballgown(dataDir=dir, samplePattern=pattern, meas=measure)
save(bg, file=file)