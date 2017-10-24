# N-CELLS
coordinates <- read.delim("enhancerAnnotation/peakAnalysis/coordinates.bed",header = FALSE)
N_ac <- read.delim(file = "enhancerAnnotation/peakAnalysis/N_ac_per_sample.tab",header = FALSE, sep = " ")
N_me <- read.delim(file = "enhancerAnnotation/peakAnalysis/N_me_per_sample.tab",header = FALSE, sep = " ")
peaksN_me_ac <- coordinates[apply(N_ac,1,all) & apply(N_me,1,all),]
peaksN_me <- coordinates[!apply(N_ac,1,all) & apply(N_me,1,all),]
peaksN_ac <- coordinates[apply(N_ac,1,all) & !apply(N_me,1,all),]
write.table(peaksN_me_ac, file = "enhancerAnnotation/specificPeaks/N_me_ac",col.names = FALSE,quote = FALSE,row.names = FALSE)
write.table(peaksN_me, file = "enhancerAnnotation/specificPeaks/N_me",col.names = FALSE,quote = FALSE,row.names = FALSE)
write.table(peaksN_ac, file = "enhancerAnnotation/specificPeaks/N_ac",col.names = FALSE,quote = FALSE,row.names = FALSE)

# T-CELLS
coordinates <- read.delim("enhancerAnnotation/peakAnalysis/coordinates.bed",header = FALSE)
T_ac <- read.delim(file = "enhancerAnnotation/peakAnalysis/T_ac_per_sample.tab",header = FALSE, sep = " ")
T_me <- read.delim(file = "enhancerAnnotation/peakAnalysis/T_me_per_sample.tab",header = FALSE, sep = " ")
peaksT_me_ac <- coordinates[apply(T_ac,1,all) & apply(T_me,1,all),]
peaksT_me <- coordinates[!apply(T_ac,1,all) & apply(T_me,1,all),]
peaksT_ac <- coordinates[apply(T_ac,1,all) & !apply(T_me,1,all),]
write.table(peaksT_me_ac, file = "enhancerAnnotation/specificPeaks/T_me_ac",col.names = FALSE,quote = FALSE,row.names = FALSE)
write.table(peaksT_me, file = "enhancerAnnotation/specificPeaks/T_me",col.names = FALSE,quote = FALSE,row.names = FALSE)
write.table(peaksT_ac, file = "enhancerAnnotation/specificPeaks/T_ac",col.names = FALSE,quote = FALSE,row.names = FALSE)