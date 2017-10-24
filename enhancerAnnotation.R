# N-CELLS
coordinates <- read.delim("enhancerAnnotation/peakAnalysis/coordinates.bed",header = FALSE)
N_ac <- read.delim(file = "enhancerAnnotation/peakAnalysis/N_ac_per_sample.tab",header = FALSE, sep = " ")
N_me <- read.delim(file = "enhancerAnnotation/peakAnalysis/N_me_per_sample.tab",header = FALSE, sep = " ")
peaksN_me_ac <- coordinates[apply(N_ac,1,all) & apply(N_me,1,all),]
peaksN_me <- coordinates[!apply(N_ac,1,all) & apply(N_me,1,all),]
peaksN_ac <- coordinates[apply(N_ac,1,all) & !apply(N_me,1,all),]

# T-CELLS
coordinates <- read.delim("enhancerAnnotation/peakAnalysis/coordinates.bed",header = FALSE)
T_ac <- read.delim(file = "enhancerAnnotation/peakAnalysis/T_ac_per_sample.tab",header = FALSE, sep = " ")
T_me <- read.delim(file = "enhancerAnnotation/peakAnalysis/T_me_per_sample.tab",header = FALSE, sep = " ")
peaksT_me_ac <- coordinates[apply(T_ac,1,all) & apply(T_me,1,all),]
peaksT_me <- coordinates[!apply(T_ac,1,all) & apply(T_me,1,all),]
peaksT_ac <- coordinates[apply(T_ac,1,all) & !apply(T_me,1,all),]

# DIFFERENCES
extractPeaks <- function(data){as.vector(apply(data,1,function(x){gsub(" ","",paste0(x[1],x[2],x[3],collapse = ""))}))}
N_me_ac <- extractPeaks(peaksN_me_ac)
T_me_ac <- extractPeaks(peaksT_me_ac)
peaksN_me_ac <- peaksN_me_ac[!N_me_ac%in%T_me_ac,]
peaksT_me_ac <- peaksT_me_ac[!T_me_ac%in%N_me_ac,]
N_me <- extractPeaks(peaksN_me)
T_me <- extractPeaks(peaksT_me)
peaksN_me <- peaksN_me[!N_me%in%T_me,]
peaksT_me <- peaksT_me[!T_me%in%N_me,]
N_ac <- extractPeaks(peaksN_ac)
T_ac <- extractPeaks(peaksT_ac)
peaksN_ac <- peaksN_ac[!N_ac%in%T_ac,]
peaksT_ac <- peaksT_ac[!T_ac%in%N_ac,]

# OUTPUT FILES
write.table(peaksT_me_ac, file = "enhancerAnnotation/specificPeaks/T_me_ac",col.names = FALSE,quote = FALSE,row.names = FALSE)
write.table(peaksT_me, file = "enhancerAnnotation/specificPeaks/T_me",col.names = FALSE,quote = FALSE,row.names = FALSE)
write.table(peaksT_ac, file = "enhancerAnnotation/specificPeaks/T_ac",col.names = FALSE,quote = FALSE,row.names = FALSE)
write.table(peaksN_me_ac, file = "enhancerAnnotation/specificPeaks/N_me_ac",col.names = FALSE,quote = FALSE,row.names = FALSE)
write.table(peaksN_me, file = "enhancerAnnotation/specificPeaks/N_me",col.names = FALSE,quote = FALSE,row.names = FALSE)
write.table(peaksN_ac, file = "enhancerAnnotation/specificPeaks/N_ac",col.names = FALSE,quote = FALSE,row.names = FALSE)
