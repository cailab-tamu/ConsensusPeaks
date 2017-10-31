files <- list.files("greatOutput/",pattern = "^100_Top_Shared_M_ac", full.names = TRUE)
out <- sapply(files,function(file){
  peakname <- gsub(".tsv_$","",gsub("100_Top_Shared_","",basename(file)))
  file <- read.delim(file,header = F, comment.char = "#", skip = 1, stringsAsFactors = FALSE)
  genes <- unique(unlist(strsplit(file[,ncol(file)],",")))
  cbind(peakname,genes)
  
})
