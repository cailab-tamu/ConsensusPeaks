# Splt files by line
for(file in list.files("consensusPeaks/", full.names = TRUE)){
  fileLines <- readLines(file)
  for (line in fileLines){
    newFileName <- paste0(basename(file),"-",gsub(" ","_",line))
    writeLines(text = line,con = paste0("greatInput/",newFileName))
  }
}

# Generate the queries
files <- list.files("greatInput/")
httpString <- sapply(files,function(file){
  paste0("curl 'http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=hg19&requestName=TAMUannotationBLUEPRINTpeaks&requestSender=Daniel+Osorio+dcosorioh_at_tamu.edu&requestURL=https://raw.githubusercontent.com/cailab-tamu/ConsensusPeaks/master/greatInput/",file,"' -o ", file,".tsv")
})
writeLines(text = httpString,con = "curlGREAT.sh", sep = "\n")

# Run the queries
system("nohup parallel -j 5 <- curlGREAT.sh &")

# Process the outputs
initFiles <- list.files("consensusPeaks/")
for (pattern in initFiles){
  files <- list.files("greatOutput/", pattern = paste0("^",pattern) , full.names = TRUE)
  annotations <- sapply(files, function(file){
    print(file)
    if(!all(grepl("^#",readLines(file)[-1]))){
      fileContent <- read.delim(file,header = FALSE, sep = "\t", comment.char = "#", skip = 1, stringsAsFactors = FALSE)
      annotations <-paste0(fileContent[,3], collapse = "\t")
    } else {
      annotations <- ""
    }
    coordinates <- gsub('.tsv','',unlist(strsplit(unlist(strsplit(basename(file),"-"))[2],"_")))
    allInfo <- paste0(c(coordinates,annotations),collapse = "\t")
    allInfo
  })
  writeLines(annotations,con = paste0("annotatedPeaks/",pattern))
}
