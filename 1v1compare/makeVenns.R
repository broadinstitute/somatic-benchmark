library(plyr)
library(venneuler)

files <- c("out1","out2","out3")

read_csv_filename <- function(filename){
  ret <- read.csv(filename, stringsAsFactors=FALSE, header=FALSE)
  ret$Source <- filename
  ret
}

import.list <- ldply(files, read_csv_filename)
v <- venneuler(import.list)
plot(v)