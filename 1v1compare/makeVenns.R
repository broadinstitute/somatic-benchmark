library(plyr)
library(venneuler)

files <- commandArgs(trailingOnly = TRUE)

read_csv_filename <- function(filename){
  ret <- read.csv(filename, stringsAsFactors=FALSE, header=FALSE)
  ret$Source <- filename
  ret
}

import.list <- ldply(files, read_csv_filename)
v <- venneuler(import.list)
plot(v)
