library(ggplot2)
library(reshape)

print("Starting up")
###################Function definitions####################################

#Merges the false positive or false negative table with the depths table to add tumor and normal depths
addDepths <- function(data, depths){
    merged <- merge(data, depths, by.x="Tumor", by.y="File")
    merged <- merge(merged, depths, by.x="Normal", by.y="File", suffix=c("Tumor","Normal"))
    merged$Tumor_Coverage <- merged$CoverageTumor
    merged$Normal_Coverage <- merged$CoverageNormal
    merged$CoverageNormal <- NULL
    merged$CoverageTumor <- NULL
    return(merged)
}

#Create the false positives graphs
graph_false_positives <- function(outputdir, depths, fpCounts) {
  fp <- read.delim(fpCounts)
  
  merged <- addDepths(fp, depths)
  qplot(Tumor_Coverage, False_Positives, data=merged,color=Tool, facets=~Normal_Coverage, geom="point" ) + theme_bw()
  ggsave(file=paste(outputdir,"/plot.png",sep=""), height=5, width = 10)
  
  
  merged <- subset(merged, Normal_Coverage > 25)
  qplot(Tumor_Coverage, False_Positives, data=merged,color=Tool, facets=~Normal_Coverage, geom="point" ) + theme_bw()
  ggsave(file=paste(outputdir,"/high_coverage_plot.png",sep=""), height=5, width = 10)
}
       
#Create the false negatives graphs
graph_false_negatives <- function(outputdir, depths, diffResults){
  fn <- read.delim(diffResults)
  
  merged <- addDepths(fn, depths)
  qplot(Tumor_Coverage, FN, data=merged,color=Tool,  facets=~Fraction, geom="point" ) + theme_bw()
  ggsave(file=paste(outputdir,"/fn_plot.png",sep=""), height=5, width = 10)
  
  merged <- subset(merged, Normal_Coverage > 25)
  qplot(Tumor_Coverage, FN, data=merged,color=Tool,  facets=~Fraction, geom="point" ) + theme_bw()
  ggsave(file=paste(outputdir,"/fn_high_coverage_plot.png",sep=""), height=5, width = 10)
}
###########################################################################


#Get command line arguments
#usage <outputdir> <fpCounts.tsv> <diffresults.tsv>

outputdir <- commandArgs(trailingOnly=TRUE)[1]
fpCounts <- commandArgs(trailingOnly=TRUE)[2]
fnCounts <- commandArgs(trailingOnly=TRUE)[3]

#Create the output directory
dir.create(outputdir,showWarnings=FALSE)

#Load the depths information
depths <- read.delim("collectedCoverage.tsv")

#Make graphs
if(file.exists(fpCounts) ) {
    print("Drawing false positive graphs")
    graph_false_positives(outputdir, depths, fpCounts)
}

if(file.exists(fnCounts) ) {
    print("Drawing false negative graphs")
    graph_false_negatives(outputdir, depths, fnCounts)
}
