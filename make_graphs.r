library(ggplot2)
library(reshape)
print("Starting up")

depths <- read.delim("collectedCoverage.tsv")
fp <- read.delim("falsePositiveCounts.tsv")

#Create the output directory
outputdir <- "graphs"
dir.create(outputdir,showWarnings=FALSE)

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

#Create the false positives graph
merged <- addDepths(fp, depths)
qplot(Tumor_Coverage, False_Positives, data=merged,color=Tool, facets=~Normal_Coverage, geom="line" ) + theme_bw()
ggsave(file="graphs/plot.png", height=5, width = 10)


merged <- subset(merged, Normal_Coverage > 25)
qplot(Tumor_Coverage, False_Positives, data=merged,color=Tool, facets=~Normal_Coverage, geom="line" ) + theme_bw()
ggsave(file="graphs/high_coverage_plot.png", height=5, width = 10)

#Create the false negatives graph
fn <- read.delim("diffResults.tsv")
merged <- addDepths(fn, depths)
qplot(Tumor_Coverage, FN, data=merged,color=Tool,  facets=~Fraction, geom="line" ) + theme_bw()
ggsave(file="graphs/fn_plot.png", height=5, width = 10)

merged <- subset(merged, Normal_Coverage > 25)
qplot(Tumor_Coverage, FN, data=merged,color=Tool,  facets=~Fraction, geom="line" ) + theme_bw()
ggsave(file="graphs/fn_high_coverage_plot.png", height=5, width = 10)
