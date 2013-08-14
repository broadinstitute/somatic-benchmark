library(ggplot2)
library(reshape)
print("Starting up")
depths <- read.delim("collectedCoverage.tsv")
fp <- read.delim("falsePositiveCounts")

addDepths <- function(data, depths){
    merged <- merge(data, depths, by.x="Tumor", by.y="File")
    merged <- merge(merged, depths, by.x="Normal", by.y="File", suffix=c("Tumor","Normal"))
    merged$Tumor_Coverage <- merged$CoverageTumor
    merged$Normal_Coverage <- merged$CoverageNormal
    merged$CoverageNormal <- NULL
    merged$CoverageTumor <- NULL
    return(merged)
}

print(merged)

merged <- addDepths(fp, depths)
qplot(Tumor_Coverage, False_Positives, data=merged,color=Tool, facets=~Normal_Coverage, geom="line" ) + theme_bw()
ggsave(file="plot.png", height=5, width = 10)

merged <- subset(merged, Normal_Coverage > 25)
qplot(Tumor_Coverage, False_Positives, data=merged,color=Tool, facets=~Normal_Coverage, geom="line" ) + theme_bw()
ggsave(file="high_coverage_plot.png", height=5, width = 10)

fn <- read.delim("diffResults.tsv")
merged <- addDepths(fn, depths)
fn 
qplot(Tumor_Coverage, FN, data=merged,color=Tool,  facets=~Fraction, geom="line" ) + theme_bw()
ggsave(file="fn_plot.png", height=5, width = 10)

merged <- subset(merged, Normal_Coverage > 25)
qplot(Tumor_Coverage, FN, data=merged,color=Tool,  facets=~Fraction, geom="line" ) + theme_bw()
ggsave(file="fn_high_coverage_plot.png", height=5, width = 10)


