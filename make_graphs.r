library(ggplot2)
library(reshape)
print("Starting up")
depths <- read.delim("collectedCoverage.tsv")
fp <- read.delim("falsePositiveCounts")

merged <- merge(fp, depths, by.x="Tumor", by.y="File")
merged <- merge(merged, depths, by.x="Normal", by.y="File", suffix=c("Tumor","Normal"))

print(merged)
merged$Tumor_Coverage <- merged$CoverageTumor
merged$Normal_Coverage <- merged$CoverageNormal
qplot(Tumor_Coverage, False_Positives, data=merged,color=Tool, facets=~Normal_Coverage, geom="line" ) + theme_bw()
ggsave(file="plot.png", height=5, width = 10)

merged <- subset(merged, Normal_Coverage > 25)
qplot(Tumor_Coverage, False_Positives, data=merged,color=Tool, facets=~Normal_Coverage, geom="line" ) + theme_bw()
ggsave(file="high_coverage_plot.png", height=5, width = 10)
