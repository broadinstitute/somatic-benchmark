library(gsalib)
library(reshape)
library(ggplot2)

d = gsa.read.gatkreport("combo.gsareport")
countVar <- d$CountVariants
qplot(AltReadFraction,nVariantLoci,data=countVar, geom="line", facets=~EvalRod)
ggsave(file="graphs/readFractionDistribution.png", height=5, width = 10)