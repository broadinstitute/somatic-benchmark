library(gsalib)
library(reshape)
library(ggplot2)

d = gsa.read.gatkreport("combo.gsareport")
countVar <- d$CountVariants
countVar$EvalRod <- countVar$EvalRod/2

vline.data <- data.frame(z = unique(countVar$EvalRod), EvalRod=unique(countVar$EvalRod))


p <- qplot(AltReadFraction,nVariantLoci,data=countVar, geom="point", facets=~EvalRod) + theme_bw()
p + geom_vline(aes(xintercept = z), vline.data)
ggsave(file="graphs/altAlleleFraction.png", height=5, width = 10)
