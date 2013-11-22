library(ggplot2)



save_with_name <- function(name,height=10, width=10) {
  name_pieces <- c(outputdir, "/", name, ".pdf")
  filename <- paste(name_pieces, collapse='')
  print(paste("Saving ",filename, sep=''))
  ggsave(file=filename, height=height, width=width, units="in", limitsize=FALSE)  
}

if( "--interactive" %in% commandArgs()){
  print("I see you're running interactively, setting default values")
  outputdir <- "."
  inputfile <- "~/Downloads/AN_TCGA_LUAD_PAIR_capture_freeze_FINAL_230.final_analysis_set.maf"
} else {
  print("reading values from the command line")
  outputdir <- commandArgs(trailingOnly=TRUE)[1]
  inputfile <- commandArgs(trailingOnly=TRUE)[2]
}


maf <- read.delim(file=inputfile,header=TRUE)

samples <- length( unique(maf$Tumor_Sample_Barcode)) 


maf$allele_fraction <- maf$t_alt_count / (maf$t_alt_count+maf$t_ref_count)
maf$Matches_COSMIC_Mutation <- ! maf$COSMIC_overlapping_mutations == ""
maf$Overlaps_DB_SNP_Site <- ! maf$dbSNP_RS ==""

 qplot(data=maf, x=allele_fraction) + theme_bw()
 save_with_name("allele_fraction_all")
# 
 qplot(data=maf, x=allele_fraction) + facet_wrap(facets=~Tumor_Sample_Barcode, ncol=4)+theme_bw() + theme(strip.text.x = element_text(size=8))
 save_with_name("allele_fraction_by_sample", height=samples/4)

#qplot(data=maf, x=allele_fraction, fill=Overlaps_DB_SNP_Site) + theme_bw()

ggplot(maf, aes(x = allele_fraction)) + geom_bar(aes(fill = Overlaps_DB_SNP_Site), position = 'fill')
ggsave(file=paste(outputdir,"/allele_fraction_all.pdf",sep=""))

qplot(data=maf, x=allele_fraction, fill=Matches_COSMIC_Mutation) + theme_bw()
ggsave(file=paste(outputdir,"/fraction_by_cosmic_overlap.pdf",sep=""))

ggplot(maf, aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Overlaps_DB_SNP_Site), position = 'fill') + coord_flip()
save_with_name("dbSnpOverlap_by_sample",height=samples/8)

ggplot(maf, aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Matches_COSMIC_Mutation), position = 'fill') + coord_flip()
save_with_name("COSMIC_overlap_by_sample",height=samples/8)

calc_length <- function( ref, tumor){
  ref <- as.character(ref)
  tumor <- as.character(tumor)
  if(ref =="-" && tumor != "-"){
    return(nchar(tumor))
  } else if( ref != "-" && tumor == "-"){
    return( - nchar(ref))
  } else return( 0 )
}


indels <- mutate(maf, Indel_Length = calc_length(Reference_Allele, Tumor_Seq_Allele1))
perc <- ddply(maf, "Tumor_Sample_Barcode", summarise, 
              percent_dbSNP = sum(Overlaps_DB_SNP_Site)/length(Overlaps_DB_SNP_Site), 
              percent_COSMIC= sum(Matches_COSMIC_Mutation)/length(Matches_COSMIC_Mutation),
              samples="all")

qplot(data=perc, samples, percent_dbSNP, geom="boxplot")
save_with_name("Overall_Cosmic_Overlap", height=5, width=3)

qplot(data=perc, samples, percent_COSMIC, geom="boxplot")
save_with_name("Overall_dbSnp_Overlap", height=5, width=3)