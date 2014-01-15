# ### Getting input file names  
if( interactive() ){
  #for developement and testing purposes
  print("I see you're running interactively, setting default values")
  outputdir <- "."
  inputfile <- "~/Downloads/An_Histo_Only.final_analysis_set.maf"
} else {
  args <- commandArgs(trailingOnly=TRUE)

  if(length(args)!=2){
    print("Usage: Rscript graphs_from_mafs.R <input.maf> <outputdirectory> ")
    quit()
  }
  inputfile <- args[1]
  print(paste("input maf =", inputfile))
  outputdir <- args[2]
  print(paste("output directory =", outputdir))
 }

if( ! file.exists(inputfile)){
    print("Input maf does not exist.  Exiting")
    quit()
} 
if( ! file.exists(outputdir) ){
    print("Output directory doesn't exist.  Creating it.")
    dir.create(outputdir)
}

###  Loading required libraries
library(ggplot2)
library(plyr)
library(gridExtra)
library(gtools)
library(gdata)


### Defining functions
NAToFalse <- function(x){
    NAToUnknown(x, unknown=FALSE, force=TRUE)
}

sort_chromosomes <- function(df){
  return(factor(df$Chromosome, mixedsort(df$Chromosome))) 
}


save_with_name <- function(name,height=10, width=10) {
  name_pieces <- c(outputdir, "/", name, ".pdf")
  filename <- paste(name_pieces, collapse='')
  print(paste("Saving ",filename, sep=''))
  ggsave(file=filename, height=height, width=width, units="in", limitsize=FALSE)  
}

cosmic_or_dbsnp <- function(is_cosmic, is_dbsnp){
  if (is_cosmic)
  {
    return("COSMIC")
  } else if(is_dbsnp){
    return("dbSNP")
  } else return("Unclassified")
}

calc_length <- function( ref, tumor, type){
  if(type == "INS"){
    return(nchar(tumor))
  } else if( type =="DEL"){
    return(nchar(ref))
  } else {
    return(0)
  }
}


### Preparing data

maf <- read.table(file=inputfile,header=TRUE, quote='', sep="\t", stringsAsFactors=FALSE)

samples <- length( unique(maf$Tumor_Sample_Barcode)) 

maf$Chromosome <- sort_chromosomes(maf)

maf$allele_fraction <- maf$t_alt_count / (maf$t_alt_count+maf$t_ref_count)

maf$Matches_COSMIC_Mutation <- NAToFalse( ! maf$COSMIC_overlapping_mutations == "")

maf$Overlaps_DB_SNP_Site <- NAToFalse(! maf$dbSNP_RS =="" )

maf$Classification <-  mapply(cosmic_or_dbsnp,maf$Matches_COSMIC_Mutation, maf$Overlaps_DB_SNP_Site)

maf$Indel_Length = mapply( calc_length, maf$Reference_Allele, maf$Tumor_Seq_Allele2, maf$Variant_Type)
perc <- ddply(maf, "Tumor_Sample_Barcode", summarise, 
              percent_dbSNP = sum(Overlaps_DB_SNP_Site)/length(Overlaps_DB_SNP_Site), 
              percent_COSMIC= sum(Matches_COSMIC_Mutation)/length(Matches_COSMIC_Mutation),
              samples="all")


maf <- mutate(maf, Tumor_Depth = t_alt_count+t_ref_count)

#### Making graphs
plot_percentage_and_count <- function(df, variable, name, outputdir){

  #sort factor by length
  df$Tumor_Sample_Barcode <- with(df, reorder(Tumor_Sample_Barcode, Tumor_Sample_Barcode, function(x) length(x)))

  percent <- ggplot(df, aes(x = Tumor_Sample_Barcode)) + geom_bar(aes_string(fill = variable ), position = 'fill') +
    coord_flip() +theme_bw(base_family='Helvetica')+ scale_fill_brewer(palette="Paired") +
    theme(legend.position = "none") +labs(y="Percent")

  counts <- ggplot(df, aes(x = Tumor_Sample_Barcode)) + geom_bar(aes_string(fill = variable)) + coord_flip() +
    theme_bw(base_family='Helvetica')+ scale_fill_brewer(palette="Paired") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank())

  g <- arrangeGrob(percent, counts, nrow=1)
  name_pieces <- c(outputdir, "/", name, ".pdf")
  filename <- paste(name_pieces, collapse='')
  print(paste("Saving ",filename, sep=''))
  ggsave(file=filename, g, height=max(samples/8,4), width=10, units="in", limitsize=FALSE)
}




ggplot(maf, aes(x=Tumor_Depth, y = allele_fraction)) + geom_point(size=.5, alpha=.1) + theme_bw() + scale_x_log10() 
save_with_name("allele_fraction_vs_Tumor_Depth")
ggplot(maf, aes(x=Tumor_Depth, y = allele_fraction, color=Classification)) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + geom_point(size=.5, alpha=.1) + geom_density2d(size=.4) + theme_bw() + scale_x_log10() 
save_with_name("allele_fraction_vs_Tumor_Depth_with_Density")


ggplot(maf, aes(x=Tumor_Depth, y = allele_fraction)) + facet_wrap(facets=~Variant_Type, drop=TRUE) + geom_point(size=.5, alpha=.1) +  theme_bw()+ scale_x_log10()
save_with_name("allele_fraction_vs_Tumor_Depth_by_VariantType")
ggplot(maf, aes(x=Tumor_Depth, y = allele_fraction, color=Classification)) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + facet_wrap(facets=~Variant_Type, drop=TRUE) + geom_point(size=.5, alpha=.1) + geom_density2d(alpha = .8) + theme_bw()+ scale_x_log10()
save_with_name("allele_fraction_vs_Tumor_Depth_by_VariantType_with_Density")


qplot(data=maf,x=allele_fraction, fill=Classification) + theme_bw()
save_with_name("allele_fraction_all_samples") 

qplot(data=maf, x=allele_fraction, fill=Classification) + facet_wrap(facets=~Tumor_Sample_Barcode, ncol=4)+theme_bw() + theme(strip.text.x = element_text(size=6))
save_with_name("allele_fraction_by_sample", height=max(samples/4,4))

qplot(data=maf, x=allele_fraction, fill=Classification) + facet_wrap(facets=~Tumor_Sample_Barcode, ncol=4, scales="free_y") +theme_bw() + theme(strip.text.x = element_text(size=6))
save_with_name("allele_fraction_by_sample_normalized", height=max(samples/4,4), width=10)

plot_percentage_and_count(maf, "Classification", "COSMIC_overlap_by_sample", outputdir)

ggplot(data=maf)+ geom_bar(subset=.(Variant_Type == "DEL"), aes(x=Indel_Length,y=-..count..,fill=Variant_Classification,stat="identity")) + geom_bar(subset=.(Variant_Type=="INS"), aes(x=Indel_Length,y=..count.., fill=Variant_Classification, stat="identity")) + ylab("Deletions - Insertions")
save_with_name("stacked_indel_lengths", height=5, width=7)

#subset to indels
indels_only <- maf[maf$Variant_Type %in% c("DEL","INS"), ]

ggplot(data=indels_only)+ geom_bar( aes(x=Indel_Length,y=..count..,fill=Variant_Classification,stat="identity")) + facet_wrap(facets=~Variant_Type, drop=TRUE)
save_with_name("indel_lengths_by_type", height=5, width=7)

plot_percentage_and_count(indels_only, "Variant_Classification", "indels_by_type", outputdir)



