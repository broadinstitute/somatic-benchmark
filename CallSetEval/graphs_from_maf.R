#### Getting input file names  
if( "--interactive" %in% commandArgs()){
  #for developement and testing purposes
  print("I see you're running interactively, setting default values")
  outputdir <- "."
  inputfile <- "~/Downloads/AN_TCGA_LUAD_PAIR_capture_freeze_FINAL_230.final_analysis_set.maf"
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

### Defining functions
sort_chromosomes <- function(df){
  return(factor(df$Chromosome, mixedsort(levels(df$Chromosome)))) 
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

calc_length <- function( ref, tumor){
  ref <- as.character(ref)
  tumor <- as.character(tumor)
  if(ref =="-" && tumor != "-"){
    return(nchar(tumor))
  } else if( ref != "-" && tumor == "-"){
    return( - nchar(ref))
  } else return( 0 )
}


### Preparing data

maf <- read.delim(file=inputfile,header=TRUE)

samples <- length( unique(maf$Tumor_Sample_Barcode)) 

maf$Chromosome <- sort_chromosomes(maf)

maf$allele_fraction <- maf$t_alt_count / (maf$t_alt_count+maf$t_ref_count)
maf$Matches_COSMIC_Mutation <- ! maf$COSMIC_overlapping_mutations == ""
maf$Overlaps_DB_SNP_Site <- ! maf$dbSNP_RS ==""

maf$Classification <-  mapply(cosmic_or_dbsnp,maf$Matches_COSMIC_Mutation, maf$Overlaps_DB_SNP_Site)

indels <- mutate(maf, Indel_Length = calc_length(Reference_Allele, Tumor_Seq_Allele1))
perc <- ddply(maf, "Tumor_Sample_Barcode", summarise, 
              percent_dbSNP = sum(Overlaps_DB_SNP_Site)/length(Overlaps_DB_SNP_Site), 
              percent_COSMIC= sum(Matches_COSMIC_Mutation)/length(Matches_COSMIC_Mutation),
              samples="all")



#### Making graphs
plot_percent_cosmic_and_dbsnp_overlap <- function(){
  percent <- ggplot(maf, aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Classification), position = 'fill') + 
    coord_flip() +theme_bw(base_family='Helvetica')+ scale_fill_brewer(palette="Paired") + 
    theme(legend.position = "none") +labs(y="Percent")
  
  counts <- ggplot(maf, aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Classification)) + coord_flip() + 
    theme_bw(base_family='Helvetica')+ scale_fill_brewer(palette="Paired") + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank())
  
  g <- arrangeGrob(percent, counts, nrow=1)
  name_pieces <- c(outputdir, "/", "COSMIC_overlap_by_sample", ".pdf")
  filename <- paste(name_pieces, collapse='')
  print(paste("Saving ",filename, sep=''))
  ggsave(file=filename, g, height=max(samples/8,4), width=10, units="in", limitsize=FALSE)  
}




qplot(data=perc, samples, percent_dbSNP, geom="boxplot") + theme_bw() 
save_with_name("Overall_Cosmic_Overlap", height=5, width=3)

qplot(data=perc, samples, percent_COSMIC, geom="boxplot") + theme_bw()
save_with_name("Overall_dbSnp_Overlap", height=5, width=3)

maf <- mutate(maf, Tumor_Depth = t_alt_count+t_ref_count)
qplot(data=maf, x=Tumor_Depth, y = allele_fraction, color = Variant_Type) + theme_bw()
save_with_name("allele_fraction_vs_Tumor_Depth")


qplot(data=maf, x=Tumor_Depth, y = allele_fraction, facets = ~Variant_Type) + theme_bw()
save_with_name("allele_fraction_vs_Tumor_Depth_by_VariantType")

qplot(data=maf,x=allele_fraction, fill=Classification) + theme_bw()
save_with_name("allele_fraction_all_samples") 

qplot(data=maf, x=allele_fraction, fill=Classification) + facet_wrap(facets=~Tumor_Sample_Barcode, ncol=4)+theme_bw() + theme(strip.text.x = element_text(size=6))
save_with_name("allele_fraction_by_sample", height=max(samples/4,4))

qplot(data=maf, x=allele_fraction, fill=Classification) + facet_wrap(facets=~Tumor_Sample_Barcode, ncol=4, scales="free_y") +theme_bw() + theme(strip.text.x = element_text(size=6))
save_with_name("allele_fraction_by_sample_normalized", height=max(samples/4,4), width=10)

plot_percent_cosmic_and_dbsnp_overlap()
