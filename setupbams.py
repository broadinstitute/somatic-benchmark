from subprocess import call
#bam = "/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.NA12878.clean.dedup.recal.bam"
#ref = "/humgen/1kg/reference/human_g1k_v37_decoy.fasta"
bam ="/xchip/cga/home/kcibul/analysis/mutect/public_data/cghub/ebdb53ae-6386-4bc4-90b1-4f249ff9fcdf/C835.HCC1143_BL.4.bam"
ref ="/humgen/gsa-hpprojects/1kg/reference/human_g1k_v37.fasta"
libraries = ["Solexa-18483","Solexa-18484","Solexa-23661"]
interval = "one_gig.interval_list"
pieces = "6"
outDir = "data_1g_wgs"

moduleDir="."




def fractureBams():
	print "Fracturing bams."
	for library in libraries:
		call(["perl",moduleDir+"/fracture_bam.pl" , bam, ref, library, interval, pieces, outDir] )
	print "Done fracturing bams."

def mergeBams():
	print "Merging fractured bams by library."
	call( ["perl", moduleDir+"/merge_bams.pl"] )
	print "Done merging fractured bams."	

fractureBams()
