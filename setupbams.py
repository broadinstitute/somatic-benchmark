from subprocess import check_call
bam = "/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.bam"
ref = "/humgen/1kg/reference/human_g1k_v37_decoy.fasta"
libraries = ["Solexa-18483","Solexa-18484","Solexa-23661"]
interval = "chr20.interval_list"
pieces = "6"
outDir = "data_1g_wgs"

moduleDir="."




def fractureBams():
	print "Fracturing bams."
	for library in libraries:
		check_call(["perl",moduleDir+"/fracture_bam.pl" , bam, ref, library, interval, pieces, outDir] )
	print "Done fracturing bams."

def mergeBams():
	print "Merging fractured bams by library."
    	check_call( ["perl", moduleDir+"/merge_bams.pl"] )
	print "Done merging fractured bams."	

def main():
    fractureBams()
    mergeBams()
