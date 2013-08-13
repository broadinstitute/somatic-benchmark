package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.function.RetryMemoryLimit

import scala.collection.immutable.Map
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.{MergeSamFiles, SortSam}
import net.sf.samtools.SAMFileHeader.SortOrder
import java.security.InvalidParameterException
import java.io.PrintWriter
import org.apache.commons.io.{FileUtils, IOUtils}
import java.util

class GenerateBenchmark extends QScript with Logging {
    qscript =>

    //TODO implement these as cmdline parameters instead of hard coding them
    val indelFile: File = new File("/humgen/1kg/DCC/ftp/technical/working/20130610_ceu_hc_trio/broad/CEU.wgs.HaplotypeCaller_bi.20130520.snps_indels.high_coverage_pcr_free.genotypes.vcf.gz")
    val referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
    val bam: File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.NA12878.bam")
    val spikeContributorBAM: File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.bam")

    val libDir: File = new File(".")

    @Input(doc = "Directory to locate output files in", shortName = "o", required = false)
    var output_dir: File = new File(libDir)

    @Argument(doc = "Run in test mode which produces a limited set of jobs", required = false)
    var is_test: Boolean = false

    @Argument(doc = "Run without generating spiked data", required = false)
    var no_spike: Boolean = false

    lazy val vcfDataDir = new File(output_dir, "vcf_data")
    lazy val spikeSitesVCF = new File(vcfDataDir, "na12878_ref_NA12891_het_chr1_high_conf.vcf")


    val intervalFile = new File(libDir, "benchmark.interval_list")

    val SAMPLE_NAME_PREFIX = "NA12878.WGS"

    //TODOAn ugly hardcoded hack.  Must eventually be replaced when the number of divisions while fracturing is allowed to be changed from 6.
    var bamNameToFileMap: Map[String, File] = null


    //Used by make_fn_data
    val alleleFractions = Set(0.04, .1, .2, .4, .8)
    val maxDepth = "123456789ABC"
    val depths = for (i <- 1 to maxDepth.length) yield maxDepth.substring(0, i)

    val PIECES = 6

    //the last library in the list is considered the "normal"
    val LIBRARIES = List("Solexa-18484", "Solexa-23661", "Solexa-18483")

    def script() = {

       //fracture bams
        val fractureOutDir = new File(output_dir, "data_1g_wgs")
        val (splitBams, fractureCmds) = FractureBams.makeFractureJobs(bam, referenceFile, LIBRARIES, PIECES, fractureOutDir)
        fractureCmds.foreach(add(_))

        qscript.bamNameToFileMap = splitBams.map((bam: File) => (bam.getName, bam)).toMap

        if( !no_spike ){
            //make vcfs
            val makeVcfs = new MakeVcfs
            makeVcfs.indelFile = qscript.indelFile
            makeVcfs.vcfOutFile = spikeSitesVCF
            add(makeVcfs)

            //use SomaticSpike to create false negative test data
            val makeFnCommands = new FalseNegativeSim(spikeSitesVCF, spikeContributorBAM)
            val (spikedBams,falseNegativeCmds) = makeFnCommands.makeFnSimCmds(alleleFractions, depths)
            falseNegativeCmds.foreach(add(_))
        }

        //merge bams
        val (mergedBams, mergers) = MergeBams.makeMergeBamsJobs(fractureOutDir)
        addAll(mergers)

        //compute depth information for generated bams
        val depthFiles = mergedBams.map(file => swapExt(file.getParent,file,"bam","depth"))

        val sounders = (mergedBams,depthFiles).zipped map( (bamFile,depthFile) => new DepthOfCoverage with GeneratorArguments {
            this.input_file:+= bamFile
            this.out = depthFile
            this.omitDepthOutputAtEachBase = true
        })

        val gatherDepths = new GatherDepths
        gatherDepths.depthFiles = depthFiles
        gatherDepths.coverageFile=new File("collectedCoverage.tsv")

        addAll(sounders)
        add(gatherDepths)


    }



    trait BaseArguments extends CommandLineFunction with RetryMemoryLimit{

    }

    trait GeneratorArguments extends CommandLineGATK with BaseArguments {
        this.reference_sequence = referenceFile
    }

    class FilterByLibrary extends PrintReads with GeneratorArguments {
        @Argument(doc = "library name")
        var library: String = _

        this.memoryLimit = 2
        this.read_filter ++= List("DuplicateRead", "FailsVendorQualityCheck", "UnmappedRead")
        this.simplifyBAM = true

        override def commandLine = super.commandLine + required("-rf", "LibraryRead") + required("--library", library)

    }

    object FractureBams {


        class SplitBam extends CommandLineFunction with BaseArguments {
            @Input(doc = "bam file to copy header from")
            var headerBam: File = _

            @Input(doc = "bam to be split, the one filtered by library")
            var nameSortedBam: File = _

            @Output(doc = "list of output files")
            var outFiles: List[File] = _

            def commandLine = required("%s/splitBam.pl".format(libDir)) +
                required(headerBam) +
                required(nameSortedBam) +
                repeat(outFiles)
        }

        def makeFractureJobs(bam: File, reference: File, libraries: Traversable[String], pieces: Int, outDir: File) = {

            def makeSingleFractureJob(libraryName: String): (List[File], List[CommandLineFunction]) = {

                def getSplitSamNames(library: String, pieces: Int): Traversable[String] = {
                    for (i <- 1 to pieces) yield getSplitFileName(library, i, "sam")
                }

                def getCoordinateSortAndConvertToBam(inputSam: File, outputBam: File): CommandLineFunction = {
                    val sort = new SortSam with BaseArguments
                    sort.memoryLimit = 2
                    sort.input :+= inputSam
                    sort.output = outputBam
                    sort.sortOrder = SortOrder.coordinate
                    sort.createIndex = true
                    sort
                }

                val libraryFiltered = new File(outDir, SAMPLE_NAME_PREFIX + ".original.regional.filtered.%s.bam".format(libraryName))
                val filter = new FilterByLibrary {
                    this.memoryLimit = 2
                    this.library = libraryName
                    this.input_file :+= bam
                    this.out = libraryFiltered
                    this.isIntermediate = true
                }

                val sortedBam = new File(outDir, SAMPLE_NAME_PREFIX + ".original.regional.namesorted.%s.bam".format(libraryName))
                val sort = new SortSam with BaseArguments {
                    this.memoryLimit = 16
                    this.maxRecordsInRam = 4000000
                    this.input :+= libraryFiltered
                    this.output = sortedBam
                    this.sortOrder = net.sf.samtools.SAMFileHeader.SortOrder.queryname
                    this.compressionLevel = 1
                    this.isIntermediate=true
                }

                val split = new SplitBam
                split.headerBam = bam
                split.nameSortedBam = sortedBam



                val splitSams: List[File] = getSplitSamNames(libraryName, pieces).map(new File(outDir, _)).toList
                split.outFiles = splitSams

                val splitBams: List[File] = splitSams.map((sam: File) => swapExt(outDir, sam, "sam", "bam"))

                val converters = (splitSams, splitBams).zipped map {
                    (samFile, outputBam) =>
                        getCoordinateSortAndConvertToBam(samFile, outputBam)
                }

                (splitBams, List(filter, sort, split) ++ converters)
            }

            val (splitBams, cmds) = (for (library <- libraries) yield makeSingleFractureJob(library)).unzip

            (splitBams.flatten, cmds.flatten)
        }
    }

    object MergeBams {
        private val outFileNameTemplate = "NA12878.somatic.simulation.merged.%s.bam"
        private lazy val BAMGROUPS = if (is_test) {
            List("123456789ABC", "DEFGHI")
        } else {
            List("123456789ABC",
                "123456789AB",
                "123456789A",
                "123456789",
                "12345678",
                "1234567",
                "123456",
                "12345",
                "1234",
                "123",
                "12",
                "1",
                "DEFGHI",
                "DEFGH",
                "DEFG",
                "DEF",
                "DE",
                "FG",
                "HI",
                "D")
        }

        def makeMergeBamsJobs(dir: File) = {
            BAMGROUPS.map {
                name =>
                    val mergedFile = new File(dir, outFileNameTemplate.format(name))
                    val inputBams = getBams(name)
                    val merge = new MergeSamFiles

                    merge.memoryLimit = 2
                    merge.input ++= inputBams
                    merge.output = mergedFile
                    merge.createIndex = true
                    merge.USE_THREADING = true
                    (mergedFile, merge)
            }.unzip
        }


    }



    class MakeVcfs extends CommandLineFunction with BaseArguments {
        @Input(doc = "vcf file containing indels to use as true indel sites") var indelFile: File = _
        @Output(doc = "dummy output for queue ordering") var vcfOutFile: File = _

        this.memoryLimit = 33

        def commandLine = "%s/make_vcfs.pl %s %s %s %s %s %s".format(libDir, indelFile, bam,
            spikeContributorBAM, intervalFile,
            referenceFile, vcfDataDir)
    }

    class FalseNegativeSim(spikeSitesVCF: File, spikeInBam: File) {
        val spikedOutputDir = new File(output_dir, "fn_data")

        def makeFnSimCmds(alleleFractions: Traversable[Double], depths: Traversable[String])= {
            val pairs = for {
                fraction <- alleleFractions
                depth <- depths
            } yield makeMixedBam(fraction, depth)
            pairs.unzip
        }

        private def makeMixedBam(alleleFraction: Double, depth: String): (File, CommandLineFunction) = {
            val tumorBams = getBams(depth)
            val outBam = new File(spikedOutputDir, deriveBamName(alleleFraction, depth))
            val outIntervals = swapExt(spikedOutputDir, outBam, "bam", "interval_list")
            val outVcf = swapExt(spikedOutputDir, outBam, "bam", "vcf")

            val spike = new SomaticSpike with GeneratorArguments
            spike.javaMemoryLimit = 4
            spike.simulation_fraction = alleleFraction
            spike.out = outBam
            spike.input_file ++= tumorBams
            spike.spiked_intervals_out = outIntervals
            spike.intervals :+= spikeSitesVCF
            spike.input_file :+= new TaggedFile(spikeContributorBAM, "spike")
            spike.variant = spikeSitesVCF
            spike.spiked_variants = outVcf
            (outBam, spike)
        }

        private def deriveBamName(alleleFraction: Double, depth: String): String = {
            val bamNameTemplate = "NA12878_%s_NA12891_%s_spikein.bam"
            bamNameTemplate.format(depth, alleleFraction)
        }
    }

    /**
     * Returns a list of bams that correspond to the encoded digitString.
     * Each digit in the string maps to specific bam.  Each library is split into multiple files and there is a unique digit
     * assigned to each file.
     * @param digitString
     * @return
     */
    def getBams(digitString: String): List[File] = {
        val bamDigitToNameMap = generateBamMap
        try {
            digitString.map(digit => bamNameToFileMap(bamDigitToNameMap(digit))).toList
        } catch {
            case e: Exception =>
                println(bamNameToFileMap)
                println(bamDigitToNameMap)
                throw e
        }
    }

    def generateBamMap: Map[Char, String] = {
        //Convert the combination of library / string into a unique character
        def calculateDigit(library: String, piece: Int) = {
            val digits = "123456789ABCDEFGHI"
            val index = LIBRARIES.indexOf(library) * PIECES + piece - 1
            digits.charAt(index)
        }


        if (PIECES > 6 || LIBRARIES.size > 3) throw new UnsupportedOperationException("Currently only supported for PIECES <= 6 and LIBRARIES.size <= 3")

        val mappings = for {
            library <- LIBRARIES
            piece <- 1 to PIECES
        } yield (calculateDigit(library, piece), getSplitFileName(library, piece, "bam"))

        mappings.toMap

    }

    def getSplitFileName(library: String, piece: Int, extension: String) = {
        val fileNameTemplate = SAMPLE_NAME_PREFIX + ".somatic.simulation.%s.%03d.%s"
        fileNameTemplate.format(library, piece, extension)
    }

    class GatherDepths extends InProcessFunction {

        @Input(doc="depth files")
        var depthFiles: Seq[File] = Nil

        @Output(doc="coverage summary file")
        var coverageFile: File = _

        def printHeaderAndContents(outputFile: File, header: String, contents: Traversable[String])={
            val writer = new PrintWriter(outputFile)
            try{
                writer.println(header)
                contents.foreach(writer.println)
            } finally {
                IOUtils.closeQuietly(writer)
            }
        }

        def extractAverageCoverage(file:File):String = {
            import scala.collection.JavaConversions._
            val COVERAGE_POSITION = 2
            val lines: util.List[String] = FileUtils.readLines(file)
            val totalLine = lines.find(_.startsWith("Total"))
            val averageCoverage = totalLine.get.split("\t")(COVERAGE_POSITION)
            averageCoverage

        }

        def extractFilename(file:File):String = {
            val splitName = swapExt(file.getParentFile(), file, ".depth.sample_summary","").getName.split('.')
            logger.debug("extractFilename: Extracting name from %s: splits as %s".format(file.getName, splitName))
            splitName(splitName.length-1)
        }

        def run() {
            val summaryFiles = depthFiles.map{ file => new File(file.getAbsolutePath+".sample_summary")}
            val filenames = summaryFiles.map(extractFilename)
            val coverage = summaryFiles.map(extractAverageCoverage )

            val header = "File\tCoverage"
            val results = (filenames, coverage).zipped map( (name, coverage) => "%s\t%s".format(name,coverage) )

            printHeaderAndContents(coverageFile, header, results)

        }
    }

}


