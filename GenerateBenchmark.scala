package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.Logging

import scala.collection.immutable.Map
import org.broadinstitute.sting.queue.extensions.gatk.{TaggedFile, SomaticSpike, CommandLineGATK, PrintReads}
import org.broadinstitute.sting.queue.extensions.picard.{MergeSamFiles, SortSam}
import net.sf.samtools.SAMFileHeader.SortOrder

class GenerateBenchmark extends QScript with Logging {
    qscript =>

    //TODO implement these as cmdline parameters instead of hard coding them
    val indelFile: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf")
    val referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
    val bam: File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.NA12878.bam")
    val spikeContributorBAM: File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.bam")

    val libDir: File = new File(".")

    @Input(doc = "Directory to locate output files in", shortName = "o", required = false)
    var output_dir: File = new File(libDir)

    lazy val vcfDataDir = new File(output_dir, "vcf_data")
    lazy val spikeSitesVCF = new File(vcfDataDir, "na12878_ref_NA12891_het_chr1_high_conf.vcf")


    val intervalFile = new File(libDir, "chr20.interval_list")


    val prefix = "chr1"

    //TODOAn ugly hardcoded hack.  Must eventually be replaced when the number of divisions while fracturing is allowed to be changed from 6.
    lazy val bamMapFile = new File(libDir, "louis_bam_1g_info.txt")
    var bamDigitToNameMap: Map[Char, String] = null
    var bamNameToFileMap: Map[String, File] = null


    //Used by make_fn_data
    val alleleFractions = Set(0.04, .1, .2, .4, .8)
    val maxDepth = "123456789ABC"
    val depths = for (i <- 1 to maxDepth.length) yield maxDepth.substring(0, i)

    val PIECES = 6
    val LIBRARIES = List("Solexa-18483", "Solexa-18484", "Solexa-23661")


    def script() = {
        qscript.bamDigitToNameMap = loadBamMap(bamMapFile)


        //make vcfs
        val makeVcfs = new MakeVcfs
        makeVcfs.indelFile = qscript.indelFile
        makeVcfs.vcfOutFile = spikeSitesVCF
        add(makeVcfs)

        //fracture bams
        val fractureOutDir = new File(output_dir, "data_1g_wgs")
        val (splitBams, fractureCmds) = FractureBams.makeFractureJobs(bam, referenceFile, LIBRARIES, intervalFile, PIECES, fractureOutDir)
        fractureCmds.foreach(add(_))

        qscript.bamNameToFileMap = splitBams.map((bam: File) => (bam.getName, bam)).toMap

        //use SomaticSpike to create false negative test data
        val makeFnCommands = new FalseNegativeSim(spikeSitesVCF, spikeContributorBAM)
        val falseNegativeCmds = makeFnCommands.makeFnSimCmds(alleleFractions, depths)
        falseNegativeCmds.foreach(add(_))

        //merge bams
        val mergers = MergeBams.makeMergeBamsJobs(fractureOutDir)
        mergers.foreach(add(_))

    }

    trait BaseArguments extends CommandLineFunction {

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
        val FILE_NAME_PREFIX = "NA12878.WGS"

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

        def makeFractureJobs(bam: File, reference: File, libraries: Traversable[String], interval: File, pieces: Int, outDir: File) = {

            def makeSingleFractureJob(libraryName: String): (List[File], List[CommandLineFunction]) = {
                def getSplitBamNames(library: String, pieces: Int): Traversable[String] = {
                    val outmask = FILE_NAME_PREFIX + ".somatic.simulation.%s.%03d.sam"
                    for (i <- 1 to pieces) yield outmask.format(library, i)
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

                val libraryFiltered = new File(outDir, FILE_NAME_PREFIX + ".original.regional.filtered.%s.bam".format(libraryName))
                val filter = new FilterByLibrary {
                    this.memoryLimit = 2
                    this.library = libraryName
                    this.input_file :+= bam
                    this.out = libraryFiltered
                    this.intervals :+= interval
                    this.isIntermediate = true
                }

                val sortedBam = new File(outDir, FILE_NAME_PREFIX + ".original.regional.namesorted.%s.bam".format(libraryName))
                val sort = new SortSam with BaseArguments {
                    this.memoryLimit = 16
                    this.maxRecordsInRam = 4000000
                    this.input :+= libraryFiltered
                    this.output = sortedBam
                    this.sortOrder = net.sf.samtools.SAMFileHeader.SortOrder.queryname
                    this.compressionLevel = 1
                }

                val split = new SplitBam
                split.headerBam = bam
                split.nameSortedBam = sortedBam

                val splitSams: List[File] = getSplitBamNames(libraryName, pieces).map(new File(outDir, _)).toList
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
        private val BAMGROUPS = List(
            "123456789ABC",
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
            "D"
        )

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
                    merge
            }
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

        def makeFnSimCmds(alleleFractions: Traversable[Double], depths: Traversable[String]): Traversable[CommandLineFunction] = {
            for {
                fraction <- alleleFractions
                depth <- depths
            } yield makeMixedBam(fraction, depth)

        }

        private def makeMixedBam(alleleFraction: Double, depth: String): CommandLineFunction = {
            val tumorBams = getBams(depth)
            val outBam = new File(spikedOutputDir, deriveBamName(alleleFraction, depth))
            val outIntervals = swapExt(spikedOutputDir, outBam, "bam", "interval_list")

            val spike = new SomaticSpike with GeneratorArguments
            spike.javaMemoryLimit = 4
            spike.simulation_fraction = alleleFraction
            spike.out = outBam
            spike.input_file ++= tumorBams
            spike.spiked_intervals_out = outIntervals
            spike.intervals :+= spikeSitesVCF
            spike.input_file :+= new TaggedFile(spikeContributorBAM, "spike")
            spike
        }

        private def deriveBamName(alleleFraction: Double, depth: String): String = {
            val bamNameTemplate = "NA12878_%s_NA12891_%s_spikein.bam"
            bamNameTemplate.format(depth, alleleFraction)
        }
    }


    def getBams(hexDigitString: String): List[File] = {
        try {
            hexDigitString.map(digit => bamNameToFileMap(bamDigitToNameMap(digit))).toList
        } catch {
            case e: Exception =>
                println(bamNameToFileMap)
                println(bamDigitToNameMap)
                throw e
        }
    }

    def loadBamMap(bamMapFile: File): Map[Char, String] = {
        println("loading file")
        def splitLine(line: String): Option[(Char, String)] = {
            try {
                val segments = line.split("\\s+")
                val char = segments(0).charAt(0)
                val name = segments(1).trim()
                Some(char, name)
            } catch {
                case e: Throwable =>
                    None
            }
        }
        val fileLines = io.Source.fromFile(bamMapFile).getLines()
        val map: Map[Char, String] = fileLines.map(splitLine).flatten.toMap

        map
    }
}


