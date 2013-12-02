import org.broadinstitute.sting.queue.extensions.gatk.VariantEval
import org.broadinstitute.sting.queue.QScript


class EvaluateCallSet extends QScript{

    @Input(doc="Somatic mutation file (currently vcf only)", shortName = "V" , fullName = "Variants")
    var variants: File = _

    @Input(doc="Reference fasta file", shortName="R", fullName="reference")
    var reference: File = _

    @Argument(shortName = "dbSNP", doc="dbSNP", required=false)
    val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_137.b37.vcf")

    @Argument(shortName = "hotSpots", doc="Vcf of cancer hotspots, defaults to Mike Lawrences curated list", required = false)
    val cancerHotSpots: File = new File("~/benchmark/data/CosmicCodingMuts_v67_20131024.vcf.gz")

    @Argument(shortName = "goldStandardIndels", doc="Path to gold standard indels", required=false)
    val goldStandardIndels: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.vcf")


    def script()= {
        val countsByAlleleFraction = new Eval(variants, "_allele_fraction", Seq("Sample", "AltReadFraction"), Nil )
        val countsByDbSNPAndFraction = new Eval(variants, "_dbSNP_and_Fraction", Seq("Sample", "AltReadFraction"), Nil){
            this.dbsnp=dbSNP
        }

        val countsByCosmicStatus = new Eval(variants, "_cosmic", Seq("Sample", "AltReadFraction"), Nil){
            this.comp :+= cancerHotSpots
        }

        add(countsByAlleleFraction)
        add(countsByDbSNPAndFraction)
        add(countsByCosmicStatus)

    }

    /*
     * Values we want:
     * Per Set/Sample
     * Allelic fraction
     * stratify by alleleic fraction
     *
     * Insertion / deletion ratio
     * Indel length distribution
     *
     * % dbSNP overlap
     * nonsense / frameshift ratio
     * % cosmic overlap
     * % panCan hotspots
     *
     * mutation rate
     *
     * Ti/TV
     *
     * tumor/normal coverage by allelic fraction
     *
     *
     * Special Splits
     *  By cosmic vs not cosmic
     *      allelic fraction
     */

    class Eval(evalVCF: File, prefix: String, extraStrats: Seq[String], extraEvals: Seq[String]) extends VariantEval {
        this.reference_sequence = reference
        this.eval :+= evalVCF
        //this.dbsnp = dbSNP
        this.doNotUseAllStandardModules = true
        //this.gold = goldStandardIndels
        this.evalModule = List("TiTvVariantEvaluator", "CountVariants", "CompOverlap", "IndelSummary", "IndelLengthHistogram") ++ extraEvals
        this.doNotUseAllStandardStratifications = true
        this.stratificationModule = Seq("EvalRod", "CompRod") ++ extraStrats
        this.memoryLimit = 8
        this.out = swapExt(evalVCF, ".vcf", prefix + ".eval")
    }
}
