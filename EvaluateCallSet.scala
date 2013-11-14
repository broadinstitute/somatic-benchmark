import org.broadinstitute.sting.queue.QScript


class EvaluateCallSet extends QScript{

    @Input(doc="Somatic mutation file (currently vcf only)", shortName = "V" , fullName = "Variants")
    var variants: File = _

    @Input(doc="Reference fasta file", shortName="R", fullName="reference")
    var reference: File = _

    def script()= {


    }
}
