#!/bin/bash -l

if [ $# -ne 4 ]
then
  echo "Usage: `basename $0`  normalBam tumorBam reference outputdir"
  exit 1
fi

libdir=/xchip/cga_home/louisb/Strelka
indiv=sample

normalBam=$1
tumorBam=$2
ref=$3
outputdir=$4
#echo "mkdir $outputdir"
#mkdir -p $outputdir

echo "Invoking ${libdir}/runStrelka.sh"
STRELKA_INSTALL_DIR=${libdir}/strelka_workflow_1.0.7

$STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow_cga.pl --normal=$normalBam --tumor=$tumorBam --ref=$ref --config=$libdir/strelka_config_bwa_cgaexome.ini --output-dir=${outputdir}

make -j 4 -C ${outputdir}
echo "Done running Strelka"

echo "Copying strelka output to final.indels.vcf"
cp ${outputdir}/results/passed.somatic.indels.vcf ${outputdir}/final.indels.vcf
cp ${outputdir}/results/passed.somatic.snvs.vcf ${outputdir}/final.snps.vcf

echo "Removing intermediates"
rm -r ${outputdir}/chromosomes

