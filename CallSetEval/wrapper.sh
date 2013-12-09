#!/bin/sh

. /broad/tools/scripts/useuse

reuse R-2.15
reuse .matlab_2013a_mcr

LIBDIR=$1
MAF=$2

echo "Making lego plots from $MAF"
$LIBDIR/lego_plots/lego_plot_wrapper $MAF legos

echo "Running graphs_from_maf.R on $MAF"
Rscript $LIBDIR/graphs_from_maf.R $MAF .


