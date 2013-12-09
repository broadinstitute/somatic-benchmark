#!/bin/sh

. /broad/tools/scripts/useuse

reuse R-2.15

LIBDIR=$1
MAF=$2

Rscript $1/graphs_from_maf.R $2 .
