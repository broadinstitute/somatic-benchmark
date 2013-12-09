%  lego_plot_wrapper.m

%  Created by Mara Rosenberg on 11/21/13.
%

function lego_plot_wrapper(maf, OUT)

%Need to make sure that this paths are added (this code should be compiled at some point so we don't have to do this - but I figured for development purposes to keep as is)
%addpath('/xchip/cga_home/mara/CancerGenomeAnalysis/trunk/matlab')
%addpath('/xchip/cga_home/mara/CancerGenomeAnalysis/trunk/matlab/seq')
%addpath('/xchip/cga_home/mara/CancerGenomeAnalysis/trunk/matlab/mike')

P.C1='exome';
P.zscale=true;
P.printrates=true

[XY,X,C1]=plotMutationSpectrumCategLegos(maf,OUT,P);

