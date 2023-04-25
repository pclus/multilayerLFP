#!/usr/bin/env -S julia --threads 16

cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
push!(LOAD_PATH, pwd())
using NeuropixelAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests

n0=226;
nf=361;

for sub_id in [8,9,10,11,12,13,15,16,18,19,20]
	subject = "suj"*string(sub_id)
	
	if !isdir("../4_outputs/pipeline/"*subject)
	    mkdir("../4_outputs/pipeline/"*subject)
	end
	mpre,mpost = export_matfile("/media/pclusella/Pandora/UPO_data/"*subject*".mat",subject)
	process_data(n0,nf;mpre,mpost);
	prepost_analysis(n0,nf;mpre=mpre,mpost=mpost,foutname=subject*"/");
	
end

