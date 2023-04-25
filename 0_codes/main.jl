#!/usr/bin/env -S julia --threads 16

cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
push!(LOAD_PATH, pwd())
using NeuropixelAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests

n0=226;
nf=361;

mpre = 9000000
mpost = 7289726

subject = "suj10"
if !isdir("../4_outputs/pipeline/"*subject)
    mkdir("../4_outputs/pipeline/"*subject)
end
mpre,mpost = export_matfile("/media/pclusella/Pandora/UPO_data/"*subject*".mat",subject)

process_data(n0,nf;mpre,mpost);

prepost_analysis(n0,nf;mpre=mpre,mpost=mpost,foutname="suj9/");

2+2

# n=size(n0:nf)[1]
# ns=n-4
# NeuropixelAnalysis.prepost_comparison("bipolar_",ns; mpre=9000000, mpost=mpost, foutname = "suj9/")
