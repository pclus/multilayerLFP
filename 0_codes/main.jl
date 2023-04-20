#!/usr/bin/env -S julia --threads 16

cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
push!(LOAD_PATH, pwd())
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
# plotlyjs()

n0=226;
nf=361;
mpost = 7289726
mpre = 9000000
process_data(n0,nf,mpre,mpost);
# prepost_analysis(n0,nf);

