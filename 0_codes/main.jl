#!/usr/bin/env -S julia --threads 16

push!(LOAD_PATH, "/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
# plotlyjs()

n0=226;
nf=361;
# process_data(n0,nf);
prepost_analysis(n0,nf);

2+2
