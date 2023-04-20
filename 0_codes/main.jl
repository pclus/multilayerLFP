#!/usr/bin/env -S julia --threads 16

cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
push!(LOAD_PATH, pwd())
using NeuropixelAnalysis
# ,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
# plotlyjs()

n0=226;
nf=361;
mpost = 7289726
mpre = 9000000
# process_data(n0,nf,mpre,mpost);
# prepost_analysis(n0,nf);

data = "bipolar_"
n=size(n0:nf)[1]
n0 = n-4
a,b,c,d,q1,q2,q3,q4 = NeuropixelAnalysis.prepost_comparison(data,n0;mpre=mpre,mpost=mpost)

id = 0
dt = 0.0004
fl = "bipolar_pre"
t,f,tfhm = timefreq(id, fl; m=mpost)