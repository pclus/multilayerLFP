#!/usr/bin/env -S julia --threads 16

cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
push!(LOAD_PATH, pwd())
using NeuropixelAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests

sub_id=ARGS[1]
subject = "suj"*string(sub_id)

# n0=0;
# nf=383;

layers = readdlm("../1_data/layers.dat",skipstart=5)
j = tryparse(Int64,sub_id)
col = findall(x -> x==j,layers[1,:])
if isempty(col)
    print("Layers not found for subject "*string(sub_id)*"\n")
    n0 = 222; nf = 363
else
    layers = Int.(384.0.-layers[2:end,col])[:,1] 
    n0 = layers[end];
    nf = layers[1]+1;
end


	
if !isdir("../4_outputs/pipeline/"*subject)
    mkdir("../4_outputs/pipeline/"*subject)
end

# Only required once if the files are saved
# mpre,mpost = export_matfile("/media/pclusella/Pandora/UPO_data/"*subject*".mat",subject)
# process_data(n0,nf;mpre,mpost);
mpre = 9000000;
mpost = 9000000;
prepost_analysis(n0,nf,subject;mpre=mpre,mpost=mpost);
