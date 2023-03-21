#!/usr/bin/env -S julia --threads 16

push!(LOAD_PATH, "/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
using FFTW

function spectral_stationarity(id,fl,Δt,mov_filter::String)

    rate = 2500;
    dt = 1/rate;
    t, y = read_channel(id, 4e-4, 900, fl);
    T  = length(y)/rate;
    ns = Int(T / Δt);

    if !isempty(mov_filter)
        dts = T/ns;
        ts = 0.5*dts:dts:T
        tr,yr = movfilter(ts,reshape(y,:,ns),mov_filter,Δt);
        y=reshape(yr,:,1)[:,1];
    end

    S = rfft(y); 
    # f = rfftfreq(length(t), 1.0/dt); 
    S0 = @. abs.(S)*exp(im*rand()*2*π);
    y0 = irfft(S0,2*length(S0)-2);

    t, f, tfhm = timefreq(y,Δt);
    smean = mean(tfhm,dims=2)[:,1];
    ts,f,tfhm0 = timefreq(y0,Δt);
    smean0 = mean(tfhm0,dims=2)[:,1];

    ns = length(ts)
    # m = [logspectral_dist(tfhm[:,i],tfhm[:,j],f) for i in 1:ns,j in 1:ns]
    # m0 = [logspectral_dist(tfhm0[:,i],tfhm0[:,j],f) for i in 1:ns,j in 1:ns]
    # heatmap(m,clim=(0.0,0.45))

    q = [logspectral_dist(tfhm[:,i],smean,f) for i in 1:ns]
    q0 = [logspectral_dist(tfhm0[:,i],smean0,f) for i in 1:ns]

    return q,q0,tfhm,tfhm0;
end

id=100; fl="cortex_pre"; Δt=10.0; mov_filter="pre";
q,q0,tfhm,tfhm0 = spectral_stationarity(id,fl,Δt,mov_filter)
qr,qr0,tfhmr,tfhmr0 = spectral_stationarity(id,fl,Δt,"")

# plots
plot([q,q0],lt=:scatter,legend=:false)
plot!([qr,qr0],lt=:scatter,legend=:false)

# tests
ApproximatePermutationTest(q0, q, mean, 1000)
ApproximateTwoSampleKSTest(q0, q)
KSampleADTest(q0,q; modified = true, nsim = 0)

# surrogate vs surrogate?
p,p0 = spectral_stationarity(id,fl,Δt,mov_filter)
plot!([q0,p0])
ApproximateTwoSampleKSTest(q0, p0)

tr, tfhmr = movfilter(ts,tfhm,"pre",Δt)
smeanr = mean(tfhmr,dims=2)[:,1]

using Gnuplot
@gp "
unset key
set lmargin 10;
set rmargin 5;
set multiplot layout 2,1;" :-
@gp :- 1 "
set log zcb;
set pale @RAINBOW;
set yrange[0:100];
set xrange[0:900];
set cbrange[1e-17:*];
set xlabel 'time [s]';
set ylabel 'freq. [Hz]'; " :-
@gp :- tfhmr "origin=(5,0) dx=10 dy=0.1 w image pixels" :-
@gp :- 2 "
set ylabel 'log-spectral dist.';
set yrange[0:*];
set size 0.8,0.5" :-
@gp :- qr "u (-5+\$0*10):1 w p ps 0.5"
@gp :- qr0 "u (-5+\$0*10):1 w p ps 0.5 pt 5"
