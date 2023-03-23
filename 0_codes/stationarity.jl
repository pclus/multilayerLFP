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

    return q,q0,tfhm,tfhm0,f;
end

id=100; fl="cortex_pre"; Δt=10.0; mov_filter="pre";
q,q0,tfhm,tfhm0,f = spectral_stationarity(id,fl,Δt,mov_filter)
qr,qr0,tfhmr,tfhmr0 = spectral_stationarity(id,fl,Δt,"")

# tests
p_perm = pvalue(ApproximatePermutationTest(q0, q, mean, 1000))
p_KS = pvalue(ApproximateTwoSampleKSTest(q0, q))
p = pvalue(KSampleADTest(q0,q; modified = true, nsim = 0))

pr_perm = pvalue(ApproximatePermutationTest(qr0, qr, mean, 1000))
pr_KS = pvalue(ApproximateTwoSampleKSTest(qr0, qr))
pr_AD = pvalue(KSampleADTest(qr0,qr; modified = true, nsim = 0))

# plots
using LaTeXStrings
dotplot([q q0], msw = 0, ms=2.0)
dotplot!([qr qr0],
    xticks=([ 1 2 3 4],["no-mov", "no-mov surr.", "mov.", "mov. surr."]),
    xtickfontsize=11,legend=:false,
    msw = 0, ms=2.0)
plot!([1:4], [mean(q), mean(q0), mean(qr), mean(qr0)],
    yerr=[std(q), std(q0), std(qr), std(qr0)],
    lt=:scatter,mc=:black,lw=3,msw=3,
    ylabel="distance to mean "*L"d(S_j,\langle S\rangle)")
annotate!(1.5,0.25,text("pₖₛ="*@sprintf("%.3e",p_KS),10))
annotate!(3.5,0.25,text("pₖₛ="*@sprintf("%.3e",pr_KS),10))


# surrogate vs surrogate?
# q ,q1 = spectral_stationarity(id,fl,Δt,mov_filter)
# plot([q1,q0])
# ApproximateTwoSampleKSTest(q1, q0)


smean = mean(tfhm,dims=2)[:,1]
smean0 = mean(tfhm0,dims=2)[:,1]
# smeanr = mean(tfhmr,dims=2)[:,1]
# smeanr0 = mean(tfhmr0,dims=2)[:,1]

sstd = std(tfhm,dims=2)[:,1]
sstd0 = std(tfhm0,dims=2)[:,1]
# sstdr = std(tfhmr,dims=2)[:,1]
# sstdr0 = std(tfhmr0,dims=2)[:,1]

plot(f,smean,ribbon=sstd,c=1,fillalpha=0.25,xlim=(0,100))
plot!(f,smean0,ribbon=sstd0,c=2,fillalpha=0.25,xlim=(0,100),
    xlabel = "Freq [Hz]", ylabel="Power [mV²]", labels=[ "no-mov." "no-mov. surr."])

# plot(f,smeanr,ribbon=sstd,c=1,fillalpha=0.25,xlim=(0,100))
# plot!(f,smeanr0,ribbon=sstd0,c=2,fillalpha=0.25,xlim=(0,100))


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
set key box center top
set ylabel 'log-spectral dist.';
set yrange[0:*];
set size 0.88,0.5" :-
@gp :- qr "u (-5+\$0*10):1 w p ps 0.5 t 'mov.'"
@gp :- qr0 "u (-5+\$0*10):1 w p ps 0.5 pt 5 t 'mov. surr.'"


# iterate overall all channels
fl="cortex_pre"; Δt=10.0; mov_filter="pre";
n=136
pvals=zeros(n)
@Threads.threads for id in 0:n-1
    q,q0,tfhm,tfhm0,f = spectral_stationarity(id,fl,Δt,mov_filter)
    pvals[id+1] = pvalue(ApproximateTwoSampleKSTest(q0, q))
end

