#!/usr/bin/env -S julia --threads 16

cd("0_codes/")
push!(LOAD_PATH, pwd())
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
using FFTW



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

n=136
lseg = [1.0, 5.0, 10.0, 20.0, 30.0, 50.0, 60.0, 90.0]
pvals = zeros(length(lseg),length(0:5:n-1))
@Threads.threads for (j,id) in collect(enumerate(0:5:n-1))
    for (i,Δt) in enumerate(lseg)
        q,q0,tfhm,tfhm0,f = spectral_stationarity(id,fl,Δt,mov_filter)
        pvals[i,j] = pvalue(ApproximateTwoSampleKSTest(q0, q))
    end
end

# plots
using LaTeXStrings,StatsPlots,Printf
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
savefig("~/Desktop/Fig2.png")

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

plot(f,smean,ribbon=sstd,c=1,fillalpha=0.25,xlim=(0,100), labels="no-mov.")
plot!(f,smean0,ribbon=sstd0,c=2,fillalpha=0.25,xlim=(0,100),
    xlabel = "Freq [Hz]", ylabel="Power [mV²]", labels="no-mov. surr.")
savefig("~/Desktop/Fig3.png")

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
savefig("~/Desktop/Fig1.png")

# iterate overall all channels
fl="cortex_pre"; Δt=10.0; mov_filter="pre";
n=136
pvals=zeros(n)
@Threads.threads for id in 0:n-1
    q,q0,tfhm,tfhm0,f = spectral_stationarity(id,fl,Δt,mov_filter)
    pvals[id+1] = pvalue(ApproximateTwoSampleKSTest(q0, q))
end

plot(pvals,xlabel="Channel", ylabel="KS p-value",yaxis=:log, lt=:scatter,legend=:false)
savefig("~/Desktop/Fig4.png")
