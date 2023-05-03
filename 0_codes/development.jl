push!(LOAD_PATH, pwd())
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
using StatsBase
using ColorSchemes,LaTeXStrings,Measures

id = 250
fl = "bipolar_pre"
dt = 0.0004
m = 9000000
t, y = read_channel(id,fl ; t0=dt, tf=m*dt, m=m)

n=2500*10
s = welch_pgram(y, 2500*10, 2500*5; onesided=eltype(y)<:Real, nfft=nextfastfft(n), fs=rate, window=hanning)

plotlyjs()
plot(s.freq*2500.0,s.power)

rate = 2500.0;
dt = 1.0 / rate;
NW = 1.0 * length(y) * dt / (2.0);  
K = 10;   
S = multispec(y, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);

plot!(S.f[1:100:end],S.S[1:100:end]*rate*pi*0.5)

using ColorSchemes
s = spectrogram(y, 2500, Int(2500/4); fs=rate, window=hanning)

plotlyjs()

@time s = mt_pgram(y[1:2500*100]; fs=rate, nw=NW, ntapers=K);
@time S = multispec(y[1:2500*100], dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);

plot(S.f,S.S./sum(S.S))
plot!(s.freq,s.power./sum(s.power))

s = mt_spectrogram(y,2500*10,0;fs=rate,nw=NW,ntapers=K)
#much faster!


Y = [y[1:2500*100] y[2500*100+1:2500*200]]
@time S = multispec(Y, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);
@time S = multispec(Y[:,1], dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);