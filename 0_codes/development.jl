push!(LOAD_PATH, pwd())
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
using StatsBase
using ColorSchemes,LaTeXStrings,Measures

id = 50
fl = "bipolar_pre"
dt = 0.0004
m = 9000000
t, y = read_channel(id,fl ; t0=dt, tf=m*dt, m=m)
t, z = read_channel(100,fl ; t0=dt, tf=m*dt, m=m,path=)



@time s = mt_spectrogram( y , 2500*10 , 0 ; fs=rate , nw=NW , ntapers=K )

# same results
plot(s.freq[1:2000],s.power[1:2000,200]*0.5)
plot!(f,tfhm[:,200])

lags = -2500*5:2500*5;
cx = crosscor(z[1:2500*10], y[1:2500*10],lags)

plot(lags,cx)