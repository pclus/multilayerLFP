#!/usr/bin/env -S julia --threads 16

push!(LOAD_PATH, "/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
# plotlyjs()

# id=220
# t,chdat = read_channel(id,200.0,201.0,"filtered_pre")
# rate = 2500.0;
# dt = 1.0 / rate;
# NW = 1.0 * length(chdat) * dt / (2.0);  
# K = 8;   
# S = multispec(chdat, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);
# p1 = plot(S.f, S.S, xlim=(0, 205), ylim=(1e-18, 1e-14), lw=2.0, yaxis=:log,
# xlabel="Freq. [Hz]",ylabel="Power [mV²]")


id=100; fl="cortex_pre"
t, f, tfhm = timefreq(id, fl);
t, f, tfhmtemp = timefreq(y);
m = [logspectral_dist(tfhm[:,i],tfhm[:,j],f) for i in 1:90,j in 1:90]
heatmap(m,clim=(0.0,0.45))

smean = mean(tfhm,dims=2)[:,1]


id=100; fl="cortex_pre"
t, y = read_channel(id, 4e-4, 900, fl)
S = rfft(y); 
f = rfftfreq(length(t), 1.0/dt); 
S0 = @. abs.(S)*exp(im*rand()*2*π);
y0 = irfft(sur_S,2*length(S0)-2)




# -------------------------------------------------------------
# Surrogates --------------------------------------------------
# -------------------------------------------------------------
using FFTW
dt = 1/2500.0
l=2000

# Fill the average spectra with points interpolating
f_end=f[end]; df=f[2]-f[1];
freqs = 0.0:df:f_end
itp = Interpolations.scale(interpolate(smean, BSpline(Linear())), freqs )
l0 = Int((10*rate + 2)/2)
df0 = f_end/(ln-1)
freqs0 =  0.0:dfn:199.9
smean0 = itp(freqs0)

# create a surrogate with random phases
sur_S = @. sqrt.(smean0)*exp(im*rand()*2*π);
# scale=(2.0*length(smean0)/dt)^0.5
scale = (2.0*l0*rate*10*df0)^0.5
sur = irfft(sur_S,2*l0-2)*scale

# interpolate back to original sampling
itp0 = extrapolate(Interpolations.scale(interpolate(S.S,BSpline(Linear())),S.f),Line())
s0 = itp0(freqs) 
plot([smean,s0,tfhm[:,1]])

logspectral_dist(smean,s0,freqs)
logspectral_dist(smean,tfhm[:,1],freqs)

# This method won't help us. There are two issues:
# (1) We are modifying the spectra resolution "ad hoc".
# This might be only a technical point, but it looks "fishy".
# (2) [The important point] The surrogates are too close to the
# average spectra. Notice that if it was not for the interpolation 
# the smean and s0 would be completely identical
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------