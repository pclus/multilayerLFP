#!/usr/bin/env -S julia --threads 16

push!(LOAD_PATH, "/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
# plotlyjs()

id=220
t,chdat = read_channel(id,200.0,201.0,"filtered_pre")
rate = 2500.0;
dt = 1.0 / rate;
NW = 1.0 * length(chdat) * dt / (2.0);  
K = 10;   
S = multispec(chdat, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);
p1 = plot(S.f, S.S, xlim=(0, 205), ylim=(1e-18, 1e-14), lw=2.0, yaxis=:log,
xlabel="Freq. [Hz]",ylabel="Power [mV²]")

#-------------------------------------------
rate = 2500.0
dt = 1.0 / rate
n = 384
Δt = 1.0      # segment duration (10 seconds)
ns = Int(900.0 / Δt)  # number of segments
m = Int(2250000 / ns) # segments of 10 seconds
# mh = m/2;       # length of the MT spectra...
l = 2000       # but...up to l is enough to get the relevant freqs
NW = 1.0 * m * dt / (2.0)
K = 8
t, chdat = read_channel(id, 4e-4, 900, "filtered_pre")    # time series starts at 4e-4

tfhm = zeros(l, ns)   # time-freq heatmap
local S    # so that S exists outside the loop
for s in 1:ns
    segdat = chdat[((s-1)*m+1):s*m]
    S = multispec(segdat, dt=dt, NW=NW, K=K)
    tfhm[:, s] = S.S[1:l]
    if s % 5 == 0
        print(s, "\r")
        flush(stdout)
    end
end
# -------------------------------------------------------

# spectra = NeuropixelAnalysis.load_precomputed() ;
# spectra_std = NeuropixelAnalysis.load_precomputed_std()
ts, f, tfhm = timefreq(id, data*"pre")


# Plain permutation test [to be incorporated in the module]-------------
id = 100
data = "cortex_"
ts, f, tfhm = timefreq(id, data*"pre")
idx, f_pre = movfilter(ts, tfhm, "pre")

t,ch = read_channel(id,1/2500.0,900,data*"pre")

heatmap(ts,f,log10.(tfhm),c=:rainbow1,clim=(-17,-15))

writedlm("../4_outputs/tfhm.dat",tfhm," ")
writedlm("../4_outputs/chdat.dat",[t ch]," ")


m = [logspectral_dist(tfhm[:,i],tfhm[:,j],f) for i in 1:90,j in 1:90]
mf = [logspectral_dist(f_pre[:,i],f_pre[:,j],f) for i in 1:82,j in 1:82]
using LaTeXStrings
heatmap(m,aspect_ratio=:equal,xlim=(0.5,90.5),ylim=(0.5,90.5),c=:afmhot,
xlabel="Time segment "*L"i",
ylabel="Time segment "*L"j",
colorbar_title="\nlog-spectral dist. "*L"d(S_i,S_j)")
heatmap(mf,aspect_ratio=:equal,xlim=(0.5,82.5),ylim=(0.5,82.5),c=:afmhot,
xlabel="Time segment "*L"i",
ylabel="Time segment "*L"j",
colorbar_title="\nlog-spectral dist. "*L"d(S_i,S_j)")

##
nspectra=copy(spectra).*0;
nspectra_std=copy(spectra_std).*0;
dists=
for k in 1:8
    for i in 1:size(spectra[k])[2]
        nspectra[k][:,i]=spectra[k][:,i]./mean(spectra[k][:,i])
        nspectra_std[k][:,i]=spectra_std[k][:,i]./mean(spectra[k][:,i])
    end
        
end
##
heatmap(((nspectra[4]-nspectra[8])'),clim=(-5,5),xlim=(0,1000),c=:bwr)

Sp=pvalue.(EqualVarianceTTest.(82, 67, nspectra[4], nspectra[8],
     nspectra_std[4].^2,nspectra_std[8].^2))

heatmap(log10.(p'),clim=(-4,0),xlim=(0,1000))

n=[ size(spectra[i])[2] for i in 1:8]
dists=[ zeros(n[i]) for i in 1:4]
for k=1:4
    for i in 1:size(spectra[k])[2]
        dists[k][i]=SpectralAnalysis.logspectral_dist(nspectra[k][:,i],nspectra[k+4][:,i],2.0)
    end
end
plot(dists[:])

id=100; fl="pre"
t, f, tfhm = timefreq_complete(id, fl);
m = [logspectral_dist(tfhm[:,i],tfhm[:,j],f) for i in 1:90,j in 1:90]


heatmap(m,clim=(0.0,0.45)) #amazing


f_idx, f_tfhm = movfilter(t, tfhm, "pre");
smean = mean(f_tfhm,dims=2)
σ = std(f_tfhm,dims=2);
rs = randomize_spectra(smean,σ,90);
mr = [logspectral_dist(rs[:,i],rs[:,j],f) for i in 1:90,j in 1:90]
heatmap(mr) # nice result, but still, maybe surrogates

m4 = [logspectral_dist(f_tfhm[:,i],f_tfhm[:,j],f) for i in 1:82,j in 1:82]
heatmap(m4) #amazing



# Surrogates
using FFTW
dt = 1/2500.0
sur_S = @. sqrt.(smean)*exp(im*rand()*2*π)
# sur_y = irfft(sur_S[1:1025],2048)
scale=(2.0*length(smean)/dt)^0.5
sur_y = irfft(sur_S,2*(12500-1))*scale; 
# plot(2.5*dt:dt:10.0,sur_y,lw=1.0,xlim=(0,0.5)) ##  dt/length(S.f)

sur1 = @. sur_y + 1e-7*randn()
# Fr = rfft(sur1); 
# wr = rfftfreq(length(sur_y), 1.0/dt); 
S = multispec(sur1[:,1], dt=dt, NW= 1.0 * length(sur1) * dt / (2.0), K=10);
plot([S.f,f,f],[S.S,tfhm[:,2],smean],yaxis=:log,lw=[2 1 0.5],xlim=(0.0,2e2))


LogSpectralDistance(mean_S[1:1025],abs2.(Fr),0.0:0.1:102.4)

plot([rs[:,1],tfhm[:,2],smean],yaxis=:log)

plot([abs2.(Fr),S.S,smean],yaxis=:log)


