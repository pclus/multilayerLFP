#!/usr/bin/env -S julia --threads 16

push!(LOAD_PATH, "/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics;
# plotlyjs()

# load_precomputed();

id=100; fl="pre"
t, f, tfhm = timefreq_complete(id, fl);
m = [logspectral_dist(tfhm[:,i],tfhm[:,j],f) for i in 1:90,j in 1:90]
heatmap(m,clim=(0.2,0.45)) #amazing


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

function timefreq_complete(id, fl)

    rate = 2500.0
    dt = 1.0 / rate
    n = 384
    Δt = 10.0      # segment duration (10 seconds)
    ns = Int(900.0 / Δt)  # number of segments
    m = Int(2250000 / ns) # segments of 10 seconds
    mh = m/2;       # length of the MT spectra...
    l = Int(mh)       # but...up to l is enough to get the relevant freqs
    NW = 1.0 * m * dt / (2.0)
    K = 8
    t, chdat = read_channel(id, 4e-4, 900, fl)    # time series starts at 4e-4


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

    freqs = collect(S.f[1:l])
    times = collect(5:10:900)

    times, freqs, tfhm
end
