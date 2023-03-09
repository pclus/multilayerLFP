#!/usr/bin/env -S julia --threads 16

push!(LOAD_PATH, "/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
using NeuropixelAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics;
# plotlyjs()

n0=226;
nf=361;
# process_data(n0,nf);


# Load time series for a specific channel from second 200 to 300:
t0 = 100.0;
tf = 200.0;
id = 100;
fl = "cortex_pre"
t, chdat = read_channel(id, t0, tf, fl);

plot(t, chdat)

using StatsBase, HypothesisTests
ν=autocor(chdat,0:10*2500)
plot(dt*(0:25000),ν,xlim=(0,1.0))
# ADFTest(chdat, :none, 500) 


# Compute PSD using Multitaper PSD ----------------------------
rate = 2500.0;
dt = 1.0 / rate;
NW = 1.0 * length(chdat) * dt / (2.0);  # real bandwith is ω*dt/2.0, and NW = N*ω*dt/2.0
K = 10;    # number of tappers (should be similar slightly less than 2*NW)
S = multispec(chdat, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);


p1 = plot(S.f, S.S, xlim=(0, 200), ylim=(1e-18, 1e-14), lw=1.0, yaxis=:log)

# Numerical check that the variance of the time signal is the overall power.
# Does not hold for the MT because of the tappering
using FFTW
Fr = rfft(chdat); 
wr = rfftfreq(length(t), 1.0/dt); 

area = 2.0*sum(abs2.(Fr))*1e-2*dt/length(chdat)
areaMT=sum(S.S)*1e-2 #df=1e-2
var(chdat)

# Compare different Power Spectral Densities

f = S.f[1:20000]
psd1 = S.S[1:20000];

using HypothesisTests

S1 = S.S[1:20000]
S2 = abs2.(Fr)[1:20000]

function KS_psd(S1,S2,f)
    df = f[2]-f[1]
    psd1 = S1/ (sum(S1)*df); # normalization
    cdf1 = [ sum(psd1[1:i])*df for i in 1:length(psd1)]

    psd2 = S2/ (sum(S2)*df); # normalization
    cdf2 = [ sum(psd2[1:i])*df for i in 1:length(psd2)]

    sampl=1000   # very sensible to this parameter
    ζ1=rand(sampl)
    fsampl1 = [f[findfirst(ζ1[i].<=cdf1)] for i in 1:sampl]

    ζ2=rand(sampl)
    fsampl2 = [f[findfirst(ζ2[i].<=cdf2)] for i in 1:sampl]

    # notice that we are using the `dev` version of the package
    # (revert with `] free HypothesisTest`)
    tKS = ApproximateTwoSampleKSTest(fsampl1,fsampl2)
    tAD = KSampleADTest(fsampl1,fsampl2)
    return tKS,tAD,cdf1,cdf2
end

t1,t2,a,b = KS_psd(S.S[1:2000],S.S[1:2000],f)



t, f, tfhm = timefreq(id, fl);
m1 = [pvalue(KS_psd(tfhm[:,i],tfhm[:,j],f)[1]) for i in 1:90,j in 1:90]
m2 = [pvalue(KS_psd(tfhm[:,i],tfhm[:,j],f)[2]) for i in 1:90,j in 1:90]
heatmap(m1,clim=(0,1e-1)) # changes A LOT depending on `sampl`
heatmap(m2,clim=(0,1e-1))
# maybe use χ2 or maybe directly test for stationarity

function LogSpectralDistance(s1,s2,f)
    y = @.  10*log10(s1/s2)^2
    integrate(f,y)^0.5 # maybe divided by f[end]^0.5
end

m3 = [LogSpectralDistance(tfhm[:,i],tfhm[:,j],f) for i in 1:90,j in 1:90]
heatmap(m3,clim=(1e1,1.5e1)) #amazing

function RandomSpectra(smean,σ,nsampl)
    n = length(smean);
    rs = zeros(n,nsampl)
    for i in 1:nsampl
        rs[:,i] = smean + σ.*randn(n);
    end
    rs[findall(rs.<0)].=1e-22; # avoid negative powers, this is problematic though
    return rs;
end

f_idx, f_tfhm = movfilter(t, tfhm, "pre");
smean = mean(f_tfhm,dims=2)
σ = std(f_tfhm,dims=2);
# σ = findmax(σ)[1]
rs = RandomSpectra(smean,σ,90);
mr = [LogSpectralDistance(rs[:,i],rs[:,j],f) for i in 1:90,j in 1:90]
heatmap(mr) # nice result, but still, maybe surrogates


# -------------------------------------------------------------
# -------------------------------------------------------------
path="/home/pclusella/Documents/Data/UPO-tACs/7_results/cortex_heatmaps/"
lfp = readdlm(path*"psd_mean_tfhm_cortex_pre.dat")
blfp = readdlm(path*"psd_mean_tfhm_bipolar_pre.dat")
csd = readdlm(path*"psd_mean_tfhm_csd_pre.dat")
kcsd = readdlm(path*"psd_mean_tfhm_kcsd_pre.dat")
freqs = 0.0:0.1:199.9

function relative_power(fhm,freqs)
    df = freqs[2]-freqs[1]
    n = size(fhm)[2]
    αband = findall(@. f>4.0 && f<22)
    γband = findall(@. f>35.0 && f<60)
    α = zeros(n)
    γ = zeros(n)
    for i in 1:n
        α[i] = mean(fhm[αband,i])*df # take the sum instead of the mean to get the area
        γ[i] = mean(fhm[γband,i])*df
    end

    return α,γ
end

α1,γ1 = relative_power(lfp,freqs)
α2,γ2 = relative_power(blfp,freqs)
α3,γ3 = relative_power(csd,freqs)
α4,γ4 = relative_power(kcsd,freqs)

dp1 = depth("lfp",size(lfp)[2])
dp2 = depth("blfp",size(blfp)[2])
dp3 = depth("csd",size(csd)[2])
dp4 = depth("kcsd",size(kcsd)[2])

plot(dp1[1:4:end],(α1./γ1)[1:4:end],lc=1)
plot!(dp2[1:4:end],(α2./γ2)[1:4:end],lc=2)
plot!(dp3[2:2:end],(α3./γ3)[2:2:end],lc=3)
plot!(dp4,α4./γ4,lc=4)


plot((α1)[1:4:end]/maximum(α1[1:4:end]),dp1[1:4:end],lc=1)
plot!((γ1)[1:4:end]/maximum(γ1[1:4:end]),dp1[1:4:end],lc=2)

plot((α2)[1:4:end]/maximum(α2[1:4:end]),dp2[1:4:end],lc=1)
plot!((γ2)[1:4:end]/maximum(γ2[1:4:end]),dp2[1:4:end],lc=2)

plot((α3)[2:2:end]/maximum(α3[2:2:end]),dp3[2:2:end],lc=1)
plot!((γ3)[2:2:end]/maximum(γ3[2:2:end]),dp3[2:2:end],lc=2)

plot(α4/maximum(α4),dp4,lc=1)
plot!(γ4/maximum(γ4),dp4,lc=2)


# -------------------------------------------------------------
# -------------------------------------------------------------
# Create a heatmap of the pre data from 200s to 300s using Multitaper:
n0=226;
nf=361;
w, lfp_cortex = heatmapMT(200.0, 300.0, "cortex_pre", 0:(nf-n0));
w, lfp = heatmapMT(200.0, 300.0, "filtered_pre", 0:383);
w, csd = heatmapMT(200.0, 300.0, "csd_pre", 1:190);
w, bip = heatmapMT(200.0, 300.0, "bipolar_pre", 1:376);
w, kcsd = heatmapMT(200.0, 300.0, "kCSD_electrodes_pre", 1:384);
w, kcsd = heatmapMT(200.0, 300.0, "kCSD_centers_pre", 1:384);

w, lfp_post = heatmapMT(200.0, 300.0, "filtered_post", 1:384);
# gr()
# heatmap(w,1:384,kcsd')

writedlm("../4_outputs/cortex_lfp.dat", lfp_cortex', " ");
writedlm("../4_outputs/lfp.dat", lfp', " ");
writedlm("../4_outputs/csd.dat", csd', " ");
writedlm("../4_outputs/bip.dat", bip', " ");
writedlm("../4_outputs/kcsd_cent.dat", kcsd', " ");
# -------------------------------------------------------------

# -------------------------------------------------------------
# Time-frequency analysis for a single channel:
fl = "cortex_pre"
id = 100
t, f, tfhm = timefreq(id, fl);
writedlm("../4_outputs/ch" * string(id) * "_tfhm_" * fl * ".dat", tfhm, " ");

# Remove the segments containing movement:
f_idx, f_tfhm = movfilter(t, tfhm, "pre");

# Compute statistics for this channel:
using Statistics, HypothesisTests
mean_psd = mean(f_tfhm, dims=2);
std_psd = std(f_tfhm, dims=2);
plot(f, mean_psd, ribbon=std_psd, c=1, fillalpha=0.25, yaxis=:log, yrange=(1e-23, 0.2e-21))

# -------------------------------------------------------------

# -------------------------------------------------------------
# Segmentation analysis for all channels and create the heatmaps for the averages
# ns=[ 384,384,190,376] # full
n=size(n0:nf)[1]
ns=[ n, n, n/2-2,n-4]    # cortex
ns=Int.(ns)
# for cond in ("pre",)
for cond in ("pre","post")
    for (i,data) in enumerate(("kcsd_","cortex_","csd_","bipolar_"))
    # for (i,data) in enumerate(("kcsd_",))
        fl=data*cond
        m, s = heatmap_segments(fl,ns[i])
        writedlm("../4_outputs/psd_mean_tfhm_" * fl * ".dat", m', ' ');
        writedlm("../4_outputs/psd_std_tfhm_" * fl * ".dat", s', ' ');
        # print(i,data,n[i],"\n")
    end
end

# heatmap(0.1:0.1:200,1:384,log.(psd_mean_tfhm))
# -------------------------------------------------------------
# n0=226;
# nf=361;
# n=size(n0:nf)[1]
# ns=[ n, n, n/2-2,n-4]    # cortex
# ns=Int.(ns)
# fl="kcsd_pre"
# m, s = heatmap_segments(fl,ns[2])
# writedlm("../4_outputs/psd_mean_tfhm_" * fl * ".dat", m', ' ');
# writedlm("../4_outputs/psd_std_tfhm_" * fl * ".dat", s', ' ');

# -------------------------------------------------------------
# Check for differences between pre and post using t-test of the
# different time windows
using Statistics, HypothesisTests

l = 2000;
n=size(n0:nf)[1]
ns=[ n, n, n/2-2,n-4]
ns=Int.(ns)

for (i,data) in enumerate(("kcsd_","cortex_","csd_","bipolar_"))
    n0=ns[i];

    psd_mean_tfhm_pre = zeros(n0, l);
    psd_std_tfhm_pre = zeros(n0, l);
    psd_mean_tfhm_post = zeros(n0, l);
    psd_std_tfhm_post = zeros(n0, l);
    pvals_tfhm = zeros(n0, l);
    pvals_uneq_tfhm = zeros(n0, l);

    state = Threads.Atomic{Int}(0);

    Threads.@threads for id in 0:n0-1

        # Prompt state
        Threads.atomic_add!(state, 1)
        print("--> ", state[], " out of ", n0, "\n")
        flush(stdout)

        # Pre
        t, f, tfhm = timefreq(id, data*"pre")
        idx, f_pre = movfilter(t, tfhm, "pre")
        psd_mean_tfhm_pre[id+1, :] = mean(f_pre, dims=2)
        psd_std_tfhm_pre[id+1, :] = std(f_pre, dims=2)

        # Post
        t, f, tfhm = timefreq(id, data*"post")
        idx, f_post = movfilter(t, tfhm, "post")
        psd_mean_tfhm_post[id+1, :] = mean(f_post, dims=2)
        psd_std_tfhm_post[id+1, :] = std(f_post, dims=2)

        # T-test. 
        for j in 1:l
            pvals_uneq_tfhm[id+1, j] = pvalue(UnequalVarianceTTest(f_pre[j, :], f_post[j, :]))
            pvals_tfhm[id+1, j] = pvalue(EqualVarianceTTest(f_pre[j, :], f_post[j, :]))# can be computed from the means
        end
    end

    writedlm("../4_outputs/psd_mean_tfhm_"*data*"pre.dat", psd_mean_tfhm_pre', ' ');
    writedlm("../4_outputs/psd_std_tfhm_"*data*"pre.dat", psd_std_tfhm_pre', ' ');
    writedlm("../4_outputs/psd_mean_tfhm_"*data*"post.dat", psd_mean_tfhm_post', ' ');
    writedlm("../4_outputs/psd_std_tfhm_"*data*"post.dat", psd_std_tfhm_post', ' ');
    writedlm("../4_outputs/psd_"*data*"pvals_eq.dat", pvals_tfhm', ' ');
    writedlm("../4_outputs/psd_"*data*"pvals.dat", pvals_uneq_tfhm', ' ');
    writedlm("../4_outputs/psd_"*data*"difference.dat", (psd_mean_tfhm_post - psd_mean_tfhm_pre)', ' ');
end



# heatmap(0.1:0.1:200,1:n,log.(psd_mean_tfhm))

#plot(f,psd_mean_tfhm_pre[id,:],ribbon=psd_std_tfhm_pre[id,:],c=1,fillalpha=0.25)
#plot!(f,psd_mean_tfhm_post[id,:],ribbon=psd_std_tfhm_post[id,:],c=2,fillalpha=0.25)
#plot!(f,(pvals_tfhm[id,:].<=1e-5).*1e-17,ylim=(1e-18,1.2e-15),lt=:scatter)

# t-test comparison
# a=randn(100);
# b=10*randn(100).+2.2;
# p=pvalue(UnequalVarianceTTest(a,b))
# p=pvalue(EqualVarianceTTest(a,b))
