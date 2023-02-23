#!/usr/bin/env -S julia --threads 16

push!(LOAD_PATH, "/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
using NeuropixelAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics;
# plotlyjs()

# Comparison between raw and filtered data----------------------


# Load time series for a specific channel from second 200 to 300:
t0 = 100;
tf = 200.0;
id = 196;
fl = "pre"
t, chdat = read_channel(id, t0, tf, fl);

# jd=Int((id-2)/2) # the formula changes if id is even
# t,chcsd=read_channel(jd,t0,tf,"csd_pre")

# Bandpass filter
bpfilter = digitalfilter(Bandpass(1.0, 300.0; fs=2500), Butterworth(3));
fil_chdat = filtfilt(bpfilter, chdat)


plot(t, chdat)
plot!(t, fil_chdat)


# Compute PSD using Multitaper PSD ----------------------------
fl = "pre"
t, chdat = read_channel(150, t0, tf, fl);
bpfilter = digitalfilter(Bandpass(1.0, 300.0; fs=2500), Butterworth(3));
fil_chdat = filtfilt(bpfilter, chdat)
rate = 2500.0;
dt = 1.0 / rate;

NW = 1.0 * length(chdat) * dt / (2.0);  # real bandwith is ω*dt/2.0, and NW = N*ω*dt/2.0
K = 10;    # number of tappers (should be similar slightly less than 2*NW)
S = multispec(chdat, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);
fS = multispec(fil_chdat, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);


p1 = plot(S.f, S.S, xlim=(0, 200), ylim=(1e-18, 1e-15), lw=1.0, yaxis=:log)
plot!(fS.f, fS.S, xlim=(0, 200), ylim=(1e-18, 1e-15), lw=1.0, yaxis=:log)
# Ftest pvalues
# plot(S.f,S.Fpval,ylim=(0.0,1e-2),lt=:scatter,xlim=(0,100))

t, kcsd = read_channel(id, t0, tf, "kCSD_electrodes_pre");
t, lcsd = read_channel(channel_idx(id), t0, tf, "csd_pre");


plot(t, kcsd)
plot!(t, lcsd)

plot(kcsd[1:10000], lcsd[1:10000], lt=:scatter, ms=2, msw=0, lw=0)
cor(kcsd[1:10000], lcsd[1:10000])

S1 = multispec(kcsd, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);
S2 = multispec(lcsd, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);

plot(S1.f, S1.S, xlim=(0, 200), ylim=(1e-18, 1e-15), lw=1.0, yaxis=:log)
plot!(S2.f, S2.S, xlim=(0, 200), ylim=(1e-18, 1e-15), lw=1.0, yaxis=:log)
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------

# -------------------------------------------------------------
# Create a heatmap of the pre data from 200s to 300s using Multitaper:
w, lfp = heatmapMT(200.0, 300.0, "filtered_pre", 1:384);
w, csd = heatmapMT(200.0, 300.0, "csd_pre", 1:190);
w, bip = heatmapMT(200.0, 300.0, "bipolar_pre", 1:376);
w, kcsd = heatmapMT(200.0, 300.0, "kCSD_electrodes_pre", 1:384);
w, kcsd = heatmapMT(200.0, 300.0, "kCSD_centers_pre", 1:384);

# gr()
# heatmap(w,1:384,kcsd')

writedlm("../4_outputs/lfp.dat", lfp', " ");
writedlm("../4_outputs/csd.dat", csd', " ");
writedlm("../4_outputs/bip.dat", bip', " ");
writedlm("../4_outputs/kcsd_cent.dat", kcsd', " ");
# -------------------------------------------------------------

# -------------------------------------------------------------
# Time-frequency analysis for a single channel:
fl = "csd_pre"
id = 140
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
n=[ 384,384,190,376]
for cond in ("pre","post")
    for (i,data) in enumerate(("kCSD_centers_","filtered_","csd_","bipolar_"))
        fl=data*cond
        m, s = heatmap_segments(fl,n[i])
        writedlm("../4_outputs/psd_mean_tfhm_" * fl * ".dat", m', ' ');
        writedlm("../4_outputs/psd_std_tfhm_" * fl * ".dat", s', ' ');
        # print(i,data,n[i],"\n")
    end
end

# heatmap(0.1:0.1:200,1:384,log.(psd_mean_tfhm))
# -------------------------------------------------------------


# -------------------------------------------------------------
# Check for differences between pre and post using t-test of the
# different time windows
using Statistics, HypothesisTests

n = 384;
l = 2000;

psd_mean_tfhm_pre = zeros(n, l);
psd_std_tfhm_pre = zeros(n, l);
psd_mean_tfhm_post = zeros(n, l);
psd_std_tfhm_post = zeros(n, l);
pvals_tfhm = zeros(n, l);
pvals_uneq_tfhm = zeros(n, l);

state = Threads.Atomic{Int}(0);

Threads.@threads for id in 1:n

    # Prompt state
    Threads.atomic_add!(state, 1)
    print("--> ", state[], " out of ", n, "\n")
    flush(stdout)

    # Pre
    t, f, tfhm = timefreq(id, "pre")
    idx, f_pre = movfilter(t, tfhm, "pre")
    psd_mean_tfhm_pre[id, :] = mean(f_pre, dims=2)
    psd_std_tfhm_pre[id, :] = std(f_pre, dims=2)

    # Post
    t, f, tfhm = timefreq(id, "post")
    idx, f_post = movfilter(t, tfhm, "post")
    psd_mean_tfhm_post[id, :] = mean(f_post, dims=2)
    psd_std_tfhm_post[id, :] = std(f_post, dims=2)

    # T-test. 
    for j in 1:l
        pvals_uneq_tfhm[id, j] = pvalue(UnequalVarianceTTest(f_pre[j, :], f_post[j, :]))
        pvals_tfhm[id, j] = pvalue(EqualVarianceTTest(f_pre[j, :], f_post[j, :]))# can be computed from the means
    end
end

writedlm("../4_outputs/psd_mean_tfhm_post.dat", psd_mean_tfhm_pre', ' ');
writedlm("../4_outputs/psd_std_tfhm_post.dat", psd_std_tfhm_pre', ' ');
writedlm("../4_outputs/psd_mean_tfhm_pre.dat", psd_mean_tfhm_post', ' ');
writedlm("../4_outputs/psd_std_tfhm_pre.dat", psd_std_tfhm_post', ' ');
writedlm("../4_outputs/psd_pvals.dat", pvals_tfhm', ' ');
writedlm("../4_outputs/psd_pvals_uneq.dat", pvals_uneq_tfhm', ' ');
writedlm("../4_outputs/psd_mean_tfhm_difference.dat", (psd_mean_tfhm_post - psd_mean_tfhm_pre)', ' ');

# heatmap(0.1:0.1:200,1:n,log.(psd_mean_tfhm))

#plot(f,psd_mean_tfhm_pre[id,:],ribbon=psd_std_tfhm_pre[id,:],c=1,fillalpha=0.25)
#plot!(f,psd_mean_tfhm_post[id,:],ribbon=psd_std_tfhm_post[id,:],c=2,fillalpha=0.25)
#plot!(f,(pvals_tfhm[id,:].<=1e-5).*1e-17,ylim=(1e-18,1.2e-15),lt=:scatter)

# t-test comparison
# a=randn(100);
# b=10*randn(100).+2.2;
# p=pvalue(UnequalVarianceTTest(a,b))
# p=pvalue(EqualVarianceTTest(a,b))
