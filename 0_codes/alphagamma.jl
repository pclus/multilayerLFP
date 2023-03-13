#!/usr/bin/env -S julia --threads 16

push!(LOAD_PATH, "/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
using NeuropixelAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics;
# plotlyjs()
gr()

path="/home/pclusella/Documents/Data/UPO-tACs/7_results/cortex_heatmaps/"
condition = "post"
p = Array{Plots.Plot{Plots.GRBackend},2}(undef,5,2)
for (i,condition) in enumerate(["pre","post"])

    lfp = readdlm(path*"psd_mean_tfhm_cortex_"*condition*".dat")
    blfp = readdlm(path*"psd_mean_tfhm_bipolar_"*condition*".dat")
    csd = readdlm(path*"psd_mean_tfhm_csd_"*condition*".dat")
    kcsd = readdlm(path*"psd_mean_tfhm_kcsd_"*condition*".dat")
    freqs = 0.0:0.1:199.9


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
    p[5,i] = plot!(dp4,α4./γ4,lc=4)

    plot((α1)[1:4:end]/maximum(α1[1:4:end]),dp1[1:4:end],lc=1)
    p[1,i] = plot!((γ1)[1:4:end]/maximum(γ1[1:4:end]),dp1[1:4:end],lc=2)

    plot((α2)[1:4:end]/maximum(α2[1:4:end]),dp2[1:4:end],lc=1)
    p[2,i] = plot!((γ2)[1:4:end]/maximum(γ2[1:4:end]),dp2[1:4:end],lc=2)

    plot((α3)[2:2:end]/maximum(α3[2:2:end]),dp3[2:2:end],lc=1)
    p[3,i] = plot!((γ3)[2:2:end]/maximum(γ3[2:2:end]),dp3[2:2:end],lc=2)

    plot(α4/maximum(α4),dp4,lc=1)
    p[4,i] = plot!(γ4/maximum(γ4),dp4,lc=2)
end

plot(p[5,1],p[5,2],layout=(1,2),ylim=(0,12.0))
plot(p[1,1],p[1,2],layout=(1,2),labels=["α" "γ"],
xlabel = "Relative power",
ylabel = "Depth [μm]",title=["pre" "post"])

