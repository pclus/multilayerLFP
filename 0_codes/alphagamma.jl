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
ns = [136 132 66 136]
α1 = zeros(ns[1],2); γ1 = zeros(ns[1],2)
α2 = zeros(ns[2],2); γ2 = zeros(ns[2],2)
α3 = zeros(ns[3],2); γ3 = zeros(ns[3],2)
α4 = zeros(ns[4],2); γ4 = zeros(ns[4],2)

for (i,condition) in enumerate(["pre","post"])

    lfp = readdlm(path*"psd_mean_tfhm_cortex_"*condition*".dat")
    blfp = readdlm(path*"psd_mean_tfhm_bipolar_"*condition*".dat")
    csd = readdlm(path*"psd_mean_tfhm_csd_"*condition*".dat")
    kcsd = readdlm(path*"psd_mean_tfhm_kcsd_"*condition*".dat")
    fr = 0.0:0.1:199.9


    α1[:,i],γ1[:,i] = relative_power(lfp,fr)
    α2[:,i],γ2[:,i] = relative_power(blfp,fr)
    α3[:,i],γ3[:,i] = relative_power(csd,fr)
    α4[:,i],γ4[:,i] = relative_power(kcsd,fr)

    global dp1 = depth("lfp",size(lfp)[2])
    global dp2 = depth("blfp",size(blfp)[2])
    global dp3 = depth("csd",size(csd)[2])
    global dp4 = depth("kcsd",size(kcsd)[2])

    # plot((α1./γ1)[1:4:end],dp1[1:4:end],lc=1)
    # plot!((α2./γ2)[1:4:end],dp2[1:4:end],lc=2)
    # plot!((α3./γ3)[2:2:end],dp3[2:2:end],lc=3)
    # p[5,i] = plot!(α4./γ4,dp4,lc=4)

    plot((α1./γ1)[1:4:end],dp1[1:4:end],lc=1)
    plot!((α2./γ2)[1:4:end],dp2[1:4:end],lc=2)
    plot!((α3./γ3)[2:2:end],dp3[2:2:end],lc=3)
    p[5,i] = plot!(α4./γ4,dp4,lc=4)

    plot((α1)[1:4:end]/maximum(α1[1:4:end]),dp1[1:4:end],lc=1)
    p[1,i] = plot!((γ1)[1:4:end]/maximum(γ1[1:4:end]),dp1[1:4:end],lc=2)

    plot((α2)[1:4:end]/maximum(α2[1:4:end]),dp2[1:4:end],lc=1)
    p[2,i] = plot!((γ2)[1:4:end]/maximum(γ2[1:4:end]),dp2[1:4:end],lc=2)

    plot((α3)[2:2:end]/maximum(α3[2:2:end]),dp3[2:2:end],lc=1)
    p[3,i] = plot!((γ3)[2:2:end]/maximum(γ3[2:2:end]),dp3[2:2:end],lc=2)

    plot(α4/maximum(α4),dp4,lc=1)
    p[4,i] = plot!(γ4/maximum(γ4),dp4,lc=2)
end

plot(p[5,1],p[5,2],layout=(1,2),xlim=(0,12.0),lw=2,
labels=["LFP" "bLFP" "CSD" "kCSD"],
xaxis = "α/γ",
yaxis = "depth [μm]",
legendfontsize=12,foreground_color_legend = nothing)


plot(p[1,1],p[1,2],layout=(1,2),labels=["α" "γ"],lw=2,
xlabel = "Relative power",
ylabel = "Depth [μm]",title=["pre" "post"],
legendfontsize=12,foreground_color_legend = nothing)

plot(p[4,1],p[4,2],layout=(1,2),labels=["α" "γ"],lw=2,
xlabel = "Relative power",
ylabel = "Depth [μm]",title=["pre" "post"],
legendfontsize=12,foreground_color_legend = nothing)

plot(α1[1:4:end,2]./α1[1:4:end,1],dp1[1:4:end],palette=:Set1_4)
plot!(α2[1:4:end,2]./α2[1:4:end,1],dp2[1:4:end])
plot!(α3[1:2:end,2]./α3[1:2:end,1],dp3[1:2:end])
np1=plot!(α4[:,2]./α4[:,1],dp4)

plot(γ1[1:4:end,2]./γ1[1:4:end,1],dp1[1:4:end],palette=:Set1_4)
plot!(γ2[1:4:end,2]./γ2[1:4:end,1],dp2[1:4:end])
plot!(γ3[1:2:end,2]./γ3[1:2:end,1],dp3[1:2:end])
np2=plot!(γ4[:,2]./γ4[:,1],dp4)

plot(np1,np2,layout=(1,2),xlim=(0.0,2.0),lw=2,
labels=["LFP" "bLFP" "CSD" "kCSD"],
xlabel = ["αₚᵣₑ/αₚₒₛₜ" "γₚᵣₑ/γₚₒₛₜ"],
ylabel = "depth [μm]",
legendfontsize=12,foreground_color_legend = nothing)




plot(4:0.1:22,lfp[40:220,70],
fillrange=1e-18.*ones(181),lc=1,fillalpha=0.25)
plot!(32:0.1:48,lfp[320:480,70],
fillrange=1e-18.*ones(161),lc=2,fillalpha=0.25)
plot!(0:0.1:199.9,lfp[:,70],
yaxis=:log,lw=2,
ylim=(0.5e-17,1e-14),xlim=(0,1e2),
xlabel="Freq. [Hz]",
ylabel="Power [mV²]",lc=3,key=:false)