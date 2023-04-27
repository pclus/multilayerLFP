#!/usr/bin/env -S julia --threads 16

cd("0_codes/")
push!(LOAD_PATH, pwd())
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
using StatsBase
using ColorSchemes,LaTeXStrings

tfhm_mean_pre=Dict()
tfhm_mean_post=Dict()
tfhm_std_pre=Dict()
tfhm_std_post=Dict()
tfhm_diff=Dict()
tfhm_pvals=Dict()

band_stats_pre=Dict()
band_stats_post=Dict()
band_pvals_alpha=Dict()
band_pvals_gamma=Dict()

Q_pre=Dict()
Q0_pre=Dict()
Q_post=Dict()
Q0_post=Dict()

for subj in [8,9,10,11,12,13,15,16,18,20]
    data = "bipolar_"
    foutname = "suj"*string(subj)*"/"
    namebase = "../4_outputs/pipeline/"*foutname*data

    tfhm_mean_pre[subj]=readdlm(namebase*"tfhm_mean_pre.dat")
    tfhm_mean_post[subj]=readdlm(namebase*"tfhm_mean_post.dat")
    tfhm_std_pre[subj]=readdlm(namebase*"tfhm_std_pre.dat")
    tfhm_std_post[subj]=readdlm(namebase*"tfhm_std_post.dat")
    tfhm_diff[subj]=readdlm(namebase*"tfhm_diff.dat")
    tfhm_pvals[subj]=readdlm(namebase*"tfhm_pvals.dat")

    band_stats_pre[subj]=readdlm(namebase*"band_stats_pre.dat")
    band_stats_post[subj]=readdlm(namebase*"band_stats_post.dat")
    band_pvals_alpha[subj]=readdlm(namebase*"band_pvals_alpha.dat")
    band_pvals_gamma[subj]=readdlm(namebase*"band_pvals_gamma.dat")

    Q_pre[subj]=readdlm(namebase*"Q_pre.dat")
    Q0_pre[subj]=readdlm(namebase*"Q0_pre.dat")
    Q_post[subj]=readdlm(namebase*"Q_post.dat")
    Q0_post[subj]=readdlm(namebase*"Q0_post.dat")
end

# Figures
function fillhoriz(data,dp)
    dp = cat(dp,[dp[end], dp[1], dp[1]],dims=1)
    ndat = cat(data,[-1, -1, data[1]],dims=1)
    return ndat,dp
end

subj=9
dp = depth("bipolar",132)



# Relative power filled with color
plot(xlabel="Rel. power",ylabel="depth [μm]")
plot!(fillhoriz(band_stats_pre[subj][:,3] + band_stats_pre[subj][:,4],dp), label = "α",
lc=1, fc=1,fill=true,fillalpha=0.25)
plot!(fillhoriz(band_stats_pre[subj][:,4],dp),label = "γ",lc=2,fc=2,
    fill=true,fillalpha=0.25,xlim=(0,1.0),ylim=(dp[1],dp[end]),framestyle=:box)
plot!(band_stats_post[subj][:,4],dp, label="",lc=2,ls=:dash)
plot!(band_stats_pre[subj][:,3]+band_stats_post[subj][:,4],dp, label="",lc=1,ls=:dash)

# Plain relative power
plot(xlabel="Rel. power",ylabel="depth [μm]")
plot!(band_stats_pre[subj][:,3],dp, label = "α", lc=1, fc=1, lw =2,
    xlim=(0,1.0),ylim=(dp[1],dp[end]),framestyle=:box)
plot!(band_stats_pre[subj][:,4],dp,label = "γ",lc=2,fc=2,lw =2)
plot!(band_stats_post[subj][:,3],dp, label="",lc=1,ls=:dash)
plot!(band_stats_post[subj][:,4],dp, label="",lc=2,ls=:dash)


# Straight power
subj=9
plot(xlabel="Rel. power",ylabel="depth [μm]")
plot!(band_stats_pre[subj][:,1],dp, label = "α", lc=1, fc=1, lw =2,xlim=:default,
    ylim=(dp[1],dp[end]),framestyle=:box)
plot!(band_stats_pre[subj][:,2],dp,label = "γ",lc=2,fc=2,lw =2)
plot!(band_stats_post[subj][:,1],dp, label="",lc=1,ls=:dash)
plot!(band_stats_post[subj][:,2],dp, label="",lc=2,ls=:dash)

# Power comparison
subj=9
plot(xlabel=L"\alpha_\textrm{post}/\alpha_\textrm{pre},\;\gamma_\textrm{post}/\gamma_\textrm{pre}")
plot!(band_stats_post[subj][:,1]./band_stats_pre[subj][:,1],dp, label = "α", lc=1, fc=1, lw =2,#xlim=(0,5e-13),
    ylim=(dp[1],dp[end]),framestyle=:box)
plot!(band_stats_post[subj][:,2]./band_stats_pre[subj][:,2],dp,label = "γ",lc=2,fc=2,lw =2)




# Straight power
plot(xlabel=L"\alpha_\textrm{post}/\alpha_\textrm{pre}",xlim=:default)
for subj in [8,10,11,12,13,15,16,18,20]

    plot!(band_stats_post[subj][:,1].-band_stats_pre[subj][:,1],dp, 
    label = "α", lc=1, fc=1, lw =1, ylim=(dp[1],dp[end]),framestyle=:box,legend=:false)
   
end
pα = plot!()

plot(xlabel=L"\gamma_\textrm{post}/\gamma_\textrm{pre}",xlim=:default)
for subj in [8,10,11,12,13,15,16,18,20]
    plot!(band_stats_post[subj][:,2].-band_stats_pre[subj][:,2],dp,
     lc=2, fc=1, lw =1, ylim=(dp[1],dp[end]),framestyle=:box,legend=:false)
end
pγ = plot!()

plot(pα,pγ,layout=(1,2))




# Heatmaps
using Plots.PlotMeasures
Plots.gr_cbar_width[]=0.01
Plots.gr_cbar_width[]=0.01
Plots.gr_colorbar_tick_size[]=0.001
Plots.gr_cbar_offsets[]=(0.005,0.02)
hm= Any[]
for subj in [8,10,11,12,13,15,16,18,20]
    p = heatmap(0.1:0.1:200,dp,log10.(tfhm_mean_pre[subj]./sum(tfhm_mean_pre[subj])),
    color=:rainbow1,xrange=(0.0,100.0),frame=:box,clim=(-6,-4.5), title="suj"*string(subj),
    rightmargin=2mm,label="suj"*string(subj),colorbar=true)
    # p = heatmap(0.1:0.1:200,dp,tfhm_mean_pre[subj],
    # color=:rainbow1,xrange=(0.0,100.0),frame=:box,colorbar_scale=:log10,
    # rightmargin=10mm)
    plot!(ylabel="Depth [μm]",xlabel="Freq. [Hz]")
    push!(hm,p)
end

plot(hm...,size=(1920,1080),bottommargin=10mm,leftmargin=10mm)



subj=8
p = heatmap(0.1:0.1:200,dp,tfhm_mean_pre[subj],
    color=:rainbow1,xrange=(0.0,100.0),frame=:box,
    rightmargin=10mm,colorbar_scale=:log10)
plot!(ylabel="Depth [μm]",xlabel="Freq. [Hz]")

a=rand(4,4)
plot(-2:1,-3:0,a,lt=:heatmap,colorbar_scale=:log10,xlim=:default)
plot(a,lt=:heatmap)
plot(log10.(a),lt=:heatmap)
heatmap(a,colorbar_scale=:log10)


