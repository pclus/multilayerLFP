#!/usr/bin/env -S julia --threads 16

cd("0_codes/")
push!(LOAD_PATH, pwd())
using NeuropixelAnalysis #,SpectralAnalysis
using DelimitedFiles, Plots, DSP, Statistics,HypothesisTests
using StatsBase
using ColorSchemes,LaTeXStrings,Measures

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

for subj in [8,10,11,12,13,15,16,18,20]
    data = "bipolar_"
    foutname = "suj"*string(subj)*"/"
    # namebase = "../4_outputs/pipeline/full/"*foutname*data
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
layers = readdlm("../1_data/layers.dat",skipstart=5)
n0 = 384-layers[end,1];
nf = 384-layers[2,1]+1;
dp = @. depth(n0:nf)
dps = @.depth(384-layers); dps = dps[2:end,:]


ss = 20
k = [dp[1:2:end] band_stats_pre[ss][1:2:end,3] band_stats_pre[ss][1:2:end,4]]
writedlm("../data_subj"*string(ss)*".dat",k,' ')

# total power, missing from the analysis...
total_pre=Dict()
for (i,tfhm) in tfhm_mean_pre
    total_pre[i] = sum(tfhm,dims=2)*0.1
end
total_post=Dict()
for (i,tfhm) in tfhm_mean_post
    total_post[i] = sum(tfhm,dims=2)*0.1
end


# Figures
function fillhoriz(data,dp)
    dp = cat(dp,[dp[end], dp[1], dp[1]],dims=1)
    ndat = cat(data,[-1, -1, data[1]],dims=1)
    return ndat,dp
end


# subj=9
# # Relative power filled with color
# plot(xlabel="Rel. power",ylabel="depth [μm]")
# plot!(fillhoriz(band_stats_pre[subj][:,3] + band_stats_pre[subj][:,4],dp), label = "α",
# lc=1, fc=1,fill=true,fillalpha=0.25)
# plot!(fillhoriz(band_stats_pre[subj][:,4],dp),label = "γ",lc=2,fc=2,
#     fill=true,fillalpha=0.25,xlim=(0,1.0),ylim=(dp[1],dp[end]),framestyle=:box)
# plot!(band_stats_post[subj][:,4],dp, label="",lc=2,ls=:dash)
# plot!(band_stats_pre[subj][:,3]+band_stats_post[subj][:,4],dp, label="",lc=1,ls=:dash)

#-----------------------------------------------------
# Relative power multipanel
#-----------------------------------------------------
# Multipanel figure with pre and post relative power (i.e., amount of α and γ within the layer).

pls = Any[]
pythonplot()
for subj in [8,10,11,12,13,15,16,18,20]
    plot(xlabel="Rel. power",ylabel="depth [μm]",title="suj"*string(subj))
    plot!(band_stats_pre[subj][:,3],dp, label = "αₚᵣₑ", lc=1, fc=1, lw =2,
        xlim=(0,1.0),ylim=(dp[1],dp[end]),framestyle=:box)
    plot!(band_stats_pre[subj][:,4],dp,label = "γₚᵣₑ",lc=2,fc=2,lw =2)
    plot!(band_stats_post[subj][:,3],dp, label="αₚₒₛₜ",lc=1,ls=:dash)
    p = plot!(band_stats_post[subj][:,4],dp, label="γₚₒₛₜ",lc=2,ls=:dash)
    push!(pls,p)
end
plot(pls...,layout=(2,5),size=(600*3,400*2),legend=:bottomright,xformatter=:auto,
bottommargin=10mm,topmargin=2mm)
savefig("../3_figures/global/Fig_RelPower.pdf")

#-----------------------------------------------------
# Relative power mean
#-----------------------------------------------------
# Same as multipanel, but with mean for each case 
# Grey lines depict individual results.
# We use strider to select only 2 out of the 4 neuropixel columns.

α= [ band_stats_pre[i][j,3] for j in 1:142 , (i,a) in band_stats_pre]
γ= [ band_stats_pre[i][j,4] for j in 1:142 , (i,a) in band_stats_pre]

mα = mean(α,dims=2)[:,1]
sα = std(α,dims=2)[:,1]
mγ = mean(γ,dims=2)[:,1]
sγ = std(γ,dims=2)[:,1]

strider = 1:2:142
plot()
for subj in [8,10,11,12,13,15,16,18,20]
    plot!(band_stats_pre[subj][strider,3],dp[strider],lw=0.1,lc=:grey,label="")
    plot!(band_stats_pre[subj][strider,4],dp[strider],lw=0.1,lc=:grey,label="")
end

plot!(mα[strider],dp[strider],lc=1,label="α")
plot!(mγ[strider],dp[strider],lc=2,label="γ")

L = @. depth(384-layers[2:end,:])[:,1]
hline!(L[1:2:end],lc=:black,lw=0.5,label="")
plot!(xlabel="Rel. power",ylabel="Depth [μm]")
plot!(size=(300,400))
savefig("../3_figures/global/mean_rel_power.pdf")

writedlm("../bipolar_relpower_9subjects.dat",[dp[strider] mα[strider] mγ[strider]],' ')

#-----------------------------------------------------
# Total power
#-----------------------------------------------------
# Figure showing the α and γ powers.
# Here the only normalization is to make sure that power accross subjects are comparable
# This means that for each subject we compute the sum of all powers across layers and normalize 
# by this quantity. In this way the changes of power across layers for a single subject 
# are mantained.
tot = Dict()
for (i,a) in total_pre 
    tot[i] = sum(total_pre[i]) 
end
α= [ band_stats_pre[i][j,1]/tot[i] for j in 1:142 , (i,a) in band_stats_pre]
γ= [ band_stats_pre[i][j,2]/tot[i] for j in 1:142 , (i,a) in band_stats_pre]

mα = mean(α,dims=2)[:,1]
sα = std(α,dims=2)[:,1]
mγ = mean(γ,dims=2)[:,1]
sγ = std(γ,dims=2)[:,1]

strider = 1:2:142

plot()
for subj in 1:9
    plot!(α[strider,subj],dp[strider],lw=0.1,lc=:grey,label="")
    plot!(γ[strider,subj],dp[strider],lw=0.1,lc=:grey,label="")
end
plot!(mα[strider],dp[strider],lc=1,label="α")
plot!(mγ[strider],dp[strider],lc=2,label="γ")

plot!(xaxis=:log)
L = @. depth(384-layers[2:end,:])[:,1]
hline!(L[1:2:end],lc=:black,lw=0.5,label="")
plot!(xlabel="Power",ylabel="Depth [μm]")
plot!(size=(300,400))
savefig("../3_figures/global/mean_power.pdf")

#-----------------------------------------------------
# Bastos relative power (mean)
#-----------------------------------------------------
# Departing from the averages obtained in the previous plot,
# we divide by the peak α and peak γ. 
# We just do this for the mean, so the individual trajectories might display
# peaks above or below 1.

maxα = maximum(mα[strider])
maxγ = maximum(mγ[strider])
plot()
for subj in 1:9
    plot!(α[strider,subj]./maxα,dp[strider],lw=0.1,lc=:grey,label="")
    plot!(γ[strider,subj]./maxγ,dp[strider],lw=0.1,lc=:grey,label="")
end
plot!(mα[strider]./maxα,dp[strider],label="α",lc=1)
plot!(mγ[strider]./maxγ,dp[strider],label="γ",lc=2)
hline!(L[1:2:end],lc=:black,lw=0.5,label="")
plot!(xlabel="Power/Max. Power",ylabel="Depth [μm]")
plot!(size=(300,400),xlim=(0,1.2))
savefig("../3_figures/global/mean_relmax_power.pdf")

#-----------------------------------------------------
# Bastos relative power (individual)
#-----------------------------------------------------
# Again, we divide each power by the peak across layers,
# but in this case we do it for each subject separately.
# Therefore, each subject and each power display a layer with rel. power equal to 1.
# But the averages are different. 

strider = 1:2:142
α= [ band_stats_pre[i][j,1]/maximum(band_stats_pre[i][strider,1]) for j in 1:142 , (i,a) in band_stats_pre]
γ= [ band_stats_pre[i][j,2]/maximum(band_stats_pre[i][strider,2]) for j in 1:142 , (i,a) in band_stats_pre]

mα = mean(α,dims=2)[:,1]
sα = std(α,dims=2)[:,1]
mγ = mean(γ,dims=2)[:,1]
sγ = std(γ,dims=2)[:,1]

plot()
for subj in 1:9
    maxα = maximum(α[strider,subj])
    maxγ = maximum(γ[strider,subj])
    plot!(α[strider,subj],dp[strider],lw=0.1,lc=:grey,label="")
    plot!(γ[strider,subj],dp[strider],lw=0.1,lc=:grey,label="")
end
plot!()

plot!(mα[strider]./maxα,dp[strider],label="α",lc=1)
plot!(mγ[strider]./maxγ,dp[strider],label="γ",lc=2)
hline!(L[1:2:end],lc=:black,lw=0.5,label="")
plot!(xlabel="Power/Max. Power",ylabel="Depth [μm]")
plot!(size=(300,400),xlim=(0,1.0))
savefig("../3_figures/global/mean_relmax_power2.pdf")

# Multipanel version:

pls = Any[]
pythonplot()
for subj in 1:9
    plot(xlabel="power/max. power",ylabel="depth [μm]") #,title="suj"*string(subj))
    maxα = maximum(α[strider,subj])
    maxγ = maximum(γ[strider,subj])
    plot!(α[strider,subj],dp[strider],label="α")
    p=plot!(γ[strider,subj],dp[strider],label="γ")

    push!(pls,p)
end
plot(pls...,layout=(2,5),size=(600*3,400*2),legend=:bottomright,xformatter=:auto,
bottommargin=10mm,topmargin=2mm)
savefig("../3_figures/global/Fig_RelPower2.pdf")



#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------


# Total power diverges a lot, so it makes sense to normalize
plot()
for subj in [8,10,11,12,13,15,16,18,20]
    plot!(total_pre[subj][:],dp, lc=3, fc=1, lw =1,#xlim=(0,5e-13),
        ylim=(dp[1],dp[end]),framestyle=:box,label="pre",title="suj"*string(subj))
end
plot!(xaxis=:log,lw=1,legend=false)

#-----------------------------------------------------
pythonplot()
plot(xlabel="Rel. power",ylabel="depth [μm]")
for subj in [8,9,10,11,12,13,15,16,18,20]
    plot!(band_stats_pre[subj][:,3],dp, label = "αₚᵣₑ suj"*string(subj), lc=1, fc=1, lw =1,
        xlim=(0,1.0),ylim=(dp[1],dp[end]),framestyle=:box)
end
plot!(lw=1)

pythonplot()
plot(xlabel="Rel. power",ylabel="depth [μm]")
for subj in [8,9,10,11,12,13,15,16,18,20]
    plot!(band_stats_pre[subj][:,4],dp, label = "γₚᵣₑ suj"*string(subj), lc=2, fc=1, lw =1,
        xlim=(0,1.0),ylim=(dp[1],dp[end]),framestyle=:box)
end
plot!(lw=1)

plot!(xlim=(0.05,0.15))

# plot(pls...,layout=(2,5),size=(600*3,400*2),legend=:bottomright,xformatter=:auto,
bottommargin=10mm,topmargin=2mm)
# savefig("../3_figures/global/Fig_RelPower.pdf")






# # Straight power
# subj=10
# plot(xlabel="bLFP power",ylabel="depth [μm]")
# plot!(band_stats_pre[subj][:,1],dp, label = L"\alpha_\textrm{pre}", lc=1, fc=1, lw =2,xlim=:default,
#     ylim=(dp[1],dp[end]),framestyle=:box)
# plot!(band_stats_pre[subj][:,2],dp,label = L"\gamma_\textrm{pre}",lc=2,fc=2,lw =2)
# plot!(band_stats_post[subj][:,1],dp, label=L"\alpha_\textrm{post}",lc=1,ls=:dash)
# plot!(band_stats_post[subj][:,2],dp, label=L"\gamma_\textrm{post}",lc=2,ls=:dash)

# Power comparison

pls = Any[]
# gr()
pythonplot()
for subj in [8,9,10,11,12,13,15,16,18,20]
    plot(xlabel="post/pre",ylabel="Depth [μm]")
    plot!(ones(size(dp)),dp,lc=:black,ls=:dash,label="",title="suj"*string(subj))
    plot!(band_stats_post[subj][:,1]./band_stats_pre[subj][:,1],dp, lc=1, fc=1, lw =2,#xlim=(0,5e-13),
        ylim=(dp[1],dp[end]),framestyle=:box,label="αₚₒₛₜ/αₚᵣₑ")
    p = plot!(band_stats_post[subj][:,2]./band_stats_pre[subj][:,2],dp,lc=2,fc=2,lw =2,
    label=L"γₚₒₛₜ/γₚᵣₑ")
    push!(pls,p)
end
plot(pls...,layout=(2,5),size=(600*3,400*2),key=true,xformatter=:auto,
bottommargin=10mm,topmargin=2mm,legend=:bottomright)
savefig("../3_figures/global/Fig_PostPre.pdf")


# Total power comparison-----------------
# plotlyjs()
# gr()
pythonplot()
pls = Any[]
for subj in [8,10,11,12,13,15,16,18,20]
    plot(xlabel="Total power",ylabel="Depth [μm]")
    plot!(total_pre[subj][:],dp, lc=3, fc=1, lw =2,#xlim=(0,5e-13),
        ylim=(dp[1],dp[end]),framestyle=:box,label="pre",title="suj"*string(subj))
    p=plot!(total_post[subj][:],dp, lc=3, fc=1, lw =1,ls=:dash,
        ylim=(dp[1],dp[end]),framestyle=:box,label="post")
        push!(pls,p)
end

using Plots.Measures
plot(pls...,layout=(2,5),size=(600*3,400*2),key=false,xformatter=:auto,
bottommargin=15mm,topmargin=2mm)
savefig("../3_figures/global/Fig_TotPower.pdf")

# Straight power
# plot(xlabel=L"\alpha_\textrm{post}/\alpha_\textrm{pre}",xlim=:default)
# for subj in [8,10,11,12,13,15,16,18,20]

#     plot!(band_stats_post[subj][:,1].-band_stats_pre[subj][:,1],dp, 
#     label = "α", lc=1, fc=1, lw =1, ylim=(dp[1],dp[end]),framestyle=:box,legend=:false)
   
# end
# pα = plot!()

# plot(xlabel=L"\gamma_\textrm{post}/\gamma_\textrm{pre}",xlim=:default)
# for subj in [8,10,11,12,13,15,16,18,20]
#     plot!(band_stats_post[subj][:,2].-band_stats_pre[subj][:,2],dp,
#      lc=2, fc=1, lw =1, ylim=(dp[1],dp[end]),framestyle=:box,legend=:false)
# end
# pγ = plot!()

# plot(pα,pγ,layout=(1,2))


# NEEDS TO ACCOUNT FOR ACTUAL LAYERS!
using Gnuplot
Gnuplot.options
@gp "
set term pngcairo enhanced size 3000,1000;
set pale @RAINBOW;
set xrange[0:100]; 
set logs zcb;
set yrange[-1545.92:-245.92];
set cbrange[1e-6:0.5e-4];
unset key;
set border lw 1;
set multiplot layout 2,5 title 'Normalized power PRE' ;
set xlabel 'Freq [Hz]'
set ylabel 'Depth [μm]'
set colorbox
set cbtics format '1e%T'" :-

for subj in [8,10,11,12,13,15,16,18,20]
    @gp :- subj "set title 'suj"*string(subj)*"'"
    @gp :- subj tfhm_mean_pre[subj]./sum(tfhm_mean_pre[subj]) "origin=(0, -1545.92) dx=0.1 dy=10 w image"
end
p1 = @gp :- ""

display(MIME("image/png"),p1)

# POST
Gnuplot.options
@gp "
set term pngcairo enhanced size 3000,1000;
set pale @RAINBOW;
set xrange[0:100]; 
set logs zcb;
set yrange[-1545.92:-245.92];
set cbrange[1e-6:0.5e-4];
unset key;
set border lw 1;
set multiplot layout 2,5 title 'Normalized power PRE' ;
set xlabel 'Freq [Hz]'
set ylabel 'Depth [μm]'
set colorbox
set cbtics format '1e%T'" :-

for subj in [8,9,10,11,12,13,15,16,18,20]
    if subj==9
        @gp :- subj "set multiplot next;" :-
    else
        @gp :- subj "set title 'suj"*string(subj)*"'"
        @gp :- subj tfhm_mean_pre[subj]./sum(tfhm_mean_pre[subj]) "origin=(0, -1545.92) dx=0.1 dy=10 w image"

    end
end
p1 = @gp :- ""

display(MIME("image/png"),p1)


#---------------------------------------------------------------
# LOAD FULL

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

for subj in [8,10,11,12,13,15,16,18,20]
    data = "bipolar_"
    foutname = "suj"*string(subj)*"/"
    namebase = "../4_outputs/pipeline/full/"*foutname*data

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
# dp = depth("bipolar",380).-2239.08

using Gnuplot
Gnuplot.options
@gp "
set term pngcairo enhanced size 3000,1000;
set pale @RAINBOW;
set xrange[0:100]; 
set logs zcb;
set yrange[-3840:0.0];
set cbrange[0.5e-6:0.1e-4];
unset key;
set border lw 1;
set multiplot layout 2,5 title 'Normalized power PRE' ;
set xlabel 'Freq [Hz]'
set ylabel 'Depth [μm]'
set colorbox
set cbtics format '1e%T'" :-

kk=1
for subj in [8,9,10,11,12,13,15,16,18,20]
    if subj==9
        @gp :- subj "set multiplot next;" :-
    else
        a = string(dps[end,kk])
        b = string(dps[1,kk])
        @gp :- subj "set object 1 rect from 0,"*a*" to 100,"*b*" front fillstyle empty lw 5"
        @gp :- subj "set title 'suj"*string(subj)*"'"
        @gp :- subj tfhm_mean_pre[subj]./sum(tfhm_mean_pre[subj]) "origin=(0, -3840.0) dx=0.1 dy=10 w image"
    end
    kk = kk+1
end
p1 = @gp :- ""

display(MIME("image/png"),p1)