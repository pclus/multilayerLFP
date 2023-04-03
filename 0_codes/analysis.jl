#!/usr/bin/env -S julia --threads 16

cd("0_codes/")
push!(LOAD_PATH, pwd())
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
using StatsBase

Δt=1
t, f, tfhm = timefreq(id, fl, Δt);
idx, tfhm = movfilter(t, tfhm, "pre")

freqs=f;

n0=226;
nf=361;

# write -------------------------------------------------------
spre,spost,pα,pγ = prepost_analysis(n0,nf)

for (i,data) in enumerate(("kcsd_","cortex_","csd_","bipolar_","filtered_"))
    writedlm("../4_outputs/"*data*"spre.dat", spre[data], ' ');
    writedlm("../4_outputs/"*data*"spost.dat", spost[data], ' ');
    writedlm("../4_outputs/"*data*"palpha.dat", pα[data], ' ');
    writedlm("../4_outputs/"*data*"pgamma.dat", pγ[data], ' ');
end
#--------------------------------------------------------------

# read --------------------------------------------------------
spre = Dict();
spost = Dict();
pα = Dict();
pγ = Dict();

for (i,data) in enumerate(("kcsd_","cortex_","csd_","bipolar_","filtered_"))
    spre[data] =readdlm("../4_outputs/"*data*"spre.dat");
    spost[data] = readdlm("../4_outputs/"*data*"spost.dat");
    pα[data] = readdlm("../4_outputs/"*data*"palpha.dat");
    pγ[data] = readdlm("../4_outputs/"*data*"pgamma.dat");
end
#--------------------------------------------------------------



dp = Dict()
for (i,data) in enumerate(("kcsd_","cortex_","csd_","bipolar_","filtered_"))
    dp[data] = depth(data[1:end-1],size(spre[data])[1])
end


# plot(spre[data][:,3],dp[data],xerr=spre[data][:,7],msc=:auto)
# p1 = plot!(spre[data][:,4],dp[data],xerr=spre[data][:,7],
# xlim=(0,1),msc=:auto,msw=0.1,
# xlabel = "Rel. power (pre)", ylabel = "depth [μm]")

data="bipolar_"

plot(spre[data][:,3],dp[data])
p1 = plot!(spre[data][:,4],dp[data],
xlim=(0,1),msc=:auto,
xlabel = "Rel. power (pre)", ylabel = "depth [μm]")

plot(spost[data][:,3],dp[data],msc=:auto)
p2 = plot!(spost[data][:,4],dp[data],
xlim=(0,1),msc=:auto,
xlabel = "Rel. power (post)", ylabel = "depth [μm]")


plot(-spre[data][:,3]+spost[data][:,3],dp[data])
p3 = plot!(-spre[data][:,4]+spost[data][:,4],dp[data],xlim=(-0.2,0.2),
xlabel = "Difference (post-pre)", ylabel = "depth [μm]")

plot(p1,p2,p3,layout=(1,3),labels=["α" "γ"],size=(800,450),lw=3)
savefig("/home/pclusella/Desktop/"*data[1:end-1]*".png")