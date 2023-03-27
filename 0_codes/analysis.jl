#!/usr/bin/env -S julia --threads 16

push!(LOAD_PATH, "/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
cd("/home/pclusella/Documents/Data/UPO-tACs/0_codes/")
using NeuropixelAnalysis,SpectralAnalysis
using DelimitedFiles, Multitaper, Plots, DSP, Statistics,HypothesisTests
# plotlyjs()


# change it to "cortex_pre"
id=100; fl="cortex_pre"; Δt=10.0
t, f, tfhm = timefreq(id, fl,Δt);
idx, tfhm = movfilter(t, tfhm, "pre")

freqs=f;
function relative_power_(freqs,tfhm)
    df = freqs[2]-freqs[1]
    αband = findall(@. freqs>4.0 && freqs<22)
    γband = findall(@. freqs>32.0 && freqs<48.0)
    α = sum(tfhm[αband,:],dims=1)[1,:]*df 
    γ = sum(tfhm[γband,:],dims=1)[1,:]*df 
    tot = sum(tfhm[:,:],dims=1)[1,:]*df 
    # mα = mean(α)
    # mγ = mean(γ)
    # sα = std(α)
    # sγ = std(γ)
    # return mα,mγ,sα,sγ
    return α,γ,tot
end

n0=226;
nf=361;
function prepost_analysis(n0,nf)
    l = 2000;
    n=size(n0:nf)[1]
    ns=[ n, n, n/2-2,n-4,384]
    ns=Int.(ns)

    stats_band_pre = Dict();
    stats_band_post = Dict();
    pvals_α = Dict();
    pvals_γ = Dict();
    for (i,data) in enumerate(("kcsd_","cortex_","csd_","bipolar_","filtered_"))
    # for (i,data) in enumerate(("kcsd_",))
        n0=ns[i];
        stats_band_pre[data], stats_band_post[data], pvals_α[data], pvals_γ[data] = prepost_comparison(data,n0)
    end

    return stats_band_pre, stats_band_post, pvals_α , pvals_γ
end

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

function prepost_comparison(data,n0)
    l = 2000

    psd_mean_tfhm_pre = zeros(n0, l);
    psd_std_tfhm_pre = zeros(n0, l);
    psd_mean_tfhm_post = zeros(n0, l);
    psd_std_tfhm_post = zeros(n0, l);
    pvals_tfhm = zeros(n0, l);
    pvals_uneq_tfhm = zeros(n0, l);
    pvals_perm_tfhm = zeros(n0, l);

    stats_band_pre = zeros(n0,8)
    stats_band_post = zeros(n0,8)
    pvals_α = zeros(n0,2)
    pvals_γ = zeros(n0,2)

    state = Threads.Atomic{Int}(0);

    Threads.@threads for id in 0:n0-1
        # Prompt state
        Threads.atomic_add!(state, 1)
        print("--> ", state[], " out of ", n0, "\n")
        flush(stdout)

        # Pre
        t, f, tfhm = timefreq(id, data*"pre")
        idx, tfhm_pre = movfilter(t, tfhm, "pre")
        psd_mean_tfhm_pre[id+1, :] = mean(tfhm_pre, dims=2)
        psd_std_tfhm_pre[id+1, :] = std(tfhm_pre, dims=2)

        # Post
        t, f, tfhm = timefreq(id, data*"post")
        idx, tfhm_post = movfilter(t, tfhm, "post")
        psd_mean_tfhm_post[id+1, :] = mean(tfhm_post, dims=2)
        psd_std_tfhm_post[id+1, :] = std(tfhm_post, dims=2)

        # band analysis 
        α_pre, γ_pre, total_pre = relative_power_(f,tfhm_pre)
        stats_band_pre[id+1,1:8] .=
                mean(α_pre),mean(γ_pre),mean(α_pre./total_pre),mean(γ_pre./total_pre),
                std(α_pre), std(γ_pre), std(α_pre./total_pre), std(γ_pre./total_pre)

        α_post, γ_post, total_post = relative_power_(f,tfhm_post)
        stats_band_post[id+1,1:8] .=
                mean(α_post),mean(γ_post),mean(α_post./total_post),mean(γ_post./total_post),
                std(α_post), std(γ_post), std(α_post./total_post), std(γ_post./total_post)

        pvals_α[id+1,1] = pvalue(ApproximatePermutationTest(α_pre, α_post, mean, 1000))
        pvals_γ[id+1,1] = pvalue(ApproximatePermutationTest(γ_pre, γ_post, mean, 1000))
        pvals_α[id+1,2] = pvalue(ApproximatePermutationTest(α_pre./total_pre, α_post./total_post, mean, 1000))
        pvals_γ[id+1,2] = pvalue(ApproximatePermutationTest(γ_pre./total_pre, γ_post./total_post, mean, 1000))

        # T-test. 
        # for j in 1:l
        #     pvals_uneq_tfhm[id+1, j] = pvalue(UnequalVarianceTTest(tfhm_pre[j, :], tfhm_post[j, :]))
        #     pvals_tfhm[id+1, j] = pvalue(EqualVarianceTTest(tfhm_pre[j, :], tfhm_post[j, :]))# can be computed from the means
        #     pvals_perm_tfhm[id+1,j] = pvalue(ApproximatePermutationTest(tfhm_pre[j, :], tfhm_post[j, :], mean, 1000))
        # end
    end

    return stats_band_pre, stats_band_post, pvals_α, pvals_γ

end

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
xlabel = "Rel. power (pre)", ylabel = "depth [μm]")


plot(spre[data][:,3]-spost[data][:,3],dp[data])
p3 = plot!(spre[data][:,4]-spost[data][:,4],dp[data],xlim=(-0.2,0.2),
xlabel = "Difference", ylabel = "depth [μm]")

plot(p1,p2,p3,layout=(1,3),labels=["α" "γ"],size=(800,450),lw=3)
savefig("/home/pclusella/Desktop/"*data[1:end-1]*".png")