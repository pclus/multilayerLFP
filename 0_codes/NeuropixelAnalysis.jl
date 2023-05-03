module NeuropixelAnalysis

using DelimitedFiles, Multitaper, Plots, DSP,Statistics,HypothesisTests,FFTW,MAT,Distributions

export DelimitedFiles, Multitaper, Plots, DSP,Statistics,HypothesisTest,FFTW,MAT,Distributions
export read_channel,channel_idx,heatmapMT,movfilter,process_data,relative_power,depth,prepost_analysis
export compute_bipolar_horizontal
export logspectral_dist,export_matfile


"""
    export_matfile(filename,varname)

Reads the original matlab files and stores pre and post data as binaries,
and the movement data in plain text.
"""
function export_matfile(filename,varname)
    file = matopen(filename);
    data = read(file,varname);
    @views write("../1_data/pre.bin",data[1,1][384:-1:1,:]');
    @views write("../1_data/post.bin",data[2,1][384:-1:1,:]');
    writedlm("../1_data/mov_pre.dat",data[1,2],' ');
    writedlm("../1_data/mov_post.dat",data[2,2],' ');
    close(file)

    mpre  = size(data[1,1])[2]
    mpost = size(data[2,1])[2]
    return mpre,mpost
end



"""
    read_channel(id, t0, tf, fl)

Obtain the time series indexed by `id` between times `t0` and `tf` 
from the binary file identified by `fl`.

# Examples
```julia-repl
julia> v158=read_channel(158,0.0004,900,"filtered_pre");
```
"""
function read_channel(id,fl; t0 = 0.0004, tf=3600.0, m=9000000,kpath="local") # Equivalent to the read_binary.c code
    id = Int(id)
    n = 384
    dt = 1.0 / 2500.0

    if kpath=="local"
        path = "../1_data/"
    else
        path = "/home/pclusella/Documents/Athena/UPO_data_processed/"*kpath*"/"
    end

    data = zeros(m)
    fin = open(path * fl * ".bin", "r")
    for i in 0:id
        read!(fin, data)
    end
    close(fin)

    m0 = Int(ceil(t0 / dt))
    mf = Int(floor(tf / dt))
    dm = mf - m0 + 1

    channel = data[m0:mf]
    t = collect(dt * (m0:mf))

    return t, channel
end

"""
    cut_cortex(fl,channels)

Read the channels listed in the array `channels` from the "fl" binary file,
filter them, and write in a 'cortex_'*fl binary file.
Already bandpass_filtered.
"""
function cut_cortex(fl,n0,nf,m=9000000) 

    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = dt*m

    fout = open("../1_data/cortex_" * fl * ".bin", "w")
    # fout = open("../1_data/cortex_" * fl * ".bin", "w")

    bpfilter = digitalfilter(Bandpass(1.0, 300.0; fs=rate), Butterworth(3))
    # bsfilter  = digitalfilter(Bandstop(49.9,50.1; fs=rate),Butterworth(1))
    # bsfilter2 = digitalfilter(Bandstop(59.9,60.1; fs=rate),Butterworth(1))

    for id in n0:nf
        t, ch = read_channel(id, fl; t0=t0, tf=tf,m=m)
        fil_ch = filtfilt(bpfilter, ch)
        # fil_ch = filtfilt(bsfilter, fil_ch)
        # fil_ch = filtfilt(bsfilter2, fil_ch)
        write(fout, fil_ch)
    end
    close(fout)
end


"""
    compute_bipolar(fl)

Store in a binary file the bipolar data obtained by computing the first 
derivative of the data using a symmetric stencil in along each neuropixel column.
"""
function compute_bipolar(flin,flout,n,m=9000000)
    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = m*dt

    fout = open("../1_data/"* flout * ".bin", "w")

    t, prev_a = read_channel(0, flin; t0=t0, tf=tf, m=m)
    t, prev_b = read_channel(1, flin; t0=t0, tf=tf, m=m)
    t, prev_c = read_channel(2, flin; t0=t0, tf=tf, m=m)
    t, prev_d = read_channel(3, flin; t0=t0, tf=tf, m=m)

    inv_dy = 1 / 20.0 # 1/μm

    nf = n - n%4 
    for id in 4:4:nf-1
        t, curr_a = read_channel(id    , flin; t0=t0, tf=tf, m=m)
        t, curr_b = read_channel(id + 1, flin; t0=t0, tf=tf, m=m)
        t, curr_c = read_channel(id + 2, flin; t0=t0, tf=tf, m=m)
        t, curr_d = read_channel(id + 3, flin; t0=t0, tf=tf, m=m)

        bip_a = - inv_dy * (prev_a - curr_a)
        bip_b = - inv_dy * (prev_b - curr_b)
        bip_c = - inv_dy * (prev_c - curr_c)
        bip_d = - inv_dy * (prev_d - curr_d)

        prev_a = curr_a
        prev_b = curr_b
        prev_c = curr_c
        prev_d = curr_d

        write(fout, bip_a)
        write(fout, bip_b)
        write(fout, bip_c)
        write(fout, bip_d)

    end
    close(fout)

end

"""
    compute_bipolar(fl)

Store in a binary file the  horizontal difference between pairs of electrodes at the same depth
"""
function compute_bipolar_horizontal(flin,flout,n,m=9000000)
    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = m*dt

    fout = open("../1_data/"* flout * ".bin", "w")

    inv_dy = 1 / 20.0 # 1/μm

    # nf = n - n%2
    nf = n
    for id in 0:2:nf-1
        t, curr_a = read_channel(id    , flin; t0=t0, tf=tf, m=m)
        t, curr_b = read_channel(id + 1, flin; t0=t0, tf=tf, m=m)

        bip_0 = - inv_dy * (curr_a - curr_b)

        write(fout, bip_0)

    end
    close(fout)

end

"""
    compute_csd(fl)

Store in a binary file the CSD obtained by computing the Laplacian 
using a 2-dimensional diagonal finite difference stencil.
"""
function compute_csd(flin,flout,n,m=9000000)
    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = dt*m

    fout = open("../1_data/" * flout * ".bin", "w")

    t, ch_5 = read_channel(0, flin; t0=t0, tf=tf, m=m)
    t, ch_6 = read_channel(1, flin; t0=t0, tf=tf, m=m)
    t, ch_7 = read_channel(2, flin; t0=t0, tf=tf, m=m)
    t, ch_8 = read_channel(3, flin; t0=t0, tf=tf, m=m)

    inv_dy2 = 1 / 25.61^2 # 1/μm
    σ = 1.0 # conductivity S/μm or 1/ (Ω μm)

    nf = n - n%4 
    for id in 4:4:nf-1

        ch_1 = ch_5
        ch_2 = ch_6
        ch_3 = ch_7
        ch_4 = ch_8

        t, ch_5 = read_channel(id + 0, flin; t0=t0, tf=tf, m=m)
        t, ch_6 = read_channel(id + 1, flin; t0=t0, tf=tf, m=m)
        t, ch_7 = read_channel(id + 2, flin; t0=t0, tf=tf, m=m)
        t, ch_8 = read_channel(id + 3, flin; t0=t0, tf=tf, m=m)

        csd_right = -σ * inv_dy2 * (-4.0 * ch_4 + ch_1 + ch_2 + ch_5 + ch_6)
        csd_left  = -σ * inv_dy2 * (-4.0 * ch_5 + ch_3 + ch_4 + ch_7 + ch_8)

        write(fout, csd_right)
        write(fout, csd_left)

    end
    close(fout)

end

"""
    channel_idx(chid)

Given the id of a channel, return its index for the csd, if it exists.
"""
function channel_idx(ch_id)
    chid=ch_id+1
    i = -1
    if chid % 4 != 0 && chid % 4 != 1
        print("Channel not available")
    else
        if chid % 4 == 0
            i = chid / 2 - 1.0
        else
            i = (chid - 1) / 2
        end
    end
    return i
end


"""
    movfilter(t, tfhm, fl)
Authomatic removal of segments including movement data.
If Q is is assigned, then it also removes segments with Q
larger than a certain threshold.
This function assumes segments of Δt=10
"""
function movfilter(t, tfhm, fl;Δt=10.0,q = zeros(1),kpath="local")

    if kpath=="local"
        path = "../1_data/"
    else
        path = "/home/pclusella/Documents/Athena/UPO_data_processed/"*kpath*"/"
    end

    mov_pre = readdlm(path * fl * ".dat", ' ')
    if q!=zeros(1)
        d=Normal(mean(q),std(q))
        ctop = quantile(d,0.975)
        tq = findall(==(1),q.>ctop)
        mov_pre = cat(mov_pre,(tq.-0.5)*Δt,dims=1)
    end
    mov_indx = floor.(mov_pre ./ Δt) * Δt .+ 0.5*Δt
    indx = t .∈ [mov_indx]
    indx = findall(==(0), indx)
    f_tfhm = tfhm[:, indx]
    indx, f_tfhm
end

"""
    process_data()
    Bandpass-filters the pre and post data, and computes the CSD 
    and bipolar binary files. 
"""
function process_data(n0,nf;mpre=9000000,mpost=9000000)
    n=length(n0:nf)
    # cut_cortex("pre",n0,nf, mpre)
    # cut_cortex("post",n0,nf, mpost)
    # compute_bipolar("cortex_pre","bipolar_pre",n,mpre)   
    # compute_bipolar("cortex_post","bipolar_post",n,mpost)
    compute_bipolar("pre","bipolar_pre",n,mpre)   
    compute_bipolar("post","bipolar_post",n,mpost)  
    # compute_csd("cortex_pre","csd_pre",n,mpre)           
    # compute_csd("cortex_post","csd_post",n,mpost)        
    compute_bipolar_horizontal("pre" ,"horizontal_pre",n,mpre)
    compute_bipolar_horizontal("post","horizontal_post",n,mpost) 
end


"""
    prepost_analysis(n0,nf)

Compute the mean spectra and perform comparison between pre and post for all datatypes.
"""
function prepost_analysis(n0,nf,subj;mpre=9000000,mpost=9000000) # ISSUE: incorporate datatypes in the function call
    l = 2000;
    n=size(n0:nf)[1]
    ns=[n-4, n, n, n/2-2,]
    ns=Int.(ns)

    stats_band_pre = Dict();
    stats_band_post = Dict();
    pvals_α = Dict();
    pvals_γ = Dict();
    Qpre = Dict();
    Q0pre = Dict();
    Qpost = Dict();
    Q0post = Dict();
    # for (i,data) in enumerate(("bipolar_","kcsd_","cortex_","csd_"))
    for (i,data) in enumerate(("bipolar_",))
        n0=ns[i];
        # stats_band_pre[data], stats_band_post[data], pvals_α[data], pvals_γ[data], 
        # Qpre[data], Q0pre[data], Qpost[data], Q0post[data] =
        prepost_comparison(data,n0,subj;mpre = mpre , mpost = mpost);
    end

    # return stats_band_pre, stats_band_post, pvals_α , pvals_γ
    print("\nAll done!\n")

    return 
end


function prepost_comparison(data,n0,subj; mpre=9000000, mpost=9000000)
    l = 2000

    psd_mean_tfhm_pre = zeros(n0, l);
    psd_std_tfhm_pre = zeros(n0, l);
    psd_mean_tfhm_post = zeros(n0, l);
    psd_std_tfhm_post = zeros(n0, l);
    # pvals_tfhm = zeros(n0, l);
    # pvals_uneq_tfhm = zeros(n0, l);
    pvals_perm_tfhm = zeros(n0, l);

    stats_band_pre = zeros(n0,10)
    stats_band_post = zeros(n0,10)
    pvals_α = zeros(n0,2)
    pvals_γ = zeros(n0,2)

    # ns_pre = NeuropixelAnalysis.numberofsegments("pre";m=mpre)
    # ns_post = NeuropixelAnalysis.numberofsegments("post";m=mpost)
    ns_pre = Int(floor(mpre/25000));
    ns_post = Int(floor(mpost/25000));
    Q_pre = zeros(n0,ns_pre)
    Q0_pre = zeros(n0,ns_pre)
    Q_post = zeros(n0,ns_post)
    Q0_post = zeros(n0,ns_post)
    # pvals_Q = zeros(n0,2)

    state = Threads.Atomic{Int}(0);

    Threads.@threads for id in 0:n0-1
        # Prompt state
        Threads.atomic_add!(state, 1)
        print("--> ", state[], " out of ", n0, "\n")
        flush(stdout)


        # Pre
        f,tfhm_pre,Q_pre[id+1,:],Q0,tr = tfhm_analysis(id, data*"pre", "pre" ; m=mpre,kpath=string(subj)); Q0_pre[id+1,tr]=Q0;
        # t, f, tfhm = timefreq(id, data*"pre"; m = mpre)
        # idx, tfhm_pre = movfilter(t, tfhm, "pre")
        psd_mean_tfhm_pre[id+1, :] = mean(tfhm_pre, dims=2)
        psd_std_tfhm_pre[id+1, :] = std(tfhm_pre, dims=2)

        # Post
        f,tfhm_post,Q_post[id+1,:],Q0,tr = tfhm_analysis(id, data*"post", "post" ; m=mpost,kpath=string(subj)); Q0_post[id+1,tr]=Q0;
        # t, f, tfhm = timefreq(id, data*"post"; m = mpost)
        # idx, tfhm_post = movfilter(t, tfhm, "post")
        psd_mean_tfhm_post[id+1, :] = mean(tfhm_post, dims=2)
        psd_std_tfhm_post[id+1, :] = std(tfhm_post, dims=2)

        # Compare Q and Q0 # DEPRECATED, we are not using surrogate analysis right now
        # pvals_Q[id+1,1] = pvalue(ApproximateTwoSampleKSTest(Q_pre[id+1,:], Q0_pre[id+1,:]))
        # pvals_Q[id+1,2] = pvalue(ApproximateTwoSampleKSTest(Q_post[id+1,:], Q0_post[id+1,:]))

        # band analysis 
        α_pre, γ_pre, total_pre = relative_power_from_segments(f,tfhm_pre)
        stats_band_pre[id+1,1:10] .=
                mean(α_pre),mean(γ_pre),mean(α_pre./total_pre),mean(γ_pre./total_pre), mean(total_pre),
                std(α_pre), std(γ_pre), std(α_pre./total_pre), std(γ_pre./total_pre), std(total_pre)

        α_post, γ_post, total_post = relative_power_from_segments(f,tfhm_post)
        stats_band_post[id+1,1:10] .=
                mean(α_post),mean(γ_post),mean(α_post./total_post),mean(γ_post./total_post), mean(total_post),
                std(α_post), std(γ_post), std(α_post./total_post), std(γ_post./total_post), std(total_post)

        pvals_α[id+1,1] = pvalue(ApproximatePermutationTest(α_pre, α_post, mean, 1000))
        pvals_γ[id+1,1] = pvalue(ApproximatePermutationTest(γ_pre, γ_post, mean, 1000))
        pvals_α[id+1,2] = pvalue(ApproximatePermutationTest(α_pre./total_pre, α_post./total_post, mean, 1000))
        pvals_γ[id+1,2] = pvalue(ApproximatePermutationTest(γ_pre./total_pre, γ_post./total_post, mean, 1000))

        # T-test. 
        for j in 1:l
            # pvals_uneq_tfhm[id+1, j] = pvalue(UnequalVarianceTTest(tfhm_pre[j, :], tfhm_post[j, :]))
            # pvals_tfhm[id+1, j] = pvalue(EqualVarianceTTest(tfhm_pre[j, :], tfhm_post[j, :]))# can be computed from the means
            pvals_perm_tfhm[id+1,j] = pvalue(ApproximatePermutationTest(tfhm_pre[j, :], tfhm_post[j, :], mean, 1000))
        end
    end

    foutname=string(subj)*"/"
    namebase = "../4_outputs/pipeline/"*foutname*data
    writedlm(namebase*"tfhm_mean_post.dat",psd_mean_tfhm_post," ")
    writedlm(namebase*"tfhm_mean_pre.dat",psd_mean_tfhm_pre," ")
    writedlm(namebase*"tfhm_std_post.dat",psd_std_tfhm_post," ")
    writedlm(namebase*"tfhm_std_pre.dat",psd_std_tfhm_pre," ")
    writedlm(namebase*"tfhm_diff.dat",psd_mean_tfhm_post-psd_mean_tfhm_pre," ")
    writedlm(namebase*"tfhm_pvals.dat",pvals_perm_tfhm," ")

    writedlm(namebase*"band_stats_pre.dat",stats_band_pre," ")
    writedlm(namebase*"band_stats_post.dat",stats_band_post," ")
    writedlm(namebase*"band_pvals_alpha.dat",pvals_α," ")
    writedlm(namebase*"band_pvals_gamma.dat",pvals_γ," ")

    writedlm(namebase*"Q_pre.dat",Q_pre," ")
    writedlm(namebase*"Q0_pre.dat",Q0_pre," ")
    writedlm(namebase*"Q_post.dat",Q_post," ")
    writedlm(namebase*"Q0_post.dat",Q0_post," ")

    # writedlm(namebase*"Q_pvals.dat",pvals_Q," ")
    print("\nDone!\n")
    # return stats_band_pre, stats_band_post, pvals_α, pvals_γ, Q_pre, Q0_pre, Q_post, Q0_post
    return 
end

function numberofsegments(fl;m=9000000,Δt=10.0)
    rate = 2500.0
    dt = 1.0 / rate
    T  = m/rate     # total length
    ns = Int(floor(T / Δt))    # number of segments
    t = collect(0.5*Δt:Δt:Δt*ns)

    mov_pre = readdlm("../1_data/mov_" * fl * ".dat", ' ')
    mov_indx = floor.(mov_pre ./ Δt) * Δt .+ 0.5*Δt
    indx = t .∈ [mov_indx]
    indx = findall(==(0), indx)

    return length(indx)
end


function tfhm_analysis(id, fl, flmov ; m=9000000, Δt=10.0,kpath="local")
    rate = 2500;
    dt = 0.0004
    t, y = read_channel(id, fl; t0=dt, tf=m*dt, m=m,kpath=kpath)
    NW = 1.0*Δt/(2.0)
    K = 8
    l = 2000;

    T  = floor(length(y)/rate);
    ns = Int(floor(T / Δt));
    y = y[1:Int(floor(m/ns)*ns)]    # length `m` might not coincide with multiple of Δt*dt

    # t, f, tfhm = timefreq(y);
    s = mt_spectrogram( y , 2500*10 , 0 ; fs=rate , nw=NW , ntapers=K );
    t = s.time; f = s.freq[1:l]; tfhm = 0.5*s.power[1:l,:];
    smean = mean(tfhm,dims=2)[:,1];
    q = [logspectral_dist(tfhm[:,i],smean,f) for i in 1:ns]
    tr, f_tfhm = movfilter(t,tfhm,flmov; Δt, q=q)

    fsmean = mean(f_tfhm,dims=2)[:,1];
    ns = length(tr)
    q0 = [logspectral_dist(f_tfhm[:,i],fsmean,f) for i in 1:ns]

    return f,f_tfhm,q,q0,tr
end



function logspectral_dist(s1,s2,f)
    # y = @.  10*log10(s1/s2)^2
    # integrate(f,y)^0.5 # maybe divided by f[end]^0.5
    y = @. log10(s1/s2)^2
    mean(y)^0.5 # maybe divided by f[end]^0.5
end



# skip: number of rows to skip; nrow = number of channels per row
ypos(y,skip,nrow,dp,dist) = -0.413*dp + (skip + floor((y-1)/nrow))*dist

function depth(type,n)
    dp = 3840
    if type=="cortex"
        return ypos.(1:n, 0, 2, dp, 20)
    elseif type=="bipolar"
        return ypos.(1:n, 2, 2, dp, 20)
    elseif type=="csd"
        return ypos.(1:n, 1, 1, dp, 20)
    elseif type=="kcsd"
        return ypos.(1:n, 0, 1, dp, 10)
    else
        print("Unkown data type "*type*"\n")
    end
end

function relative_power_from_mean(fhm,freqs)
    df = freqs[2]-freqs[1]
    n = size(fhm)[2]
    αband = findall(@. freqs>4.0 && freqs<22)
    γband = findall(@. freqs>32.0 && freqs<48.0)
    α = zeros(n)
    γ = zeros(n)
    for i in 1:n
        α[i] = mean(fhm[αband,i])*df # take the sum instead of the mean to get the area
        γ[i] = mean(fhm[γband,i])*df
    end

    return α,γ
end

function relative_power_from_segments(freqs,tfhm)
    df = freqs[2]-freqs[1]
    αband = findall(@. freqs>4.0 && freqs<22)
    γband = findall(@. freqs>32.0 && freqs<48.0)
    α = sum(tfhm[αband,:],dims=1)[1,:]*df 
    γ = sum(tfhm[γband,:],dims=1)[1,:]*df 
    tot = sum(tfhm[:,:],dims=1)[1,:]*df 
    return α,γ,tot
end


function load_precomputed()
    path="/home/pclusella/Documents/Data/UPO-tACs/7_results/cortex_heatmaps/"
    condition = "pre"
    lfp_pre = readdlm(path*"psd_mean_tfhm_cortex_"*condition*".dat")
    blfp_pre = readdlm(path*"psd_mean_tfhm_bipolar_"*condition*".dat")
    csd_pre = readdlm(path*"psd_mean_tfhm_csd_"*condition*".dat")
    kcsd_pre = readdlm(path*"psd_mean_tfhm_kcsd_"*condition*".dat")

    condition = "post"
    lfp_post = readdlm(path*"psd_mean_tfhm_cortex_"*condition*".dat")
    blfp_post = readdlm(path*"psd_mean_tfhm_bipolar_"*condition*".dat")
    csd_post = readdlm(path*"psd_mean_tfhm_csd_"*condition*".dat")
    kcsd_post = readdlm(path*"psd_mean_tfhm_kcsd_"*condition*".dat")

    return [lfp_pre,blfp_pre,csd_pre,kcsd_pre,lfp_post,blfp_post,csd_post,kcsd_post]
end

# not exported
function load_precomputed_std()
    path="/home/pclusella/Documents/Data/UPO-tACs/7_results/cortex_heatmaps/"
    condition = "pre"
    lfp_pre = readdlm(path*"psd_std_tfhm_cortex_"*condition*".dat")
    blfp_pre = readdlm(path*"psd_std_tfhm_bipolar_"*condition*".dat")
    csd_pre = readdlm(path*"psd_std_tfhm_csd_"*condition*".dat")
    kcsd_pre = readdlm(path*"psd_std_tfhm_kcsd_"*condition*".dat")

    condition = "post"
    lfp_post = readdlm(path*"psd_std_tfhm_cortex_"*condition*".dat")
    blfp_post = readdlm(path*"psd_std_tfhm_bipolar_"*condition*".dat")
    csd_post = readdlm(path*"psd_std_tfhm_csd_"*condition*".dat")
    kcsd_post = readdlm(path*"psd_std_tfhm_kcsd_"*condition*".dat")

    return [lfp_pre,blfp_pre,csd_pre,kcsd_pre,lfp_post,blfp_post,csd_post,kcsd_post]
end

end # module
