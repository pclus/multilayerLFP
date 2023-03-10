module NeuropixelAnalysis

using DelimitedFiles, Multitaper, Plots, DSP,Statistics;

export DelimitedFiles, Multitaper, Plots, DSP,Statistics
export read_channel,channel_idx,heatmapMT,timefreq,movfilter,heatmap_segments, process_data,relative_power,depth

"""
    read_channel(id, t0, tf, fl)

Obtain the time series indexed by `id` between times `t0` and `tf` 
from the binary file identified by `fl`.

# Examples
```julia-repl
julia> v158=read_channel(158,0.0004,900,"filtered_pre");
```
"""
function read_channel(id, t0, tf, fl) # Equivalent to the read_binary.c code
    id = Int(id)
    n = 384
    m = 2250000
    dt = 1.0 / 2500.0


    data = zeros(m)
    fin = open("../1_data/" * fl * ".bin", "r")
    for i in 0:id
        read!(fin, data)
    end
    close(fin)

    m0 = Int(ceil(t0 / dt))
    mf = Int(floor(tf / dt))
    dm = mf - m0 + 1

    channel = data[m0:mf]
    t = collect(dt * (m0:mf))

    # writedlm("../temp.dat",[t channel]," ");

    return t, channel
end

"""
    cut_cortex(fl,channels)

Read the channels listed in the array `channels` from the "fl" binary file,
filter them, and write in a 'cortex_'*fl binary file.
Already bandpass_filtered.
"""
function cut_cortex(fl,n0,nf) 

    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = 900.0

    bpfilter = digitalfilter(Bandpass(1.0, 300.0; fs=2500), Butterworth(3))
    fout = open("../1_data/cortex_" * fl * ".bin", "w")

    for id in n0:nf
        t, ch = read_channel(id, t0, tf, fl)
        fil_ch = filtfilt(bpfilter, ch)
        write(fout, fil_ch)
    end
    close(fout)
end


"""
    bandpass_filter(fl)

Bandpass filter all data from binary file referenced by `fl` between 1 and 300Hz.
Uses a Butterworth filter of order 3 and it stores it in a new binary file 
named `../1_data/filtered_...`.
"""
function bandpass_filter(fl) # possibly deprecated
    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = 900.0
    n = 384

    bpfilter = digitalfilter(Bandpass(1.0, 300.0; fs=rate), Butterworth(3))
    # bsfilter  = digitalfilter(Bandstop(49.95,50.05; fs=rate),Butterworth(1))
    # bsfilter2 = digitalfilter(Bandstop(59.95,60.05; fs=rate),Butterworth(1))

    fout = open("../1_data/filtered_" * fl * ".bin", "w")

    for id in 0:n-1
        t, ch = read_channel(id, t0, tf, fl)
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
SHOULD BE CHANGED TO TAKE next-current INSTEAD OF next-previous
"""
function compute_bipolar(flin,flout,n)
    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = 900.0

    fout = open("../1_data/"* flout * ".bin", "w")

    t, prev_a = read_channel(0, t0, tf, flin)
    t, prev_b = read_channel(1, t0, tf, flin)
    t, prev_c = read_channel(2, t0, tf, flin)
    t, prev_d = read_channel(3, t0, tf, flin)

    inv_dy = 1 / 20.0 # 1/μm

    nf = n - n%4 
    for id in 4:4:nf-1
        t, curr_a = read_channel(id    , t0, tf, flin)
        t, curr_b = read_channel(id + 1, t0, tf, flin)
        t, curr_c = read_channel(id + 2, t0, tf, flin)
        t, curr_d = read_channel(id + 3, t0, tf, flin)

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
    compute_csd(fl)

Store in a binary file the CSD obtained by computing the Laplacian 
using a 2-dimensional diagonal finite difference stencil.
"""
function compute_csd(flin,flout,n)
    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = 900.0

    fout = open("../1_data/" * flout * ".bin", "w")

    t, ch_5 = read_channel(0, t0, tf, flin)
    t, ch_6 = read_channel(1, t0, tf, flin)
    t, ch_7 = read_channel(2, t0, tf, flin)
    t, ch_8 = read_channel(3, t0, tf, flin)

    inv_dy2 = 1 / 25.61^2 # 1/μm
    σ = 1.0 # conductivity S/μm or 1/ (Ω μm)

    nf = n - n%4 
    for id in 4:4:nf-1

        ch_1 = ch_5
        ch_2 = ch_6
        ch_3 = ch_7
        ch_4 = ch_8

        t, ch_5 = read_channel(id + 0, t0, tf, flin)
        t, ch_6 = read_channel(id + 1, t0, tf, flin)
        t, ch_7 = read_channel(id + 2, t0, tf, flin)
        t, ch_8 = read_channel(id + 3, t0, tf, flin)

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
    heatmapMT(t0, tf, fl, channels)
    
Computes the spectral heatmap with Multitaper PSD for specific channels and a time frame
"""
function heatmapMT(t0, tf, fl, channels)

    rate = 2500.0
    dt = 1.0 / rate
    n = size(channels, 1)
    m = 250001
    NW = 1.0 * m * dt / (2.0)  #bandwith is W*dt/2.0, and NW = N*W*dt/2.0
    K = 8    # number of tappers (should be less than 2*NW)

    l = 30002 # WARNING: This assumes relevant frequencies below entry l of the spectra
    hm = zeros(n, l)

    local S
    Threads.@threads for id in channels
        t, chdat = read_channel(id, t0, tf, fl)
        S = multispec(chdat, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=true)
        hm[id+1, :] = S.S[1:l]
        if id % 5 == 0
            print(id, "\r")
            flush(stdout)
        end
    end
    print(stdout, "\n")
    freqs = collect(S.f[1:l])

    freqs, hm
end

"""
    timefreq(id, fl)

Computes the time-frequency heatmap with MT-PSD using segments of 10s
"""
function timefreq(id, fl)

    rate = 2500.0
    dt = 1.0 / rate
    n = 384
    Δt = 10.0      # segment duration (10 seconds)
    ns = Int(900.0 / Δt)  # number of segments
    m = Int(2250000 / ns) # segments of 10 seconds
    # mh = m/2;       # length of the MT spectra...
    l = 2000       # but...up to l is enough to get the relevant freqs
    NW = 1.0 * m * dt / (2.0)
    K = 8
    t, chdat = read_channel(id, 4e-4, 900, fl)    # time series starts at 4e-4

    #--------------------------------------------------------------
    # Non-overlapping windows <default>
    #--------------------------------------------------------------
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
    #--------------------------------------------------------------
    # Overlapping windows <just for display>
    #--------------------------------------------------------------
    # tfhm = zeros(l,900);   
    # sec = 2250000/rate# for every second
    # local S;    # so that S exists outside the loop
    # for s in 5:895
    #     segdat = chdat[Int(1+(s-5)*rate):Int((s+5)*rate)]
    #     S  = multispec(segdat, dt=dt, NW=NW, K=K);
    #     tfhm[:,s] = S.S[1:l];
    #     if s%5==0
    #         print(s,"\r");
    #         flush(stdout);
    #     end
    # end
    #--------------------------------------------------------------

    freqs = collect(S.f[1:l])
    times = collect(5:10:900)

    times, freqs, tfhm
end

"""
    movfilter(t, tfhm, fl)
Authomatic removal of segments including movement data.
This function assumes segments of Δt=10
"""
function movfilter(t, tfhm, fl)
    mov_pre = readdlm("../1_data/mov_" * fl * ".dat", ' ')
    Δt = 10.0
    mov_indx = Int.(floor.(mov_pre ./ Δt) * Δt .+ 5)
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
function process_data(n0,nf)
    n=length(n0:nf)
    cut_cortex("pre",n0,nf)
    cut_cortex("post",n0,nf)
    compute_bipolar("cortex_pre","bipolar_pre",n)   
    compute_bipolar("cortex_post","bipolar_post",n) 
    compute_csd("cortex_pre","csd_pre",n)           
    compute_csd("cortex_post","csd_post",n)         
end

"""
    heatmap_segments(fl)

Create the segments of 10s, compute their MT spectra and provide their mean and standard deviation.
This is done for all channels, removing the time segments that contain movement.
"""
function heatmap_segments(fl,n)
    # n = 384;
    l = 2000;

    flmv=findall("_pre",fl)
    if isempty(flmv)
        flmv="post"
    else
        flmv="pre"
    end

    psd_mean_tfhm = zeros(n,l);
    psd_std_tfhm  = zeros(n,l);

    state = Threads.Atomic{Int}(0);
    
    Threads.@threads for id in 0:n-1
        Threads.atomic_add!(state, 1)
        print("--> ",state[]," out of ",n,"\n");
        flush(stdout);
        t,f,tfhm = timefreq(id,fl);
        f_idx,f_tfhm = movfilter(t,tfhm,flmv);
        psd_mean_tfhm[id+1,:] = mean(f_tfhm,dims=2);  
        psd_std_tfhm[id+1,:] = std(f_tfhm,dims=2);
    end
    return psd_mean_tfhm,psd_std_tfhm
end

# """
#     ttest_eqvar(fmean1,fmean2,fstd1,fstd2)

# T test of equal variances using pre-saved heatmaps (should have same size)
# """
# function ttest_eqvar(fmean1,fmean2,fstd1,fstd2)
#     A=dlmread(fmean1)
#     B=dlmread(fmean2)
#     sA=dlmread(fstd1)
#     sB=dlmread(fstd2)
#     nx= ;
#     ny= ;
#     test = EqualVarianceTTest.(nx,ny,A,B,sA,sB)
# end

"""
    prepost_analysis(n0,nf)

Compute the mean spectra and perform comparison between pre and post for all datatypes.
"""
function prepost_analysis(n0,nf)
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
end


# skip: number of rows to skip; nrow = number of channels per row
ypos(y,skip,nrow,dp,dist) = -0.413*dp + (skip + floor((y-1)/nrow))*dist

function depth(type,n)
    dp = 3840
    if type=="lfp"
        return ypos.(1:n, 0, 2, dp, 20)
    elseif type=="blfp"
        return ypos.(1:n, 2, 2, dp, 20)
    elseif type=="csd"
        return ypos.(1:n, 1, 1, dp, 20)
    elseif type=="kcsd"
        return ypos.(1:n, 0, 1, dp, 10)
    else
        print("Unkown data type "*type*"\n")
    end
end

function relative_power(fhm,freqs)
    df = freqs[2]-freqs[1]
    n = size(fhm)[2]
    αband = findall(@. freqs>4.0 && freqs<22)
    γband = findall(@. freqs>35.0 && freqs<58)
    α = zeros(n)
    γ = zeros(n)
    for i in 1:n
        α[i] = mean(fhm[αband,i])*df # take the sum instead of the mean to get the area
        γ[i] = mean(fhm[γband,i])*df
    end

    return α,γ
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