module NeuropixelAnalysis

using DelimitedFiles, Multitaper, Plots, DSP,Statistics;

export DelimitedFiles, Multitaper, Plots, DSP,Statistics
export read_channel,channel_idx,heatmapMT,timefreq,movfilter,heatmap_segments

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
    for i in 1:id
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
    bandpass_filter(fl)

Bandpass filter all data from binary file referenced by `fl` between 1 and 300Hz.
Uses a Butterworth filter of order 3 and it stores it in a new binary file 
named `../1_data/filtered_...`.
"""
function bandpass_filter(fl)
    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = 900.0
    n = 384

    bpfilter = digitalfilter(Bandpass(1.0, 300.0; fs=2500), Butterworth(3))
    fout = open("../1_data/filtered_" * fl * ".bin", "w")

    for id in 1:384
        t, ch = read_channel(id, t0, tf, fl)
        fil_ch = filtfilt(bpfilter, ch)
        write(fout, fil_ch)
    end
    close(fout)

end

"""
    compute_bipolar(fl)

Store in a binary file the bipolar data obtained by computing the first 
derivative of the data using a symmetric stencil in along each neuropixel column.
"""
function compute_bipolar(fl)
    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = 900.0
    n = 384

    fout = open("../1_data/bipolar_" * fl * ".bin", "w")

    fl = "filtered_" * fl
    t, prev_a = read_channel(1, t0, tf, fl)
    t, prev_b = read_channel(2, t0, tf, fl)
    t, prev_c = read_channel(3, t0, tf, fl)
    t, prev_d = read_channel(4, t0, tf, fl)

    t, curr_a = read_channel(5, t0, tf, fl)
    t, curr_b = read_channel(6, t0, tf, fl)
    t, curr_c = read_channel(7, t0, tf, fl)
    t, curr_d = read_channel(8, t0, tf, fl)

    inv_dy = 1 / 20.0 # 1/μm

    for id in 9:4:n
        t, next_a = read_channel(id, t0, tf, fl)
        t, next_b = read_channel(id + 1, t0, tf, fl)
        t, next_c = read_channel(id + 2, t0, tf, fl)
        t, next_d = read_channel(id + 3, t0, tf, fl)

        bip_a = 0.5 * inv_dy * (prev_a - next_a)
        bip_b = 0.5 * inv_dy * (prev_b - next_b)
        bip_c = 0.5 * inv_dy * (prev_c - next_c)
        bip_d = 0.5 * inv_dy * (prev_d - next_d)

        prev_a = curr_a
        prev_b = curr_b
        prev_c = curr_c
        prev_d = curr_d

        curr_a = next_a
        curr_b = next_b
        curr_c = next_c
        curr_d = next_d

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
function compute_csd(fl)
    rate = 2500.0
    dt = 1 / rate
    t0 = dt
    tf = 900.0
    n = 384

    fout = open("../1_data/csd_" * fl * ".bin", "w")
    fl = "filtered_" * fl

    t, ch_5 = read_channel(1, t0, tf, fl)
    t, ch_6 = read_channel(2, t0, tf, fl)
    t, ch_7 = read_channel(3, t0, tf, fl)
    t, ch_8 = read_channel(4, t0, tf, fl)

    inv_dy2 = 1 / 25.61^2 # 1/μm
    σ = 1.0 # conductivity S/μm or 1/ (Ω μm)

    for id in 5:4:n

        ch_1 = ch_5
        ch_2 = ch_6
        ch_3 = ch_7
        ch_4 = ch_8

        t, ch_5 = read_channel(id + 0, t0, tf, fl)
        t, ch_6 = read_channel(id + 1, t0, tf, fl)
        t, ch_7 = read_channel(id + 2, t0, tf, fl)
        t, ch_8 = read_channel(id + 3, t0, tf, fl)

        csd_right = -σ * inv_dy2 * (-4.0 * ch_4 + ch_1 + ch_2 + ch_5 + ch_6)
        csd_left = -σ * inv_dy2 * (-4.0 * ch_5 + ch_3 + ch_4 + ch_7 + ch_8)

        write(fout, csd_right)
        write(fout, csd_left)

    end
    close(fout)

end

"""
    channel_idx(chid)

Given the id of a channel, return its index for the csd, if it exists
"""
function channel_idx(chid)
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
        hm[id, :] = S.S[1:l]
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
function process_data()
    bandpass_filter("pre")
    bandpass_filter("post")
    compute_bipolar("pre")
    compute_bipolar("post")
    compute_csd("pre")
    compute_csd("post")
end

"""
    heatmap_segments(fl)

Create the segments of 10s, compute their MT spectra and provide their mean and standard deviation.
This is done for all channels, removing the time segments that contain movement.
"""
function heatmap_segments(fl,n)
    # n = 384;
    l = 2000;

    flmv=findall("pre",fl)
    if isempty(flmv)
        flmv="post"
    else
        flmv="pre"
    end

    psd_mean_tfhm = zeros(n,l);
    psd_std_tfhm  = zeros(n,l);

    state = Threads.Atomic{Int}(0);
    
    Threads.@threads for id in 1:n
        Threads.atomic_add!(state, 1)
        print("--> ",state[]," out of ",n,"\n");
        flush(stdout);
        t,f,tfhm = timefreq(id,fl);
        f_idx,f_tfhm = movfilter(t,tfhm,flmv);
        psd_mean_tfhm[id,:] = mean(f_tfhm,dims=2);  
        psd_std_tfhm[id,:] = std(f_tfhm,dims=2);
    end
    return psd_mean_tfhm,psd_std_tfhm
end

end # module