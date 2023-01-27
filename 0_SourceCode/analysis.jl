#!/usr/bin/env -S julia --threads 16
using DelimitedFiles,Multitaper,Plots;

# -------------------------------------------------------------
# -------------------------------------------------------------
# FUNCTIONS:
# -------------------------------------------------------------
# Reads the binary data for a specific channel, time window, and file
# -------------------------------------------------------------
function read_channel(id,t0,tf,fl) # Equivalent to the read_binary.c code
    id=Int(id);
    n = 384; 
    m = 2250000;
    dt=1.0/2500.0;

    
    data = zeros(m);
    fin=open("../1_Raw/"*fl*".bin","r");
    for i in 1:id
        read!(fin,data)
    end
    fclose(fin);

    m0=Int(ceil(t0/dt));
    mf=Int(floor(tf/dt));
    dm=mf-m0+1

    channel=data[m0:mf]
    t=collect(dt*(m0:mf));

    # writedlm("../temp.dat",[t channel]," ");
    
    return t, channel;
end

# -------------------------------------------------------------
# Creates the CSD binary file by computing the Laplacian [UNDER DEVELOPMENT]
# -------------------------------------------------------------
function csd(fl)
    rate = 2500.0;
    dt = 1/rate;
    t0 = dt;
    tf = 900.0; 
    n = 384
    read_channel(id,t0,tf,fl)

    column = 1; #odd
    column = 2; #even
    for i in 1:n

    end


end

# -------------------------------------------------------------
# Creates the bipolar potential binary file by computing a 
# 2-sided derivative for each of the 4 columns of Neuropixel
# -------------------------------------------------------------
function bipolar(fl)
    rate = 2500.0;
    dt = 1/rate;
    t0 = dt;
    tf = 900.0; 
    n = 384;

    fout = open("../1_Raw/bipolar_"*fl*".bin","a");

    t,prev_a = read_channel(1,t0,tf,fl)
    t,prev_b = read_channel(2,t0,tf,fl)
    t,prev_c = read_channel(3,t0,tf,fl)
    t,prev_d = read_channel(4,t0,tf,fl)

    t,curr_a = read_channel(5,t0,tf,fl)
    t,curr_b = read_channel(6,t0,tf,fl)
    t,curr_c = read_channel(7,t0,tf,fl)
    t,curr_d = read_channel(8,t0,tf,fl)

    inv_dy=1/20.0; # 1/μm

    for id in 9:4:n
        t,next_a = read_channel(id  ,t0,tf,fl)
        t,next_b = read_channel(id+1,t0,tf,fl)
        t,next_c = read_channel(id+2,t0,tf,fl)
        t,next_d = read_channel(id+3,t0,tf,fl)

        bip_a = 0.5*inv_dy*( prev_a - next_a);
        bip_b = 0.5*inv_dy*( prev_b - next_b);
        bip_c = 0.5*inv_dy*( prev_c - next_c);
        bip_d = 0.5*inv_dy*( prev_d - next_d);

        prev_a=curr_a;
        prev_b=curr_b;
        prev_c=curr_c;
        prev_d=curr_d;

        curr_a=next_a;
        curr_b=next_b;
        curr_c=next_c;
        curr_d=next_d;

        write(fout,bip_a);
        write(fout,bip_b);
        write(fout,bip_c);
        write(fout,bip_d);

    end
    close(fout)

end


#--------------------------------------------------------------
# Spectral heatmap with Multitaper PSD for all channels and a time frame
#--------------------------------------------------------------
function heatmapMT(t0,tf,fl)
    
    rate = 2500.0;
    dt = 1.0/rate;
    n = 384;
    m = 250001;
    NW = 1.0*m*dt/(2.0) ;  #bandwith is W*dt/2.0, and NW = N*W*dt/2.0
    K = 8;    # number of tappers (should be less than 2*NW)


    l = 30002; # WARNING: This assumes relevant frequencies below entry l of the spectra
    hm = zeros(n,l);

    local S;  
    Threads.@threads for id in 1:n
        t,chdat=read_channel(id,t0,tf,fl);
        S  = multispec(chdat, dt=dt, NW=NW, K=K,jk=true,Ftest=true, a_weight=true);
        hm[id,:] = S.S[1:l] 
        if id%5==0
            print(id,"\r");
            flush(stdout);
        end
    end
    print(stdout,"\n");
    freqs=collect(S.f[1:l]);

    freqs,hm;
end

#--------------------------------------------------------------
# Time-frequency heatmap with MT-PSD using segments of 10s
#--------------------------------------------------------------
function timefreq(id,fl)

    rate = 2500.0;
    dt = 1.0/rate;
    n = 384;
    Δt = 10.0;      # segment duration (10 seconds)
    ns = Int(900.0/Δt);  # number of segments
    m = Int(2250000/ns); # segments of 10 seconds
    # mh = m/2;       # length of the MT spectra...
    l = 2000;       # but...up to l is enough to get the relevant freqs
    NW = 1.0*m*dt/(2.0) ;  
    K = 8;
    t,chdat = read_channel(id,4e-4,900,fl);    # time series starts at 4e-4

    #--------------------------------------------------------------
    # Non-overlapping windows <default>
    #--------------------------------------------------------------
    tfhm = zeros(l,ns);   # time-freq heatmap
    local S;    # so that S exists outside the loop
    for s in 1:ns
        segdat = chdat[((s-1)*m+1):s*m]
        S  = multispec(segdat, dt=dt, NW=NW, K=K);
        tfhm[:,s] = S.S[1:l];
        if s%5==0
            print(s,"\r");
            flush(stdout);
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
    
    freqs = collect(S.f[1:l]);
    times = collect(5:10:900);

    times,freqs,tfhm;
end

#--------------------------------------------------------------
# Authomatic removal of segments including (pre-) movement data
# notice that this function assumes Δt=10
#--------------------------------------------------------------
function movfilter(t,tfhm,fl)
    mov_pre=readdlm("../1_Raw/mov_"*fl*".dat",' ');
    Δt=10.0;
    mov_indx=Int.(floor.(mov_pre./Δt)*Δt.+5);
    indx = t.∈ [mov_indx];
    indx = findall(==(0),indx)
    f_tfhm = tfhm[:,indx];
    indx,f_tfhm
end
#--------------------------------------------------------------
# -------------------------------------------------------------

# --------------------------------------------------------------
# --------------------------------------------------------------
# WORKFLOW:
# --------------------------------------------------------------
# --------------------------------------------------------------

# Load time series for a specific channel from second 200 to 300:
t0=200;
tf=210;
id=5;
fl="pre"

t,chdat=read_channel(id,t0,tf,fl);

# Comapring bipolar and raw data:
#--------------------------------------------------------------
# Load time series for a specific channel from second 200 to 300:
t0=1/2500.0;
tf=300;

fl="pre"
t,ch1=read_channel(1,t0,tf,fl);
t,ch5=read_channel(5,t0,tf,fl);
t,ch9=read_channel(9,t0,tf,fl);

fl="bipolar_pre"
t,bi5=read_channel(1,t0,tf,fl);


plot(ch1[1:1000])
plot!(ch5[1:1000])
# plot!(ch9[1:1000])
plot!(data[1:1000])
plot!((ch1[1:1000]-ch9[1:1000]))


# Compute PSD using Multitaper PSD ----------------------------
rate=2500.0;
dt=1.0/rate;

NW = 1.0*length(chdat)*dt/(2.0) ;  # real bandwith is ω*dt/2.0, and NW = N*ω*dt/2.0
K  = 10;    # number of tappers (should be similar slightly less than 2*NW)
@time S  = multispec(chdat, dt=dt, NW=NW, K=K,jk=true,Ftest=true, a_weight=false);

plot(S,xlim=(0,2000),ylim=(1e-18,5e-15))


p1=plot(S.f,S.S,xlim=(0,200),ylim=(1e-18,1e-15),lw=1.0,yaxis=:log)
# Ftest pvalues
# plot(S.f,S.Fpval,ylim=(0.0,1e-2),lt=:scatter,xlim=(0,100))



# # Compute the PSD using standard FFT:--------------------------
# using FFTW
# rate=2500.0;
# dt=1.0/rate;
# Fr = rfft(chdat); 
# wr = rfftfreq(length(t), rate); 
# plot(wr, abs.(Fr),xlim=(0,100)) # Logscale with yaxis=:log

# plot(wr,abs2.(Fr),xlim=(0,100),ylim=(0,1e-7),lw=0.5,yaxis=:log)
# plot!(S.f,S.S.*length(S.f)/dt,yaxis=:log,lw=2,ylim=(1e-10,5e-7))
# #--------------------------------------------------------------




# --- some notes about  multitaper and standard fft------------
#
# plot(S,xlim=(0,200),ylim=(1e-20,5e-16)) # authomatic MT plot
#
# Comparison between multitaper and standard fft:
# Notice that multitaper package normalizes by dt/N, similar to what I am doing in C.
# plot(S.f,S.S.*length(S.f)/dt)
# plot!(wr,abs2.(Fr),xlim=(0,100),ylim=(0,1e-7),lw=0.1,yaxis=:log)
#
# Plot including 95% confidence intervals, as in the authomatic version...
# does not look right...it seems that the authomatic recipe plots the confidence
# for the log-spectra, and not the true spectra
# using StatsFuns
# z=norminvcdf(0,1,0.975);
# plot(S.f,S.S,xlim=(0,200),ylim=(0,2e-16))
# plot!(S.f,S.S.*exp.(-z*sqrt.(S.jkvar)),fillrange=S.S.*exp.(2*z*sqrt.(S.jkvar)),c=1,alpha=0.35)
# plot!(yaxis=:log,ylim=(1e-20,5e-16))
#
# precompute tappers # or does not work, or it does not improve performance
# @time dpss=dpss_tapers(length(chdat),NW,K);
# @time Z  = multispec(chdat, dt=dt, NW=NW, K=K,jk=true,Ftest=true, a_weight=true, dpVec=dpss);
# plot!(Z.f,Z.S,xlim=(0,200),ylim=(1e-18,1e-15),lw=0.2,yaxis=:log)
# -------------------------------------------------------------

# -------------------------------------------------------------
# Create a heatmap of the pre data from 200s to 300s using Multitaper:
a,b=heatmapMT(200.0,300.0,fl)
# writedlm("../4_outputs/psd_mthm.dat",b," ");
# -------------------------------------------------------------

# -------------------------------------------------------------
# Time-frequency analysis for a single channel:
fl="pre"
id=150
t,f,tfhm = timefreq(id,fl);
writedlm("../4_outputs/ch"*string(id)*"_tfhm_HD.dat",tfhm," ");

# Remove the segments containing movement:
f_idx,f_tfhm = movfilter(t,tfhm,fl);

# Compute statistics for this channel:
using Statistics, HypothesisTests
mean_psd=mean(f_tfhm,dims=2);
std_psd =std(f_tfhm,dims=2);
plot(f,mean_psd,ribbon=std_psd,c=1,fillalpha=0.25)

# -------------------------------------------------------------

# -------------------------------------------------------------
# Segmentation analysis for all channels and create the heatmaps for the averages
#using Statistics
#n = 384;
#l = 2000;
#
#psd_mean_tfhm = zeros(n,l);
#psd_std_tfhm  = zeros(n,l);
#
#state = Threads.Atomic{Int}(0);
#
#Threads.@threads for id in 1:n
#    Threads.atomic_add!(state, 1)
#    print("--> ",state[]," out of ",n,"\n");
#    flush(stdout);
#    t,f,tfhm = timefreq(id,fl);
#    f_idx,f_tfhm = movfilter(t,tfhm,fl);
#    psd_mean_tfhm[id,:] = mean(f_tfhm,dims=2);  
#    psd_std_tfhm[id,:] = std(f_tfhm,dims=2);
#end
#
#writedlm("../4_outputs/psd_mean_tfhm_"*fl*".dat",psd_mean_tfhm',' ');
#writedlm("../4_outputs/psd_std_tfhm_"*fl*".dat",psd_std_tfhm',' ');
#
#heatmap(0.1:0.1:200,1:n,log.(psd_mean_tfhm))
# -------------------------------------------------------------


# -------------------------------------------------------------
# Check for differences between pre and post using t-test of the
# different time windows
using Statistics,HypothesisTests

n = 384;
l = 2000;

psd_mean_tfhm_pre = zeros(n,l);
psd_std_tfhm_pre  = zeros(n,l);
psd_mean_tfhm_post = zeros(n,l);
psd_std_tfhm_post  = zeros(n,l);
pvals_tfhm = zeros(n,l);
pvals_uneq_tfhm = zeros(n,l);

state = Threads.Atomic{Int}(0);

Threads.@threads for id in 1:n

    # Prompt state
    Threads.atomic_add!(state, 1)              
    print("--> ",state[]," out of ",n,"\n");
    flush(stdout);

    # Pre
    t,f,tfhm = timefreq(id,"pre");
    idx,f_pre = movfilter(t,tfhm,"pre");
    psd_mean_tfhm_pre[id,:] = mean(f_pre,dims=2);  
    psd_std_tfhm_pre[id,:] = std(f_pre,dims=2);

    # Post
    t,f,tfhm = timefreq(id,"post");
    idx,f_post = movfilter(t,tfhm,"post");
    psd_mean_tfhm_post[id,:] = mean(f_post,dims=2);  
    psd_std_tfhm_post[id,:] = std(f_post,dims=2);

    # T-test. 
    for j in 1:l
        pvals_uneq_tfhm[id,j]=pvalue(UnequalVarianceTTest(f_pre[j,:],f_post[j,:]));
        pvals_tfhm[id,j]=pvalue(EqualVarianceTTest(f_pre[j,:],f_post[j,:]));# can be computed from the means
    end
end

writedlm("../4_outputs/psd_mean_tfhm_post.dat",psd_mean_tfhm_pre',' ');
writedlm("../4_outputs/psd_std_tfhm_post.dat",psd_std_tfhm_pre',' ');
writedlm("../4_outputs/psd_mean_tfhm_pre.dat",psd_mean_tfhm_post',' ');
writedlm("../4_outputs/psd_std_tfhm_pre.dat",psd_std_tfhm_post',' ');
writedlm("../4_outputs/psd_pvals.dat",pvals_tfhm',' ');
writedlm("../4_outputs/psd_pvals_uneq.dat",pvals_uneq_tfhm',' ');
writedlm("../4_outputs/psd_mean_tfhm_difference.dat",(psd_mean_tfhm_post-psd_mean_tfhm_pre)',' ');

# heatmap(0.1:0.1:200,1:n,log.(psd_mean_tfhm))

#plot(f,psd_mean_tfhm_pre[id,:],ribbon=psd_std_tfhm_pre[id,:],c=1,fillalpha=0.25)
#plot!(f,psd_mean_tfhm_post[id,:],ribbon=psd_std_tfhm_post[id,:],c=2,fillalpha=0.25)
#plot!(f,(pvals_tfhm[id,:].<=1e-5).*1e-17,ylim=(1e-18,1.2e-15),lt=:scatter)

# t-test comparison
# a=randn(100);
# b=10*randn(100).+2.2;
# p=pvalue(UnequalVarianceTTest(a,b))
# p=pvalue(EqualVarianceTTest(a,b))
