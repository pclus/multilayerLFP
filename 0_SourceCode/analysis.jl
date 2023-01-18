#!/usr/bin/env -S julia --threads 16
using DelimitedFiles,Multitaper,Plots;

# -------------------------------------------------------------
# -------------------------------------------------------------
# FUNCTIONS:
# -------------------------------------------------------------
# Reads the "pre" data for a specific channel and time window:
# -------------------------------------------------------------
function read_channel(id,t0,tf,fl) # Equivalent to the read_binary.c code
    id=Int(id);
    # fin=open("../1_Raw/pre.bin","r");
    fin=open("../1_Raw/"*fl*".bin","r");
    n = 384; 
    m = 2250000;
    dt=1.0/2500.0;

    m0=Int(ceil(t0/dt));
    mf=Int(floor(tf/dt));

    data=zeros(n);

    for i in 1:(m0-1)
        read!(fin,data);
    end
    dm=mf-m0+1;
    channel=zeros(dm);
    for i in 1:dm
        if !eof(fin)
            read!(fin,data);
            channel[i]=data[id];
        else
            print("End-of-file at ",i," out of ",dm);
        end
    end
    close(fin);

    t=collect(dt*(m0:mf));
    writedlm("../temp.dat",[t channel]," ",);
    return t, channel;
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

    # This should be parallelized
    for id in 1:n
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
    tfhm = zeros(l,ns);                     # time-freq heatmap
    
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
tf=300;
id=150;
fl="post"

t,chdat=read_channel(id,t0,tf,fl);


# Compute PSD using Multitaper PSD ----------------------------
rate=2500.0;
dt=1.0/rate;

NW = 1.0*length(chdat)*dt/(2.0) ;  # real bandwith is ω*dt/2.0, and NW = N*ω*dt/2.0
K  = 40;    # number of tappers (should be similar slightly less than 2*NW)
@time S  = multispec(chdat, dt=dt, NW=NW, K=K,jk=true,Ftest=true, a_weight=false);

plot(S,xlim=(0,200),ylim=(1e-18,5e-15))


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
# a,b=heatmapMT(200.0,300.0)
# writedlm("../4_outputs/psd_mthm.dat",b," ");
# -------------------------------------------------------------

# -------------------------------------------------------------
# Time-frequency analysis for a single channel:
fl="pre"
id=150
t,f,tfhm = timefreq(id,fl);
# writedlm("../4_outputs/ch"*string(id)*"_tfhm.dat",tfhm," ");

# Remove the segments containing movement:
f_idx,f_tfhm = movfilter(t,tfhm,fl);

# Compute statistics for this channel:
using Statistics
mean_psd=mean(f_tfhm,dims=2);
std_psd =std(f_tfhm,dims=2);
plot(f,mean_psd,ribbon=std_psd,c=1,fillalpha=0.25)

# Classify segments using kmeans
using Clustering
km = kmeans(f_tfhm,2,tol=1e-50)
ikm = assignments(km)
ckm = km.centers;
plot(f_tfhm[:,:],c=ikm[:]',lw=0.5,legend=false,ylim=(5e-18,5e-15),alpha=0.2)
plot!(ckm[:,:],c=[1 2],lw=3)

# T-test between the two groups
m=size(f_tfhm,1);
pvals=zeros(m)
for j in 1:m
    pvals[j]=pvalue( UnequalVarianceTTest(f_tfhm[j,findall(ikm.==1)] , f_tfhm[j,findall(ikm.==2)]) );
end
plot(pvals.<1e-3,lt=:bar)

t,chdat=read_channel(id,dt,900,fl);
writedlm("../4_outputs/ch"*string(id)*".dat",[t chdat]," ");
writedlm("../4_outputs/ch"*string(id)*"_tfhm.dat",tfhm," ");
writedlm("../4_outputs/ch"*string(id)*"_tfhm.dat",tfhm," ");
writedlm("../4_outputs/km_ch"*string(id)*"_tfhm_k1.dat",f_tfhm[:,findall(ikm.==1)]);
writedlm("../4_outputs/km_ch"*string(id)*"_tfhm_k2.dat",f_tfhm[:,findall(ikm.==2)]);
writedlm("../4_outputs/km_ch"*string(id)*"_idx.dat",[t[f_idx] ikm]);
writedlm("../4_outputs/km_ch"*string(id)*"_centers.dat",ckm);
writedlm("../4_outputs/km_ch"*string(id)*"_pvals.dat",pvals);

# Silhouette validation of k-means
# using LinearAlgebra
# s=size(f_tfhm,2)
# dists=zeros(s,s);
# for i in 1:s
#     for j in 1:s
#         dists[i,j]=norm(f_tfhm[:,j]-f_tfhm[:,i]);
#     end
# end
# sil=silhouettes(km,dists)
# mean(sil)
#
# Repeat for all channels [BUG: paralleization finishes with unfinished work, unknown reason]
n = 384;
l = 2000;
pvals=zeros(n,l)
state = Threads.Atomic{Int}(0);
Threads.@threads for id in failed
    Threads.atomic_add!(state, 1)
    print("--> ",state[]," out of ",n,"\n");
    flush(stdout);

    t,f,tfhm = timefreq(id,fl);
    f_idx,f_tfhm = movfilter(t,tfhm,fl);

    km = kmeans(f_tfhm,2,tol=1e-50)
    ikm = assignments(km)

    for j in 1:l
        pvals[id,j]=pvalue( UnequalVarianceTTest(f_tfhm[j,findall(ikm.==1)] , f_tfhm[j,findall(ikm.==2)]) );
    end
end

writedlm("../4_outputs/km_pvals.dat",pvals," ");
#
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
# Time-frequency analysis: check α and γ variations and correlations

# # Single channel:
# id=200;
# t,f,tfhm = timefreq(id,fl);

# # # Remove the segments containing movement:
# idx,f_tfhm = movfilter(t,tfhm,fl);

# # # Compute statistics for this channel:
# mean_psd=mean(f_tfhm,dims=2);
# std_psd =std(f_tfhm,dims=2);
# plot(f,mean_psd,ribbon=std_psd,c=1,fillalpha=0.25,xlim=(0,70))

# # Separate α and γ powers:
# len = size(f_tfhm)[2]
# α = zeros(len,2);
# γ = zeros(len,2);
# αr = 81:201;    # check f[αr]
# γr = 341:461;   # check f[γr]
# # αr = 21:201;    # check f[αr]
# # γr = 301:481;   # check f[γr]
# α[:,1] = mean(f_tfhm[αr,:],dims=1);
# α[:,2] = std(f_tfhm[αr,:],dims=1);
# γ[:,1] = mean(f_tfhm[γr,:],dims=1);
# γ[:,2] = std(f_tfhm[γr,:],dims=1);

# cor(α[:,1],γ[:,1]) # Pearson correlation

# # plot results
# l1="α ("*string(f[αr][1])*"-"*string(f[αr][end])*" Hz)";
# l2="γ ("*string(f[γr][1])*"-"*string(f[γr][end])*" Hz)";
# plot(1:len,α[:,1],ribbon=α[:,2],fillalpha=0.3,labels=l1)
# p1=plot!(1:len,γ[:,1],ribbon=γ[:,2],fillalpha=0.3, 
# xlabel="Time [s]", ylabel="Power [V²]",labels=l2)

# p2 = plot(α[:,1],γ[:,1],xerr=α[:,2],yerr=γ[:,2], key=false,
#  lt=:scatter,xlabel=l1,ylabel=l2,xticks=2e-16:4e-16:2e-15,yticks=5e-17:5e-17:2e-16)
# # plot(log.(α[:,1]),log.(γ[:,1]),lt=:scatter)

# # Pearson correlations for all channels
# n = 384;
# l = 2000;

# pearson = zeros(n,1);

# state = Threads.Atomic{Int}(0);

# Threads.@threads for id in 1:n
#     Threads.atomic_add!(state, 1)
#     print("--> ",state[]," out of ",n,"\n");
#     flush(stdout);
#     t,f,tfhm = timefreq(id,fl);
#     idx,f_tfhm = movfilter(t,tfhm,fl);    
#     αr = 81:201;    # check f[αr]
#     γr = 341:461;   # check f[γr]
#     # αr = 21:201;    # check f[αr]
#     # γr = 301:481;   # check f[γr]
#     pearson[id] = cor( mean(f_tfhm[αr,:],dims=1)', mean(f_tfhm[γr,:],dims=1)')[1]
# end


# p3=plot(pearson,1:384, ylabel="Channel",xlabel="R",legend=false,c=3)

# l=@layout [a ; [b c]]
# pl=plot(p1,p2,p3,layout=l,size=(800,800))
# savefig(pl,"../3_figures/figure3.pdf")
# # writedlm("../4_outputs/pearson_alpha_gamma_wider.dat",pearson,' ');



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
