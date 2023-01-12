using DelimitedFiles,Multitaper,Plots;

# -------------------------------------------------------------
# -------------------------------------------------------------
# FUNCTIONS:
# -------------------------------------------------------------
# Reads the "pre" data for a specific channel and time window:
# -------------------------------------------------------------
function read_channel(id,t0,tf) # Equivalent to the read_binary.c code
    id=Int(id);
    fin=open("../1_Raw/pre.bin","r");
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
function heatmapMT(t0,tf)
    
    rate = 2500.0;
    dt = 1.0/rate;
    n = 384;
    m = 250001;
    NW = 1.0*m*dt/(2.0) ;  #bandwith is W*dt/2.0, and NW = N*W*dt/2.0
    K = 10;    # number of tappers (should be less than 2*NW)


    l = 30002; # WARNING: This assumes relevant frequencies below entry l of the spectra
    hm = zeros(n,l);

    # This should be parallelized
    for id in 1:n
        t,chdat=read_channel(id,t0,tf);
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
function timefreq(id)

    rate = 2500.0;
    dt = 1.0/rate;
    n = 384;
    Δt = 10.0;      # segment duration (10 seconds)
    ns = Int(900.0/Δt);  # number of segments
    m = Int(2250000/ns); # segments of 10 seconds
    # mh = m/2;       # length of the MT spectra...
    l = 2000;       # but...up to l is enough to get the relevant freqs
    NW = 1.0*m*dt/(2.0) ;  
    K = 10;
    
    t,chdat = read_channel(id,4e-4,900);    # time series starts at 4e-4
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
function movfilter(t,tfhm)
    mov_pre=readdlm("../1_Raw/mov_pre.dat",' ');
    Δt=10.0;
    mov_indx=Int.(floor.(mov_pre./Δt)*Δt.+5);
    indx = t.∈ [mov_indx];
    indx = findall(==(0),indx)
    f_tfhm = tfhm[:,indx];
    f_tfhm
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
id=350;

t,chdat=read_channel(id,t0,tf);

# Compute the PSD using standard FFT:--------------------------
# using FFTW
# rate=2500.0;
# dt=1.0/rate;
# Fr = rfft(chdat); 
# wr = rfftfreq(length(t), rate); 
# plot(wr, abs.(Fr),xlim=(0,100)) # Logscale with yaxis=:log
#--------------------------------------------------------------


# Compute PSD using Multitaper PSD ----------------------------
rate=2500.0;
dt=1.0/rate;

NW = 1.0*length(chdat)*dt/(2.0) ;  #bandwith is W*dt/2.0, and NW = N*W*dt/2.0
K  = 10;    # number of tappers (should be less than 2*NW)
S  = multispec(chdat, dt=dt, NW=NW, K=K,jk=true,Ftest=true, a_weight=true);

plot(S.f,S.S,xlim=(0,200),ylim=(1e-18,1e-15),lw=1.0,yaxis=:log)

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
# -------------------------------------------------------------

# -------------------------------------------------------------
# Create a heatmap of the pre data from 200s to 300s using Multitaper:
# a,b=heatmapMT(200.0,300.0)
# writedlm("../4_outputs/psd_mthm.dat",b," ");
# -------------------------------------------------------------

# -------------------------------------------------------------
# Time-frequency analysis for a single channel:
# t,f,tfhm = timefreq(150);
# # writedlm("../4_outputs/tfhm"*string(id)*".dat",tfhm," ");

# # Remove the segments containing movement:
# f_tfhm = movfilter(tfhm);

# # Compute statistics for this channel:
# using Statistics
# mean_psd=mean(f_tfhm,dims=2);
# std_psd =std(f_tfhm,dims=2);
# plot(f,mean_psd,ribbon=std_psd,c=1,fillalpha=0.25)
# -------------------------------------------------------------

# -------------------------------------------------------------
# Repeat the analysis for all channels and create the heatmaps
using Statistics
n = 384;
l = 2000;

psd_mean_tfhm = zeros(n,l);
psd_std_tfhm  = zeros(n,l);

state = Threads.Atomic{Int}(0);

Threads.@threads for id=1:n
    Threads.atomic_add!(state, 1)
    print("--> ",state[]," out of ",n,"\n");
    flush(stdout);
    t,f,tfhm = timefreq(id);
    f_tfhm = movfilter(t,tfhm);
    psd_mean_tfhm[id,:] = mean(f_tfhm,dims=2);  
    psd_std_tfhm[id,:] = std(f_tfhm,dims=2);
end

# writedlm("../4_outputs/psd_mean_tfhm.dat",psd_mean_tfhm',' ');
# writedlm("../4_outputs/psd_std_tfhm.dat",psd_std_tfhm',' ');

heatmap(0.1:0.1:200,1:n,log.(psd_mean_tfhm))
# -------------------------------------------------------------

# -------------------------------------------------------------
# Time-frequency analysis: check α and γ variations and correlations

# Single channel:
id=200;
t,f,tfhm = timefreq(id);

# # Remove the segments containing movement:
f_tfhm = movfilter(t,tfhm);

# # Compute statistics for this channel:
mean_psd=mean(f_tfhm,dims=2);
std_psd =std(f_tfhm,dims=2);
plot(f,mean_psd,ribbon=std_psd,c=1,fillalpha=0.25,xlim=(0,70))

# Separate α and γ powers:
len = size(f_tfhm)[2]
α = zeros(len,2);
γ = zeros(len,2);
αr = 81:201;    # check f[αr]
γr = 341:461;   # check f[γr]
# αr = 21:201;    # check f[αr]
# γr = 301:481;   # check f[γr]
α[:,1] = mean(f_tfhm[αr,:],dims=1);
α[:,2] = std(f_tfhm[αr,:],dims=1);
γ[:,1] = mean(f_tfhm[γr,:],dims=1);
γ[:,2] = std(f_tfhm[γr,:],dims=1);

cor(α[:,1],γ[:,1]) # Pearson correlation

# plot results
l1="α ("*string(f[αr][1])*"-"*string(f[αr][end])*" Hz)";
l2="γ ("*string(f[γr][1])*"-"*string(f[γr][end])*" Hz)";
plot(1:len,α[:,1],ribbon=α[:,2],fillalpha=0.3,labels=l1)
p1=plot!(1:len,γ[:,1],ribbon=γ[:,2],fillalpha=0.3, 
xlabel="Time [s]", ylabel="Power [V²]",labels=l2)

p2 = plot(α[:,1],γ[:,1],xerr=α[:,2],yerr=γ[:,2], key=false,
 lt=:scatter,xlabel=l1,ylabel=l2,xticks=2e-16:4e-16:2e-15,yticks=5e-17:5e-17:2e-16)
# plot(log.(α[:,1]),log.(γ[:,1]),lt=:scatter)

# Pearson correlations for all channels
n = 384;
l = 2000;

pearson = zeros(n,1);

state = Threads.Atomic{Int}(0);

Threads.@threads for id=1:n
    Threads.atomic_add!(state, 1)
    print("--> ",state[]," out of ",n,"\n");
    flush(stdout);
    t,f,tfhm = timefreq(id);
    f_tfhm = movfilter(t,tfhm);    
    αr = 81:201;    # check f[αr]
    γr = 341:461;   # check f[γr]
    # αr = 21:201;    # check f[αr]
    # γr = 301:481;   # check f[γr]
    pearson[id] = cor( mean(f_tfhm[αr,:],dims=1)', mean(f_tfhm[γr,:],dims=1)')[1]
end


p3=plot(pearson,1:384, ylabel="Channel",xlabel="R",legend=false,c=3)

l=@layout [a ; [b c]]
pl=plot(p1,p2,p3,layout=l,size=(800,800))
savefig(pl,"../3_figures/figure3.pdf")
# writedlm("../4_outputs/pearson_alpha_gamma_wider.dat",pearson,' ');



# -------------------------------------------------------------