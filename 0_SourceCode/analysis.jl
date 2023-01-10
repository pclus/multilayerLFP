using DelimitedFiles

t0=200;
tf=300;
id=350;

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

t,chdat=read_channel(id,t0,tf);

# # Real FFT ------------------------------------------------
# using FFTW,Plots

# rate=2500.0;
# dt=1.0/rate;

# Fr = rfft(chdat); 
# wr = rfftfreq(length(t), rate); 
# plot(wr, abs.(Fr),xlim=(0,100)) # Logscale with yaxis=:log
# #----------------------------------------------------------

#------------------------------------------------------------
# Multitaper PSW
#------------------------------------------------------------
using Plots,Multitaper,StatsFuns ;

rate=2500.0;
dt=1.0/rate;

NW = 1.0*length(chdat)*dt/(2.0) ;  #bandwith is W*dt/2.0, and NW = N*W*dt/2.0
K  = 10;    # number of tappers (should be less than 2*NW)
S  = multispec(chdat, dt=dt, NW=NW, K=K,jk=true,Ftest=true, a_weight=true);

plot(S.f,S.S,xlim=(0,200),ylim=(1e-18,1e-15),lw=1.0,yaxis=:log)

# --- some notes about  multitaper:
# plot(S,xlim=(0,200),ylim=(1e-20,5e-16)) # authomatic plot
#
# comparison between multitaper and standard fft.
# Notice that multitaper package normalizes by dt/N, similar to what I am doing in C.
# plot(S.f,S.S.*length(S.f)/dt)
# plot!(wr,abs2.(Fr),xlim=(0,100),ylim=(0,1e-7),lw=0.1,yaxis=:log)
#
# Plot including 95% confidence intervals, as in the authomatic version...
# does not look right...it seems that the authomatic recipe plots the confidence
# for the log-spectra, and not the true spectra
# z=norminvcdf(0,1,0.975);
# plot(S.f,S.S,xlim=(0,200),ylim=(0,2e-16))
# plot!(S.f,S.S.*exp.(-z*sqrt.(S.jkvar)),fillrange=S.S.*exp.(2*z*sqrt.(S.jkvar)),c=1,alpha=0.35)
# plot!(yaxis=:log,ylim=(1e-20,5e-16))
# ---

#------------------------------------------------------------
# Spectral heatmap with Multitaper PSD

#------------------------------------------------------------
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


a,b=heatmapMT(200.0,300.0)
writedlm("../4_outputs/psd_mthm.dat","b"," ");

#------------------------------------------------------------
# Time-frequency analysis
#------------------------------------------------------------

id=100;
t0=200;
tf=210;
t,chdat=read_channel(id,t0,tf);

using Plots,Multitaper,StatsFuns ;

rate=2500.0;
dt=1.0/rate;

NW = 1*length(chdat)*dt/(2.0) ;  #bandwith is W*dt/2.0, and NW = N*W*dt/2.0
K  = 10;    # number of tappers (should be less than 2*NW)
S  = multispec(chdat, dt=dt, NW=NW, K=K);

plot(S.f,S.S,xlim=(0,200),ylim=(1e-20,1e-13),lw=1.0,yaxis=:log)
plot(S,xlim=(0,200));

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

t,f,tfhm = timefreq(150);
writedlm("../4_outputs/tfhm"*string(id)*".dat",tfhm," ");