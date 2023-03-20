# Check bipolars ----------------------------------------------
 t0=100
 tf=110
 t,bip=read_channel(100,t0,tf,"bipolar_pre");
 t,ch1=read_channel(100,t0,tf,"pre");
 t,ch2=read_channel(104,t0,tf,"pre");
 t,ch3=read_channel(108,t0,tf,"pre");

 plot(bip)
 plot!((ch1-ch3)/(2*20.0),lw=0.5)

 sum(abs.(bip-(ch1-ch3)/(2*20.0)))

 # Check CSD ---------------------------------------------------
  # Check CSD ---------------------------------------------------
 t,cs=read_channel(100,t0,tf,"csd_pre");
 t,ch5=read_channel(201,t0,tf,"pre");
 t,ch3=read_channel(199,t0,tf,"pre");
 t,ch4=read_channel(200,t0,tf,"pre");
 t,ch7=read_channel(203,t0,tf,"pre");
 t,ch8=read_channel(204,t0,tf,"pre");

 plot(cs,lw=2)
 plot!((4*ch5-ch3-ch4-ch7-ch8)/(25.61^2),lw=1,xlim=(0,100))

 sum(abs.(cs-(4*ch5-ch3-ch4-ch7-ch8)/(25.61^2)))

 # distance analysis, to validate distance between sites -------
 diff=zeros(96,6);
 for i in 1:4:384
     t,ch1=read_channel(i+0,t0,tf,fl);
     t,ch2=read_channel(i+1,t0,tf,fl);
     t,ch3=read_channel(i+2,t0,tf,fl);
     t,ch4=read_channel(i+3,t0,tf,fl);
     diff[Int(floor(i/4))+1,1]=norm(ch1-ch2)
     diff[Int(floor(i/4))+1,2]=norm(ch1-ch3)
     diff[Int(floor(i/4))+1,3]=norm(ch1-ch4)
     diff[Int(floor(i/4))+1,4]=norm(ch2-ch3)
     diff[Int(floor(i/4))+1,5]=norm(ch2-ch4)
     diff[Int(floor(i/4))+1,6]=norm(ch3-ch4)
 end

# bandstop filter validation ----------------------------------
# -------------------------------------------------------------
rate = 2500.0
dt = 1 / rate
t0 = dt
tf = 900.0
id = 100

fl="cortex_pre"
t, chdat = read_channel(id, t0, tf, fl);

bsfilter = digitalfilter(Bandstop(49.95,50.05; fs=rate),Butterworth(1))
fil = filtfilt(bsfilter, chdat)
bsfilter2 = digitalfilter(Bandstop(59.95,60.05; fs=rate),Butterworth(1))
fil = filtfilt(bsfilter2, fil)

NW = 1.0 * length(chdat) * dt / (2.0);  # real bandwith is ω*dt/2.0, and NW = N*ω*dt/2.0
K = 10;    # number of tappers (should be similar slightly less than 2*NW)
S = multispec(chdat, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);
S2 = multispec(fil, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);
p1 = plot(S.f, S.S, xlim=(40, 70), ylim=(1e-18, 1e-14), lw=1.0, yaxis=:log)
plot!(S2.f, S2.S, xlim=(40, 70), ylim=(1e-18, 1e-14), lw=0.5, yaxis=:log)
# -------------------------------------------------------------



# -------------------------------------------------------------
# Numerical check that the variance of the time signal is the overall power.
# Does not hold for the MT because of the tappering
# -------------------------------------------------------------
rate = 2500.0;
dt = 1.0 / rate;
NW = 1.0 * length(chdat) * dt / (2.0);  
K = 10;   
S = multispec(chdat, dt=dt, NW=NW, K=K, jk=true, Ftest=true, a_weight=false);


p1 = plot(S.f, S.S, xlim=(0, 200), ylim=(1e-18, 1e-14), lw=1.0, yaxis=:log)


using FFTW
Fr = rfft(chdat); 
wr = rfftfreq(length(t), 1.0/dt); 

area = 2.0*sum(abs2.(Fr))*1e-2*dt/length(chdat)
areaMT=sum(S.S)*1e-2 #df=1e-2
var(chdat)
# -------------------------------------------------------------


# autocorrelations --------------------------------------------
using StatsBase, HypothesisTests

rate = 2500.0
dt = 1 / rate
t0 = dt
tf = 900.0
id = 100

fl="cortex_pre"
t, chdat = read_channel(id, t0, tf, fl);

ν=autocor(chdat,0:10*2500)
plot(dt*(0:25000),ν,xlim=(0,1.0))
#--------------------------------------------------------------

#--------------------------------------------------------------
# check spectral differences with KS and AD tests (perform poorer than logspectral_distance)
id=100; fl="pre"
t, f, tfhm = timefreq(id, fl);

m1 = [pvalue(KS_psd(tfhm[:,i],tfhm[:,j],f)[1]) for i in 1:90,j in 1:90]
m2 = [pvalue(KS_psd(tfhm[:,i],tfhm[:,j],f)[2]) for i in 1:90,j in 1:90]
heatmap(m1,clim=(0,1e-1)) # changes A LOT depending on `sampl`
heatmap(m2,clim=(0,1e-1))
#--------------------------------------------------------------




# -------------------------------------------------------------
# -------------------------------------------------------------
# Create a heatmap of the pre data from 200s to 300s using Multitaper:
n0=226;
nf=361;
w, lfp_cortex = heatmapMT(200.0, 300.0, "cortex_pre", 0:(nf-n0));
w, lfp = heatmapMT(200.0, 300.0, "filtered_pre", 0:383);
w, csd = heatmapMT(200.0, 300.0, "csd_pre", 1:190);
w, bip = heatmapMT(200.0, 300.0, "bipolar_pre", 1:376);
w, kcsd = heatmapMT(200.0, 300.0, "kCSD_electrodes_pre", 1:384);
w, kcsd = heatmapMT(200.0, 300.0, "kCSD_centers_pre", 1:384);

w, lfp_post = heatmapMT(200.0, 300.0, "filtered_post", 1:384);
# gr()
# heatmap(w,1:384,kcsd')

writedlm("../4_outputs/cortex_lfp.dat", lfp_cortex', " ");
writedlm("../4_outputs/lfp.dat", lfp', " ");
writedlm("../4_outputs/csd.dat", csd', " ");
writedlm("../4_outputs/bip.dat", bip', " ");
writedlm("../4_outputs/kcsd_cent.dat", kcsd', " ");
# -------------------------------------------------------------
# -------------------------------------------------------------

# -------------------------------------------------------------
# Time-frequency analysis for a single channel:
fl = "cortex_pre"
id = 100
t, f, tfhm = timefreq(id, fl);
writedlm("../4_outputs/ch" * string(id) * "_tfhm_" * fl * ".dat", tfhm, " ");

# Remove the segments containing movement:
f_idx, f_tfhm = movfilter(t, tfhm, "pre");

# Compute statistics for this channel:
using Statistics, HypothesisTests
mean_psd = mean(f_tfhm, dims=2);
std_psd = std(f_tfhm, dims=2);
plot(f, mean_psd, ribbon=std_psd, c=1, fillalpha=0.25, yaxis=:log, yrange=(1e-23, 0.2e-21))
#
# -------------------------------------------------------------

# -------------------------------------------------------------
# Frequency for 1 second resolution
rate = 2500.0
dt = 1.0 / rate
n = 384
Δt = 1.0      # segment duration (10 seconds)
ns = Int(900.0 / Δt)  # number of segments
m = Int(2250000 / ns) # segments of 10 seconds
# mh = m/2;       # length of the MT spectra...
l = 1000       # but...up to l is enough to get the relevant freqs
NW = 1.0 * m * dt / (2.0)
K = 8
t, chdat = read_channel(id, 4e-4, 900, "filtered_pre")    # time series starts at 4e-4

tfhm = zeros(l,9000);   
for s in 5:8995
    segdat = chdat[Int(round(1+(s-5)*0.1*rate)):Int(round((s+5)*0.1*rate))]
    S  = multispec(segdat, dt=dt, NW=NW, K=K);
    tfhm[:,s] = S.S[1:l];
    if s%5==0
        print(s,"\r");
        flush(stdout);
    end
end
k=log10.(tfhm[:,1000:2000])
heatmap(k[:,600:800],ylim=(0,100),clim=(-17,-15))
# -------------------------------------------------------------


# -------------------------------------------------------------
id=100; fl="pre"
t, f, tfhm = timefreq_complete(id, fl);
m = [logspectral_dist(tfhm[:,i],tfhm[:,j],f) for i in 1:90,j in 1:90]
heatmap(m,clim=(0.0,0.45))

f_idx, f_tfhm = movfilter(t, tfhm, "pre");
smean = mean(f_tfhm,dims=2)
σ = std(f_tfhm,dims=2);
rs = randomize_spectra(smean,σ,90);
mr = [logspectral_dist(rs[:,i],rs[:,j],f) for i in 1:90,j in 1:90]
heatmap(mr) # nice result, but still, maybe surrogates

m4 = [logspectral_dist(f_tfhm[:,i],f_tfhm[:,j],f) for i in 1:82,j in 1:82]
heatmap(m4) #amazing
# -------------------------------------------------------------



plotlyjs()
plot([blfp_pre[:,100],blfp_post[:,100]],labels=["pre" "post"],yaxis=:log)