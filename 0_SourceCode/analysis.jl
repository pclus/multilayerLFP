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
        read!(fin,data);
        channel[i]=data[id];
    end
    close(fin);

    t=collect(dt*(m0:mf));
    writedlm("../temp.dat",[t channel]," ",);
    return t, channel;
end

t,chdat=read_channel(id,t0,tf);

using FFTW,Plots

rate=2500.0;
dt=1.0/rate;

# freqs=fftfreq(length(t),2500)
psd=abs2.(rfft(chdat))
len=length(psd)
freqs=collect(1:len)