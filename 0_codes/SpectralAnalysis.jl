module SpectralAnalysis

using NeuropixelAnalysis, DelimitedFiles, Multitaper, Plots, DSP, Statistics;
using NumericalIntegration,HypothesisTests,FFTW

export KS_psd,randomize_spectra, spectral_stationarity
export logspectral_dist

function KS_psd(S1,S2,f)
    df = f[2]-f[1]
    psd1 = S1/ (sum(S1)*df); # normalization
    cdf1 = [ sum(psd1[1:i])*df for i in 1:length(psd1)]

    psd2 = S2/ (sum(S2)*df); # normalization
    cdf2 = [ sum(psd2[1:i])*df for i in 1:length(psd2)]

    sampl=1000   # very sensible to this parameter
    ζ1=rand(sampl)
    fsampl1 = [f[findfirst(ζ1[i].<=cdf1)] for i in 1:sampl]

    ζ2=rand(sampl)
    fsampl2 = [f[findfirst(ζ2[i].<=cdf2)] for i in 1:sampl]

    # notice that we are using the `dev` version of the Hyp.Test. package
    # (revert with `] free HypothesisTest`)
    tKS = ApproximateTwoSampleKSTest(fsampl1,fsampl2)
    tAD = KSampleADTest(fsampl1,fsampl2)
    return tKS,tAD,cdf1,cdf2
end


"""
    spectral_stationarity(id,fl,Δt,mov_filter::String)

Given a channel index id and a file to test, computes the tfhm and its surrogate.
This function should be incorporated in the main analysis with modifications to 
avoid computing the tfhm twice.
"""
function spectral_stationarity(id,fl,Δt,mov_filter::String)

    rate = 2500;
    dt = 1/rate;
    t, y = read_channel(id, 4e-4, 900, fl);
    T  = length(y)/rate;
    ns = Int(T / Δt);

    if !isempty(mov_filter)
        dts = T/ns;
        ts = 0.5*dts:dts:T
        tr,yr = movfilter(ts,reshape(y,:,ns),mov_filter,Δt);
        y=reshape(yr,:,1)[:,1];
    end

    S = rfft(y); 
    # f = rfftfreq(length(t), 1.0/dt); 
    S0 = @. abs.(S)*exp(im*rand()*2*π);
    y0 = irfft(S0,2*length(S0)-2);

    t, f, tfhm = timefreq(y,Δt);
    smean = mean(tfhm,dims=2)[:,1];
    ts,f,tfhm0 = timefreq(y0,Δt);
    smean0 = mean(tfhm0,dims=2)[:,1];

    ns = length(ts)
    
    q = [logspectral_dist(tfhm[:,i],smean,f) for i in 1:ns]
    q0 = [logspectral_dist(tfhm0[:,i],smean0,f) for i in 1:ns]

    return q,q0,tfhm,tfhm0,f;
end


function spectral_stationarity(id,fl,Δt,mov_filter::String)

    rate = 2500;
    dt = 1/rate;
    t, y = read_channel(id, 4e-4, 900, fl);
    T  = length(y)/rate;
    ns = Int(T / Δt);

    if !isempty(mov_filter)
        dts = T/ns;
        ts = 0.5*dts:dts:T
        tr,yr = movfilter(ts,reshape(y,:,ns),mov_filter,Δt);
        y=reshape(yr,:,1)[:,1];
    end

    S = rfft(y); 
    # f = rfftfreq(length(t), 1.0/dt); 
    S0 = @. abs.(S)*exp(im*rand()*2*π);
    y0 = irfft(S0,2*length(S0)-2);

    t, f, tfhm = timefreq(y,Δt);
    smean = mean(tfhm,dims=2)[:,1];
    ts,f,tfhm0 = timefreq(y0,Δt);
    smean0 = mean(tfhm0,dims=2)[:,1];

    ns = length(ts)
    
    q = [logspectral_dist(tfhm[:,i],smean,f) for i in 1:ns]
    q0 = [logspectral_dist(tfhm0[:,i],smean0,f) for i in 1:ns]

    return q,q0,tfhm,tfhm0,f;
end


end # module