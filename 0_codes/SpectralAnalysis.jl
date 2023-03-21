module SpectralAnalysis

using NeuropixelAnalysis, DelimitedFiles, Multitaper, Plots, DSP, Statistics;
using NumericalIntegration,HypothesisTests

export KS_psd,randomize_spectra
export logspectral_dist, kolmogorov_smirnov_dist
export kullback_leibler_divergence, itakura_saito_divergence

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

function logspectral_dist(s1,s2,f)
    # y = @.  10*log10(s1/s2)^2
    # integrate(f,y)^0.5 # maybe divided by f[end]^0.5
    y = @. log10(s1/s2)^2
    mean(y)^0.5 # maybe divided by f[end]^0.5
end

function randomize_spectra(smean,σ,nsampl)
    n = length(smean);
    rs = zeros(n,nsampl)
    for i in 1:nsampl
        rs[:,i] = smean + σ.*randn(n);
    end
    rs[findall(rs.<0)].=1e-22; # avoid negative powers, this is problematic though
    return rs;
end

function kolmogorov_smirnov_dist(S1,S2,f)
    df = f[2]-f[1]
    psd1 = S1/ (sum(S1)*df); # normalization
    cdf1 = [ sum(psd1[1:i])*df for i in 1:length(psd1)]

    psd2 = S2/ (sum(S2)*df); # normalization
    cdf2 = [ sum(psd2[1:i])*df for i in 1:length(psd2)]

    return maximum(cdf1-cdf2)
end

function kullback_leibler_divergence(S1,S2,f)
    df = f[2]-f[1]
    psd1 = S1/ (sum(S1)*df); # normalization
    psd2 = S2/ (sum(S2)*df); # normalization
    # psd1 = S1/ mean(S1)
    # psd2 = S2/ mean(S2)

    return sum(@. psd1*log(psd1/psd2))
end

function itakura_saito_divergence(S1,S2,f)
    df = f[2]-f[1]
    psd1 = S1/ (sum(S1)*df); # normalization
    psd2 = S2/ (sum(S2)*df); # normalization

    return sum(@. psd1/psd2 - log(psd1/psd2) - 1.0)
end

end # module