# Analysis of LFP from UPO-tACS Dataset

## Index
<!-- vim-markdown-toc GFM -->

* [Processing of the `.mat` files](#processing-of-the-mat-files)
	* [Local Field Potential files](#local-field-potential-files)
	* [Reading binary files with C](#reading-binary-files-with-c)
	* [Movement data](#movement-data)
* [Data Overview](#data-overview)
* [Bandpass filter, bipolar, and CSD](#bandpass-filter-bipolar-and-csd)
* [Frequency analysis](#frequency-analysis)
* [Postprocessing and statistical analysis](#postprocessing-and-statistical-analysis)
* [Next steps:](#next-steps)

<!-- vim-markdown-toc -->


## Processing of the `.mat` files

### Local Field Potential files

The "pre" and "post" data in the original Matlab file have been transcribed to binary files for convenience.
The "pre" data has been stored in `1_data/pre.bin`, and the "post" data has been stored in `1_data/pre.bin`.
These files have been created by running the code `tobinary.m` in Matlab (remember it is possible to run matlab from a terminal with `matlab -nodesktop`).
This can be done from the command line with:

```
FN=\'/media/pclusella/Pandora/UPO_data/suj9.mat\'
matlab -nodesktop -nosplash -nodisplay -r "tobinary($FN); exit()"
```

Notice that, in the original Matlab matrices, the order of the channels has been inverted (i.e., `Suj{1}(1,:)` contains the data for channel 384).
In the binary files we revert back so that each channel ID matches its row in the file.
Additionally, we transpose the matrix so that each column is the data from one channel.
This is very convinient for later reading of the files, because Matlab stores matrices in column major.

We have 12 subjects, with ids: 8,9,10,11,12,13,15,16,18,19,20,*27*

Some comments about the data:
	1. *Suj9* does not have the entire hour post-stimulations, it goes up to 7289726 time recordings instead of 9e6
	2. *Suj19* does not have movement because it did not register.
	3. *Suj27* is the one under anesthesia

### Reading binary files with C

The binary files contain a matrix with dimensions 2250000x384.
The Julia function `read_channel()` in `NeuropixelAnalysis.jl` takes care of this (see below).
However we include an additional C program, `readbin.c`, that
can be used without need of Julia.
Source code is stored in the `0_codes` directory and can be compiled
with `gcc` from that directory using `gcc read_binary.c -o ../readbin -O3`.
The usage of the function is then:

```
./readbin <infile> <outfile> <id> <t0> <tf>`
```
As a result, the time series of channel `<id>` (from 0 to 383) stored in the binary `<infile>` 
from time `<t0>` to `<tf>` (from 0.0004 to 900) will be copied in the newly created file `<outfile>` (will override existing files).
For example:

```
./readbin 1_data/pre.bin ts.dat 350 200 300
```

stores the time series of channel 350 from 200s to 300s to the `ts.dat` file.


### Movement data

The movement can be converted to `.dat` file in octave:
```
load('raw/mov.mat');
dlmwrite('mov_pre.dat',mov{1},' ')
dlmwrite('mov_dur.dat',mov{2},' ')
dlmwrite('mov_post.dat',mov{3},' ')
```
This is done only once, and the results are stored in the `1_data` folder.


## Data Overview


If `pre1.dat` contains all the 900s of a channel, then in gnuplot one can use:

```
plot 'mov_pre.dat' u 1:(-1e-5):(0):(2e-5) w vectors nohead lc 'gray' , 'pre1.dat' ev 10 w l lc 1 lw 1
```

A quick way to overview the data from different channels. First we generate the data from the binaries:

```
for i in {55,56,57,58,59};
do 
	./readbin 1_data/pre.bin t${i}.dat ${i} 0.0004 900.0; 
done
```

Then we plot in gnuplot:

```
j=1;
set multiplot layout 3,2;
do for[k=55:59]{
	plot 't'.k.'.dat' ev 10 w l lw 1 lc j t 't'.k;
	j=j+1
};
unset multiplot;
```

Also, we can generate a statistical summary in gnuplot:

```
set print 'stats_filtered_pre.dat'
do for[k=1:384]{
	cmd="./readbin 1_data/filtered_pre.bin temp.dat ".k." 0.0004 900.0"
	system(cmd)
	stats 'temp.dat' nooutput
	print k,STATS_mean_y,STATS_stddev_y
}
set print
```

## Bandpass filter, bipolar, and CSD

The Julia code `analysis.jl` contains the functions `bandpass_filter()`, `compute_bipolar()`, and `compute_csd()`
that read and filter the original LFP data and compute, respectively, the bipolar potential
differences and the CSD:

- `bandpass_filter()` applies a Butterworth filter of order 3 to bandpass the original LFP data between 1-300Hz. This does not affect most of the frequency analysis (since we work in frequencies below 200Hz), but significantly reduces the noisy fluctuations.

- `compute_bipolar()` computes the first derivative along each column of the prove of the filtered data. 
Therefore, in the `bipolar_<flag>.bin` binary files the 4 first and the 4 last channels
are missing. So there are 376 channels in the binaries, and `read_channel(id,t0,tf,"bipolar_pre")` reads the bipolar data for channel `id+4`.

- `compute_csd()` computes the CSD using Poisson's equation on the filtered data. It can only be calculated within the inner columns of the prove.
Therefore, in `csd_<flag>.bin` only 384/2-2=190 channels are computed (the two inner columns minus
the first and last row, i.e., channels 1 and 384). The relation between what is passed to `read_channel` and the actual
channel is given in the following table:

<div align="center">

|`id`| 1 | 2 | 3 | 4 |... | i | ... | 188 | 189 | 190 |
|:-: |:-:|:-:|:-:|:-:|:-:| :-: | :-: | :-: | :-:|  :-:|
| `channel` | 4 | 5 | 8 | 9 | ...  | 2i+1+i%2 | ... | 377 | 380 | 381 |

</div>

Alternatively we also use the [kernel-CSD](https://doi.org/10.1162/neco_a_00236).
For this we use the [Python library](http://biorxiv.org/lookup/doi/10.1101/708511) created by the same group. This method has several advantages (and some constrains):
- It works with arbitrary positions of the electrodes, thus we can also exclude broken electrodes.
- It is more robust to measurement noise and broken electrodes.
- Assumes arbitrary number of sources with Gaussian profiles.
- On the downside, has many parameters to tune, which might cause artifacts if not set properly.

We are still generating, validating, and comparing the results, but it produces better results than finite differences.
The comparison with bipolar, standard CSD, and kCSD should appear in a separate document.

Cross-validation of kCSD method gives the following optimal parameters to estimate
CSD using only cortex electrodes: R=30.0, h=10, lambda=1e-10.
The cross-validation results do not seem to change greatly for R between 10 and 40,
 lambda 1e-12 to 1e-8, and h between 1 and 20.

## Frequency analysis

Use script `freq_image.sh` (Tools folder).
It can be used to generate psd for each channel and gather 
them to produce a spectral heatmap:

```
set logs zcb;
set cbrange[1e-18:1e-15]
df=0.019073486328125
plot 'psd_rawhm.dat' matrix u (df*($2-1)):1:3 w ima
```

The julia script `analysis.jl` also produces a heatmap `psd_mthm.dat` using a multitaper power spectra.
You can plot the result using, for instance, the gnuplot calls:

```
df=0.00999996000016
set logs zcb
set cbrange[1e-18:1e-15]
plot 'psd_mthm.dat' u (df*$1):2:3 matrix w ima
```

The same script contains the function `timefreq` that can be used to perform a time-frequency heatmap
for specific channels.
In the case of segments of 10 seconds, the results outputted by `writedlm` can be plot in gnuplot with

```
set cbrange[1e-17:5e-15]
set log zcb
df=0.0999960001599936
plot 'tfhm100.dat' mat u (5+10*$1):(df*$2):3 w ima, '../../1_data/mov_pre.dat' u 1:(150)
```

Then, we can use the function `movfilter` to remove the data segments that
contain movement recordings, and average then all the resulting PSDs to
generate a global average heatmap.

These can be plot for specific channels, including the standard deviations, in gnuplot:

```
df=0.1;
id=150 # channel id
comm="awk -v k=".id." 'NR==FNR{a[FNR]=$k} NR>FNR{b[FNR]=$k} END {for(i in a){print a[i],b[i]}}' psd_mean_tfhm.dat psd_std_tfhm.dat"
plot "< ".comm u ($0*df):($1-$2):($1+$2) w filledcu fs solid 0.5, '' u ($0*df):1 w l lc 1
```
In order to generate the files to plot externally use, in a terminal,

```
CH=100 # or any other id
awk -v k=$CH 'NR==FNR{a[FNR]=$k} NR>FNR{b[FNR]=$k} END {for(i in a){print a[i],b[i]}}' psd_mean_tfhm.dat psd_std_tfhm.dat > psd_mean_tf_ch${CH}.dat
```

To plot the total heatmap:

```
df=0.1
set cbrange[1e-18:1e-15]
plot 'psd_mean_tfhm.dat' matrix u (df*$2):1:3 w ima
```

## Postprocessing and statistical analysis

1. Postprocessing:
	- Bandpass filter (Butterworth 3rd order at 300Hz). 
	- Bandstop (notch) filter at 50 and 60Hz (Butterworth 1st order, 0.1Hz bandwidth, i.e., 49.95-50.05Hz.
	- Compute bLFP, sCSD, kCSD.
	- Consider only signals from electrodes within the cortex.
2. For each datatype (LFP,bLFP,sCSD,kCSD) test for stationarity using Witt et al. (1997). Should see that signals are stationarity if movement is excluded (segments of 10s).
3. Consider each segment as an independent recording, at compare pre and post datasets (t-test, permutation test, cluster-based permutation test [spatial and spectral clusters]).


## Next steps:


0. Code cleaning:
	-[ ] Clear the pipelene
1. Documentation: create three separate documentation blocks, with their repos
	-[ ] Pipeline (main repo)
	-[ ] kCSD
	-[ ] Surrogate analysis
2. Further research:
	-[ ] Correlations accross channels
	-[ ] Artifacts in time series (aside from movement). In particular `suj9`
	-[ ] Anesthesized animal
