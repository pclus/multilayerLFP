# Analysis of LFP from UPO-tACS Dataset

## Index
<!-- vim-markdown-toc GFM -->

* [Processing of the `.mat` files](#processing-of-the-mat-files)
	* [Local Field Potential files](#local-field-potential-files)
	* [Reading binary files with C](#reading-binary-files-with-c)
	* [Movement data](#movement-data)
* [Data Overview](#data-overview)
* [Frequency analysis:](#frequency-analysis)

<!-- vim-markdown-toc -->


## Processing of the `.mat` files

### Local Field Potential files

The "pre" and "post" data in the Matlab file have been transcribed to binary files for convenience.
The "pre" data has been stored in `1_Raw/pre.bin`, and the "post" data has been stored in `1_Raw/pre.bin`.
These files have been created by running this code in Matlab:

```
mf=matfile('1_Raw/Suj9.mat');

pre=mf.Suj9(1,1);
fid=fopen("pre.bin","w"); fwrite(fid,pre{1}(384:-1:1,:)','double'); fclose(fid)

post=mf.Suj9(2,1);
fid=fopen("post.bin","w"); fwrite(fid,post{1}(384:-1:1,:)','double'); fclose(fid)
```

Notice that, in the original Matlab matrices, the order of the channels has been inverted (i.e., `Suj{1}(1,:)` contains the data for channel 384).
In the binary files we revert back so that each channel ID matches its row in the file.
Additionally, we transpose the matrix so that each column is the data from one channel.
This is very convinient for later reading of the files, because Matlab stores matrices in column major.

### Reading binary files with C

To read these binary files notice that the matrix dimensions of each, `pre` and `post`, are 2250000x384.
The Julia function `read_channel()` in `analysis.jl` takes care of this (see later).
However we include an additional C function that can be used directly

The function `readbin` can be used to read such binary files
Source code is stored in the `0_SourceCode`directory.
The usage of the function is 

```
./readbin <infile> <outfile> <id> <t0> <tf>`
```
and outputs in `<outfile>` the time series of channel `<id>` (from 1 to 384)
from time `<t0>` to `<tf>` (from 0.0004 to 900).
The data must be stored in the binary `<infile>`.
For example:

```
./readbin 1_Raw/pre.bin t2.dat 350 200 300
```

stores the time series of channel 350 from 200s to 300s to the t2.dat file.


### Movement data

The movement can be converted to `.dat` file in octave:
```
load('mov.mat');
dlmwrite('mov_pre.dat',mov{1},' ')
dlmwrite('mov_dur.dat',mov{2},' ')
dlmwrite('mov_post.dat',mov{3},' ')
```
This is done only once, and the results are stored in the `1_Raw` folder.


## Data Overview


If `pre1.dat` contains all the 900s of a channel, then in gnuplot one can use:

```
plot 'mov_pre.dat' u 1:(-1e-5):(0):(2e-5) w vectors nohead lc 'gray' , 'pre1.dat' ev 10 w l lc 1 lw 1
```
## Frequency analysis:

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
plot 'tfhm100.dat' mat u (5+10*$1):(df*$2):3 w ima, '../../1_Raw/mov_pre.dat' u 1:(150)
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

Next steps:

1. Look for differences in time-frequency accross channels: Kmeans might be missleading (cloud points), should use different tools:
	- [] Moving window time-freq heatmap to visualize better the differences.
	- [] Clusting algorithm with authomatic cluster detection.
	- ...
2. Band-pass filter or freq. bands comparison 
	- [] Repeat pre/post analysis with single-long chunk of time-series.
3. Ask UPO:
	- [] Are the movement artifacts because of the movement, or the change in "brain state"? Ask UPO
	- [] What happens with the 50Hz and 60 Hz artifacts?
4. Other:
	- [] P-values heatmaps, use 4 colors (3 + white), one for each decade.
	- [] T-test assumes normalitiy, should check there are no long tails at least.

