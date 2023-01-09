# Initial processing of the files:

## LFP Data:

The "pre" data in the matlab file has been transcribed to a binary file using
For convenience, a binary file with all the "pre" data has been stored in `Raw/pre.bin`.
This file has been created by running this code in Matlab:

```
mf=matfile('Raw/Suj9.mat');
pre=mf.Suj9(1,1);
fid=fopen("pre.bin","w"); fwrite(fid,pre{1}(:,:),'double'); fclose(fid)
```
To read this binary file notice that the matrix dimensions of `pre` are 384x2250000
and that Matlab stores in column major.

The function `readbin` can be used to read such binary file.
Source code is stored in the `Tools`directory.
The usage of the function is 

```
./readbin <infile> <outfile> <id> <t0> <t0>`
```
and outputs in `<outfile>` the time series of channel `<id>` (from 1 to 384)
from time `<t0>` to `<tf>` (from 0 to 900).
The data must be stored in the binary `<infile>`.
For example:

```
./readbin Raw/pre.bin t2.dat 350 200 300
```

stores the time series of channel 350 from 200s to 300s to the t2.dat file.

### Other options:

If you do not want to use the binary file you can work directly in Matlab,
outputting single channels:

```
mf=matfile('Raw/Suj9.mat');
pre=mf.Suj9(1,1);

i=1;
dt=1.0/2500.0;
l=numel(pre{1}(i,:));
data=zeros(l,2);
data(:,1)=(1:l).*dt;
data(:,2)=pre{1}(i,:);
dlmwrite('pre1.dat',data,'delimiter',' ','precision','%.8g');

```
Alternatively (not recommended), one can simply use:
```
mf=matfile('Raw/Suj9.mat');
pre=mf.Suj9(1,1);
dlmwrite('pre.dat',pre{1}(1,:),'delimiter','\n','precision','%.8g');
```
and only one column is outputted. Remember then to use double format to output
with awk if further processing is done:
```
awk '{printf "%.16g %.16g\n",NR/2500.0,$1}' pre.dat > ppre1.dat
```

These short scripts output the data from channel `i` to a `.dat` file.
The true time series can then be plot using gnuplot:

```
plot 'pre.dat' u 1:2 
```

or `u ($0/2500.0):1`, depending on your format,
since time step is dt=1/2500.0).
The time series can be manually cut using:

```
awk '$1>=200.0 && $1<=300{print $0}' pre1.dat > cut_pre1.dat
```

## Movement data:

The movement can be converted to dat file in octave:
```
load('mov.mat');
dlmwrite('mov_pre.dat',mov{1},' ')
dlmwrite('mov_dur.dat',mov{2},' ')
dlmwrite('mov_post.dat',mov{3},' ')
```
This is done only once, and the results are stored in the "Raw" folder.


## Data Overview:


If `pre1.dat` contains all the 900s of a channel, then in gnuplot one can use:

```
plot 'mov_pre.dat' u 1:(-1e-5):(0):(2e-5) w vectors nohead lc 'gray' , 'pre1.dat' ev 10 w l lc 1 lw 1
```
## Frequency analysis:

Use script `freq_image.sh` (Tools folder).
It can be used to generate psd for each channel and gather 
them to produce a spectral heatmap:

```
df=0.019073486328125
plot 'freqs.dat' matrix u (df*($2-1)):1:3 w ima
```

The julia script `analysis.jl` also produces a heatmap `psd_mthm.dat` using a multitaper power spectra.
You can plot the result using, for instance, the gnuplot calls:

```
df=0.00999996000016
set logs zcb
set cbrange[1e-18:1e-15]
plot 'psd_mthm.dat' u (df*$1):2:3 matrix w ima
```

Next steps:

1. Parallelize the julia code
2. Time-freq analysis
3. Post-data
