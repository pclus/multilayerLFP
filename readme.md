# Initial processing of the files:

## LFP Data:

In matlab:

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

(or `u ($0/2500.0):1`, depending on your format,
since time step is dt=1/2500.0).

Finally, for convenience, a binary file with all the data is stored in `Raw/pre.bin`:

```
mf=matfile('Raw/Suj9.mat');
pre=mf.Suj9(1,1);
fid=fopen("pre.bin","w"); fwrite(fid,pre{1}(:,:),'double'); fclose(fid)
```
To read this binary file notice that the matrix dimensions of `pre` are 384x2250000
and that Matlab stores in column major.

## Movement data:

The movement can be converted to dt file in octave:
```
load('mov.mat');
dlmwrite('mov_pre.dat',mov{1},' ')
dlmwrite('mov_dur.dat',mov{2},' ')
dlmwrite('mov_post.dat',mov{3},' ')
```
This is done only once, and the results are stored in the "Raw" folder.


## Data Overview:

In gnuplot: one can use:


```
plot 'mov_pre.dat' u 1:(-1e-5):(0):(2e-5) w vectors nohead lc 'gray' , 'pre1.dat' ev 10 w l lc 1 lw 1
```
# Frequency analysis:

Cut the time series using

```
awk '$1>=200.0 && $1<=300{print $0}' pre1.dat > cut_pre1.dat
```
