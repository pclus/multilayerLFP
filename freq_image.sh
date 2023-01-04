#!/bin/bash

# read binary, create psd's and create image
touch freqs.dat
for i in {1..384}; do
	./readbin Raw/pre.bin temp.dat $i 200 300
	./psd temp.dat
	awk '$1<=250.0{print $0}' psd_temp.dat > 4_outputs/psd_t${i}.dat

	awk '$1<=250.0{print $2}' psd_temp.dat > temp2.dat
	paste -d ' ' freqs.dat temp2.dat > temp3.dat;
	mv temp3.dat freqs.dat
done
rm -f temp.dat
rm -f temp2.dat
rm -f psd_temp.dat

# create image using existing psd's using the filtered version:
touch freqs_filtered.dat
for i in {1..384};
do
        awk '$1<=250.0{print $3}' 4_outputs/psd_t${i}.dat > temp2.dat
        paste -d ' ' freqs_filtered.dat temp2.dat > temp3.dat;
        mv temp3.dat freqs_filtered.dat
done
rm -f temp2.dat

