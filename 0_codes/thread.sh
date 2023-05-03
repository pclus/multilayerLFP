#!/bin/bash

#for i in {8,9,10,11,12,13,15,16,18,20,19};
for i in {10,};
do
	julia main.jl $i
	DEST="~/Documents/Athena/UPO_data_processed/suj${i}"
	mkdir -p ${DEST}
	mv ../1_data/*bin ${DEST}/
	mv ../1_data/mov*dat ${DEST}/ 
done

