#!/bin/bash

for i in {20,};
do
	julia main.jl $i
done

#DEST="/media/pclusella/Seagate\ Hub/Athena/UPO_data_processed/Suj${i}"
#mkdir -p ${DEST}
#mv ../*bin ${DEST}/
#mv ../mov*dat ${DEST}/ 
