#!/bin/bash

for i in {8,10,11,12,13,15,16,18,20};
do
	julia main.jl $i
done
