#!/bin/bash

for i in {16,18,19,20};
do
	julia main.jl $i
done
