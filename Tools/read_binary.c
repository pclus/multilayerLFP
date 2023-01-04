/* Date: Tue 03 Jan 2023 07:41:00 PM CET
 * Author: Pau Clusella 
 * Email: pau.clusella@upf.edu
 * Github: pclus
 *
 * Reads binary file with pre or post data of the tACs experiment provided by UPO.
 * Outputs time window for a specific channel.
 * 
 * Compile with:
 * gcc read_binary.c -o ../readbin -O3 -lm
 *
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv)
{
	if(argc!=6)
	{
		fprintf(stderr,"USAGE: ./readbin \
				<file_in> <file_out> <channel> <time init> <time end>\n");
	}

	int n = 384; 
	int m = 2250000; 
	double dt=1.0/2500.0;

	FILE *fin=fopen(argv[1],"rb");
	FILE *fout=fopen(argv[2],"w");
	int id=atoi(argv[3]);
	double t0=atof(argv[4]);
	double tf=atof(argv[5]);
	int m0=ceil(t0/dt);
	int mf=floor(tf/dt);
	
	//-------------------------
	fprintf(stderr,"Assuming %d channels %d long sampled at %lf Hz\n",n,m,1/dt);
	fprintf(stderr,"Output channel: %d\n",id);
	fprintf(stderr,"From position %d to %d\n",m0,mf);
	//-------------------------
	id=id-1; 	// We read channels starting from 1

	double *data=malloc(sizeof(double)*n);
	size_t s;
	for(long unsigned i=1;i<m0;i++)
	{
		s=fread(data,sizeof(double),n,fin);
	}
	for(long unsigned i=m0;i<=mf;i++)
	{
		s=fread(data,sizeof(double),n,fin);
		fprintf(fout,"%lf %.16g\n",i*dt,data[id]);
	}
	fclose(fin);
	fclose(fout);
	
	free(data);
	return 1;
}
