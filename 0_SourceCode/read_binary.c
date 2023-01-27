/* Date: Fri 27 ene 2023 14:45:35 CET
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

/* New version (binary files are 2250000x384)*/
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
	fprintf(stderr,"Assuming %d channels %d units long sampled at %lf Hz\n",n,m,1/dt);
	fprintf(stderr,"Output channel: %d\n",id);
	fprintf(stderr,"From position %d to %d\n",m0,mf);
	//-------------------------

	double *data=malloc(sizeof(double)*m);
	size_t s;
	for(long unsigned i=0;i<id;i++)
	{
		s=fread(data,sizeof(double),m,fin);
	}
	for(long unsigned i=m0;i<=mf;i++) 
	{
		fprintf(fout,"%lf %.16g\n",i*dt,data[i-1]); // m0 is never 0
	}
	fclose(fin);
	fclose(fout);
	
	free(data);
	return 1;
}


/* Old version (binary files are 384x2250000)*/
/*
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
	fprintf(stderr,"Assuming %d channels %d units long sampled at %lf Hz\n",n,m,1/dt);
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
*/
