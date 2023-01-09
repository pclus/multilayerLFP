/*
*  
*  The data is saved in a <file in> and you want to compute the data
*  from the entry <k0> to <kf>. The time step between to points in <file in> is <dt>.
*  The program won't analyse exactly the part you ask, but shorter so that is a power of 2.
*  This value is the 'length'
*  
*  COMPILE USING
*  	gcc psd_fast.c -o psd -lgsl -lgslcblas -lm
*  	gcc psd_fast.c -o psd -L/usr/local/lib -lgsl -lgslcblas -lm
*  PLOT PSD IN GNUPLOT USING 
* 	plot 'psd.dat' w l
* 
* 16th of January, 2017
* << Updated at March 2021 to read binary files >>
* 
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_filter.h>

#define RAND gsl_rng_uniform(r)//((double)rand())/RAND_MAX 
#define RANDI(n) ((int) n*RAND) 
#define RANDIFF (-1+2.*RAND)

///////////// GSL RNG ////////////////////
const gsl_rng_type * T;
gsl_rng * r;
////////////////////////////////////
int count_lines(char *name);
int gaussian_filter(double *power, int nn);


// arguments <file in name> <k initial> <k end>
int main(int argc, char **argv ){
//------------------------------------------------------------------
	/* Initialization */
//------------------------------------------------------------------
	if(argc!=2){
		fprintf(stderr,"USAGE: <file in name>\n");
		return -1;
	}
	fprintf(stderr,"Initializing\n");
	int i,j=0, flag;
	int k0=0, kf=count_lines(argv[1]); fprintf(stderr,"\n %d lines detected\n",kf);
	//int k0=atoi(argv[2]), kf=atoi(argv[3]);
	int length=kf-k0;
	length=(int) pow(2,floor(log2(length)));
	fprintf(stderr,"Time series length = %d\n",length);
        double a,b, dt; //=atof(argv[4]);

	char *name=argv[1];
	char *nameout=malloc(sizeof(char)*120);
	sprintf(nameout,"psd_%s",argv[1]);


	
        double *x=malloc(sizeof(double)*length*2);
	if(x==NULL){
		fprintf(stderr,"ERROR: Not enough memory.\n");
		fprintf(stderr,"ERROR: length=%d\n",length);
		return -1;
	}
        gsl_complex_packed_array v=x;
	
  ///////////// GSL RNG ////////////////////
  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  srand(time(0));
  gsl_rng_set(r,time(0));
  ////////////////////////////////////////// 

//-------------------------------------------------------------------
	/* Reading input file */
//-------------------------------------------------------------------
	
	fprintf(stderr,"Start reading\n");
	FILE *fin;
	if((fin=fopen(argv[1],"r"))==NULL){
		fprintf(stderr,"ERROR: File %s does not exist.\n", argv[1]);
		return -1;
	}
	double mean=0, t, t0; dt=0;
        for(i=0;(i<length && flag!=EOF);i++){
		flag=fscanf(fin,"%lf %lf\n",&t, x+i);
		if(i>0) dt+=t-t0;
		t0=t;
		mean+=x[i]; 
        }

	fprintf(stderr,"Time step = %lf\n",dt);
	dt=dt/(1.0*length-1);
	fprintf(stderr,"Time step = %lf\n",dt);
	fclose(fin);
	if(flag==EOF){
		fprintf(stderr,"ERROR: Length too large, length=%d\n", length);
		fprintf(stderr,"ERROR DATA:\n k0=%d\nkend=%d\ni=%d\n", k0, kf, i);
		return -1;
	}

	mean=mean/(1.0*length);
	for(i=0;i<length;i++) x[i]=x[i]-mean;

	fprintf(stderr,"Reading finished\n");

//-------------------------------------------------------------------
	/* Computing PSD and writting output */
//-------------------------------------------------------------------
// The power being multiplyed by dt and divided by the length.
// This might not be what you want in some cases (this is basically a problem of units).

	fprintf(stderr,"Start computing PSD\n");

        gsl_fft_real_radix2_transform (v, 1, length);
	double *power=malloc(sizeof(double)*(length/2+1));
	double *gaussian=malloc(sizeof(double)*(length/2+1));

	gaussian[0]=power[0]=(x[0]*x[0])*0.5*dt/(length*M_PI);
        for(i=1;i<length/2-1;i++){
		gaussian[i]=power[i]=(x[i]*x[i]+x[length-i]*x[length-i])*0.5*dt/(length*M_PI);
	}
	gaussian[i]=power[length/2]=(x[length/2]*x[length/2])*0.5*dt/(length*M_PI);

	/* GAUSSIAN FILTER */
	gaussian_filter(gaussian, (length/2+1));

	/* PRINT RESULTS */
        FILE *fout=fopen(nameout,"w");
        fprintf(fout,"%.16g %.16g %.16g\n",0.0/(length*dt),power[0],gaussian[0]);
        for(i=1;i<length/2-1;i++){
                fprintf(fout,"%.16g %.16g %.16g\n",i/(length*dt),power[i],gaussian[i]);
	}
        fprintf(fout,"%.16g %.16g %.16g\n",(0.5*length-1)/(length*dt),power[length/2], gaussian[length/2]);
	fclose(fout);

	/* FINISH */
	free(gaussian);
	free(power);
	free(nameout);
	free(v);
	fprintf(stderr,"Finished\n");
        return -1;
}

int gaussian_filter(double *power, int nn){
	int K_window=50;
	double alpha=5;
	int i;
	gsl_vector *input=gsl_vector_alloc(nn); for(i=0;i<nn;i++) gsl_vector_set(input,i,power[i]);
	gsl_vector *output=gsl_vector_alloc(nn); 

	// GSL_FILTER_END_PADZERO or GSL_FILTER_END_PADVALUE or  GSL_FILTER_END_TRUNCATE
	gsl_filter_gaussian_workspace * ws=gsl_filter_gaussian_alloc(K_window);
	gsl_filter_gaussian(GSL_FILTER_END_PADZERO, alpha, 0, input, output, ws);
	gsl_filter_gaussian_free(ws);
	
	for(i=0;i<nn;i++){
		power[i]=gsl_vector_get(output,i);
	}
	
	gsl_vector_free(input);
	gsl_vector_free(output);

	return 1;
}

int count_lines(char *name){
	FILE *stream;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	
	stream = fopen(name, "r");
	if (stream == NULL){
		fprintf(stderr,"ERROR: File %s does not exist.\n", name);
		return -1;
	}

	int count=0;
	while ((read = getline(&line, &len, stream)) != -1) {
		count++;
	}
	
	free(line);
	fclose(stream);
	return count;
}
