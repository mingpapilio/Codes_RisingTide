/*
 * Continual generation model with good/ bad year
 * Focal point is founder effect and the selection of better survival/ better birth species
 * First day of building: 19 July, 2017
 * Transformed to short-long term variability model from 23 Feb, 2018
 * File instruction:
 * 1. Put "gen_beta.h" and "gen_beta.c" in the folder containing this file
 * 2. Put the dsfmt folder and the folder containing this file into the same folder
 ************************************************************************************************

Execution: (g++ and gcc both work)
gcc cdmc_ts.c gen_beta.h gen_beta.c -lm -stdlib=libstdc++
./a.out

Plot with gnuplot:
gnuplot
plot 'summary.txt' using 3:1 title 'rising-tide' with lines lc rgb 'orange',\
'summary.txt' using 3:2 title 'bet-hedging' with lines lc rgb 'skyblue'

 ************************************************************************************************
 * Key parameters
 * mean_env: The average environmental condition
 * shape_long_var, shape_short_var: The shape parameter of Beta distribution for generating environemntal conditions, both long-term and short-term. Larger velues indicate narrower distribution. 
 * T_long, T_short: The relative durations of shorter- and longer-term conditions (T_long= T_short*m)  
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>	
#include "gen_beta.h"

// Functions
    void message_error(char error_text[]);              // printing errors
    double *d_vector(long size);                        // creating a vector
    double **d_matrix(long size_row, long size_column); // creating a matrix
    void free_d_vector(double *x);                      // erasing a vector
    void free_d_matrix(double **x);                     // erasing a matrix
    void rk4(double p[], double k1[],int n, double t, double h, double pout[],void(*diff)(double,double[],double[]));
    void differential(double time,double in[],double out[]);
    double dx_dt(double time, double vr[]);             //vr is a vector of variable
    double dy_dt(double time, double vr[]);
    double normal_dist_BM (double mean, double sd, double u1, double u2);
    double factorial (double x); 
    double mean (double x[], int length);

// Global variables
    // Index of variables
        int i_spcl= 1;          // Index of specialist (rising-tide strategy)
        int i_genl= 2;          // Index of generalist (bet-hedging strategy)
    // differential equation paramters (K, b, d will change according to environmental conditions)
        double b_spcl;          // Intrinsic growth rate of specialists
        double b_genl;          // Intrinsic growth rate of generallists
        double d_spcl;          // Environment-dependent mortality of specialists
        double d_genl;          // Environment-dependent mortality of generalists 
        double K_spcl;          // Carrying capacity of specialists 
        double K_genl;          // Carrying capacity of generalists 
        double a_sg= 0.4;       // Intensity of inter-strategic interactions (specialists to generalists)
        double a_gs= 0.4;       // Intensity of inter-strategic interactions (generalists to specialists)

// Main function
int main (void)
{
    // Switches
        int s4= 1;				        // Switch of terminating tiny populaitons
        int s6= 2;                      // Specify the type of distribution in environmental conditions
    // Temporal variables
        int N_gen= 1000;                // Number of generations (calculations)
        int N_event= 20;                // Number of events within a generation (calculation), equivilent to m in the text
    // Basic variables
		int i,j,c1,c2,c3;               // For loop counters
		double pp=             0.0;	    // Temp for random number
        double ext_thr=        0.5;     // Threshold of population for terminating (related to s4)
        int rep=               100.0;   // Number of repetition
    // Performance curve parameters
                double scale_K=     250.0;                      // Scaling coefficient of carrying capacity
                double scale_b=     0.5;                        // Scaling coefficient of intrinsic growth rate
                double scale_d=     0.01;                       // Scaling coefficient of mortality rate
                double range_bio=   100.0;                      // The responsive environmental conditions (0-range_bio) of the biological performance curves
            // Beta dist
                double shape_spcl=	     50.0;                  // Shape of specialist's performance curve in carrying capacity (The larger, the narrower)
                double shape_genl=	     5.0;                   // Shape of generalist's performance curve in carrying capacity
                double skew_spcl=        1.0;                   // Skewness of specialist's performance curve, alpha= skew*shape*2
                double skew_genl=        1.0;                   // Skewness of generalist's performance curve
                double p_mean=           1/(1+skew_spcl);       // Input of the maximum of the Beta shape function, for inversing the shape in mortality dunction (i.e. d_max)
                double x_beta1, x_beta2;
                double B_beta,tmp_beta1, tmp_beta2, alpha_beta, beta_beta;
                // Calculating the maximum mortality for inversing the performance curves (assuming specialists are higher)
                    double d_max;
                        alpha_beta= shape_spcl;
                        beta_beta= shape_spcl*skew_spcl;
                        tmp_beta1= alpha_beta-1;
                        tmp_beta2= beta_beta-1;
                        B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);                    
                    d_max= pow(p_mean,tmp_beta1)*pow(p_mean,tmp_beta2)/B_beta;
    // Temperaure distribution (environmental factors)
        double tmp_env;                                 // temporal storage the sampled long-term variation
        double curr_env;                                // the current environmental condition
        // normal dist (switch: s6)
            double mean_tmp= 	25.0;
            double sd_tmp=	 	5.0;
            double u1, u2;
        // uniform dist (switch: s6)
            double lower= 		5.0;
            double upper= 		45.0;
            double base= 		lower;
            double range= 		upper- lower;
        // Beta dist (switch: s6)
            // scale of environmental variation (longer-term)
            double range_env=           100.0;          // specify the range of environmental conditions to [0:range_env]
            // scale of environmental variation (shorter-term)
            double range_short_env=     60.0;           // specify the range of environmental conditions to [0:range_short_env]
            // average environmental condition
            double mean_env=            50.0;
            //******** environment= mean+ sampled_value -1/2*sample_range ********//
            double shape_long_var=	    100.0;          // Shape of environmental distribution (longer-term)
            double shape_short_var=     100.0;          // Shape of environmental distribution (shorter-term)
            gen_beta_param beta_tmp;	           	    // create the random beta type variable [0:1]
            gen_beta_param beta_short_tmp;
            gen_beta_initialize(&beta_tmp, shape_long_var, shape_long_var); // assuming no skewness
            gen_beta_initialize(&beta_short_tmp, shape_short_var, shape_short_var); 

	// Creating temporal space
		double *p= d_vector(2);
        double *count= d_vector(3);
        double abf_spcl, abf_genl;
        double *K_spcl_within= d_vector(N_event);
        double *K_genl_within= d_vector(N_event);
        double *b_spcl_within= d_vector(N_event);
        double *b_genl_within= d_vector(N_event);
        double *d_spcl_within= d_vector(N_event);
        double *d_genl_within= d_vector(N_event);
        double K_1, K_2, b_1, b_2, d_1, d_2, tmp_x, tmp_y;
	// Output
		FILE *out;
		out= fopen("summary.txt","w");
        fprintf(out,"shape_spcl\tshape_genl\tlong_var\tshort_var\tdom_1\tdom_2\tcoexist\n");
    // Starting to list all param space    
    for (c1=1; c1<=5; c1++){
        shape_long_var=              1000*pow(10, -0.5*(c1-1));
        for (c2=1; c2<= 5; c2++){
            shape_short_var=        1000*pow(10, -0.5*c2);
            gen_beta_initialize(&beta_tmp, shape_long_var, shape_long_var); // assuming symmetric distribution (a=b)
            gen_beta_initialize(&beta_short_tmp, shape_short_var, shape_short_var); // a=b        
            for (c3=1; c3<=3; c3++) count[c3]=0;
            // start
            for (c3=1; c3<= rep; c3++){
            // Initialization
                p[i_spcl]= 250;// initial population of specialists
                p[i_genl]= 250;// initial population of generalists

            // Main loop of time series
            for(j=1; j<=N_gen; j++){   // avoiding misscounting or other errors
                tmp_env= mean_env+ gen_beta(&beta_tmp)*range_env- range_env/2;
                for (i=1; i<=N_event; i++){
                    // Starting a short-term variation
                        curr_env= tmp_env+ gen_beta(&beta_short_tmp)*range_short_env -range_short_env/2;
                        // Calculating parameters from Beta probability distribution function
                            x_beta1= curr_env/range_bio;
                            x_beta2= 1-x_beta1;
                            // Specialist
                                // Parameters of the Bata function
                                alpha_beta= shape_spcl;
                                beta_beta= shape_spcl*skew_spcl;
                                tmp_beta1= alpha_beta-1;
                                tmp_beta2= beta_beta-1;
                                B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);
                                // Calculation
                                K_spcl= scale_K* pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta;
                                b_spcl= scale_b* pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta;
                                d_spcl= scale_d* (d_max-pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta);
                                if(K_spcl!=K_spcl) K_spcl=0;  // hadling the conditions out of range
                                if(b_spcl!=b_spcl) b_spcl=0;  // hadling the conditions out of range
                                if(d_spcl!=d_spcl) d_spcl=0;  // hadling the conditions out of range
                            // Generalist
                                // Parameters of the Beta function
                                alpha_beta= shape_genl;
                                beta_beta= shape_genl*skew_genl;
                                tmp_beta1= alpha_beta-1;
                                tmp_beta2= beta_beta-1;
                                B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);
                                // Calculation
                                K_genl= scale_K* pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta;
                                b_genl= scale_b* pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta; 
                                d_genl= scale_d* (d_max-pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta);
                                if(K_genl!=K_genl) K_genl=0;  // hadling the conditions out of range
                                if(b_genl!=b_genl) b_genl=0;  // hadling the conditions out of range 
                                if(d_genl!=d_genl) d_genl=0;  // hadling the conditions out of range
                    // Define the performance where environmental condition is out of range
                        if (curr_env> range_bio || curr_env<0){
                            K_spcl= K_genl= 0;
                            b_spcl= b_genl= 0;
                            d_spcl= d_genl= scale_d*d_max;
                        }

                    K_spcl_within[i]= K_spcl;
                    K_genl_within[i]= K_genl;
                    b_spcl_within[i]= b_spcl;
                    b_genl_within[i]= b_genl;
                    d_spcl_within[i]= d_spcl;
                    d_genl_within[i]= d_genl;
                }
                K_1= mean(K_spcl_within, N_event);
                K_2= mean(K_genl_within, N_event);
                b_1= mean(b_spcl_within, N_event);
                b_2= mean(b_genl_within, N_event);
                d_1= mean(d_spcl_within, N_event);
                d_2= mean(d_genl_within, N_event);

                tmp_x= p[i_spcl];
                tmp_y= p[i_genl];
                if(K_1>0&& b_1>0) abf_spcl= (1+ b_1*(1- tmp_x/K_1- a_gs*tmp_y/K_1)- d_1);
                else abf_spcl= (1- d_1);
                if(K_2>0&& b_2>0) abf_genl= (1+ b_2*(1- tmp_y/K_2- a_sg*tmp_x/K_2)- d_2);
                else abf_genl= (1- d_2);
                p[i_spcl]= tmp_x*abf_spcl;
                p[i_genl]= tmp_y*abf_genl;
                if(p[i_spcl]<0) p[i_spcl]=0;
                if(p[i_genl]<0) p[i_genl]=0;

                // Coexisting
                if(j/N_gen>0.99) {
                    j= N_gen+1;
                    count[3]+= 1;
                }
                // Specialist wins
                if(p[i_genl]==0 && p[i_spcl]>0){
                    j= N_gen+1;
                    count[1]+= 1;
                }
                // Generalist wins
                if(p[i_spcl]==0 && p[i_genl]>0){
                    j= N_gen+1;
                    count[2]+= 1;
                }
                // Double Extinction
                if(p[i_genl]==0 && p[i_spcl]==0){
                    j= N_gen+1;
                }
            }}
        //print
        fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",shape_spcl,shape_genl,shape_long_var,shape_short_var,count[1],count[2],count[3]);
    }}

	free_d_vector(p);
    free_d_vector(count);
    free_d_vector(K_spcl_within);
    free_d_vector(K_genl_within);
    free_d_vector(b_spcl_within);
    free_d_vector(b_genl_within);
    free_d_vector(d_spcl_within);
    free_d_vector(d_genl_within);
	fclose(out);
	return 0;
}

void message_error(char error_text[]) //standard error handler
{
	printf("There are some errors...\n");
	printf("%s\n",error_text);
	printf("...now existing to system...\n");
	exit(1);
}

double dx_dt(double time, double vr[])	// Solitary
{
    if(K_spcl>0 && b_spcl>0) return (b_spcl- b_spcl/K_spcl*vr[i_spcl]- b_spcl/K_spcl*a_gs*vr[i_genl]- d_spcl)*vr[i_spcl];
    else return -1*d_spcl*vr[i_spcl];
}
double dy_dt(double time, double vr[])	// Cooperative
{
    if(K_genl>0 && b_genl>0) return (b_genl- b_genl/K_genl*vr[i_genl]- b_genl/K_genl*a_sg*vr[i_spcl]- d_genl)*vr[i_genl];
    else return -1*d_genl*vr[i_genl];
}
void differential(double time, double in[], double out[])
{
	out[i_spcl]= dx_dt(time,in);
	out[i_genl]= dy_dt(time,in);
}
void rk4(double p[], double k1[],int n, double t, double h, double pout[],void(*diff)(double,double[],double[]))
{
	int i;
	double tt,*k2,*k3,*k4,*pp;

	k2= d_vector(n);
	k3= d_vector(n);
	k4= d_vector(n);
	pp= d_vector(n);

	for (i=1;i<=n;i++) pp[i]= p[i]+ k1[i]*h/2;
	(*diff)(t+h/2,pp,k2);
	for (i=1;i<=n;i++) pp[i]= p[i]+ k2[i]*h/2;
	(*diff)(t+h/2,pp,k3);
	for (i=1;i<=n;i++) pp[i]= p[i]+ k3[i]*h;
	(*diff)(t+h,pp,k4);
	
	for(i=1;i<=n;i++) pout[i]= p[i]+ (k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*h/6;

	free_d_vector(k2);
	free_d_vector(k3);
	free_d_vector(k4);
	free_d_vector(pp);
}
double *d_vector(long size) 
{
	double *x;
	x= (double *) malloc((size_t)((size+1)*sizeof(double)));
	if(x==NULL) message_error((char *)"Allocation failure in d_vector()");
	return x;
}
double **d_matrix(long size_row, long size_column)
{
	double **x;
	long i;
	long size_row_P= size_row+1;
	long size_column_P= size_column+1;

	x= (double **) malloc((size_t)(size_row_P*sizeof(double *)));               //first dimension
	if (x==NULL) message_error((char *)"Allocation failure in d_vector()");
	x[0]= (double *) malloc((size_t)(size_row_P*size_column_P*sizeof(double))); //second dimension
	if (x[0]==NULL) message_error((char *)"Allocation failure in d_vector()");
	for(i=1;i<size_row_P;i++) x[i]= x[0]+ i*size_column_P;
	return x;
}
void free_d_matrix(double **x)
{
	free(x[0]);
	free(x);
}
void free_d_vector(double *x) {	free(x);}
double normal_dist_BM (double mean, double sd, double u1, double u2)
{
    /* 
     * Using Box-Muller method to generate pseudo-normal distributed numbers in [0,1]
     * Constructed in Feb, 2018
     */
    double z1;
	z1= sqrt(-2* log(u1))* cos(2* M_PI* u2);
	return z1*sd+ mean;
}
double factorial (double x)
{
	int i;
	double tmp, out;
	if(x< 0) return 0;	// not sure
	if(x<= 1) return 1;
	else{
		tmp= x;
		out= 1;
		for(i=1; i<= x; i++){
			out= out*tmp;
			tmp= tmp-1;
		}
		return out;
	}
}
double mean (double x[], int length)
{
	int i;
	double tmp=0.0;
	for (i=1; i<= length; i++){
		tmp+= x[i];
	}
	return tmp/length;
}

