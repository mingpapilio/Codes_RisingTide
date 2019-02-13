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
        int s4= 1;				    // Switch of terminating tiny populaitons
        int s6= 2;                  // Specify the type of distribution in environmental conditions
        int s9= 0;                  // Switch of fixing K, this overwrites any calculated K
    // Temporal variables
        int T=          20000; 	    // Duration of simulation
        double deltat=  0.005; 	    // Length of time step
        int T_short=    1;          // The length of short-term variation 
        double m=       20.0;       // The ratio of spans in sampling new long-term variation to short-term variation
        int T_long=     T_short*m;  // The length of long-term vatiation 
    // Basic variables
		int i,j;                        // For loop counters
		double t=              0.0;     // Time logs
        double T_remain_long=  0.0;     // Counter of remaining time to resample long-term variations
        double T_remain_short= 0.0;     // Counter of remaining time to resample short-term variations
		double pp=             0.0;	    // Temp for random number
        double ext_thr=        0.5;     // Threshold of population for terminating (related to s4)
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
    // Checking the simulation
        int restart= 1;
        int redo_counter= 0;
        int redo_limit= 100;
	// Creating temporal space
		double *p= d_vector(2);
		double *dfdt= d_vector(2);
        double abf_spcl, abf_genl;
	// Output
		FILE *out;
		out= fopen("summary.txt","w");
	// Initialization
		p[i_spcl]= 250;// initial population of specialists
		p[i_genl]= 250;// initial population of generalists
		fprintf(out,"Nr\tNb\ttime\tenvironment\tabf_r\tabf_b\n");
		fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",p[i_spcl],p[i_genl],t,0.0,0.0,0.0);

    // Main loop of time series
    while(restart==1&& redo_counter< redo_limit){   // avoiding misscounting or other errors
        for (i=1; i<=T/deltat; i++){
            // Starting a long-term variation
                if (T_remain_long<=0){
                    T_remain_long= T_long;
                    if (s6==2){                     // Beta distribution
                        tmp_env= mean_env+ gen_beta(&beta_tmp)*range_env- range_env/2;
                        curr_env= tmp_env+ gen_beta(&beta_short_tmp)*range_short_env- range_short_env/2;
                        T_remain_short= T_short;
                    }
                    if (s6==1){                     // Normal distribution
                        u1= gen_unif(&beta_tmp);
                        u2= gen_unif(&beta_tmp);
                        curr_env= normal_dist_BM (mean_tmp, sd_tmp, u1, u2);
                    }
                    if (s6==0){                     // Uniform distribution
                        curr_env= base+ range*gen_unif(&beta_tmp);
                    }
                }
            // Starting a short-term variation
                if(T_remain_short<= 0){
                    curr_env= tmp_env+ gen_beta(&beta_short_tmp)*range_short_env -range_short_env/2;
                    T_remain_short= T_short;
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
                }
            // Define the performance where environmental condition is out of range
                if (curr_env> range_bio || curr_env<0){
                    K_spcl= K_genl= 0;
                    b_spcl= b_genl=0;
                    d_spcl= scale_d*d_max;
                    d_genl= scale_d*d_max;
                }
                if(s9==1){
                    K_spcl= 1000.0;
                    K_genl= 1000.0;
                }
            // Calculating per capita growth rate (fitness)
                if(i%50==0) {
                    if(K_spcl>0 && b_spcl>0) abf_spcl= (b_spcl- b_spcl/K_spcl*p[i_spcl]- b_spcl/K_spcl*a_gs*p[i_genl]- d_spcl);
                    else abf_spcl= -1*d_spcl;
                    if(K_genl>0 && b_genl>0) abf_genl= (b_genl- b_genl/K_genl*p[i_genl]- b_genl/K_genl*a_sg*p[i_spcl]- d_genl);
                    else abf_genl= -1*d_genl;
                }
            // Calculating the change in population size with time change "deltat"
                differential(t,p,dfdt);
                rk4(p, dfdt, 2, t, deltat, p, differential);
            // Processing time
                T_remain_long-= deltat;
                T_remain_short-= deltat;
                t+= deltat;
            // Record current population
                if(i%50==0) fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",p[i_spcl],p[i_genl],t, curr_env, abf_spcl, abf_genl);
            // extinction of population
                if(s4==1){
                    if(p[i_spcl]<ext_thr) p[i_spcl]=0;
                    if(p[i_genl]<ext_thr) p[i_genl]=0;
                }
            // Warning message
                if(p[i_spcl]<0||p[i_genl]<0) {
                    printf("Dynamics crashed.\n");
                    restart=1;
                    redo_counter+= 1;                
                }
        }
        if(t/T>0.99) restart= 0;
        else {
            printf("Simulation is not completed, t is %lf.\n",t);
            restart=1;
            redo_counter+= 1;
        }
    }
    if (redo_counter>= redo_limit) printf("Serious termination problem after %d attempts.\n", redo_limit);

	free_d_vector(p);
	free_d_vector(dfdt);
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
