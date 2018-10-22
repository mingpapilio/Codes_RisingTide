/*
 * Continual generation model with good/ bad year
 * Focal point is founder effect and the selection of better survival/ better birth species
 * First day of building: 19 July, 2017
 * Transformed to short-long term variability model from 23 Feb, 2018
 * File instruction:
 * 1. Put "gen_beta.h" and "gen_beta.c" in the folder containing this file
 * 2. Put the dsfmt folder and the folder containing this file into the same folder
 ************************************************************************************************
Execution:
gcc stoc_cenv_cg.c gen_beta.h gen_beta.c -lm -stdlib=libstdc++
./a.out
 ************************************************************************************************
 * Key parameters
 * T_episode, T_phase: The relative durations of shorter- and longer-term conditions, and the duration of environmental conditions  
 * shape_tmp, shape_short_tmp: The shape parameter of Beta distribution for generating environemntal conditions, both long-term (i.e. shape_tmp) and short-term (i.e. shape_short_tmp). Larger velues indicate narrower distribution. 
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
    void differential_g(double time,double in[],double out[]);
    double dx_dt_g(double time, double vr[]);           //vr is a vector of variable
    double dy_dt_g(double time, double vr[]);
    void differential_b(double time,double in[],double out[]);
    double dx_dt_b(double time, double vr[]);
    double dy_dt_b(double time, double vr[]);
    double normal_dist_BM (double mean, double sd, double u1, double u2);
    double factorial (double x); 

// Global variables
    // Index of variables
        int i_x= 1; //index of x (N1)<- specialist (rising-tide strategy)
        int i_y= 2; //index of y (N2)<- generalist (bet-hedging strategy)
    // differential equation paramters (K, b, d will change according to environmental conditions)
        double b_g1= 0.05;          // Intrinsic growth rate of specialists in good years (here bad years are not used)
        double b_g2= 0.05;          // Intrinsic growth rate of generallists in good years
        double b_b1= 0.05;          // Intrinsic growth rate of specialists in bad years
        double b_b2= 0.05;          // Intrinsic growth rate of generallists in bad years
        double d_g1= 0.01;          // Mortality of specialists in good years
        double d_g2= 0.01;          // Mortality of generalists in good years
        double d_b1= 0.01;          // Mortality of specialists in bad years
        double d_b2= 0.01;          // Mortality of generalists in bad years
        double K_g1= 1000;          // Carrying capacity of specialists in good years
        double K_g2= 1000;          // Carrying capacity of generalists in good years
        double K_b1= 1000;          // Carrying capacity of specialists in bad years
        double K_b2= 1000;          // Carrying capacity of generalists in bad years
        double a_g12= 0.4;          // Intensity of inter-strategic interactions in good years (specialists to generalists)
        double a_g21= 0.4;          // Intensity of inter-strategic interactions in good years (generalists to specialists)
        double a_b12= 1;            // Intensity of inter-strategic interactions in bad years (specialists to generalists)
        double a_b21= 1;            // Intensity of inter-strategic interactions in bad years (generalists to specialists)

// Main function
int main (void)
{
    // Switches and crucial variables
        // Switches
            int s1= 0;				// Switch of good/ bad year *** not used in this model ***
            int s3= 1;				// Switch of temporal stochasticity
            int s4= 1;				// Switch of terminating tiny populaitons
            int s6= 2;              // Specify the type of distribution in environmental conditions
            int s7= 1;              // Specify the type of distribution in biological performance curve
            int s8= 1;              // Switch of making birth and death dependent to environmental conditions
        // Temporal variables
            int T=          1000; 	// Duration of simulation
            double deltat=  0.005; 	// Length of time step
            double r1=      2.0;   	// Rate of getting good season *** not used in this model ***
            // (keep it larger than 1 to avoid using the equations of bad years)
            int T_episode=  20;     // Good season length ***Note the difference of double and int
            int T_phase=    1;      // The length of short-term variation 

    // Basic variables
		int i,j;                    // For loop counters
		double t=              0.0; // Time logs
        double T_remain_long=  0.0; // Counter of remaining season length (entire episode)
        double T_remain_short= 0.0; // Counter of remaining time within one short-term condition
		double pp=             0.0;	// Temp for random number
        double ext_thr=        0.5; // Threshold of population for terminating (refer to s4)
    // performance curve parameters
		// performance curve (how carrying capacity changes in resonpse to environmental condition)
                double scale_K=     250;                        // Scaling coefficient of carrying capacity
                double scale_b=     0.5;                        // Scaling coefficient of intrinsic growth rate
                double scale_d=     1E-6;                       // Scaling coefficient of mortality rate
            // normal dist
                double mean_spcl=	25;                         // Optimal condition of specialist strategy
                double sd_spcl=     5;                          // S.D. of the performance curve
                double mean_genl=   25;                         // Optimal condition of generalist strategy
                double sd_genl=     50;                         // S.D. of the performance curve
            // beta dist
                double shape_spcl=	     30.0;                  // Shape of specialist's performance curve in carrying capacity (The larger, the narrower)
                double shape_genl=	     4.5;                   // Shape of generalist's performance curve in carrying capacity
                double shape_birth_spcl= 30.0;                  // Shape of specialist's performance curve in intrinsic growth rate
                double shape_birth_genl= 4.5;                   // Shape of generalist's performance curve in intrinsic growth rate
                double shape_death_spcl= 30.0;                  // Shape of specialist's performance curve in mortality
                double shape_death_genl= 4.5;                   // Shape of generalist's performance curve in mortality
                double skew_spcl=        0.333;                 // Skewness of specialist's performance curve, alpha= skew*shape*2
                double skew_genl=        0.333;                 // Skewness of generalist's performance curve
                double shift_death=      0.0;                   // Tuing the difference in mortality between two strategies
                double p_mean=           1/(1+skew_spcl);       // Input of the maximum of the Beta shape function, for inversing the shape in mortality dunction (i.e. d_max)
                double x_beta1, x_beta2;
                double B_beta,tmp_beta1, tmp_beta2, alpha_beta, beta_beta;
                // Calculating the maximum mortality for inversing the performance curves (assuming specialists are higher)
                    double d_max;
                        alpha_beta= shape_death_spcl;
                        beta_beta= shape_death_spcl*skew_spcl;
                        tmp_beta1= alpha_beta-1;
                        tmp_beta2= beta_beta-1;
                        B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);                    
                    d_max= pow(p_mean,tmp_beta1)*pow(p_mean,tmp_beta2)/B_beta;
        // temperaure distribution (environmental factors)
            double tmp_tmp;                                 // temporal storage the mean environmental condition of an episode
            double curr_tmp;                                // the current environmental condition of a phase
            // normal dist
                double mean_tmp= 	25;
                double sd_tmp=	 	5;
                double u1, u2;
            // uniform dist
                double lower= 		5;
                double upper= 		45;
                double base= 		lower;
                double range= 		upper- lower;
			// beta dist
                // scale of environmental variation (episode level, longer-term)
                double range_env=           100.0;          // specify the range of environmental conditions to [0:range_env]
                // scale of environmental variation (phase level, shorter-term)
                double range_short_beta=    60.0;           // specify the range of environmental conditions to [0:range_short_beta]
                // average environmental condition
                double mean_beta=           75.0;
                //******** K= mean+ sampled_value -1/2*sample_range ********//
                double shape_tmp=	        100.0;          // Shape of environmental distribution (longer-term)
                double shape_short_tmp=     1000.0;         // Shape of environmental distribution (shorter-term)
                gen_beta_param beta_tmp;	           	    // create the random beta type variable [0:1]
                gen_beta_param beta_short_tmp;
				gen_beta_initialize(&beta_tmp, shape_tmp, shape_tmp); // assuming symmetric distribution (a=b)
                gen_beta_initialize(&beta_short_tmp, shape_short_tmp, shape_short_tmp); // a=b
    // Checking the simulation
        int restart= 1;
        int redo_counter= 0;
        int redo_limit= 100;
	// Creating temporal space
		double *p= d_vector(2);
		double *dfdt= d_vector(2);
	// Output
		FILE *out;
		out= fopen("out.txt","w");
	// Initialization
		p[i_x]= 250;// initial population of specialists
		p[i_y]= 250;// initial population of generalists
		//fprintf(out,"r1 is %lf, r2 is %lf, db is %lf, duration is %d\n",b_g1,b_g2,d_b2,T);
		fprintf(out,"N1\t\tN2\t\ttime\t\tenvironment\n");
		fprintf(out,"%lf\t%lf\t%lf\t%lf\n",p[i_x],p[i_y],t,0.0);

    // Main loop of time series
    while(restart==1&& redo_counter< redo_limit){   // avoiding misscounting or other errors
        for (i=1; i<=T/deltat; i++){
            // Starting an episode
            if (T_remain_long<=0){
                if(s3==1){                          // Stochastic model
                    s1= 1;		                    // starting of an episode
                    T_remain_long= T_episode;
                    // performance curve initialization
                    // Calculating mean environmental state of the episode (i.e. tmp_tmp)
                    // Calculating current environmental state (i.e. curr_tmp)
                    if (s6==2){                     // Beta distribution
                        tmp_tmp= mean_beta+ gen_beta(&beta_tmp)*range_env- range_env/2;
                        curr_tmp= tmp_tmp+ gen_beta(&beta_short_tmp)*range_short_beta- range_short_beta/2;
                        T_remain_short= T_phase;
                    }
                    if (s6==1){                     // Normal distribution
                        u1= gen_unif(&beta_tmp);
                        u2= gen_unif(&beta_tmp);
                        curr_tmp= normal_dist_BM (mean_tmp, sd_tmp, u1, u2);
                    }
                    if (s6==0){                     // Uniform distribution
                        curr_tmp= base+ range*gen_unif(&beta_tmp);
                    }
                }
                if(s3==0){                          // Deterministic model
                    if(s1== 0||s1== 2){
                        s1= 1;		                // Good year starts
                        T_remain_long= T_episode*2*r1;}
                    else{
                        s1= 2;		                // Bad year starts
                        T_remain_long= T_episode*(2-2*r1);}
                }
            }
            // Starting a phase
            if(T_remain_short<= 0){
                curr_tmp= tmp_tmp+ gen_beta(&beta_short_tmp)*range_short_beta -range_short_beta/2;
                T_remain_short= T_phase;
                if (s7==1){                     // beta distribution
                    x_beta1= curr_tmp/range_env;
                    x_beta2= 1-x_beta1;
                        alpha_beta= shape_spcl;
                        beta_beta= shape_spcl*skew_spcl;
                        tmp_beta1= alpha_beta-1;
                        tmp_beta2= beta_beta-1;
                        B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);
                    K_g1= scale_K* pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta;
                    if(K_g1!=K_g1) K_g1=0;      // hadling the conditions out of range
                    if(s8==1){
                            alpha_beta= shape_birth_spcl;
                            beta_beta= shape_birth_spcl*skew_spcl;
                            tmp_beta1= alpha_beta-1;
                            tmp_beta2= beta_beta-1;
                            B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);
                        b_g1= scale_b* pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta;
                        if(b_g1!=b_g1) b_g1=0;  // hadling the conditions out of range
                            alpha_beta= shape_death_spcl;
                            beta_beta= shape_death_spcl*skew_spcl;
                            tmp_beta1= alpha_beta-1;
                            tmp_beta2= beta_beta-1;
                            B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);
                        d_g1= scale_d* (d_max-pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta)+ shift_death;
                        if(d_g1!=d_g1) d_g1=0;  // hadling the conditions out of range
                    }
                        alpha_beta= shape_genl;
                        beta_beta= shape_genl*skew_genl;
                        tmp_beta1= alpha_beta-1;
                        tmp_beta2= beta_beta-1;
                        B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);
                    K_g2= scale_K* pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta;
                    if(K_g2!=K_g2) K_g2=0;      // hadling the conditions out of range
                    if(s8==1){
                            alpha_beta= shape_birth_genl;
                            beta_beta= shape_birth_genl*skew_genl;
                            tmp_beta1= alpha_beta-1;
                            tmp_beta2= beta_beta-1;
                            B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);               
                        b_g2= scale_b* pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta;
                        if(b_g2!=b_g2) b_g2=0;  // hadling the conditions out of range
                            alpha_beta= shape_death_genl;
                            beta_beta= shape_death_genl*skew_genl;
                            tmp_beta1= alpha_beta-1;
                            tmp_beta2= beta_beta-1;
                            B_beta= tgamma(alpha_beta)*tgamma(beta_beta)/tgamma(alpha_beta+beta_beta);    
                        d_g2= scale_d* (d_max-pow(x_beta1,tmp_beta1)*pow(x_beta2,tmp_beta2)/B_beta);
                        if(d_g2!=d_g2) d_g2=0;  // hadling the conditions out of range
                    }
                }
            }
            // Define the performance where environmental condition is out of range
            if (curr_tmp> range_env || curr_tmp<0){
                if(s7==1){
                    K_g1= K_g2= 0;
                    if(s8==1){
                        b_g1= b_g2=0;
                        d_g1= scale_d*d_max+ shift_death;
                        d_g2= scale_d*d_max;
                    }
                }
            }
            // Calculating the change in population size with time change "deltat"
            if(s1==1){
                differential_g(t,p,dfdt);
                rk4(p, dfdt, 2, t, deltat, p, differential_g);}
            if(s1==2){
                differential_b(t,p,dfdt);
                rk4(p, dfdt, 2, t, deltat, p, differential_b);}
            T_remain_long-= deltat;
            T_remain_short-= deltat;
            t+= deltat;
            // Record current population
            if(i%50==0) fprintf(out,"%lf\t%lf\t%lf\t%lf\n",p[i_x],p[i_y],t, curr_tmp);
            // extinction of population
            if(s4==1){
                if(p[i_x]<ext_thr) p[i_x]=0;
                if(p[i_y]<ext_thr) p[i_y]=0;
            }
            // Warning message
            if(p[i_x]<0||p[i_y]<0) {
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
    if (redo_counter>= redo_limit) printf("Serious termination problem after 100 attempts.\n");
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

double dx_dt_g(double time, double vr[])	// Solitary
{
    if(K_g1>0 && b_g1>0) return (b_g1- b_g1/K_g1*vr[i_x]- b_g1/K_g1*a_g21*vr[i_y]- d_g1)*vr[i_x];
    else return -1*d_g1*vr[i_x];
}
double dy_dt_g(double time, double vr[])	// Cooperative
{
    if(K_g2>0 && b_g2>0) return (b_g2- b_g2/K_g2*vr[i_y]- b_g2/K_g2*a_g12*vr[i_x]- d_g2)*vr[i_y];
    else return -1*d_g2*vr[i_y];
}
void differential_g(double time, double in[], double out[])
{
	out[i_x]= dx_dt_g(time,in);
	out[i_y]= dy_dt_g(time,in);
}
double dx_dt_b(double time, double vr[])
{return (b_b1- b_b1/K_b1*vr[i_x]- b_b1/K_b1*a_b21*vr[i_y]- d_b1)*vr[i_x];}
double dy_dt_b(double time, double vr[])
{return (b_b2- b_b2/K_b2*vr[i_y]- b_b2/K_b2*a_b12*vr[i_x]- d_b2)*vr[i_y];}
void differential_b(double time, double in[], double out[])
{
	out[i_x]= dx_dt_b(time,in);
	out[i_y]= dy_dt_b(time,in);
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
	if(x==NULL) message_error("Allocation failure in d_vector()");
	return x;
}
double **d_matrix(long size_row, long size_column)
{
	double **x;
	long i;
	long size_row_P= size_row+1;
	long size_column_P= size_column+1;

	x= (double **) malloc((size_t)(size_row_P*sizeof(double *)));               //first dimension
	if (x==NULL) message_error("Allocation failure in d_vector()");
	x[0]= (double *) malloc((size_t)(size_row_P*size_column_P*sizeof(double))); //second dimension
	if (x[0]==NULL) message_error("Allocation failure in d_vector()");
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
