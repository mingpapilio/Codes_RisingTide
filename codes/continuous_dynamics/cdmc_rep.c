/*
 * Continuous populaiton dynamics model
 * This file counts the outcome of repeated simulations for the contour plots
 * File instruction:
 * 1. Put "gen_beta.h" and "gen_beta.c" in the folder containing this file
 * 2. Put the dsfmt folder and the folder containing this file into the same folder
 ************************************************************************************************
 # Execution (if gcc doesn't work, try "g++")
gcc cdmc_rep.c gen_beta.h gen_beta.c -lm
./a.out
 # Plot with gnuplot
gnuplot
plot 'summary.txt' using 3:1 title 'rising-tide' with lines lc rgb 'orange',\
'summary.txt' using 3:2 title 'bet-hedging' with lines lc rgb 'skyblue'

 ************************************************************************************************
 * Key parameters
 * s1, s2: Determine the type of distribution function in biological response and environmental variation
 * s3: Determine whether intrinsic growth rate and carrying capacity respond to environment differently
 * s5: Determine whether the average fitness of bet-hedging strategy is lower than that of the rising-tide strategy
 * mean_env: The average environmental condition
 * shape_long_var, shape_short_var: The shape parameter of Beta distribution for generating environemntal conditions, both long-term and short-term. Larger velues indicate narrower distribution. 
 * shape_ris, shape_bet: the shape of performance curves, higher values bring narrower curves (Beta function, s1==2)
 * T_long, T_short: The relative durations of short- and long-term conditions (T_long= T_short*m)  
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
    double normal_dist_BM (double mean, double sd, double u1, double u2);           // Generates a random number from a normal distribution
    double Beta_function (double x_beta1, double alpha_beta, double beta_beta);     // Gives the density of Beta distribution

// Global variables
    // Switches
        int s1= 2;                  // Specify the type of performance curve (1: normal distribution, 2: beta distribution)
        int s2= 2;                  // Specify the type of distribution in environmental conditions  (1: normal distribution, 2: beta distribution)
        int s3= 1;                  // Whether intrinsic growth rate, b, and carrying capacity, K, respond to environment differently, only applicable to beta distribution performance curves (1: both b and K follow rising-tide/ bet-hedging shape settings, 2: only b follows, 3: only K follows, requires s1==2)
        int s4= 1;				    // Switch of terminating tiny populaitons (see ext_thr)
        int s5= 1;                  // Switch of smaller average fitness for bet-hedging strategy (see gamma_bet)
    // Index of variables
        int i_ris= 1;               // Index of the rising-tide strategy
        int i_bet= 2;               // Index of the bet-hedging strategy
    // differential equation paramters (K, b, d will change according to environmental conditions)
        double b_ris;               // Intrinsic growth rate of rising-tide strategy
        double b_bet;               // Intrinsic growth rate of bet-hedging strategy
        double d_ris;               // Environment-dependent mortality of rising-tide strategy
        double d_bet;               // Environment-dependent mortality of bet-hedging strategy
        double K_ris;               // Carrying capacity of rising-tide strategy
        double K_bet;               // Carrying capacity of bet-hedging strategy
        double a_rb=        0.4;    // Intensity of inter-strategic interactions (rising-tide to bet-hedging)
        double a_br=        0.4;    // Intensity of inter-strategic interactions (bet-hedging to rising-tide)
        double gamma_ris=   1.0;    // Scaling coefficient of average fitness (area of performance curve)
        double gamma_bet=   0.9;    // Scaling coefficient of average fitness (area of performance curve)
        double scale_K=     250.0;  // Scaling coefficient of carrying capacity
        double scale_b=     0.5;    // Scaling coefficient of intrinsic growth rate
        double scale_d=     0.01;   // Scaling coefficient of mortality rate

// Main function
int main (void)
{
    // Temporal variables
        int T=          20000;          // Duration of simulation
        double deltat=  0.005; 	        // Length of time step
        int T_short=    1;              // The length of short-term variation 
        double m=       20.0;           // The ratio of spans in sampling new long-term variation to short-term variation
        int T_long=     T_short*m;      // The length of long-term vatiation 
    // Basic variables
		int i,j,k;                      // For loop counters
		double t=              0.0;     // Time logs
        double T_remain_long=  0.0;     // Counter of remaining time to resample long-term variations
        double T_remain_short= 0.0;     // Counter of remaining time to resample short-term variations
		double pp=             0.0;	    // Temp for random number
        double ext_thr=        0.5;     // Threshold of population for terminating (related to s4)
        int rep=               100; // Number of repetition
    // Performance curve parameters
        // Normal distribution (s1==1)
            double mean_ris=    50.0;           // Optimal environmental condition of the rising-tide strategy (average of the distribution)
            double mean_bet=    50.0;           // Optimal environmental condition of the bet-hedging strategy (average of the distribution)
            double sd_ris=      1.0;            // Width parameter of the performance of the rising-tide strategy (sd of the distribution)
            double sd_bet=      10.0;           // Width parameter of the performance of the bet-hedging strategy (sd of the distribution)
            // Calculating the maximum mortality for inversing the performance curves (assuming specialists are higher)
            double d_max;
            double env_norm= exp(0)/sqrt(2*M_PI*sd_ris*sd_ris);
            d_max= env_norm;
        // Beta distribution (s1==2)
            double range_bio=   100.0;          // The responsive environmental conditions (0-range_bio) of the biological performance curves
            double shape_ris=   50.0;           // Shape of rising-tide strategy's performance curve in carrying capacity (The larger, the narrower)
            double shape_bet=   5.0;            // Shape of bet-hedging strategy's performance curve in carrying capacity
            double skew_ris=    1.0;            // Skewness of rising-tide strategy's performance curve, alpha= skew*shape*2
            double skew_bet=    1.0;            // Skewness of bet-hedging strategy's performance curve
            // skew=1 means no skewness of the performance curve
            double p_mean=      1/(1+skew_ris); // Input of the maximum of the Beta shape function, for inversing the shape in mortality dunction (i.e. d_max)
            double x_beta1, x_beta2;
            double B_beta,tmp_beta1, tmp_beta2, alpha_beta, beta_beta, env_beta;
            // Calculating the maximum mortality for inversing the performance curves (assuming specialists are higher)
                alpha_beta=     shape_ris;
                beta_beta=      shape_ris*skew_ris;
            d_max= Beta_function(p_mean, alpha_beta, beta_beta);
    // Environmental paramters
        double tmp_env;                         // Temporal storage the sampled long-term variation
        double curr_env;                        // The current environmental condition
        // Normal distribution (s2==1)
            double mean_tmp= 	        50.0;   // Average environmental condition
            double sd_tmp=	 	        5.0;    // Size of long-term variation (standard deviation of the distribution)
            double sd_short_tmp=        5.0;    // Size of short-term variation
            double u1, u2;
        // Beta distribution (s2==2)
            // average environmental condition
            double mean_env=            50.0;   // Average environmental condition
            // range of environmental variation (long-term)
            double range_env=           100.0;  // Specify the range of environmental conditions to [0:range_env]
            // range of environmental variation (short-term)
            double range_short_env=     60.0;   // Specify the range of environmental conditions to [0:range_short_env]
            //******** environment= mean+ sampled_value -1/2*sample_range ********//
            double shape_long_var=	    100.0;  // Shape of environmental distribution (long-term)
            double shape_short_var=     100.0;  // Shape of environmental distribution (short-term)
            gen_beta_param beta_tmp;            // Create the beta type variable
            gen_beta_param beta_short_tmp;      // Create the beta type variable
            gen_beta_initialize(&beta_tmp, shape_long_var, shape_long_var); // assuming environmental distribution has no skewness
            gen_beta_initialize(&beta_short_tmp, shape_short_var, shape_short_var); 
    // Checking the simulation
        int restart= 1;
        int redo_counter= 0;
        int redo_limit= 5;
	// Creating temporal space
		double *p= d_vector(2);
		double *dfdt= d_vector(2);
        double *count= d_vector(3);
	// Output
		FILE *out;
		out= fopen("summary.txt","w");
	// Initialization
		p[i_ris]= 250;                          // initial population of rising-tide strategy
		p[i_bet]= 250;                          // initial population of bet-hedging strategy
        if (s1==2) fprintf(out,"shape_ris\tshape_bet\tlong_var\tshort_var\tdom_1\tdom_2\tcoexist\n");
        if (s1==1) fprintf(out,"sd_ris\tsd_bet\tlong_var\tshort_var\tdom_1\tdom_2\tcoexist\n");
    // Start 
    for(i=1; i<=5; i++){
            if(s2==2) shape_long_var=       1000*pow(10, -0.5*(i-1));
            if(s2==1) sd_tmp=               0.5+ 5*(i-1);
        for(j=1; j<=5; j++){
                if(s2==2) shape_short_var=  1000*pow(10, -0.5*j);
                if(s2==1) sd_short_tmp=     0.5+ 5*(j-1);
                // reset
                gen_beta_initialize(&beta_tmp, shape_long_var, shape_long_var);       // assuming environmental distribution has no skewness
                gen_beta_initialize(&beta_short_tmp, shape_short_var, shape_short_var); 
            for (k=1; k<=3; k++) count[k]= 0;
            for (k=1; k<= rep; k++){
                    restart= 1;
                // Main loop of time series
                while(restart==1&& redo_counter< redo_limit){               // avoiding misscounting or other errors
                    // Reset
                    t=0;
                    p[i_ris]= 250;
                    p[i_bet]= 250;
                    while (t<= T){
                        // Starting a long-term variation
                            if (T_remain_long<=0){
                                if (s2==2){                     
                                    // Beta distribution
                                    tmp_env= mean_env+ gen_beta(&beta_tmp)*range_env- range_env/2;
                                }
                                if (s2==1){                     
                                    // Normal distribution
                                    u1= gen_unif(&beta_tmp);
                                    u2= gen_unif(&beta_tmp);
                                    tmp_env= normal_dist_BM (mean_tmp, sd_tmp, u1, u2);
                                }
                                T_remain_long= T_long;
                            }
                        // Starting a short-term variation
                            if(T_remain_short<= 0){
                                if(s2==2){
                                    curr_env= tmp_env+ gen_beta(&beta_short_tmp)*range_short_env -range_short_env/2;
                                }
                                if(s2==1){
                                    u1= gen_unif(&beta_tmp);
                                    u2= gen_unif(&beta_tmp);
                                    curr_env= normal_dist_BM (tmp_env, sd_short_tmp, u1, u2);                                    
                                }
                                T_remain_short= T_short;
                                // Updating the biological parameter values
                                    if(s1==2){
                                    // Beta distribution performance curve
                                    x_beta1= curr_env/range_bio;
                                    x_beta2= 1-x_beta1;
                                    // Rising-tide strategy
                                        // Parameters of the Bata function
                                        alpha_beta= shape_ris;
                                        beta_beta= shape_ris*skew_ris;
                                        env_beta= Beta_function(x_beta1, alpha_beta, beta_beta);
                                        // Calculation
                                        K_ris= scale_K* env_beta;
                                        b_ris= scale_b* env_beta;
                                        d_ris= scale_d* (d_max-env_beta);
                                        if(K_ris!=K_ris) K_ris=0;  // hadling the values out of the range of Beta function
                                        if(b_ris!=b_ris) b_ris=0;  // hadling the values out of the range of Beta function
                                        if(d_ris!=d_ris) d_ris=0;  // hadling the values out of the range of Beta function
                                    // Bet-hedging strategy
                                        // Parameters of the Beta function
                                        alpha_beta= shape_bet;
                                        beta_beta= shape_bet*skew_bet;
                                        env_beta= Beta_function(x_beta1, alpha_beta, beta_beta);
                                        // Calculation
                                        K_bet= scale_K* env_beta;
                                        b_bet= scale_b* env_beta; 
                                        d_bet= scale_d* (d_max-env_beta);
                                        if(K_bet!=K_bet) K_bet=0;  // hadling the values out of the range of Beta function
                                        if(b_bet!=b_bet) b_bet=0;  // hadling the values out of the range of Beta function
                                        if(d_bet!=d_bet) d_bet=0;  // hadling the values out of the range of Beta function
                                    // Different mean fitness
                                        if(s5==1){
                                            K_bet= scale_K* env_beta* gamma_bet;
                                            b_bet= scale_b* env_beta* gamma_bet;
                                            d_bet= scale_d* (d_max-env_beta* gamma_bet);
                                        }
                                    // The two strategy only differ in reproduction response, making the carrying capacities equal
                                        if(s3==2){
                                            // Parameters of the Bata function
                                            alpha_beta= shape_bet;
                                            beta_beta= shape_bet*skew_bet;
                                            env_beta= Beta_function(x_beta1, alpha_beta, beta_beta);
                                            // Calculation
                                            K_ris= scale_K* env_beta;
                                            K_bet= K_ris;
                                        }
                                    // The two strategy only differ in survival response, making the birth rates equal
                                        if(s3==3){
                                            // Parameters of the Bata function
                                            alpha_beta= shape_bet;
                                            beta_beta= shape_bet*skew_bet;
                                            env_beta= Beta_function(x_beta1, alpha_beta, beta_beta);
                                            // Calculation
                                            b_ris= scale_b* env_beta;
                                            b_bet= b_ris;
                                        }
                                    }
                                    if(s1==1){
                                    // Normal distributon performance curve
                                        // Rising-tide strategy
                                        env_norm= exp(-1*(curr_env-mean_ris)*(curr_env-mean_ris)/2/sd_ris/sd_ris)/sqrt(2*M_PI*sd_ris*sd_ris);
                                        K_ris= scale_K*env_norm;
                                        b_ris= scale_b*env_norm;
                                        d_ris= scale_d*(d_max-env_norm);
                                        // Bet-hedging strategy
                                        env_norm= exp(-1*(curr_env-mean_bet)*(curr_env-mean_bet)/2/sd_bet/sd_bet)/sqrt(2*M_PI*sd_bet*sd_bet);
                                        K_bet= scale_K* env_norm;
                                        b_bet= scale_b* env_norm; 
                                        d_bet= scale_d* (d_max-env_norm);
                                        if(s5==1){
                                            K_bet= scale_K* env_norm* gamma_bet;
                                            b_bet= scale_b* env_norm* gamma_bet;
                                            d_bet= scale_d* (d_max- env_norm* gamma_bet);
                                        }
                                    }
                            }
                        // Define the performance where environmental condition is out of range
                            if(s1==2){
                                if (curr_env> range_bio || curr_env<0){
                                    K_ris= K_bet= 0;
                                    b_ris= b_bet=0;
                                    d_ris= scale_d*d_max;
                                    d_bet= scale_d*d_max;
                                }
                            }
                        // Calculating the change in population size with time change "deltat"
                            differential(t,p,dfdt);
                            rk4(p, dfdt, 2, t, deltat, p, differential);
                        // Processing time
                            T_remain_long-= deltat;
                            T_remain_short-= deltat;
                            t+= deltat;
                        // extinction of population
                            if(s4==1){
                                if(p[i_ris]<ext_thr) p[i_ris]=0;
                                if(p[i_bet]<ext_thr) p[i_bet]=0;
                            }
                        // Warning message
                            if(p[i_ris]<0||p[i_bet]<0) {
                                printf("Dynamics crashed.\n");
                                restart=1;
                                redo_counter+= 1;                
                            }
                        // Rising-tide strategy wins
                            if(p[i_bet]==0 && p[i_ris]>0){
                                restart= 0;
                                t= T+0.1;
                                count[1]+= 1;
                            }
                        // Bet-hedging strategy wins
                            if(p[i_ris]==0 && p[i_bet]>0){
                                restart= 0;
                                t= T+ 0.1;
                                count[2]+= 1;
                            }
                        // Double Extinction
                            if(p[i_bet]==0 && p[i_ris]==0){
                                restart= 0;
                                t= T+ 0.1;
                            }
                        // Coexisting (polymorphism)
                            if(t/T>0.99 && restart==1) {
                                restart= 0;
                                t= T+ 0.1;
                                count[3]+= 1;
                            }
                    }
                    // Checking
                    if(t/T<= 0.99){
                        printf("Simulation is not completed, t is %lf.\n",t);
                        restart=1;
                        redo_counter+= 1;
                    }
                }
                if (redo_counter>= redo_limit) printf("Serious termination problem after %d attempts.\n", redo_limit);
            }
            // Print outcome
            if (s1==2) fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",shape_ris,shape_bet,shape_long_var,shape_short_var,count[1],count[2],count[3]);
            if (s1==1) fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",sd_ris,sd_bet,sd_tmp,sd_short_tmp,count[1],count[2],count[3]);
    }}

	free_d_vector(p);
	free_d_vector(dfdt);
    free_d_vector(count);
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

double dx_dt(double time, double vr[])	// Rising-tide strategy
{
    if(s3==1) return (b_ris- scale_b/scale_K*vr[i_ris]- scale_b/scale_K*a_br*vr[i_bet]- d_ris)*vr[i_ris];
    else{
        if(K_ris>0) return (b_ris- b_ris/K_ris*vr[i_ris]- b_ris/K_ris*a_br*vr[i_bet]- d_ris)*vr[i_ris];
        else return (b_ris-d_ris)*vr[i_ris];
    }
}
double dy_dt(double time, double vr[])	// Bet-hedging strategy
{
    if(s3==1) return (b_bet- scale_b/scale_K*vr[i_bet]- scale_b/scale_K*a_rb*vr[i_ris]- d_bet)*vr[i_bet];
    else{
        if(K_bet>0) return (b_bet- b_bet/K_bet*vr[i_bet]- b_bet/K_bet*a_rb*vr[i_ris]- d_bet)*vr[i_bet];
        else return (b_bet-d_bet)*vr[i_bet];
    }
}
void differential(double time, double in[], double out[])
{
	out[i_ris]= dx_dt(time,in);
	out[i_bet]= dy_dt(time,in);
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
    // Using Box-Muller method to generate pseudo-normal distributed numbers in [0,1]
    double z1;
	z1= sqrt(-2* log(u1))* cos(2* M_PI* u2);
	return z1*sd+ mean;
}
double Beta_function (double x_beta1, double alpha_beta, double beta_beta)
{
    double x_beta2, tmp_beta1, tmp_beta2, logAB_beta, logA_beta, logB_beta;
    x_beta2= 1- x_beta1;
    tmp_beta1= alpha_beta-1;
    tmp_beta2= beta_beta-1;
    logAB_beta= lgamma(alpha_beta)+lgamma(beta_beta)-lgamma(alpha_beta+beta_beta);
    logA_beta= tmp_beta1*log(x_beta1);
    logB_beta= tmp_beta2*log(x_beta2);
    return exp(logA_beta+logB_beta-logAB_beta);
}