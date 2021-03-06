/*
 * Continual generation model with good/ bad year
 * Focal point is founder effect and the selection of better survival/ better birth species
 * First day of building: 19 July, 2017
 * Transformed to short-long term variability model from 23 Feb, 2018
 * File instruction:
 * 1. Put "gen_beta.h" and "gen_beta.c" in the folder containing this file
 * 2. Put the dsfmt folder and the folder containing this file into the same folder
 ************************************************************************************************

# Execution: (if gcc doesn't work, try "g++")
gcc cdmc_climate.c 
./a.out

# Plot with gnuplot:
gnuplot
plot 'summary.txt' using 3:1 title 'rising-tide' with lines lc rgb 'orange',\
'summary.txt' using 3:2 title 'bet-hedging' with lines lc rgb 'skyblue'

 ************************************************************************************************
 * Key parameters
 * env_in: Input environment file
 * range_bio: The range of environmental condition where biological performance curves are responsive
 * shift_env: Dealing with negative environmental values (e.g. temperature)
 * s5: Determine whether the average fitness of bet-hedging strategy is lower than that of the rising-tide strategy
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>	

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
    double factorial (double x); 

// Global variables
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
        double a_rb= 0.4;           // Intensity of inter-strategic interactions (rising-tide to bet-hedging)
        double a_br= 0.4;           // Intensity of inter-strategic interactions (bet-hedging to rising-tide)
        double gamma_ris=   1.0;    // Scaling coefficient of average fitness (area of performance curve)
        double gamma_bet=   0.9;    // Scaling coefficient of average fitness (area of performance curve)

// Main function
int main (void)
{
    // Switches
        int s4= 1;				                // Switch of terminating tiny populaitons (see ext_thr)
        int s5= 1;                              // Switch of smaller average fitness for bet-hedging strategy (see gamma_bet)
    // Temporal variables
        int T=          1308;                   // Duration of simulation (length of the climate data)
        double deltat=  0.005;                  // Length of time step
        int T_short=    1;                      // The length of short-term variation 
        double m=       20.0;                   // The ratio of spans in sampling new long-term variation to short-term variation
        int T_long=     T_short*m;              // The length of long-term vatiation 
    // Basic variables
		int i,j;                                // For loop counters
		double t=               0.0;            // Time logs
        double T_remain_long=   0.0;            // Counter of remaining time to resample long-term variations
        double T_remain_short=  0.0;            // Counter of remaining time to resample short-term variations
		double pp=              0.0;            // Temp for random number
        double ext_thr=         0.5;            // Threshold of population for terminating (related to s4)
    // Performance curve parameters
            double scale_K=     250.0;          // Scaling coefficient of carrying capacity
            double scale_b=     0.5;            // Scaling coefficient of intrinsic growth rate
            double scale_d=     0.01;           // Scaling coefficient of mortality rate
        // Beta dist
            double range_bio=   50.0;           // The responsive environmental conditions (0-range_bio) of the biological performance curves
            double shape_ris=	50.0;           // Shape of rising-tide strategy's performance curve in carrying capacity (The larger, the narrower)
            double shape_bet=	5.0;            // Shape of bet-hedging strategy's performance curve in carrying capacity
            double skew_ris=    1.0;            // Skewness of rising-tide strategy's performance curve, alpha= skew*shape*2
            double skew_bet=    1.0;            // Skewness of bet-hedging strategy's performance curve
            // skew=1 means no skewness of the performance curve
            double p_mean=      1/(1+skew_ris); // Input of the maximum of the Beta shape function, for inversing the shape in mortality dunction (i.e. d_max)
            double x_beta1, x_beta2;
            double B_beta,tmp_beta1, tmp_beta2, alpha_beta, beta_beta, env_beta;
            // Calculating the maximum mortality for inversing the performance curves (assuming rising-tide strategy is higher)
            double d_max;
                alpha_beta= shape_ris;
                beta_beta= shape_ris*skew_ris;
            d_max= Beta_function(p_mean, alpha_beta, beta_beta);
    // Environmental paramters
        double tmp_env;                         // Temporal storage the sampled long-term variation
        double curr_env;                        // The current environmental condition
        double shift_env=       0.0;            // Shifting the environmental conditions, for negative temperature data
    // Checking the simulation
        int restart= 1;
        int redo_counter= 0;
        int redo_limit= 5;
    // Initialization of real data
        FILE *env_in;
        // env_in= fopen("../climate_files/Toolik_pre.txt","r");
        // env_in= fopen("../climate_files/MtHehuan_tmp.txt","r");
        env_in= fopen("../climate_files/MpalaStation_tmp.txt","r");
        // env_in= fopen("../climate_files/AS_pre.txt","r");
        int length_env= T;
        int seq= 1;
        double *env= d_vector(length_env);
        for(i=1; i<= length_env; i++) fscanf(env_in, "%lf\n", &env[i]);
	// Creating temporal space
		double *p= d_vector(2);
		double *dfdt= d_vector(2);
        double abf_ris, abf_bet;
	// Output
		FILE *out;
		out= fopen("summary.txt","w");
	// Initialization
		p[i_ris]= 250;                          // initial population of rising-tide strategy
		p[i_bet]= 250;                          // initial population of bet-hedging strategy
		fprintf(out,"Nr\tNb\ttime\tenvironment\tabf_r\tabf_b\n");
		fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",p[i_ris],p[i_bet],t,0.0,0.0,0.0);

    // Main loop of time series
    while(restart==1&& redo_counter< redo_limit){   // avoiding misscounting or other errors
        for (i=1; i<=T/deltat; i++){
            // Starting a short-term variation
                if(T_remain_short<= 0){
                    curr_env= env[seq]+ shift_env;
                    seq+= 1;
                    T_remain_short= T_short;
                    // Calculating parameters from Beta probability distribution function
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
                        if(K_ris!=K_ris) K_ris=0;  // hadling the conditions out of range
                        if(b_ris!=b_ris) b_ris=0;  // hadling the conditions out of range
                        if(d_ris!=d_ris) d_ris=0;  // hadling the conditions out of range
                    // Bet-hedging strategy
                        // Parameters of the Beta function
                        alpha_beta= shape_bet;
                        beta_beta= shape_bet*skew_bet;
                        env_beta= Beta_function(x_beta1, alpha_beta, beta_beta);
                        // Calculation
                        K_bet= scale_K* env_beta;
                        b_bet= scale_b* env_beta; 
                        d_bet= scale_d* (d_max-env_beta);
                        if(K_bet!=K_bet) K_bet=0;  // hadling the conditions out of range
                        if(b_bet!=b_bet) b_bet=0;  // hadling the conditions out of range 
                        if(d_bet!=d_bet) d_bet=0;  // hadling the conditions out of range
                    // Different carrying capacity curve
                        if(s5==1){
                            // Parameters of the Beta function
                            alpha_beta= shape_bet;
                            beta_beta= shape_bet*skew_bet;
                            env_beta= Beta_function(x_beta1, alpha_beta, beta_beta);
                            // Calculation
                            K_bet= scale_K* env_beta* 0.9;
                            b_bet= scale_b* env_beta* 0.9; 
                            d_bet= scale_d* (d_max-env_beta* 0.9);
                        }
                }
            // Define the performance where environmental condition is out of range
                if (curr_env> range_bio || curr_env<0){
                    K_ris= K_bet= 0;
                    b_ris= b_bet=0;
                    d_ris= scale_d*d_max;
                    d_bet= scale_d*d_max;
                }
            // Calculating per capita growth rate (fitness)
                if(i%50==0) {
                    if(K_ris>0 && b_ris>0) abf_ris= (b_ris- b_ris/K_ris*p[i_ris]- b_ris/K_ris*a_br*p[i_bet]- d_ris);
                    else abf_ris= -1*d_ris;
                    if(K_bet>0 && b_bet>0) abf_bet= (b_bet- b_bet/K_bet*p[i_bet]- b_bet/K_bet*a_rb*p[i_ris]- d_bet);
                    else abf_bet= -1*d_bet;
                }
            // Calculating the change in population size with time change "deltat"
                differential(t,p,dfdt);
                rk4(p, dfdt, 2, t, deltat, p, differential);
            // Processing time
                T_remain_long-= deltat;
                T_remain_short-= deltat;
                t+= deltat;
            // Record current population
                if(i%50==0) fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",p[i_ris],p[i_bet],t, curr_env, abf_ris, abf_bet);
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

double dx_dt(double time, double vr[])	// Rising-tide strategy
{
    if(K_ris>0 && b_ris>0) return (b_ris- b_ris/K_ris*vr[i_ris]- b_ris/K_ris*a_br*vr[i_bet]- d_ris)*vr[i_ris];
    else return -1*d_ris*vr[i_ris];
}
double dy_dt(double time, double vr[])	// Bet-hedging strategy
{
    if(K_bet>0 && b_bet>0) return (b_bet- b_bet/K_bet*vr[i_bet]- b_bet/K_bet*a_rb*vr[i_ris]- d_bet)*vr[i_bet];
    else return -1*d_bet*vr[i_bet];
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
    /* 
     * Using Box-Muller method to generate pseudo-normal distributed numbers in [0,1]
     * Constructed in Feb, 2018
     */
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
