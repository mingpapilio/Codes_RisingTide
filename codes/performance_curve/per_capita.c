/*
 * Per capita growth rate of rising-tide and bet-hedging strategy across the environmental gradient
 * This file calculates the per capita growth rates from the environment-dependent parameters and the population sizes
 * File instruction:
 * 1. Put "gen_beta.h" and "gen_beta.c" in the folder containing this file
 * 2. Put the dsfmt folder and the folder containing this file into the same folder
 ************************************************************************************************
 # Execution (g++ and gcc both work) #
gcc per_capita.c
./a.out
 # Plot with gnuplot #
gnuplot
plot 'summary.txt' using 1:2 title 'rising-tide' with lines lc rgb 'orange',\
'summary.txt' using 1:3 title 'bet-hedging' with lines lc rgb 'skyblue'

 ************************************************************************************************
 * Key parameters
 * s1: Determine the type of distribution function in biological response
 * s3: Determine whether intrinsic growth rate and carrying capacity respond to environment differently
 * s5: Determine whether the average fitness of bet-hedging strategy is lower than that of the rising-tide strategy
 * mean_ris, mean_bet, sd_ris, sd_bet: Shape parameters (i.e. mean and standard deviation) of the two strategies (Gaussian function, s1==1)
 * shape_ris, shape_bet: the shape of performance curves, higher values bring narrower curves (Beta function, s1==2)
 * p[i_ris, i_bet]: The hypothetical population structure
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
    void free_d_vector(double *x);                      // erasing a vector
    double Beta_function (double x_beta1, double alpha_beta, double beta_beta);     // Gives the density of Beta distribution

int main (void){
    // Switches
        int s1= 2;                  // Specify the type of performance curve (1: normal distribution, 2: beta distribution)
        int s3= 1;                  // Whether intrinsic growth rate, b, and carrying capacity, K, respond to environment differently, only applicable to beta distribution performance curves (1: both b and K follow rising-tide/ bet-hedging shape settings, 2: only b follows, 3: only K follows, requires s1==2)
        int s5= 1;                  // Switch of smaller average fitness for bet-hedging strategy (see gamma_bet)
    // Index of variables
        int i_ris= 1;               // Index of the rising-tide strategy
        int i_bet= 2;               // Index of the bet-hedging strategy
    // Creating temporal space
		double *p= d_vector(2);
        int i, j;
        double curr_env, abf_ris, abf_bet;
        double b_ris, b_bet, d_ris, d_bet, K_ris, K_bet, a_rb, a_br;
        a_rb= a_br= 0.4;
        double gamma_ris=   1.0;    // Scaling coefficient of average fitness (area of performance curve)
        double gamma_bet=   0.9;    // Scaling coefficient of average fitness (area of performance curve)
    // Environmental parameters
        double env_limit=   100.0;
        double env_span=    0.1;
        double rep_limit= env_limit/ env_span;
    // Performance curve parameters
        double scale_K=     250.0;  // Scaling coefficient of carrying capacity
        double scale_b=     0.5;    // Scaling coefficient of intrinsic growth rate
        double scale_d=     0.01;   // Scaling coefficient of mortality rate
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
    // Output
		FILE *out;
		out= fopen("summary.txt","w");
	// Initialization
		p[i_ris]= 100;                          // Assumed population of rising-tide strategy
		p[i_bet]= 100;                          // Assumed population of bet-hedging strategy
		fprintf(out,"Environment\tabf_r\tabf_b\n");
    // Main loop
        for (i=1; i< rep_limit; i++){
            curr_env= i*0.1;
            x_beta1= curr_env/range_bio;
            x_beta2= 1-x_beta1;
            if(s1==2){
            // Rising-tide
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
            // Bet-hedging
                // Parameters of the Beta function
                alpha_beta= shape_bet;
                beta_beta= shape_bet*skew_bet;
                env_beta= Beta_function(x_beta1, alpha_beta, beta_beta);
                // Calculation
                K_bet= scale_K* env_beta;
                b_bet= scale_b* env_beta; 
                d_bet= scale_d* (d_max-env_beta);
                if(K_ris!=K_ris) K_ris=0;  // hadling the values out of the range of Beta function
                if(b_ris!=b_ris) b_ris=0;  // hadling the values out of the range of Beta function
                if(d_ris!=d_ris) d_ris=0;  // hadling the values out of the range of Beta function
            // Different carrying capacity curve
                if(s5==1){
                    K_bet= scale_K* env_beta* gamma_bet;
                    b_bet= scale_b* env_beta* gamma_bet; 
                    d_bet= scale_d* (d_max-env_beta* gamma_bet);
                }
                // The two strategy only differ in reproduction response
                if(s3==2){
                    // Parameters of the Bata function
                    alpha_beta= shape_bet;
                    beta_beta= shape_bet*skew_bet;
                    env_beta= Beta_function(x_beta1, alpha_beta, beta_beta);
                    // Calculation
                    K_ris= scale_K* env_beta;
                    K_bet= scale_K* env_beta;
                }
                // The two strategy only differ in survival response
                if(s3==3){
                    // Parameters of the Bata function
                    alpha_beta= shape_bet;
                    beta_beta= shape_bet*skew_bet;
                    env_beta= Beta_function(x_beta1, alpha_beta, beta_beta);
                    // Calculation
                    b_ris= scale_b* env_beta;
                    b_bet= scale_b* env_beta;
                }
            }
        // Normal dist
            if(s1==1){
                // Rising-tide
                env_norm= exp(-1*(curr_env-mean_ris)*(curr_env-mean_ris)/2/sd_ris/sd_ris)/sqrt(2*M_PI*sd_ris*sd_ris);
                K_ris= scale_K*env_norm;
                b_ris= scale_b*env_norm;
                d_ris= scale_d*(d_max-env_norm);
                // if(curr_env>88.5) printf("%lf\t", log(K_ris));
                // Bet-hedging
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
        // Print
                if(s3==1) abf_ris= (b_ris- scale_b/scale_K*p[i_ris]- scale_b/scale_K*a_br*p[i_bet]- d_ris);
                else{
                    if(K_ris>0) abf_ris= (b_ris- b_ris/K_ris*p[i_ris]- b_ris/K_ris*a_br*p[i_bet]- d_ris);
                    else abf_ris= b_ris-d_ris;
                }
                if(s3==1) abf_bet= (b_bet- scale_b/scale_K*p[i_bet]- scale_b/scale_K*a_rb*p[i_ris]- d_bet);
                else{
                    if(K_bet>0) abf_bet= (b_bet- b_bet/K_bet*p[i_bet]- b_bet/K_bet*a_rb*p[i_ris]- d_bet);
                    else abf_bet= b_bet-d_bet;
                }
                fprintf(out, "%lf\t%lf\t%lf\n", curr_env, abf_ris, abf_bet);
        }
	free_d_vector(p);
    fclose(out);
    return 0;
}


// Functions
void message_error(char error_text[]) //standard error handler
{
	printf("There are some errors...\n");
	printf("%s\n",error_text);
	printf("...now existing to system...\n");
	exit(1);
}
double *d_vector(long size) 
{
	double *x;
	x= (double *) malloc((size_t)((size+1)*sizeof(double)));
	if(x==NULL) message_error((char *)"Allocation failure in d_vector()");
	return x;
}
void free_d_vector(double *x) {	free(x);}
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
