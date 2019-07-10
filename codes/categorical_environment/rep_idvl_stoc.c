/*
 * Continuous growth and individual-based model with discrete environemnts (high/ low quality environments)
 * This file is designed for repeating simulations and calculating proportions of fixation
 * First day of building: July 19, 2017
 * Modified from competitive Lotka-Volterra differential equations
 * Last modified: Jul 10, 2019
 **************************************************************************************
Execution:
gcc rep_idvl_idvl_stoc.c -o a.out
./a.out
 ***************************************************************************************
 * Key parameters
 * The global parameters of differential equations (e.g. r_highqual_ris)
 * r1: probability of getting high quality environment
 * s3: switch of deterministic/ stochastic model
 * T_h, T_l: duration of a year in the simulation, correlated to growth rate
 * p[rising-tide], p[bet_hedging]: initial popualtion size (or density)
 */
#include <stdio.h>
#include <stdlib.h>
#include "../dSFMT-src-2.2.3/dSFMT.c"
#include <math.h>
#include <time.h>

// Functions
    void message_error(char error_text[]);              // printing errors
    double *d_vector(long size);                        // creating a vector
    double **d_matrix(long size_row, long size_column); // creating a matrix
    void free_d_vector(double *x);                      // erasing a vector
    void free_d_matrix(double **x);                     // erasing a matrix
    void rk4(double p[], double k1[],int n, double t, double h, double pout[],void(*diff)(double,double[],double[]));
    void differential_h(double time,double in[],double out[]);
    double dx_dt_h(double time, double vr[]);           //vr is a vector of variable
    double dy_dt_h(double time, double vr[]);
    void differential_l(double time,double in[],double out[]);
    double dx_dt_l(double time, double vr[]);
    double dy_dt_l(double time, double vr[]);
    double mean (double x[], int length);
    double var (double x[], int length);

// Global variables
    // Index of variables
        int rising_tide= 1; 		            // Index of x (N1)<- rising_tide
        int bet_hedging= 2; 	                // Index of y (N2)<- bet_hedging
    // Differential equation paramters
        double r_highqual_ris=	0.1;		// Intrinsic growth rate of N1 (high quality environment) // no age classes
        double r_highqual_bet=	0.1; 		// Intrinsic growth rate of N2 (high quality environment)
        double r_lowqual_ris=	0.05;		// Intrinsic growth rate of N1 (low quality environment)
        double r_lowqual_bet=	0.05;		// Intrinsic growth rate of N2 (low quality environment)
        double K_highqual_ris= 850.0;		// Carrying capacity of N1 (high quality environment)
        double K_highqual_bet= 590.0;		// Carrying capacity of N2 (high quality environment)
        double K_lowqual_ris=	125.0;		// Carrying capacity of N1 (low quality environment)
        double K_lowqual_bet=	150.0;		// Carrying capacity of N2 (low quality environment)
        double alpha_highqual_rb=	1.0;	// Intensity of inter-strategic interactions in high quality environments (rising_tide to bet_hedging)
        double alpha_highqual_br=	1.0;	// Intensity of inter-strategic interactions in high quality environments (bet_hedging to rising_tide)
        double alpha_lowqual_rb=	1.0;	// Intensity of inter-strategic interactions in low quality environments (rising_tide to bet_hedging)
        double alpha_lowqual_br=	1.0;	// Intensity of inter-strategic interactions in low quality environments (bet_hedging to rising_tide)
    // Type II reponse to population size for bet_hedging (cancelled in current parameter values)
        double shift_c=			0.0;		// y-axis shift
        double shape_c=			0.0;		// shape control
        double scale_c= 		1.0;	    // scaling for type two functional response

// Main function
int main (void)
{
    // Switches and crucial variables
        // Switches
		int s1= 0;				// Switch of high/ low quality environment
		int s3= 1;				// Switch of stochasticity
		int s4= 0;				// Switch of terminating tiny populaitons
	// Basic variables
		int i,j,k;				// For loop counters
		double r1= 0.7; 		// Rate of getting high environment
		double pp= 0.0;			// Temp for random number
		double pp2= 0.0;		// Another temp for random number
		double ext_thr= 1.0;	// Threshold of population for terminating (refer to s4)
        // temporal space for determining the probability of each type of events
		double total_prob, p1, p2, p3, S;
		int rid;				// integer for random selection
		int idc; 				// interger for id counter (surviving individual)
		int total_pop;			// temp of total population
		double ops= 0.0;		// temp of opportunity of selection
		int rep= 1000;
   	// Random number genertor
		int seed;
		dsfmt_t dsfmt;
		seed= time(NULL);
		if(seed==0)seed= 1;
		dsfmt_init_gen_rand(&dsfmt,seed);
	// Matrix and log settings
		int pop_limit=  1E5;    // Upper limit of the number of created individuals (each strategy)
		int pop_count=  0;      // Counting the total population size 
		int id_ris=    0.0;    // ID of the rising_tide individuals
		int id_bet=    0.0;    // ID of the bet_hedging individuals
	// Temporal settings
		double t=       0.0;	// Time logs
		int T=          15000; 	// Duration of simulation (unit is different with the deterministic model)
		double deltat=  0.005; 	// Length of time step
		int T_h=        200; 	// high environment length ***Note the difference of double and int
		int T_l=        200; 	// low environment length ***Note the difference of double and int
		double cc=      0.0;	// Counter of remaining environment length
	// Print temporal spaces
		int o1, o2, o3, o4;
	// Creating temporal space
		double *pop_size= d_vector(2);
		double *dfdt= d_vector(2);
		//double **ris_pop= d_matrix(pop_limit, 6);
		// 1= ID, 2= live/dead, 3= birth, 4= dead, 5= parent, 6=?
		//double **bet_pop= d_matrix(pop_limit, 6);
			int ID= 1;
			int state= 2;
			int birth= 3;
			int death= 4;
			int parent= 5;
			int span= 6;
		//double **offspring= d_matrix(total_pop,2);
			int oct= 0;			// offspring counter
	// Output
		FILE *out, *ris_log, *bet_log;
		out= fopen("idv_bet_summary.txt","w");

	// Initialization
		pop_size[rising_tide]= 500;					// initial population of rising_tide strategy
		pop_size[bet_hedging]= 500;					// initial population of bet_hedging strategy
		total_pop= pop_size[rising_tide]+pop_size[bet_hedging];
		fprintf(out,"rep\toutcome\ttime\n");

	for(k=1;k<=rep;k++){
		t=cc= 0.0;
		pop_size[rising_tide]= 500;					// initial population of rising_tide strategy
		pop_size[bet_hedging]= 500;					// initial population of bet_hedging strategy
		total_pop= pop_size[rising_tide]+pop_size[bet_hedging];

        // Main loop
        while (t< T){
            if (cc<=0){				// high/ low quality environment determination
                if(s3==1){
                    s1= 0;
                    pp= dsfmt_genrand_open_open(&dsfmt);
                    if(pp< r1){
                        s1= 1;		// high quality environment starts
                        cc= T_h;}
                    else{
                        s1= 2;		// low quality environment starts
                        cc= T_l;}
                    }
                if(s3==0){
                    if(s1== 0||s1== 2){
                        s1= 1;		// high quality environment starts
                        cc= T_h*2*r1;}
                    else{
                        s1= 2;		// low quality environment starts
                        cc= T_l*(2-2*r1);}
                    }
                }
            if (s1== 1){			// sum of probability in high quality environment
                if(pop_size[bet_hedging]>0){
                        total_prob= r_highqual_ris+ r_highqual_bet*scale_c*(pop_size[bet_hedging]/(shape_c+pop_size[bet_hedging])+shift_c)+ r_highqual_ris/ K_highqual_ris* pop_size[rising_tide]+ alpha_highqual_br* r_highqual_ris/ K_highqual_ris* pop_size[bet_hedging]+ r_highqual_bet/ K_highqual_bet* pop_size[bet_hedging]+ alpha_highqual_rb* r_highqual_bet/ K_highqual_bet* pop_size[rising_tide];
                    }
                else{               // dealing with zero division problem
                        total_prob= r_highqual_ris+ r_highqual_bet*scale_c*shift_c+ r_highqual_ris/ K_highqual_ris* pop_size[rising_tide]+ alpha_highqual_rb* r_highqual_bet/ K_highqual_bet* pop_size[rising_tide];
                    }
                    //*(pop_size[bet_hedging]/allee_coe-1)
                p1= r_highqual_ris/ total_prob;
                if(pop_size[bet_hedging]>0){
                        p2= p1+ r_highqual_bet*scale_c*(pop_size[bet_hedging]/(shape_c+pop_size[bet_hedging])+shift_c)/ total_prob;
                    }
                else p2= p1;
                p3= p2+ (r_highqual_ris/ K_highqual_ris* pop_size[rising_tide]+ alpha_highqual_br* r_highqual_ris/ K_highqual_ris* pop_size[bet_hedging])/ total_prob;
            }
            if (s1== 2){			// sum of probability in low quality environment
                if(pop_size[bet_hedging]>0){
                    total_prob= r_lowqual_ris+ r_lowqual_bet*scale_c*(pop_size[bet_hedging]/(shape_c+pop_size[bet_hedging])+shift_c)+ r_lowqual_ris/ K_lowqual_ris* pop_size[rising_tide]+ alpha_lowqual_br* r_lowqual_ris/ K_lowqual_ris* pop_size[bet_hedging]+ r_lowqual_bet/ K_lowqual_bet* pop_size[bet_hedging]+ alpha_lowqual_rb* r_lowqual_bet/ K_lowqual_bet* pop_size[rising_tide];
                    }
                else{
                    total_prob= r_lowqual_ris+ r_lowqual_bet*scale_c*shift_c+ r_lowqual_ris/ K_lowqual_ris* pop_size[rising_tide]+ alpha_lowqual_rb* r_lowqual_bet/ K_lowqual_bet* pop_size[rising_tide];
                    }
                p1= r_lowqual_ris/ total_prob;
                if(pop_size[bet_hedging]>0){
                    p2= p1+ r_lowqual_bet*scale_c*(pop_size[bet_hedging]/(shape_c+pop_size[bet_hedging])+shift_c)/ total_prob;
                    }
                else p2= p1;
                p3= p2+ (r_lowqual_ris/ K_lowqual_ris* pop_size[rising_tide]+ alpha_lowqual_br* r_lowqual_ris/ K_lowqual_ris* pop_size[bet_hedging])/ total_prob;
            }
                
            pp= dsfmt_genrand_open_open(&dsfmt);			// random number to decide the length of this step
            S= -log(pp)/total_prob/10;						// length of time advanced

        // individual level events
            // log matrix intialization
            pp= dsfmt_genrand_open_open(&dsfmt);			// random number to decide the event
            if(pp<p1&& pop_size[rising_tide]>0) {				// Birth of rising_tide individual
                pop_size[rising_tide]= pop_size[rising_tide]+1;
                }
            else{
                if(pp<p2&& pop_size[bet_hedging]>0) {		// Birth of bet_hedging individual
                    pop_size[bet_hedging]= pop_size[bet_hedging]+1;
                }
                else{
                    if(pp<p3&& pop_size[rising_tide]>0) { 		// Death of rising_tide individual
                        pop_size[rising_tide]= pop_size[rising_tide]-1;
                        }
                    else {									// Death of bet_hedging individual
                        if(pop_size[bet_hedging]>0){
                            pop_size[bet_hedging]= pop_size[bet_hedging]-1;
                        }
                        }
                }}

            t+= S;
            cc-= S;
            if (t>= T){
                fprintf(out,"%d\t%d\t%lf\n",k,0,t);
            }
            if (pop_size[bet_hedging]==0){
                fprintf(out,"%d\t%d\t%lf\n",k,-1,t);
                t= T+ 0.1;
            }
            if (pop_size[rising_tide]==0){
                fprintf(out,"%d\t%d\t%lf\n",k,1,t);
                t= T+ 0.1;
            }
	}}

	free_d_vector(pop_size);
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
double dx_dt_h(double time, double vr[])	// high quality environment rising_tide
{return (r_highqual_ris- r_highqual_ris/K_highqual_ris*vr[rising_tide]- r_highqual_ris/K_highqual_ris*alpha_highqual_br*vr[bet_hedging])*vr[rising_tide];}
double dy_dt_h(double time, double vr[])	// high quality environment bet_hedging
{return (r_highqual_bet- r_highqual_bet/K_highqual_bet*vr[bet_hedging]- r_highqual_bet/K_highqual_bet*alpha_highqual_rb*vr[rising_tide])*vr[bet_hedging];}
void differential_h(double time, double in[], double out[])
{
	out[rising_tide]= dx_dt_h(time,in);
	out[bet_hedging]= dy_dt_h(time,in);
}
double dx_dt_l(double time, double vr[])	// low quality environment rising_tide
{return (r_lowqual_ris- r_lowqual_ris/K_lowqual_ris*vr[rising_tide]- r_lowqual_ris/K_lowqual_ris*alpha_lowqual_br*vr[bet_hedging])*vr[rising_tide];}
double dy_dt_l(double time, double vr[])	// low quality environment bet_hedging
{return (r_lowqual_bet- r_lowqual_bet/K_lowqual_bet*vr[bet_hedging]- r_lowqual_bet/K_lowqual_bet*alpha_lowqual_rb*vr[rising_tide])*vr[bet_hedging];}
void differential_l(double time, double in[], double out[])
{
	out[rising_tide]= dx_dt_l(time,in);
	out[bet_hedging]= dy_dt_l(time,in);
}
void rk4(double pop_size[], double k1[],int n, double t, double h, double pout[],void(*diff)(double,double[],double[]))
{
	int i;
	double tt,*k2,*k3,*k4,*pp;

	k2= d_vector(n);
	k3= d_vector(n);
	k4= d_vector(n);
	pp= d_vector(n);

	for (i=1;i<=n;i++) pp[i]= pop_size[i]+ k1[i]*h/2;
	(*diff)(t+h/2,pp,k2);
	for (i=1;i<=n;i++) pp[i]= pop_size[i]+ k2[i]*h/2;
	(*diff)(t+h/2,pp,k3);
	for (i=1;i<=n;i++) pp[i]= pop_size[i]+ k3[i]*h;
	(*diff)(t+h,pp,k4);
	
	for(i=1;i<=n;i++) pout[i]= pop_size[i]+ (k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*h/6;

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

	x= (double **) malloc((size_t)(size_row_P*sizeof(double *))); //first dimension
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
double mean (double x[], int length)
{
	int i;
	double tmp=0.0;
	for (i=1; i<= length; i++){
		tmp+= x[i];
	}
	return tmp/length;
}
double var (double x[], int length)
{
	int i;
	double mean, t1, t2;
	mean= t1= t2= 0;
	for(i=1; i<= length; i++){
		t1+= x[i];
	}
	mean= t1/length;
	for(i=1; i<= length; i++){
		t2+= (x[i]-mean)*(x[i]-mean);
	}
	return t2/(length-1); // variance of samples, not matrix
}
