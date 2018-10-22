/*
 * Continual generation model with good/ bad year
 * Focal point is founder effect and the selection of better survival/ better birth species
 * First day of building: July 19, 2017
 * Modified from competitive Lotka-Volterra differential equations
 * This file is for alternative mechanism scenario
 * Allee effect is included in here as type (II) functional response
 * Last modified: Sep 12, 2017
 * Replication settings designed: Sep 20, 2017
***************************************************************************************
Execution:
gcc rep_idvl_bet_1011.c -o a.out
./a.out
***************************************************************************************
 * Key parameters
 * The global parameters of differential equations (e.g. r_goodyear_spcl)
 * r1: probability of getting good year
 * s3: switch of deterministic/ stochastic model
 * T_g, T_b: duration of a year in the simulation, correlated to growth rate
 * p[specialsit], p[generalist]: initial popualtion size (or density)
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
    void differential_g(double time,double in[],double out[]);
    double dx_dt_g(double time, double vr[]);           //vr is a vector of variable
    double dy_dt_g(double time, double vr[]);
    void differential_b(double time,double in[],double out[]);
    double dx_dt_b(double time, double vr[]);
    double dy_dt_b(double time, double vr[]);
    double mean (double x[], int length);
    double var (double x[], int length);

// Global variables
    // Index of variables
        int specialist= 1; 		            // Index of x (N1)<- specialist
        int generalist= 2; 	                // Index of y (N2)<- generalist
    // Differential equation paramters
        double r_goodyear_spcl=	0.1;		// Intrinsic growth rate of N1 (good year) // no age classes
        double r_goodyear_genl=	0.1; 		// Intrinsic growth rate of N2 (good year)
        double r_badyear_spcl=	0.05;		// Intrinsic growth rate of N1 (bad year)
        double r_badyear_genl=	0.05;		// Intrinsic growth rate of N2 (bad year)
        double K_goodyear_spcl= 850.0;		// Carrying capacity of N1 (good year)
        double K_goodyear_genl= 590.14;		// Carrying capacity of N2 (good year)
        double K_badyear_spcl=	125.0;		// Carrying capacity of N1 (bad year)
        double K_badyear_genl=	150.0;		// Carrying capacity of N2 (bad year)
        double alpha_goodyear_sg=	1.0;	// Intensity of inter-strategic interactions in good years (specialists to generalists)
        double alpha_goodyear_gs=	1.0;	// Intensity of inter-strategic interactions in good years (generalists to specialists)
        double alpha_badyear_sg=	1.0;	// Intensity of inter-strategic interactions in bad years (specialists to generalists)
        double alpha_badyear_gs=	1.0;	// Intensity of inter-strategic interactions in bad years (generalists to specialists)
    // Type II reponse to population size for generalists (cancelled in current parameter values)
        double shift_c=			0.0;		// y-axis shift
        double shape_c=			0.0;		// shape control
        double scale_c= 		1.0;	    // scaling for type two functional response

// Main function
int main (void)
{
    // Switches and crucial variables
        // Switches
		int s1= 0;				// Switch of good/ bad year
		int s3= 1;				// Switch of stochasticity
		int s4= 0;				// Switch of terminating tiny populaitons
	// Basic variables
		int i,j,k;				// For loop counters
		double r1= 0.7; 		// Rate of getting good season
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
		int id_spcl=    0.0;    // ID of the specialist individuals
		int id_genl=    0.0;    // ID of the generalist individuals
	// Temporal settings
		double t=       0.0;	// Time logs
		int T=          15000; 	// Duration of simulation (unit is different with the deterministic model)
		double deltat=  0.005; 	// Length of time step
		int T_g=        200; 	// Good season length ***Note the difference of double and int
		int T_b=        200; 	// Bad season length ***Note the difference of double and int
		double cc=      0.0;	// Counter of remaining season length
	// Print temporal spaces
		int o1, o2, o3, o4;
	// Creating temporal space
		double *pop_size= d_vector(2);
		double *dfdt= d_vector(2);
		//double **spcl_pop= d_matrix(pop_limit, 6);
		// 1= ID, 2= live/dead, 3= birth, 4= dead, 5= parent, 6=?
		//double **genl_pop= d_matrix(pop_limit, 6);
			int ID= 1;
			int state= 2;
			int birth= 3;
			int death= 4;
			int parent= 5;
			int span= 6;
		//double **offspring= d_matrix(total_pop,2);
			int oct= 0;			// offspring counter
	// Output
		FILE *out, *spcl_log, *genl_log;
		out= fopen("idv_bet_summary.txt","w");

	// Initialization
		pop_size[specialist]= 500;					// initial population of specialist strategy
		pop_size[generalist]= 500;					// initial population of generalist strategy
		total_pop= pop_size[specialist]+pop_size[generalist];
		fprintf(out,"rep\toutcome\ttime\n");

	for(k=1;k<=rep;k++){
		t=cc= 0.0;
		pop_size[specialist]= 500;					// initial population of specialist strategy
		pop_size[generalist]= 500;					// initial population of generalist strategy
		total_pop= pop_size[specialist]+pop_size[generalist];

        // Main loop
        while (t< T){
            if (cc<=0){				// good/ bad year determination
                if(s3==1){
                    s1= 0;
                    pp= dsfmt_genrand_open_open(&dsfmt);
                    if(pp< r1){
                        s1= 1;		// good year starts
                        cc= T_g;}
                    else{
                        s1= 2;		// bad year starts
                        cc= T_b;}
                    }
                if(s3==0){
                    if(s1== 0||s1== 2){
                        s1= 1;		// Good year starts
                        cc= T_g*2*r1;}
                    else{
                        s1= 2;		// Bad year starts
                        cc= T_g*(2-2*r1);}
                    }
                }
            if (s1== 1){			// sum of probability in good year
                if(pop_size[generalist]>0){
                        total_prob= r_goodyear_spcl+ r_goodyear_genl*scale_c*(pop_size[generalist]/(shape_c+pop_size[generalist])+shift_c)+ r_goodyear_spcl/ K_goodyear_spcl* pop_size[specialist]+ alpha_goodyear_gs* r_goodyear_spcl/ K_goodyear_spcl* pop_size[generalist]+ r_goodyear_genl/ K_goodyear_genl* pop_size[generalist]+ alpha_goodyear_sg* r_goodyear_genl/ K_goodyear_genl* pop_size[specialist];
                    }
                else{               // dealing with zero division problem
                        total_prob= r_goodyear_spcl+ r_goodyear_genl*scale_c*shift_c+ r_goodyear_spcl/ K_goodyear_spcl* pop_size[specialist]+ alpha_goodyear_sg* r_goodyear_genl/ K_goodyear_genl* pop_size[specialist];
                    }
                    //*(pop_size[generalist]/allee_coe-1)
                p1= r_goodyear_spcl/ total_prob;
                if(pop_size[generalist]>0){
                        p2= p1+ r_goodyear_genl*scale_c*(pop_size[generalist]/(shape_c+pop_size[generalist])+shift_c)/ total_prob;
                    }
                else p2= p1;
                p3= p2+ (r_goodyear_spcl/ K_goodyear_spcl* pop_size[specialist]+ alpha_goodyear_gs* r_goodyear_spcl/ K_goodyear_spcl* pop_size[generalist])/ total_prob;
            }
            if (s1== 2){			// sum of probability in bad year
                if(pop_size[generalist]>0){
                    total_prob= r_badyear_spcl+ r_badyear_genl*scale_c*(pop_size[generalist]/(shape_c+pop_size[generalist])+shift_c)+ r_badyear_spcl/ K_badyear_spcl* pop_size[specialist]+ alpha_badyear_gs* r_badyear_spcl/ K_badyear_spcl* pop_size[generalist]+ r_badyear_genl/ K_badyear_genl* pop_size[generalist]+ alpha_badyear_sg* r_badyear_genl/ K_badyear_genl* pop_size[specialist];
                    }
                else{
                    total_prob= r_badyear_spcl+ r_badyear_genl*scale_c*shift_c+ r_badyear_spcl/ K_badyear_spcl* pop_size[specialist]+ alpha_badyear_sg* r_badyear_genl/ K_badyear_genl* pop_size[specialist];
                    }
                p1= r_badyear_spcl/ total_prob;
                if(pop_size[generalist]>0){
                    p2= p1+ r_badyear_genl*scale_c*(pop_size[generalist]/(shape_c+pop_size[generalist])+shift_c)/ total_prob;
                    }
                else p2= p1;
                p3= p2+ (r_badyear_spcl/ K_badyear_spcl* pop_size[specialist]+ alpha_badyear_gs* r_badyear_spcl/ K_badyear_spcl* pop_size[generalist])/ total_prob;
            }
                
            pp= dsfmt_genrand_open_open(&dsfmt);			// random number to decide the length of this step
            S= -log(pp)/total_prob/10;						// length of time advanced

        // individual level events
            // log matrix intialization
            pp= dsfmt_genrand_open_open(&dsfmt);			// random number to decide the event
            if(pp<p1&& pop_size[specialist]>0) {				// Birth of specialist individual
                pop_size[specialist]= pop_size[specialist]+1;
                }
            else{
                if(pp<p2&& pop_size[generalist]>0) {		// Birth of generalist individual
                    pop_size[generalist]= pop_size[generalist]+1;
                }
                else{
                    if(pp<p3&& pop_size[specialist]>0) { 		// Death of specialist individual
                        pop_size[specialist]= pop_size[specialist]-1;
                        }
                    else {									// Death of generalist individual
                        if(pop_size[generalist]>0){
                            pop_size[generalist]= pop_size[generalist]-1;
                        }
                        }
                }}

            t+= S;
            cc-= S;
            if (t>= T){
                fprintf(out,"%d\t%d\t%lf\n",k,0,t);
            }
            if (pop_size[generalist]==0){
                fprintf(out,"%d\t%d\t%lf\n",k,-1,t);
                t= T+ 0.1;
            }
            if (pop_size[specialist]==0){
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
double dx_dt_g(double time, double vr[])	// good year specialist
{return (r_goodyear_spcl- r_goodyear_spcl/K_goodyear_spcl*vr[specialist]- r_goodyear_spcl/K_goodyear_spcl*alpha_goodyear_gs*vr[generalist])*vr[specialist];}
double dy_dt_g(double time, double vr[])	// good year generalist
{return (r_goodyear_genl- r_goodyear_genl/K_goodyear_genl*vr[generalist]- r_goodyear_genl/K_goodyear_genl*alpha_goodyear_sg*vr[specialist])*vr[generalist];}
void differential_g(double time, double in[], double out[])
{
	out[specialist]= dx_dt_g(time,in);
	out[generalist]= dy_dt_g(time,in);
}
double dx_dt_b(double time, double vr[])	// bad year specialist
{return (r_badyear_spcl- r_badyear_spcl/K_badyear_spcl*vr[specialist]- r_badyear_spcl/K_badyear_spcl*alpha_badyear_gs*vr[generalist])*vr[specialist];}
double dy_dt_b(double time, double vr[])	// bad year generalist
{return (r_badyear_genl- r_badyear_genl/K_badyear_genl*vr[generalist]- r_badyear_genl/K_badyear_genl*alpha_badyear_sg*vr[specialist])*vr[generalist];}
void differential_b(double time, double in[], double out[])
{
	out[specialist]= dx_dt_b(time,in);
	out[generalist]= dy_dt_b(time,in);
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
