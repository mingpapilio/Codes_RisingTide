/*
 * Continuous growth and deterministic model with discrete environemnts (good/ bad years)
 * First day of building: July 19, 2017
 * Modified from competitive Lotka-Volterra differential equations
 * Last modified: Oct 22, 2018
 ************************************************************************************************
Execution:
gcc det_denv_cg.c -o a.out
./a.out
 ************************************************************************************************
 * Key parameters
 * The global parameters of differential equations (e.g. r_g1)
 * r1: probability of getting good year
 * s3: switch of temporal stochasticity
 * T_g, T_b: duration of a year in the simulation, correlated to growth rate
 * p[i_x], p[i_y]: initial popualtion size (or density)
 */
#include <stdio.h>
#include <stdlib.h>
#include "../dSFMT-src-2.2.3/dSFMT.c"
#include <time.h>
#include <math.h>

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

// Global variables
    // Index of variables
        int i_x= 1;             //index of x (N1)<- specialist (rising-tide strategy)
        int i_y= 2;             //index of y (N2)<- generalist (bet-hedging strategy)
    // Differential equation paramters
        double r_g1= 0.1;		// Intrinsic growth rate of N1 (good year) // no age classes
        double r_g2= 0.1; 		// Intrinsic growth rate of N2 (good year)
        double r_b1= 0.05;		// Intrinsic growth rate of N1 (bad year)
        double r_b2= 0.05;		// Intrinsic growth rate of N2 (bad year)
        double K_g1= 850.0;	    // Carrying capacity of N1 (good year)
        double K_g2= 590.14;	// Carrying capacity of N2 (good year)
        double K_b1= 125.0;		// Carrying capacity of N1 (bad year)
        double K_b2= 150.0;		// Carrying capacity of N2 (bad year)
        double a_g12= 1.0;		// Intensity of inter-strategic interactions in good years (N1 to N2)
        double a_g21= 1.0;		// Intensity of inter-strategic interactions in good years (N2 to N1)
        double a_b12= 1.0;		// Intensity of inter-strategic interactions in bad years (N1 to N2)
        double a_b21= 1.0;		// Intensity of inter-strategic interactions in bad years (N2 to N1)

// Main function
int main (void)
{
	// Switches and crucial variables
        // Seitches
		int s1= 0;				// Switch of good/ bad year (determines the type of year, do NOT change)
		int s3= 0;				// Switch of stochasticity (deterministic or stochastic model)
		int s4= 0;				// Switch of terminating tiny populaitons
    
    // Basic variables
		int i,j;				// For loop counters
		double t= 0.0;			// Time logs
		int T= 2000; 			// Duration of simulation
		double deltat= 0.005; 	// Length of time step
		int T_g= 20; 			// Length of a good year ***Note the difference of double and int
		int T_b= 20; 			// Length of a bad year ***Note the difference of double and int
		double r1= 0.3; 		// Chance of getting a good year
		double cc= 0.0;			// Counter of remaining time of a year
		double pp= 0.0;			// Temp for random number
		double ext_thr= 0.5;	// Threshold of population for terminating (refer to s4)

   	// Initialization of random number genertor
		int seed;
		dsfmt_t dsfmt;
		seed= time(NULL);
		if(seed==0)seed= 1;
		dsfmt_init_gen_rand(&dsfmt,seed);
	// Creating temporal space
		double *p= d_vector(2);
		double *dfdt= d_vector(2);
	// Output
		FILE *out;
		out= fopen("out.txt","w");
	// Initialization
		p[i_x]= 500;            // initial population of N1
		p[i_y]= 500;            // initial population of N2
		//fprintf(out,"r1 is %lf, r2 is %lf, db is %lf, duration is %d\n",b_g1,b_g2,d_b2,T);
		fprintf(out,"N1\t\tN2\t\ttime\n");
		fprintf(out,"%lf\t%lf\t%lf\n",p[i_x],p[i_y],t);

	// Main loop of time series
	for (i=1; i<=T/deltat; i++){
		if (cc<=0){			    // Another year
			if(s3==1){          // Stochastic model
				s1= 0;
                // Draw a random number
				pp= dsfmt_genrand_open_open(&dsfmt);
				if(pp< r1){
					s1= 1;		// Good year starts
					cc= T_g;}
				else{
					s1= 2;		// Bad year starts
					cc= T_b;}
			}
			if(s3==0){          // Deterministic model
				if(s1== 0||s1== 2){
					s1= 1;		// Good year starts
					cc= T_g*2*r1;}
				else{
					s1= 2;		// Bad year starts
					cc= T_g*(2-2*r1);}
			}
		}
		if(s1==1){              // Good year
			differential_g(t,p,dfdt);
			rk4(p, dfdt, 2, t, deltat, p, differential_g);}
		if(s1==2){              // Bad year
			differential_b(t,p,dfdt);
			rk4(p, dfdt, 2, t, deltat, p, differential_b);}
		cc-=deltat;
		t+=deltat;
		// record current population
		if(i%100==0) fprintf(out,"%lf\t%lf\t%lf\n",p[i_x],p[i_y],t);
		// extinction of population
		if(s4==1){
			if(p[i_x]<ext_thr) p[i_x]=0;
			if(p[i_y]<ext_thr) p[i_y]=0;
		}
		// warning message
		if(p[i_x]<0||p[i_y]<0) {
			printf("dynamics crashed");
			exit(2);}
	}

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
{return (r_g1- r_g1/K_g1*vr[i_x]- r_g1/K_g1*a_g21*vr[i_y])*vr[i_x];}
double dy_dt_g(double time, double vr[])	// Cooperative
{return (r_g2- r_g2/K_g2*vr[i_y]- r_g2/K_g2*a_g12*vr[i_x])*vr[i_y];}
void differential_g(double time, double in[], double out[])
{
	out[i_x]= dx_dt_g(time,in);
	out[i_y]= dy_dt_g(time,in);
}
double dx_dt_b(double time, double vr[])
{return (r_b1- r_b1/K_b1*vr[i_x]- r_b1/K_b1*a_b21*vr[i_y])*vr[i_x];}
double dy_dt_b(double time, double vr[])
{return (r_b2- r_b2/K_b2*vr[i_y]- r_b2/K_b2*a_b12*vr[i_x])*vr[i_y];}
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
