/*
 * Continuous population dynamics model with categorical environment (high quality/ low quality)
 * This file generates the repeated simulations under certain probability of encountering high quality environments (r1) 
 */
#include <stdio.h>
#include <stdlib.h>
#include "../dSFMT-src-2.2.3/dSFMT.c"
#include <time.h>
#include <math.h>

// Functions
void message_error(char error_text[]);
double *d_vector(long size); //vector creator
double **d_matrix(long size_row, long size_column); //matrix creater
void free_d_vector(double *x);
void free_d_matrix(double **x);
void rk4(double p[], double k1[],int n, double t, double h, double pout[],void(*diff)(double,double[],double[]));
void differential_h(double time,double in[],double out[]);
double dx_dt_h(double time, double vr[]); //vr is variable
double dy_dt_h(double time, double vr[]);
void differential_l(double time,double in[],double out[]);
double dx_dt_l(double time, double vr[]); //vr is variable
double dy_dt_l(double time, double vr[]);

// Global variables
// Index of variables
	int i_x= 1; //index of x (N1)<- rising-tide strategy
	int i_y= 2; //index of y (N2)<- bet-hedging strategy
// differential equation paramters
	double r_h1= 0.1;		// Growth rate of N1 (high quality environment) // no age classes
	double r_h2= 0.1; 		// Growth rate of N2 (high quality environment)
	double r_l1= 0.05;		// Growth rate of N1 (low quality environment)
	double r_l2= 0.05;		// Growth rate of N2 (low quality environment)
	double K_h1= 850.0;	    // Carrying capacity of N1 (high quality environment)
	double K_h2= 590.0;	    // Carrying capacity of N2 (high quality environment)
	double K_l1= 125.0;		// Carrying capacity of N1 (low quality environment)
	double K_l2= 150.0;		// Carrying capacity of N2 (low quality environment)
	double a_h12= 1.0;		// level of inter-strategic competition relation
	double a_h21= 1.0;		// level of inter-strategic competition relation
	double a_l12= 1.0;		// level of inter-strategic competition relation
	double a_l21= 1.0;		// level of inter-strategic competition relation
    // Allee effect coefficients, do not change them (no Allee effect considered now)
	double shift_c= 0.0;	// y-axis shifting for type II response
	double shape_c= 0.0;	// shape coefficient for type II functional response
	double scale_c= 1.0;	// scaling coefficient

// Main function
int main (void)
{
    // Switches
		int s3= 1;				// Switch of stochasticity
		int s4= 1;				// Switch of extinction (go to zero if it is close to)
	// Local variables
		int i,j,k;				// For loop counters
		double t= 0.0;			// Time logs
		int T= 2000; 			// Duration of simulation
		double deltat= 0.005; 	// Length of time step
		int T_h= 20; 			// high quality season length ***Note the difference of double and int
		int T_l= 20; 			// low quality season length ***Note the difference of double and int
		double r1= 0.5; 		// Rate of getting high quality season
		double cc= 0.0;			// Counter of remaining season length
		double pp= 0.0;			// Temp for random number
		double et= 0.1;			// Threshold of extinction
		int s1= 0;				// Determine the type of environment
		int rep= 1000;			// Number of repetitions
   	// Random number genertor
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
		out= fopen("summary.txt","w");
	// Initialization
		//fprintf(out,"r1 is %lf, r2 is %lf, db is %lf, duration is %d\n",b_h1,b_h2,d_l2,T);
		fprintf(out,"rep\toutcome\ttime\n");
	// Starting of repetition, and reseting parameter values
		for (k=1; k<= rep; k++){
			p[i_x]= 500;// initial population of x (rising-tide)
			p[i_y]= 500;// initial population of y (bet-hedging)
			cc= 0.0;
			t=0.0;

	// Main loop of time series
		for (i=1; i<=T/deltat; i++){
			if (cc<=0){			// another season
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
						cc= T_h*r1/(1-r1);}
					else{
						s1= 2;		// low quality environment starts
						cc= T_l;}
				}
			}
			if(s1==1){
				differential_h(t,p,dfdt);
				rk4(p, dfdt, 2, t, deltat, p, differential_h);}
			if(s1==2){
				differential_l(t,p,dfdt);
				rk4(p, dfdt, 2, t, deltat, p, differential_l);}
			cc-=deltat;
			t+=deltat;

			// if there is no fixation
			if(i==T/deltat){
				fprintf(out,"%d\t%d\t%lf\n",k,0,t);
				/////////////////////////////////// 0= no fixation (coexistence or polymorphism)
			}
			// extinction of population (fixation)
			if(s4==1){
				if(p[i_x]<et) {
					p[i_x]=0;
					fprintf(out,"%d\t%d\t%lf\n",k,1,t);
					/////////////////////////////////// 1= fixation of bet-hedging strategy
					i= T/deltat;
				}
				if(p[i_y]<et) {
					p[i_y]=0;
					fprintf(out,"%d\t%d\t%lf\n",k,-1,t);
					//////////////////////////////////// -1= fixation of rising-tide strategy
					i= T/deltat;
				}
			}
			// warning message
			if(p[i_x]<0||p[i_y]<0) {
				printf("dynamics crashed");
				exit(1);}
		}}

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

double dx_dt_h(double time, double vr[])	// rising-tide in high quality environment
{return (r_h1- r_h1/K_h1*vr[i_x]- r_h1/K_h1*a_h21*vr[i_y])*vr[i_x];}
double dy_dt_h(double time, double vr[])	// bet-hedging
{
    if(vr[i_y]>0) return (r_h2*scale_c*(shift_c+vr[i_y]/(shape_c+vr[i_y]))- r_h2/K_h2*vr[i_y]- r_h2/K_h2*a_h12*vr[i_x])*vr[i_y];
    else return 0;
}
void differential_h(double time, double in[], double out[])
{
	out[i_x]= dx_dt_h(time,in);
	out[i_y]= dy_dt_h(time,in);
}
double dx_dt_l(double time, double vr[])	// rising-tide in low quality environment
{return (r_l1- r_l1/K_l1*vr[i_x]- r_l1/K_l1*a_l21*vr[i_y])*vr[i_x];}
double dy_dt_l(double time, double vr[])
{
    if(vr[i_y]>0) return (r_l2*scale_c*(shift_c+vr[i_y]/(shape_c+vr[i_y]))- r_l2/K_l2*vr[i_y]- r_l2/K_l2*a_l12*vr[i_x])*vr[i_y];
    else return 0;
}
void differential_l(double time, double in[], double out[])
{
	out[i_x]= dx_dt_l(time,in);
	out[i_y]= dy_dt_l(time,in);
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
