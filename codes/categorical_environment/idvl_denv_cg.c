/*
 * Continuous growth and individual-based model with discrete environemnts (high/ low quality environments)
 * First day of building: July 19, 2017
 * Modified from competitive Lotka-Volterra differential equations
 * Last modified: Jul 10, 2019
 ***************************************************************************************
 # Execution #
gcc idvl_denv_cg.c
./a.out
 # Plot with gnuplot # 
gnuplot
plot 'pop_dmc.txt' using 1:3 title 'rising-tide' with lines lc rgb 'orange',\
'pop_dmc.txt' using 1:4 title 'bet-hedging' with lines lc rgb 'skyblue'
 ***************************************************************************************
 * Key parameters
 * The global parameters of differential equations (e.g. r_high_ris)
 * r1: Probability of getting high quality environment
 * s3: Switch of temporal stochasticity
 * T_h, T_l: Duration of a period of single environment
 * pop_size[i_ris], pop_size[i_bet]: Initial popualtion size
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
        int i_ris= 1; 		            // Index of rising-tide (rising-tide strategy)
        int i_bet= 2; 	                // Index of bet-hedging (bet-hedging strategy)
    // Differential equation paramters
        double r_high_ris=	0.1;		// Intrinsic growth rate of N1 (high quality environment) // no age classes
        double r_high_bet=	0.1; 		// Intrinsic growth rate of N2 (high quality environment)
        double r_low_ris=	0.05;		// Intrinsic growth rate of N1 (low quality environment)
        double r_low_bet=	0.05;		// Intrinsic growth rate of N2 (low quality environment)
        double K_high_ris= 850.0;		// Carrying capacity of N1 (high quality environment)
        double K_high_bet= 590.14;		// Carrying capacity of N2 (high quality environment)
        double K_low_ris=	125.0;		// Carrying capacity of N1 (low quality environment)
        double K_low_bet=	150.0;		// Carrying capacity of N2 (low quality environment)
        double alpha_high_rb=	1.0;	// Intensity of inter-strategic interactions in high quality environments (rising-tide to bet-hedging)
        double alpha_high_br=	1.0;	// Intensity of inter-strategic interactions in high quality environments (bet-hedging to rising-tide)
        double alpha_low_rb=	1.0;	// Intensity of inter-strategic interactions in low quality environments (rising-tide to bet-hedging)
        double alpha_low_br=	1.0;	// Intensity of inter-strategic interactions in low quality environments (bet-hedging to rising-tide)

// Main function
int main (void)
{
    // Switches and crucial variables
        // Switches
		int s3= 1;				// Switch of stochasticity
		int s4= 0;				// Switch of terminating tiny populaitons
	// Basic variables
		int i,j;				// For loop counters
		double r1=      0.3; 	// Chance of getting a high quality environment
		int env_type=	0;		// Specify the environment type
		double pp=      0.0;	// Temp for random number
		double pp2=     0.0;	// Another temp for random number
		double ext_thr= 1.0;	// Threshold of population for terminating (refer to s4)
        // temporal space for determining the probability of each type of events
        double total_prob, p1, p2, p3, p6, S;
        int p4, p5;             // temporal space for opportunity for selection
		int rid;				// integer for random selection
		int idc; 				// interger for id counter (surviving individual)
		int total_pop;			// temp of total population
   	// Random number genertor
		int seed;
		dsfmt_t dsfmt;
		seed= time(NULL);
		if(seed==0)seed= 1;
		dsfmt_init_gen_rand(&dsfmt,seed);
	// Matrix and log settings
		int pop_limit=  1E5;    // Upper limit of the number of created individuals (each strategy)
		int pop_count=  0;      // Counting the total population size 
		int id_ris=    0.0;    // ID of the rising-tide individuals
		int id_bet=    0.0;    // ID of the bet-hedging individuals
	// Temporal settings
		double t=       0.0;	// Time logs
		int T=          15000; 	// Duration of simulation (unit is different with the deterministic model)
		double deltat=  0.005; 	// Length of time step
		int T_h=        200; 	// Length of a period of high quality environment ***If s3==0, the length change proportionally with r1
		int T_l=        200; 	// Length of a period of low quality environment 
        double cc=      0.0;	// Counter of remaining time of a period
        int scale_time= 10;     // Control the pace of each dynamic step between events
    // Opportunity for selection settings
		double ops=     0.0;	// Temp of opportunity of selection
		double abf_ris=   0.0;	// Absolute fitness (rising-tide)
        double abf_bet=   0.0;	// Absolute fitness (bet-hedging)
        double ops_span= 200.0;     // Duration for updating opportunity for selection
        double ops_count= ops_span; // Time counter
	// Print temporal spaces
		int o1, o2, o3, o4;
	// Creating temporal space
		double *pop_size= d_vector(2);
		double *dfdt= d_vector(2);
		double **ris_pop= d_matrix(pop_limit, 6);
		// 1= ID, 2= live/dead, 3= birth time, 4= death time, 5= parent (haploid), 6=life span
		double **bet_pop= d_matrix(pop_limit, 6);
			int ID= 1;
			int state= 2;
			int birth= 3;
			int death= 4;
			int parent= 5;
			int span= 6;
		//double **offspring= d_matrix(total_pop,2);
            int oct= 0;			// offspring counter (for calculating fitness)
	// Output
		FILE *out, *ris_log, *bet_log;
		out= fopen("pop_dmc.txt","w");
		ris_log= fopen("idv_ris_log.txt","w");
		bet_log= fopen("idv_bet_log.txt","w");
	// Initialization
		pop_size[i_ris]= 500;					// initial population of rising-tide strategy
        pop_size[i_bet]= 500;					// initial population of bet-hedging strategy
            int tp_size= pop_size[i_ris]+pop_size[i_bet];
            total_pop= pop_size[i_ris]+pop_size[i_bet];
            double *offspring= d_vector(total_pop);
            double *refitness= d_vector(total_pop);
            p4= pop_size[i_ris];
		fprintf(out,"time\t\tenv\tris\tbet\tTotal\topp for sel\tabf_ris\tabf_bet\n");
		fprintf(out,"%lf\t%d\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\n",t,env_type,pop_size[i_ris],pop_size[i_bet],total_pop,ops,abf_ris,abf_bet);
		for (i=1; i<= pop_size[i_ris]; ++i){	// initialize individuals of the rising-tide strategy
			id_ris+= 1;
			ris_pop[i][ID]= id_ris;
			ris_pop[i][state]= 1;
			ris_pop[i][birth]= t;
			ris_pop[i][parent]= 0;
			}
		for (i=1; i<= pop_size[i_bet]; ++i){	// initialize individuals of the bet-hedging strategy
			id_bet+= 1;
			bet_pop[i][ID]= id_bet;
			bet_pop[i][state]= 1;
			bet_pop[i][birth]= t; 
			bet_pop[i][parent]= 0;
			}

	// Main loop
	while (t< T){
        // high/ low quality environment determination
        if (cc<=0){				
			if(s3==1){
				env_type= 0;
				pp= dsfmt_genrand_open_open(&dsfmt);
				if(pp< r1){
					env_type= 1;		// high quality environment starts
					cc= T_h;}
				else{
					env_type= 2;		// low quality environment starts
					cc= T_l;}
				}
			if(s3==0){
				if(env_type== 0||env_type== 2){
					env_type= 1;		// high quality environment starts
					cc= T_h*2*r1;}
				else{
					env_type= 2;		// low quality environment starts
					cc= T_l*(2-2*r1);}
				}
            }
        // sum of probability in high quality environment
        if (env_type== 1){
			total_prob= r_high_ris+ r_high_bet+ r_high_ris/ K_high_ris* pop_size[i_ris]+ alpha_high_br* r_high_ris/ K_high_ris* pop_size[i_bet]+ r_high_bet/ K_high_bet* pop_size[i_bet]+ alpha_high_rb* r_high_bet/ K_high_bet* pop_size[i_ris];
			if(pop_size[i_ris]>0){
				p1= r_high_ris/ total_prob;
			}
            else p1= 0;
            if(pop_size[i_bet]>0){
                p2= p1+ r_high_bet/ total_prob;
            }
            else p2= p1;
			if(pop_size[i_ris]>0){
				p3= p2+ (r_high_ris/ K_high_ris* pop_size[i_ris]+ alpha_high_br* r_high_ris/ K_high_ris* pop_size[i_bet])/ total_prob;
			}
			else p3= p2;
            // Calculating absolute fitness
            abf_ris= r_high_ris- r_high_ris/ K_high_ris* pop_size[i_ris]- alpha_high_br* r_high_ris/ K_high_ris* pop_size[i_bet];
            abf_bet= r_high_bet- r_high_bet/ K_high_bet* pop_size[i_bet]- alpha_high_rb* r_high_bet/ K_high_bet* pop_size[i_ris];
		}
        // sum of probability in low quality environment
        if (env_type== 2){
            total_prob= r_low_ris+ r_low_bet+ r_low_ris/ K_low_ris* pop_size[i_ris]+ alpha_low_br* r_low_ris/ K_low_ris* pop_size[i_bet]+ r_low_bet/ K_low_bet* pop_size[i_bet]+ alpha_low_rb* r_low_bet/ K_low_bet* pop_size[i_ris];
			if(pop_size[i_ris]>0){
	            p1= r_low_ris/ total_prob;
			}
			else p1= 0;
            if(pop_size[i_bet]>0){
                p2= p1+ r_low_bet/ total_prob;
            }
            else p2= p1;
			if(pop_size[i_ris]>0){
				p3= p2+ (r_low_ris/ K_low_ris* pop_size[i_ris]+ alpha_low_br* r_low_ris/ K_low_ris* pop_size[i_bet])/ total_prob;
			}
			else p3= p2;
            // Calculating  absolute fitness
            abf_ris= r_low_ris- r_low_ris/ K_low_ris* pop_size[i_ris]- alpha_low_br* r_low_ris/ K_low_ris* pop_size[i_bet];
            abf_bet= r_low_bet- r_low_bet/ K_low_bet* pop_size[i_bet]- alpha_low_rb* r_low_bet/ K_low_bet* pop_size[i_ris];
		}

		pp= dsfmt_genrand_open_close(&dsfmt);			    // random number
        S= -1*log(pp)/total_prob/scale_time;				// duration of time advanced in this event (only one event happens per calculation)
        // Error handling
        if (S>0) ;
        else printf("%lf\n",S);
        while (S!=S){
			S= -1*log(pp)/total_prob/scale_time;			// resample for S
			printf("Encountering NAN in S.\n");
        }
    
	// individual level events
		// log matrix intialization
		pp= dsfmt_genrand_open_open(&dsfmt);			    // random number to decide which event to occur
		if(pp<p1&& pop_size[i_ris]>0) {				// --- Birth of rising-tide individual ---
			pop_size[i_ris]= pop_size[i_ris]+1;
			id_ris+= 1;
			ris_pop[id_ris][ID]= id_ris;
			ris_pop[id_ris][state]= 1;
			ris_pop[id_ris][birth]= t;
			// decide the parent
				pp2= pop_size[i_ris]*dsfmt_genrand_open_open(&dsfmt);
				rid= round(pp2);
				idc= 0;
				for(i=1; i<= id_ris;++i){
					if(ris_pop[i][state]> 0) idc+= 1;      // count the living individual
					if(idc== rid) {
                        ris_pop[id_ris][parent]= ris_pop[i][ID];
                        if (idc<= p4) offspring[idc]+= 1;
						i= id_ris;
					}
				}
			}
		else{
			if(pp<p2&& pop_size[i_bet]>0) {		    // --- Birth of bet-hedging individual ---
				pop_size[i_bet]= pop_size[i_bet]+1;
				id_bet+= 1;
				bet_pop[id_bet][ID]= id_bet;
				bet_pop[id_bet][state]= 1;
				bet_pop[id_bet][birth]= t;
				// decide the parent
					pp2= pop_size[i_bet]*dsfmt_genrand_open_open(&dsfmt);
					rid= round(pp2);
					idc= oct= 0;
					for(i=1; i<= id_bet;++i){
						if(bet_pop[i][state]>0) idc+= 1;
						if(idc==rid) {
							bet_pop[id_bet][parent]= bet_pop[i][ID];
                            oct= idc+ p4;
                            if(oct<= total_pop) offspring[oct]+= 1;
							i= id_bet;
						}
					}
			    }
			else{
				if(pp<p3&& pop_size[i_ris]>0) { 		// --- Death of rising-tide individual ---
					pp2= pop_size[i_ris]*dsfmt_genrand_open_open(&dsfmt);
					rid= round(pp2);
					idc= 0;
					for(i=1; i<= id_ris; ++i){
						if(ris_pop[i][state]> 0) idc+= 1;
						if(idc== rid){
							ris_pop[i][state]= 0;
							ris_pop[i][death]= t;
							ris_pop[i][span]= t- ris_pop[i][birth];
							i= id_ris;
						}
					}
					pop_size[i_ris]= pop_size[i_ris]-1;
					}
				else {									    // --- Death of bet-hedging individual ---
					if(pop_size[i_bet]>0){
						pp2= pop_size[i_bet]*dsfmt_genrand_open_open(&dsfmt);
						rid= round(pp2);
						idc= 0;
						for(i=1; i<= id_bet; ++i){
							if(bet_pop[i][state]> 0) idc+= 1;
							if(idc== rid){
								bet_pop[i][state]= 0;
								bet_pop[i][death]= t;
								bet_pop[i][span]= t- bet_pop[i][birth];
								i= id_bet;
							}
						}
						pop_size[i_bet]= pop_size[i_bet]-1;
					}
					}
			}}
        // Updating time and population size logger
		t+= S;
        cc-= S;
        ops_count-= S;
        int tp_size= pop_size[i_ris]+pop_size[i_bet];
        // Calculating opportunity for selection
        /*
        - Define the starting point and calaulate the population size
        - temporally store the number of offspring produced of each individual
        - determine whether it is the end of the period for calculating ops
        - if yes, calculate and fix ops to there (it will be printed in each dynamic step), and recalculate the population size
        */
        if(ops_count<=0){                                   // updating opportunity for selection
            p6= mean(offspring,total_pop);                  // average of fitness
            for(i=1; i<= total_pop; i++){
                refitness[i]= offspring[i]/p6;
            }
            ops= var(refitness,total_pop);
            // reset the space
            total_pop= pop_size[i_ris]+ pop_size[i_bet];	// update population size
            offspring= (double *) realloc(offspring,(size_t)((total_pop+1)*sizeof(double)));
            refitness= (double *) realloc(refitness,(size_t)((total_pop+1)*sizeof(double)));
            for(i=1; i<= total_pop; i++) offspring[i]= 0;	// fills with zero
            p4= pop_size[i_ris];                       // int double comparison?
            ops_count= ops_span;
        }
        fprintf(out,"%lf\t%d\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\n",t,env_type,pop_size[i_ris],pop_size[i_bet],tp_size,ops,abf_ris,abf_bet);

	}
	// Print
		fprintf(ris_log,"ID\tstate\tbirth\t\tdeath\t\tparent\tlife span\n");
		fprintf(bet_log,"ID\tstate\tbirth\t\tdeath\t\tparent\tlife span\n");
		for (i=1; i<=id_ris;++i){
			o1= ris_pop[i][ID];
			o2= ris_pop[i][state];
			o3= ris_pop[i][parent];
			fprintf(ris_log,"%d\t%d\t%lf\t%lf\t%d\t%lf\n",o1,o2,ris_pop[i][birth],ris_pop[i][death],o3,ris_pop[i][span]);
		}
		for (i=1; i<=id_bet;++i){
			o1= bet_pop[i][ID];
			o2= bet_pop[i][state];
			o3= bet_pop[i][parent];
			fprintf(bet_log,"%d\t%d\t%lf\t%lf\t%d\t%lf\n",o1,o2,bet_pop[i][birth],bet_pop[i][death],o3,bet_pop[i][span]);
		}

	free_d_vector(pop_size);
	free_d_vector(dfdt);
	free_d_matrix(bet_pop);
    free_d_matrix(ris_pop);
    free_d_vector(offspring);
    free_d_vector(refitness);
	fclose(out);
	fclose(ris_log);
	fclose(bet_log);
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
{return (r_high_ris- r_high_ris/K_high_ris*vr[i_ris]- r_high_ris/K_high_ris*alpha_high_br*vr[i_bet])*vr[i_ris];}
double dy_dt_h(double time, double vr[])	// bet-hedging in high quality environment
{return (r_high_bet- r_high_bet/K_high_bet*vr[i_bet]- r_high_bet/K_high_bet*alpha_high_rb*vr[i_ris])*vr[i_bet];}
void differential_h(double time, double in[], double out[])
{
	out[i_ris]= dx_dt_h(time,in);
	out[i_bet]= dy_dt_h(time,in);
}
double dx_dt_l(double time, double vr[])	// rising-tide in low quality environment
{return (r_low_ris- r_low_ris/K_low_ris*vr[i_ris]- r_low_ris/K_low_ris*alpha_low_br*vr[i_bet])*vr[i_ris];}
double dy_dt_l(double time, double vr[])	// bet-hedging in low quality environment
{return (r_low_bet- r_low_bet/K_low_bet*vr[i_bet]- r_low_bet/K_low_bet*alpha_low_rb*vr[i_ris])*vr[i_bet];}
void differential_l(double time, double in[], double out[])
{
	out[i_ris]= dx_dt_l(time,in);
	out[i_bet]= dy_dt_l(time,in);
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
void free_d_vector(double *x) {	
    free(x);
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
	return t2/(length-1);
}