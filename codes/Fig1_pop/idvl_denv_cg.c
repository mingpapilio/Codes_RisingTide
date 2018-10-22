/*
 * Individual-based model with continuous growth in discrete environemnts (good/ bad years)
 * First day of building: July 19, 2017
 * Modified from competitive Lotka-Volterra differential equations
 * Last modified: Oct 22, 2018
 ***************************************************************************************
Execution:
gcc idvl_bet2_1010.c -o a.out
./a.out
 ***************************************************************************************
 * Key parameters
 * The global parameters of differential equations (e.g. r_goodyear_spcl)
 * r1: probability of getting good year
 * s3: switch of temporal stochasticity
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
		int i,j;				// For loop counters
		double r1=      0.3; 	// Chance of getting a good year
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
		int id_spcl=    0.0;    // ID of the specialist individuals
		int id_genl=    0.0;    // ID of the generalist individuals
	// Temporal settings
		double t=       0.0;	// Time logs
		int T=          15000; 	// Duration of simulation (unit is different with the deterministic model)
		double deltat=  0.005; 	// Length of time step
		int T_g=        200; 	// Length of a good year ***Note the difference of double and int
		int T_b=        200; 	// Length of a bad year ***Note the difference of double and int
        double cc=      0.0;	// Counter of remaining time of a year
        int scale_time= 10;     // Control the pace of each dynamic step between events (a coefficient)
    // Opportunity for selection settings
		double ops=     0.0;	// Temp of opportunity of selection
		double abf_spcl=   0.0;	// Absolute fitness (specialist)
        double abf_genl=   0.0;	// Absolute fitness (generalist)
        double ops_span= 200.0;     // Duration for updating opportunity for selection
        double ops_count= ops_span; // Time counter
	// Print temporal spaces
		int o1, o2, o3, o4;
	// Creating temporal space
		double *pop_size= d_vector(2);
		double *dfdt= d_vector(2);
		double **spcl_pop= d_matrix(pop_limit, 6);
		// 1= ID, 2= live/dead, 3= birth time, 4= death time, 5= parent (haploid), 6=life span
		double **genl_pop= d_matrix(pop_limit, 6);
			int ID= 1;
			int state= 2;
			int birth= 3;
			int death= 4;
			int parent= 5;
			int span= 6;
		//double **offspring= d_matrix(total_pop,2);
            int oct= 0;			// offspring counter (for calculating fitness)
	// Output
		FILE *out, *spcl_log, *genl_log;
		out= fopen("pop_dmc.txt","w");
		spcl_log= fopen("idv_spcl_log.txt","w");
		genl_log= fopen("idv_genl_log.txt","w");
	// Initialization
		pop_size[specialist]= 500;					// initial population of specialist strategy
        pop_size[generalist]= 500;					// initial population of generalist strategy
            int tp_size= pop_size[specialist]+pop_size[generalist];
            total_pop= pop_size[specialist]+pop_size[generalist];
            double *offspring= d_vector(total_pop);
            double *refitness= d_vector(total_pop);
            p4= pop_size[specialist];
		//fprintf(out,"r1 is %lf, r2 is %lf, db is %lf, duration is %d\n",b_g1,b_g2,d_b2,T);
		fprintf(out,"time\t\tseason\tspecialist\tgeneralist\tTotal\tOpp of sel\tabf_spcl\tabf_genl\n");
		fprintf(out,"%lf\t%d\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\n",t,s1,pop_size[specialist],pop_size[generalist],total_pop,ops,abf_spcl,abf_genl);
		for (i=1; i<= pop_size[specialist]; ++i){	// initialize individuals of the specialist strategy
			id_spcl+= 1;
			spcl_pop[i][ID]= id_spcl;
			spcl_pop[i][state]= 1;
			spcl_pop[i][birth]= t;
			spcl_pop[i][parent]= 0;
			}
		for (i=1; i<= pop_size[generalist]; ++i){	// initialize individuals of the generalist strategy
			id_genl+= 1;
			genl_pop[i][ID]= id_genl;
			genl_pop[i][state]= 1;
			genl_pop[i][birth]= t; 
			genl_pop[i][parent]= 0;
			}

	// Main loop
	while (t< T){
        // good/ bad year determination
        if (cc<=0){				
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
        // sum of probability in good year
        if (s1== 1){
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
            // Calculating absolute fitness
            abf_spcl= r_goodyear_spcl- r_goodyear_spcl/ K_goodyear_spcl* pop_size[specialist]- alpha_goodyear_gs* r_goodyear_spcl/ K_goodyear_spcl* pop_size[generalist];
            if(pop_size[generalist]>0){
                abf_genl= r_goodyear_genl*scale_c*(pop_size[generalist]/(shape_c+pop_size[generalist])+shift_c)- r_goodyear_genl/ K_goodyear_genl* pop_size[generalist]- alpha_goodyear_sg* r_goodyear_genl/ K_goodyear_genl* pop_size[specialist];
            }
            else abf_genl= r_goodyear_genl*scale_c*(shift_c)- alpha_goodyear_sg* r_goodyear_genl/ K_goodyear_genl* pop_size[specialist];
            }
        // sum of probability in bad year
        if (s1== 2){
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
            // Calculating  absolute fitness
            abf_spcl= r_badyear_spcl- r_badyear_spcl/ K_badyear_spcl* pop_size[specialist]- alpha_badyear_gs* r_badyear_spcl/ K_badyear_spcl* pop_size[generalist];
            if(pop_size[generalist]>0){
                abf_genl= r_badyear_genl*scale_c*(pop_size[generalist]/(shape_c+pop_size[generalist])+shift_c)- r_badyear_genl/ K_badyear_genl* pop_size[generalist]- alpha_badyear_sg* r_badyear_genl/ K_badyear_genl* pop_size[specialist];
            }
            else abf_genl= r_badyear_genl*scale_c*(shift_c)- alpha_badyear_sg* r_badyear_genl/ K_badyear_genl* pop_size[specialist];
			}

		pp= dsfmt_genrand_close_open(&dsfmt);			    // random number
        S= -1*log(pp)/total_prob/scale_time;				// duration of time advanced in this event (only one event happens per calculation)
        // Error habdling
        if (S>0) ;
        else printf("%lf\n",S);
        while (S!=S){
			S= -1*log(pp)/total_prob/scale_time;			// resample for S
			printf("Encountering NAN in S.\n");
        }
    
	// individual level events
		// log matrix intialization
		pp= dsfmt_genrand_open_open(&dsfmt);			    // random number to decide which event to occur
		if(pp<p1&& pop_size[specialist]>0) {				// --- Birth of specialist individual ---
			pop_size[specialist]= pop_size[specialist]+1;
			id_spcl+= 1;
			spcl_pop[id_spcl][ID]= id_spcl;
			spcl_pop[id_spcl][state]= 1;
			spcl_pop[id_spcl][birth]= t;
			// decide the parent
				pp2= pop_size[specialist]*dsfmt_genrand_open_open(&dsfmt);
				rid= round(pp2);
				idc= 0;
				for(i=1; i<= id_spcl;++i){
					if(spcl_pop[i][state]> 0) idc+= 1;      // count the living individual
					if(idc== rid) {
                        spcl_pop[id_spcl][parent]= spcl_pop[i][ID];
                        if (idc<= p4) offspring[idc]+= 1;
						i= id_spcl;
					}
				}
			}
		else{
			if(pp<p2&& pop_size[generalist]>0) {		    // --- Birth of generalist individual ---
				pop_size[generalist]= pop_size[generalist]+1;
				id_genl+= 1;
				genl_pop[id_genl][ID]= id_genl;
				genl_pop[id_genl][state]= 1;
				genl_pop[id_genl][birth]= t;
				// decide the parent
					pp2= pop_size[generalist]*dsfmt_genrand_open_open(&dsfmt);
					rid= round(pp2);
					idc= oct= 0;
					for(i=1; i<= id_genl;++i){
						if(genl_pop[i][state]>0) idc+= 1;
						if(idc==rid) {
							genl_pop[id_genl][parent]= genl_pop[i][ID];
                            oct= idc+ p4;
                            if(oct<= total_pop) offspring[oct]+= 1;
							i= id_genl;
						}
					}
			    }
			else{
				if(pp<p3&& pop_size[specialist]>0) { 		// --- Death of specialist individual ---
					pp2= pop_size[specialist]*dsfmt_genrand_open_open(&dsfmt);
					rid= round(pp2);
					idc= 0;
					for(i=1; i<= id_spcl; ++i){
						if(spcl_pop[i][state]> 0) idc+= 1;
						if(idc== rid){
							spcl_pop[i][state]= 0;
							spcl_pop[i][death]= t;
							spcl_pop[i][span]= t- spcl_pop[i][birth];
							i= id_spcl;
						}
					}
					pop_size[specialist]= pop_size[specialist]-1;
					}
				else {									    // --- Death of generalist individual ---
					if(pop_size[generalist]>0){
						pp2= pop_size[generalist]*dsfmt_genrand_open_open(&dsfmt);
						rid= round(pp2);
						idc= 0;
						for(i=1; i<= id_genl; ++i){
							if(genl_pop[i][state]> 0) idc+= 1;
							if(idc== rid){
								genl_pop[i][state]= 0;
								genl_pop[i][death]= t;
								genl_pop[i][span]= t- genl_pop[i][birth];
								i= id_genl;
							}
						}
						pop_size[generalist]= pop_size[generalist]-1;
					}
					}
			}}
        /*
        - Define the starting point and calaulate the population size
        - temporally store the number of offspring produced of each individual
        - determine whether it is the end of the period for calculating ops
        - if yes, calculate and fix ops to there (it will be printed in each dynamic step), and recalculate the population size
        */
        // Updating time and population size logger
		t+= S;
        cc-= S;
        ops_count-= S;
        int tp_size= pop_size[specialist]+pop_size[generalist];
        // Calculating opportunity for selection
        if(ops_count<=0){                                   // updating opportunity for selection
            p6= mean(offspring,total_pop);                  // average of fitness
            for(i=1; i<= total_pop; i++){
                refitness[i]= offspring[i]/p6;
            }
            ops= var(refitness,total_pop);
            // reset the space
            total_pop= pop_size[specialist]+ pop_size[generalist];	// update population size
            offspring= (double *) realloc(offspring,(size_t)((total_pop+1)*sizeof(double)));
            refitness= (double *) realloc(refitness,(size_t)((total_pop+1)*sizeof(double)));
            for(i=1; i<= total_pop; i++) offspring[i]= 0;	// fills with zero
            p4= pop_size[specialist];                       // int double comparison?
            ops_count= ops_span;
        }
        fprintf(out,"%lf\t%d\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\n",t,s1,pop_size[specialist],pop_size[generalist],tp_size,ops,abf_spcl,abf_genl);

	}
	// Print
		fprintf(spcl_log,"ID\tstate\tbirth\t\tdeath\t\tparent\tlife span\n");
		fprintf(genl_log,"ID\tstate\tbirth\t\tdeath\t\tparent\tlife span\n");
		for (i=1; i<=id_spcl;++i){
			o1= spcl_pop[i][ID];
			o2= spcl_pop[i][state];
			o3= spcl_pop[i][parent];
			fprintf(spcl_log,"%d\t%d\t%lf\t%lf\t%d\t%lf\n",o1,o2,spcl_pop[i][birth],spcl_pop[i][death],o3,spcl_pop[i][span]);
		}
		for (i=1; i<=id_genl;++i){
			o1= genl_pop[i][ID];
			o2= genl_pop[i][state];
			o3= genl_pop[i][parent];
			fprintf(genl_log,"%d\t%d\t%lf\t%lf\t%d\t%lf\n",o1,o2,genl_pop[i][birth],genl_pop[i][death],o3,genl_pop[i][span]);
		}

	free_d_vector(pop_size);
	free_d_vector(dfdt);
	free_d_matrix(genl_pop);
    free_d_matrix(spcl_pop);
    free_d_vector(offspring);
    free_d_vector(refitness);
	fclose(out);
	fclose(spcl_log);
	fclose(genl_log);
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