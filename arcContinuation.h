#include <stdlib.h>
#include <math.h>
#define LAPACKDOUBLE
//#define PRINTNEWROOT
#define Pi M_PI
#define LOGFILE

#ifdef LAPACKDOUBLE
typedef double val_type;
#elif LAPACKCOMPLEX
typedef double complex val_type;
#endif


typedef double time_type;

val_type Power(val_type v, double exp);
val_type Sqrt(val_type v);
val_type Tan(val_type v);
val_type Cot(val_type v);


double Norm(int dim, val_type *v);


#ifdef MATLAB
  #include <mex.h>
  typedef ptrdiff_t lapack_type;
#else
  typedef int lapack_type;
#endif

typedef struct
    {
    int maxnewtoniter;
    int maxinnersteps;
    int maxstepattempts; /* this value depends on maxinnersteps, but at first steps it is a higher value */
    int prec;       // when Newton-method does not get the accuracy required in maxstepnumber iterations
                // we see if the norm of the increment is less than pow(10,prec) in which case 
                // we accept the solution as new root 
    int acc;
    int iter_opt;   // optimun number of iterations in the newton method aplication. 
		// The deltaT will be adecuated depending in this value. The goal is to get the deltaT 
		// that produces convergence in iter_opt iterations.
    int maxnumstepswithsamejac;         // 1 --> the same jacobian will be used in next step (2 steps with same jacobian)
				    // 2--> it will try to use two more times the same jacobian
                                    // 0 --> allways will get the new jacobian
    int numstepswithsamejac;  // var to be modified at each Newton iteration step
    int maxdeltacondfails;
    int approxjacobian;
    int numjac;
    double stepsize;  // The user can set the initial stepsize. Starting from t0 the process goes to t1
		  // and the initial stepsize can be set in the .conf file
/*
# the mecanism to control the minimun stepsize:
#    If the average_stepsize goes beyond min_av_stepsize the process will stop.
#    average_stepsize = \Sum_i^s stepsize_{-i}/s
#    where s is the number of steps taken into account for the estimation. (if s == 1 then 
#    the stepsize can not be beyond min_av_stepsize) s can be set in the configuration file:
#          num_steps_for_stepsize_control: is the number of steps taken into account to get the
#                average_stepsize
#          min_av_stepsize: The minimun average_stepsize acepted. Format #.#e#
*/
    double min_av_stepsize;
    int num_steps_for_stepsize_control;
    int extra_eqs;
    double iterationdecreasefactor; // the factor indicating the minimum decrease of the delta at each Newton iteration (0.8)
    int max_steps;  // the maximum number of steps to be needed to avance from tinit until tend (400) 
    double closure_tolerance;
    } option_struct;


