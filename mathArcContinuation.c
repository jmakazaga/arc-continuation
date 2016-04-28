#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <homotopia.h>
#include <mathlink.h>
typedef int lapack_type;

#ifdef LOGFILE
extern FILE *loga;
#endif

void return_values(int dim, double *vptr, double *bptr, double sa, double sb, int info)
{
if (!MLPutFunction(stdlink, "List",5))
    {
    MLErrorMessage(stdlink);
    };
if (!MLPutReal64List( stdlink, vptr, dim)) 
    {
    MLErrorMessage(stdlink);
    }
if (!MLPutReal64List( stdlink, bptr, dim)) 
    {
    MLErrorMessage(stdlink);
    }
if (!MLPutReal64( stdlink,sa)) 
    {
    MLErrorMessage(stdlink);
    }
if (!MLPutReal64( stdlink,sb)) 
    {
    MLErrorMessage(stdlink);
    }
if (!MLPutInteger( stdlink,info)) 
    {
    MLErrorMessage(stdlink);
    }
if(! MLEndPacket(stdlink))
    {
    MLErrorMessage(stdlink);
    }
if(! MLFlush(stdlink))
    {
    MLErrorMessage(stdlink);
    }
}

void mathArcContinuation( double *var0, long n, double *b0, long nb, double init, double end, int action, double *preal, long numpreal, int *pint, long numpint)
{
/* needed at output */ 
double *var1;
double *b1;
double ta;
double tb;
int  info;
/* local vars */
int i, dim;

info = 1;
/* test n == nb */
if (n != nb) info = 0;
  else dim = n;
/* allocate memory */
if (info && ((var1 = (double *)malloc ((dim)*sizeof(double)))==NULL))  // Vector dim n
    {
    info = 0;
#ifdef LOGFILE
    fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
    }
if (info && ((b1 = (double *)malloc ((dim)*sizeof(double)))==NULL))  // Vector dim n
    {
    info = 0;
#ifdef LOGFILE
    fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
    }
/********************/
/* call to function */
/********************/
arcContinuation(dim, var0, b0, init, end, action, preal, numpreal, pint, numpint, var1, b1, &ta, &tb, &info);

/* return values */
return_values(dim, var1, b1, ta,tb, info);

/* free memory */
free(var1);
free(b1);
}

