#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <arcContinuation.h>
#include <mathlink.h>

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



void return_values_2(val_type *vptr, int n, double norma3, double *dptr, long dnum, int *iptr, long inum)
{
if (!MLPutFunction(stdlink, "List",2))
    {
    MLErrorMessage(stdlink);
    };
if (!MLPutReal64List( stdlink, vptr, n)) 
    {
    MLErrorMessage(stdlink);
    }
if (!MLPutReal64( stdlink,norma3)) 
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

void return_matrix(int dim1, int dim2, val_type *m)
{   
int dims[2];

dims[0] = dim1;
dims[1] = dim2;

if(! MLPutReal64Array(stdlink, (double *)m, (int *)dims, (const char **)0, 2))
        { /* unable to send the double array to lp */ }
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
    fprintf(loga,"malloc problem! tfinal = %lf\n",init);
#endif
    }
if (info && ((b1 = (double *)malloc ((dim)*sizeof(double)))==NULL))  // Vector dim n
    {
    info = 0;
#ifdef LOGFILE
    fprintf(loga,"malloc problem! tfinal = %lf\n",init);
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





void realRest(double *v, long vlen, double t, double *preal, long prlen, int *pint, long pilen)
{
val_type *bal1;
double norma3;
int dim;
int extra_eqs;

extra_eqs=extra_len();
dim = vlen;
if ((bal1 = (val_type *)malloc ((dim+extra_eqs)*sizeof(val_type)))==NULL)  // Vector dim n+extra_len
        {
        return;
        }
rest(dim, (val_type *)v, t, bal1, preal, prlen, pint, pilen);             // bal1 = F(t,v)
norma3 = Norm(dim,bal1);
return_values_2(bal1,dim+extra_eqs,norma3,preal,prlen,pint,pilen);
free(bal1);
}

void realJacobian(double *v, long vlen, double t, double *preal, long prlen, int *pint, long pilen)
{
val_type *jac;
int dim;

#ifdef LOGFILE
loga = fopen("./realJacobian.log","w");
fprintf(loga,"realJacobian-en hasiera...\n");
#endif
dim = vlen;
if ((jac = (val_type *)malloc ((dim*dim)*sizeof(val_type)))==NULL)  // Vector dim n
        {
#ifdef LOGFILE
fprintf(loga,"end of realJacobian! (no malloc)\n");
fclose(loga);
#endif
        return;
        }
jacobian(dim,(val_type *)v,t,jac, preal, prlen, pint, pilen);
return_matrix(dim,dim,jac);
free(jac);
#ifdef LOGFILE
fprintf(loga,"end of realJacobian!\n");
fclose(loga);
#endif
}

void realFindroot(double *v, long vlen, double t, double *preal, long prlen, int *pint, long pilen, double mathacc, double mathprec, int mathapproxjac)
{
val_type *jac;
val_type *bal1;
val_type *delta;
int *ipiv;
int dim;
int rootfound;
double deltafactor;
int extra_eqs;
double last_value;
option_struct options;

extra_eqs = extra_len(); /* extra_len returns "Length[rest] - Length[vars]" because it is posible to have more equations than variables */
dim = vlen; 
options.acc = mathacc;
options.prec = mathprec;
#ifdef LOGFILE
loga = fopen("./findroot.log","w");
fprintf(loga,"math realFindroot-en hasiera...\n");
#endif
options.approxjacobian = mathapproxjac;
if ((bal1 = (val_type *)malloc ((dim)*sizeof(val_type)))==NULL)  // Vector dim n 
        {
#ifdef LOGFILE
fprintf(loga,"end of realFindroot! (no malloc)\n");
fclose(loga);
#endif
        return;
        }
if ((jac = (val_type *)malloc ((dim*dim)*sizeof(val_type)))==NULL)  // Matrix dim n x n
        {
#ifdef LOGFILE
fprintf(loga,"end of realFindroot! (no malloc)\n");
fclose(loga);
#endif
        free(bal1);
        return;
        }

if ((delta = (val_type *)malloc ((dim+extra_eqs)*sizeof(val_type)))==NULL)  // Vector dim n+ m (m is the number of extra conditions)
        {
#ifdef LOGFILE
fprintf(loga,"end of realFindroot! (no malloc)\n");
fclose(loga);
#endif
        free(bal1);
	free(jac);
        return;
        }
if ((ipiv = (int *)malloc (dim*sizeof(int)))==NULL)  // Vector dim n
        {
#ifdef LOGFILE
fprintf(loga,"end of realFindroot! (no malloc)\n");
fclose(loga);
#endif
        free(bal1);
	free(jac);
        free(delta);
        return;
        }
options.maxnumstepswithsamejac = 0;
options.numstepswithsamejac = 0;
options.maxdeltacondfails = 2;
options.numjac = 0;
options.maxnewtoniter = 20;
/*
int findrootbysnewton(int dim, double *jac0, lapack_type *ipiv, double t, double *bal0, double *bal1, double *delta, double *fnptr, double *ba, double *xa, double ta,
		      double *preal, long numpreal, int *pint, long numpint, double *factorptr, option_struct *optr)
*/
rootfound = findrootbysnewton(dim, jac,ipiv, t, v, &(bal1[0]),&(delta[0]),&last_value,(double *)0,(double *)0,(double)0.0, preal,prlen,pint,pilen, &deltafactor, &options);
return_values_2(bal1,dim,(double) rootfound,preal,prlen,pint,pilen);
free(jac);
free(ipiv);
free(bal1);
free(delta);
#ifdef LOGFILE
fprintf(loga,"end of realFindroot!\n");
fclose(loga);
#endif
}





int main(int argc, char *argv[])
{
return MLMain(argc, argv);
}



