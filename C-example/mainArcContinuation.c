#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double val_type;
typedef double time_type;

extern void arcContinuation(/* input */int dim, double *var0, double *b0, double init, double end, int action, double *preal, long numpreal, int *pint, long numpint,
		     /* output */ double *var1, double *b1, double *taptr, double *tbptr, int * infoptr);


extern double Power(double x, double exp);
/*
{
return(pow(x,exp));
}
*/

/*
res = {x^2 + y^2 + z^2 - 1, x^3 + y^3 + 1/2 z^3};
glist = {x^5 + y^5 + 1/2 z^5};
*/

void rest(int dim, val_type *var, time_type t, val_type *res, double *preal, long prelem, int *pint, long pielem)
{

res[0]=-1 + Power(var[0],2) + Power(var[1],2) + Power(var[2],2);
res[1]=Power(var[0],3) + Power(var[1],3) + Power(var[2],3)/2.;
res[2]=Power(var[0],5) + Power(var[1],5) + Power(var[2],5)/2.;
}

void jacobian(int dim, val_type *var, time_type t, val_type *jac, double *preal, long prelem, int *pint, long pielem)
{

//jac[0] = D[res[0],var[0]]
jac[0] = 2*var[0];
//jac[1] = D[res[1],var[0]]
jac[1] = 3*Power(var[0],2);
//jac[2] = D[res[2],var[0]]
jac[2] = 5*Power(var[0],4);
//jac[3] = D[res[0],var[1]]
jac[3] = 2*var[1];
//jac[4] = D[res[1],var[1]]
jac[4] = 3*Power(var[1],2);
//jac[5] = D[res[2],var[1]]
jac[5] = 5*Power(var[1],4);
//jac[6] = D[res[0],var[2]]
jac[6] = 2*var[2];
//jac[7] = D[res[1],var[2]]
jac[7] = (3*Power(var[2],2))/2.;
//jac[8] = D[res[2],var[2]]
jac[8] = (5*Power(var[2],4))/2.;
}

int main(int argc, char *argv[])
{
int dim,i;
double *xptr,*nxptr, *f;
double *betaptr,*nbetaptr;
double tinit,tend, ntinit, ntend;
int action;  // 0 -> go to tend,  1-> find x value that minimizes f_n(x) 2-> find x that makes f_n(x) = 0
long prealkop, pintkop;
double *preal;
int *pint;
int info;

dim = 3;
xptr = (double *)malloc(dim*sizeof(double));
betaptr = (double *)malloc(dim*sizeof(double));
nxptr = (double *)malloc(dim*sizeof(double));
nbetaptr = (double *)malloc(dim*sizeof(double));
f = (double *)malloc(dim*sizeof(double));

// action 1 -> looks for x that minimize last f(x)
// action 2 -> looks for x that makes f(x) = 0
// action 0 -> just moves on the implicit curve
action = 2; 
// initial values of the parameters:
xptr[0]= 0.470369;
xptr[1]= 0.470369;
xptr[2]= -0.746664;
// different initial values (no min)
//xptr[0]= -0.621391; xptr[1]=  -0.0298134;xptr[2]=  0.782904;
printf("This example shows how to use the code (function ArcContinuation) to move on an implicit curve\n");
printf("We have three variables and two equations:\n");
printf(" the equations are:\n f1(x) = x^2 + y^2 + z^2 - 1\n f2(x)= x^3 + y^3 + 1/2 z^3\n");

printf("We have a root of the functions. The values of the initial root (values of x) are:\n");
for (i = 0; i<dim; i++)
    printf("x[%d] = %lf\n", i,xptr[i]);
printf("we have an extra function,(f3(x) = x^5 + y^5 + 1/2 z^5) and we can move on the implicit curve while testing the value of the extra function\n");
// call to function that evaluates f1, f2 and the extra condition
rest(dim, xptr, tinit, f, (double *)0, 0L, (int *)0, 0L);
printf("\nthe values of the functions with initial x are:\n");
for (i = 0; i<dim-1; i++)
    printf("  f%d(x) = %.10lf\n", i,f[i]);
    
printf("The value for the extra function for this x vector is\n  f%d(x)%.10lf\n\n", dim-1,f[dim-1]);
printf("\nWe are interested on points on the implicit curve that makes 0 the extra condition, so we call arcContinuation with action 2\n");
betaptr[0]=1.0;
betaptr[1]=0.0;
betaptr[2]=0.0;
tinit = 0.0;
tend = 10.0;

arcContinuation(dim, xptr, betaptr, tinit, tend, action, 
                preal, prealkop, pint, pintkop,
		nxptr, nbetaptr, &ntinit, &ntend, &info);
		
printf("return value = %d\n",info);
		
printf("\nAfter the continuation with action = %d x values are:\n",action);
for (i = 0; i<dim; i++)
    printf("x[%d] = %lf ,", i,nxptr[i]);
printf("\n    and the values of the functions with new x are:\n");
rest(dim, nxptr, tinit, f, (double *)0, 0L, (int *)0, 0L);
for (i = 0; i<dim-1; i++)
    printf("   f%d(x) = %.10lf ,", i,f[i]);
printf("\n");
    
printf("The value for the extra function for this x vector is\n  f%d(x)%.10lf\n", dim-1,f[dim-1]);

}

