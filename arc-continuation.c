#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <arcContinuation.h>


#ifdef SOLVER_GSL
  #include <gsl/gsl_linalg.h>
#endif

#ifdef LOGFILE
FILE *loga;
#endif

    
double Norm(int dim, val_type *v);



//copy the vector pointed by oriptr into destptr
void vectorcp(int dim, val_type *oriptr, val_type *destptr);

//Normalizes the vector
void vectornormalize(int dim, val_type *vptr)
{
double length;

length = Norm(dim,vptr);

while( dim > 0 ) 
  {
  dim--;
  vptr[dim] = vptr[dim]/length;
  }
}
//copy the vector pointed by oriptr into destptr, but it normalizes the vector
void vectornormalizecp(int dim, val_type *oriptr, val_type *destptr)
{
double length;

length = Norm(dim,oriptr);

while( dim > 0 ) 
  {
  dim--;
  destptr[dim] = oriptr[dim]/length;
  }
}


#ifdef LOGFILE
void printmat(int dim, double *jac)
{
int i,j;

fprintf(loga,"\n");
for (i =0; i<dim; i++)
  {
  for (j=0; j<dim; j++) fprintf(loga,"%lg ",jac[i*dim+j]);
  fprintf(loga,"\n");
  }
}
#endif

void jacobian(int dim, val_type *varvals, time_type t, val_type *jac, double *preal, long prelem, int *pint, long pielem);
void rest(int dim, val_type *varvals, time_type t, val_type *res, double *preal, long prelem, int *pint, long pielem);
void projection(int dim,val_type *varvals, time_type t, double *preal, long prelem, int *pint, long pielem);

//This functions evaluates the Jacobian and solves the linear function: J b1 = (0,0,0...0,1)^T
// the last equation of F(x) shoud be the equation of the hiperplane:
//             \Sum_i b0_i (x_i - X_i^0) = t
// So, the last row of the jacobian should have at each element the partial derivative of the hiperplane with 
// respect to x_i, that is b0_i.
// But the last function of F(x) is another function so we have to change these values after having computed the jacobian.
// After the change we will solve the ecuation J b1  =  (0,0,0...0,1)^T and its solution, b1, is such that makes b0 \dot b1 = 1
// the solution b1 is orthonormal to the hiperplane defined by \Sum_i b0_i (x_i - X_i^0) = t
// after having solved the linear system the solution is placed in b1, 
// 
void hiperplaneortho(int dim, double *varvals, double  *ba, double *b1, double t, double *xa, double *taptr, double *preal, long prelem, int *pint, long pielem,option_struct *optr)
{
int i;
val_type * jac0;
#ifdef SOLVER_GSL
    /* GSL */
    gsl_matrix_view m; 
    gsl_vector_view b;
    gsl_vector *x;
    gsl_permutation * p;
    int s,info;

    p = gsl_permutation_alloc (dim);
#else
    char trans;
    lapack_type nsys, info, lda, ldb, dimlapack;
    lapack_type *ipiv;
    /* LAPACK */
    dimlapack = dim;
    lda = dim;
    ldb = dim;
    trans = 'N';
    nsys = 1;   // number of systems to be solved each time by zgetrs_
    if ((ipiv = (lapack_type *)malloc (dim*sizeof(lapack_type)))==NULL)  // Vector dim n
      {
#ifdef LOGFILE
      fprintf(loga,"malloc problem! hiperplaneortho ipiv\n");
#endif
      return;
      }
#endif

    if ((jac0 = (val_type *)malloc (dim*dim*sizeof(val_type)))==NULL) // n x n Matrix
      {
#ifdef LOGFILE
      fprintf(loga,"malloc problem! hiperplaneortho matrix\n");
#endif
      #ifdef SOLVER_GSL
      gsl_permutation_free (p);
      #endif
      free(ipiv);
      return;
      }

    // b1 = (0,0,0,...0,0,1)^t
    for (i=0;i<dim;i++) b1[i]=0;
    b1[dim-1]=1.0;
#ifdef SOLVER_GSL
    b = gsl_vector_view_array (b1, dim);
    m  = gsl_matrix_view_array (jac0, dim, dim);
#endif

/**** Now we have to get the jacobian but with the derivatives of the hiperplane in their places (tha matrix comes transposed!!!)
    the first dim elements of jac0 are the derivatives of f_i respect x_1, but we have to change the last one: the (dim-1)th element.
    the second group of dim values are derivatives of f_i respect x_2. We have to change the last (dim + dim-1)th  and put b0_2
    ...
    the jth group are derivatives respect x_j -> change the last of them which is the ((j-1)*dim + dim-1)th element of jac0 by b0_j
 */
  // i = Newjacobian(dim, varvals, t, jac0, preal, prelem, pint, pielem);
  i = get_advancingjacobian( dim, varvals, xa, ba, t, *taptr, jac0, preal, prelem, pint, pielem,optr);

  if ( i==0 ) // out of memory in New Jacobian computation
    {
    #ifdef SOLVER_GSL
    gsl_permutation_free (p);
    #endif
    free(ipiv);
    free(jac0);
    return;
    }
  #ifdef SOLVER_GSL
  irauli(dim,jac0);
  #endif

  #ifdef LOGFILE
  fprintf(loga,"LU factorization in hiperplaneortho...\n");
  #endif
  #ifdef LAPACKDOUBLE
    #ifdef SOLVER_GSL
    //gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int * signum)
    info = gsl_linalg_LU_decomp (&m.matrix, p, &s);
    #else
    dgetrf_( &dimlapack,&dimlapack, jac0, &dimlapack, ipiv,&info);
    #endif
    if (info !=0)
      {
#ifdef LOGFILE
      fprintf(loga,"malloc problem! hiperplaneortho LU factorization\n");
#endif
      #ifdef SOLVER_GSL
      gsl_permutation_free (p);
      #endif
      free(jac0);
      free(ipiv);
      return;
      }
  #endif
#ifdef SOLVER_GSL
  x = gsl_vector_alloc (dim);
#endif
#ifdef LAPACKDOUBLE
    #ifdef SOLVER_GSL
    gsl_linalg_LU_solve (&m.matrix,  p,  &b.vector, x);   // solve A. x = b
    vectornormalizecp(dim,gsl_vector_ptr(x,0),b1);       // copy x values to b1
    #else
    dgetrs_( &trans, &dimlapack, &nsys, jac0, &lda, ipiv, b1, &ldb, &info ); // Solve: Jac . delta = b1
    vectornormalize(dim,b1);
    #endif
#endif
  *taptr = t;
  vectorcp(dim,varvals,xa);
  #ifdef SOLVER_GSL
  gsl_permutation_free (p);
  gsl_vector_free (x);
  #endif
  free(jac0);
  free(ipiv);
}


int get_advancingjacobian(int dim, double *bal0, double *xa, double *ba, double t, double ta, double *jac0, double *preal, long numpreal, int *pint, long numpint,option_struct *optr)
{
int retval;
int i;

retval = Newjacobian(dim, bal0, t, jac0, preal, numpreal, pint, numpint,optr); // remember jacobian comes transposed
if (retval == 1)
    {
    for (i=0;i<dim; i++) jac0[i*dim + dim-1] = ba[i];  // we change the dderivative df_n(x)/d(x_i) by derivatives of hiperplane h(x): dh(x)/dx_i = ba_i 
    }
return(retval);
}

int Newjacobian(int dim, val_type *varvals, time_type t, val_type *jac, double *preal, long prelem, int *pint, long pielem,option_struct *optr)
{
int i,j;
double h;  // it can be val_type, but for analitic functions the derivative is independent of the directión of the change
val_type *fxminus, *fxplus,ori;
int extra_eqs;

if (!(optr->approxjacobian))
    {
    jacobian(dim,varvals,t,jac, preal, prelem, pint, pielem);
    #ifdef LOGFILE
    fprintf(loga,"new jacobian!\n");
    #ifdef LONGLOG
    fprintf(loga,"variables: \n");
    for (i=0;i<dim;i++) 
    #ifdef LAPACKDOUBLE
    	fprintf(loga," %lg,",varvals[i]);
    #elif LAPACKCOMPLEX
    	fprintf(loga," %lg %lg I,",creal(varvals[i]),cimag(varvals[i]));
    #endif
    fprintf(loga,"\nand jacobian\n");
    printmat(dim,jac);
    #endif
    #endif
    }
  else
  {
  h = 0.00001;
  extra_eqs = 0; //extra_len();
  fxminus = (val_type *)malloc((dim+extra_eqs) *sizeof(val_type));
  fxplus = (val_type *)malloc((dim+extra_eqs) *sizeof(val_type));
  if ((fxplus == NULL)||(fxminus ==NULL))
    {
    #ifdef LOGFILE
    fprintf(loga,"NewJacobian: out of memory!\n");
    #endif
    return(0);
    }
  for (i=0; i<dim; i++)
    {
    ori = varvals[i];
    varvals[i]-=h;
    rest(dim,varvals,t,fxminus, preal, prelem, pint, pielem);
    varvals[i] = ori+h;
    rest(dim,varvals,t,fxplus, preal, prelem, pint, pielem);
    varvals[i] = ori;
    for (j=0; j<dim; j++) 
        // be carefull!! we need traspossed matrix...
        jac [i*dim + j] = (fxplus[j] -fxminus[j])/(2*h);
    }
  #ifdef LOGFILE
  fprintf(loga,"new approximated jacobian!\n");
  #ifdef LONGLOG
  fprintf(loga,"variables: \n");
  for (i=0;i<dim;i++) 
    #ifdef LAPACKDOUBLE
    	fprintf(loga," %lg,",varvals[i]);
    #elif LAPACKCOMPLEX
    	fprintf(loga," %lg %lg I,",creal(varvals[i]),cimag(varvals[i]));
    #endif 
  fprintf(loga,"\nand approximated jacobian\n");
  printmat(dim,jac);
  #endif
  #endif
  free(fxminus);
  free(fxplus);
  } 
return(1);
}

val_type Sqrt(val_type v)
{
return(sqrt(v));
}

val_type Cot(val_type v)
{
val_type res;

res = tan(v);
return(1/res);
}

val_type Tan(val_type v)
{
val_type res;

res = tan(v);
return(res);
}

val_type Power(val_type v, double exp)
{
val_type res;

res = pow(v,exp);
return(res);
}

double Norm(int dim, val_type *v)
{
int i;
double emaitza;

emaitza = 0;
for (i = 0; i<dim; i++)
    {
    #ifdef LAPACKDOUBLE
    emaitza += v[i] * v[i];
    #elif LAPACKCOMPLEX
    auskalo_zer();
    #endif
    }
#ifdef LAPACKDOUBLE
emaitza  = sqrt(emaitza);
#elif LAPACKCOMPLEX
auskalo_zer_moduloa();
#endif
//dznrm2_(&dim,v,&emaitza);
return(emaitza);

}


void vectordiff(int n, val_type * v1, val_type * v2, val_type * v3) // v3 = v1-v2
{
while( n > 0 ) 
  {
  n--;
  v3[n] = v1[n] - v2[n];
  }
}
void vectorcp(int n, val_type * v1, val_type * v2) // v2 = v1
{
while( n > 0 ) 
  {
  n--;
  v2[n] = v1[n];
  }
}


int  hasnan(int dim, val_type *jac0)
{
int res;
int i;

dim = dim*dim;
res=0;
for (i=0; (i<dim) && !res; i++)
    if (isnan(jac0[i]))
             res = 1; 
return(res);
}

int  vhasnan(int dim, val_type *v)
{
int res;
int i;

res=0;
for (i=0; (i<dim) && !res; i++)
    if (isnan(v[i]))
             res = 1; 
return(res);
}


void irauli(int d, val_type *m)
{
int i,j;
val_type lag;

for (i=1; i<d; i++)
  for (j=0; j<i; j++)
    {
    lag = m[i*d+j];
    m[i*d+j]= m[j*d+i];
    m[j*d+i] = lag;
    }
}




void evaluate_function(int dim, double *var0, double *vara, double *ba, double t0, double ta, double *res, double *fnptr, double *preal, long numpreal, int * pint, long numpint)
{
int i;

rest(dim, var0, t0, res, preal, numpreal, pint, numpint);    // res = F(var0,t0)
*fnptr = res[dim-1];  /* we save the las value returned by the user function */
for (i=0, res[dim-1]=0.; i<dim; i++)
    res[dim-1] += ba[i]*(var0[i]-vara[i]);
res[dim-1] -= t0-ta;
}




int findrootbysnewton(int dim, double *jac0, lapack_type *ipiv, double t, double *bal0, double *bal1, double *delta, double *fnptr, double *ba, double *xa, double ta,
		      double *preal, long numpreal, int *pint, long numpint, double *factorptr, option_struct *optr)
/*==========================================================================================================
 *
 *      Mikel 24-10-2012 (emaitzak zehaztu ditut)       
 *
 *      Return values:
 *          = 1. Root has been found
 *          = 2. nanfound problem
 *          = 3  could not factorize jacobian or could not solve linear system
 *          = 4  out of memory
 *          = 0. maximal number of iterations reached.
 *               or (normdelta < 0.8*normdelta_prev)
 *
 *         AND
 *          bal1 points to new root
 *          delta points to f(bal1)
 *
 *==========================================================================================================*/
{                
int i,j;
int nanfound;
int deltacondfails;
val_type *newbal;
//val_type *resvec;
 double normdelta,normdelta0,normdelta_prev;
double auxnorm;
 int deltacond,acccond;
int iters;
#ifdef SOLVER_GSL
    /* GSL */
    gsl_matrix_view m; 
    gsl_vector_view b;
    gsl_vector *x;
    gsl_permutation * p;
    int s,info;

    p = gsl_permutation_alloc (dim);
    m  = gsl_matrix_view_array (jac0, dim, dim);
#else
    char trans;
    lapack_type nsys, info, lda, ldb, dimlapack;
    /* LAPACK */
    dimlapack = dim;
    lda = dim;
    ldb = dim;
    trans = 'N';
    nsys = 1;   // number of systems to be solved each time by zgetrs_
#endif
if ((newbal = (val_type *)malloc (dim*sizeof(val_type))) == NULL)   // Vector dim n
        {
	perror("error allocating memory for v2 of newton iteration");
        #ifdef LOGFILE
        fprintf(loga,"COMPUTING %d JACOBIAN.............\n",optr->numjac);
        #endif
        return(4);
        }
info = 0;
vectorcp(dim,bal0,newbal);					    // newbal = bal0
if (ba != (double *)0)  // we are moving in the curve, so we take into account dim-1 equations and the hiperplane
    evaluate_function(dim, bal0, xa, ba, t, ta, delta, fnptr, preal, numpreal, pint, numpint); // delta = F(bal0)
  else   // we take into account the dim equations
    rest(dim,bal0,t,delta, preal, numpreal, pint, numpint);

#ifdef LOGFILE
    fprintf(loga,"starting findroot:\n");
#ifdef LONGLOG
    fprintf(loga,"initial x:\n");
    for (i = 0; i< dim; i++) 
	{
		fprintf(loga," %lg,",bal0[i]);
	}
    fprintf(loga,"\n");
    fprintf(loga,"and initial f(x):\n");
    for (i = 0; i< dim; i++) 
	{
		fprintf(loga," %lg,",delta[i]);
	}
    fprintf(loga,"\n");
    fprintf(loga,"and initial real parameter list:\n");
    for (i = 0; i< prelem; i++) 
	{
		fprintf(loga," %lg,",preal[i]);
	}
    fprintf(loga,"\n");
#endif
#endif
#ifdef SOLVER_GSL
b = gsl_vector_view_array (delta, dim);
#endif
if (optr->maxnumstepswithsamejac == optr->numstepswithsamejac)
  {
  optr->numjac ++;
  optr->numstepswithsamejac = 0;
  #ifdef LOGFILE
  fprintf(loga,"COMPUTING %d JACOBIAN.............\n",optr->numjac);
  #endif
  if (ba != (double *)0) // the jacobian has as last row derivatives of the equation of the hiperplane (and we change the values by hand)
      i = get_advancingjacobian(dim, bal0, xa, ba, t, ta, jac0, preal, numpreal, pint, numpint,optr);
    else  // the jacobian uses the dim equations (the last is the equatiation we are trying to fulfill)
      i = Newjacobian(dim, bal0, t, jac0, preal, numpreal, pint, numpint,optr);
  if ( i==0 ) // out of memory in New Jacobian computation
    {
    free(newbal);
    return(4);
    }
#ifdef SOLVER_GSL
  irauli(dim,jac0);
#endif
  #ifdef LOGFILE
  fprintf(loga,"Jacobien computed\n");
  #endif
  // The LU decomposition of the jacobian is done only at the beginning
  #ifdef LOGFILE
  fprintf(loga,"LU factorization...\n");
  #endif
  #ifdef LAPACKDOUBLE
    #ifdef SOLVER_GSL
    //gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int * signum)
    info = gsl_linalg_LU_decomp (&m.matrix, p, &s);
    #else
    dgetrf_( &dimlapack,&dimlapack, jac0, &dimlapack, ipiv,&info);
    #endif
  #elif LAPACKCOMPLEX
  zgetrf_( &dimlapack,&dimlapack, jac0, &dimlapack, ipiv,&info);
  #endif
  #ifdef LOGFILE
  fprintf(loga,"     done! info = %d (should be 0)",info);
  #endif
  }
  else
  {
  #ifdef LOGFILE
  fprintf(loga,"we are going to reuse the jacobian.\n");
  #endif
  optr->numstepswithsamejac ++;
  }
if (info !=0)
    {
    free(newbal);
    return(3);
    }
#ifdef SOLVER_GSL
x = gsl_vector_alloc (dim);
#endif
iters = 0;
nanfound = 0;
deltacond =1;
deltacondfails=0;
normdelta = 0;
acccond = 1;
while ( (nanfound == 0) &&
        (info == 0) && 
        (deltacond == 1) && //(normdelta < 0.8*normdelta_prev) &&
        (acccond == 1) && //(normdelta/normdelta0 > pow(10.0,-acc)) &&
        (iters < optr->maxnewtoniter))
    {
#ifdef LOGFILE
    fprintf(loga,"in the iterations... %d iteration. Matrix:\n",iters);
    printmat(dim,jac0);
    fprintf(loga,"and vector b:\n");
    for (i = 0; i< dim; i++) 
	{
		fprintf(loga," %lg,",delta[i]);
	}
    fprintf(loga,"\n");
#endif
    iters ++;
#ifdef LAPACKDOUBLE
    #ifdef SOLVER_GSL
    gsl_linalg_LU_solve (&m.matrix,  p,  &b.vector, x);   // solve A. x = b
    vectorcp(dim,gsl_vector_ptr(x,0),delta);       // copy x to delta: delta = x (element by element)
    #else
    dgetrs_( &trans, &dimlapack, &nsys, jac0, &lda, ipiv, delta, &ldb, &info ); // Solve: Jac . delta = F(bal0)
    #endif
#elif LAPACKCOMPLEX
    zgetrs_( &trans, &dimlapack, &nsys, jac0, &lda, ipiv, delta, &ldb, &info ); // Solve: Jac . delta = F(bal0)
#endif
    nanfound =vhasnan(dim,delta); 
    if (nanfound) 
        {
        #ifdef LOGFILE
	fprintf(loga,"NAN found!! %d iteration\n",iters);
	#ifdef LONGLOG
        for (i = 0; i< dim; i++) 
#ifdef LAPACKDOUBLE
		fprintf(loga," %lg,",delta[i]);
#elif LAPACKCOMPLEX
		fprintf(loga," %lg %lg I,",creal(delta[i]),cimag(delta[i]));
#endif
        fprintf(loga,"\n");
        #endif
        #endif
	}
    normdelta_prev = normdelta;
    normdelta = Norm(dim,delta);   //normdelta  = ||delta||
#ifdef LOGFILE
    fprintf(loga,"   ||delta|| = %lg \n",normdelta);
    #ifdef LONGLOG
    fprintf(loga,"   delta: ");
    for (i=0;i<dim;i++) 
#ifdef LAPACKDOUBLE
	fprintf(loga," %lg,",delta[i]);
#elif LAPACKCOMPLEX
	fprintf(loga," %lg %lg I,",creal(delta[i]),cimag(delta[i]));
#endif
    fprintf(loga,"\n");
    #endif
#endif 
    if (iters == 1) normdelta0 = normdelta;
    if (iters > 1) 
      {
	if (normdelta > (optr->iterationdecreasefactor * normdelta_prev))  // after second iteration we control 
                                                // wether the change decrements at each iteration. If do not, we will stop
	    {
	        //if (normdelta > normdelta_prev) 
			deltacondfails ++;
                #ifdef LOGFILE
		fprintf(loga,"    deltacondfails = %d (%lg > %lg)",deltacondfails,normdelta,0.8 * normdelta_prev);
                #endif 
		if (deltacondfails > optr->maxdeltacondfails) deltacond = 0;
	    }
	if (normdelta/normdelta0 < pow(10.0,-optr->acc)) {acccond = 0; deltacond=1;} // if we reach accuracy deltacond does not mind
        #ifdef LOGFILE
	fprintf(loga,"    accuracy cond = is %lg < %lg \n",normdelta/normdelta0, pow(10.0,-optr->acc));
        #endif 
      }
    vectordiff(dim,newbal,delta, &(bal1[0]));                     // bal1 = newbal - delta
    vectorcp(dim,bal1,newbal);                                    // newbal = bal1;  newbal is the new approximation of the root
    if (ba != (double *)0)
        evaluate_function(dim, newbal, xa, ba, t, ta, &(delta[0]), fnptr, preal, numpreal, pint, numpint);   // delta = F(newbal)   
      else
	rest(dim,newbal,t,&(delta[0]), preal, numpreal, pint, numpint);
    } //END OF NEWTON ITERATIONS


#ifdef LOGFILE
if (nanfound != 0) fprintf(loga,"end of the newton iterations because NaN on solving equationd by lapack\n");
if (info != 0) fprintf(loga,"end of the newton iterations by info != 0\n");
if (deltacond == 0) fprintf (loga,"end newton iterations because difference between approximations do not decrease after %d iterations\n",iters);
if (iters >= optr->maxnewtoniter) fprintf (loga,"end newton iterations by cause 4 (too many iterations...)\n");
#endif
 if ((info == 0) && (nanfound ==0)&& ((acccond == 0) || (normdelta < pow(10.0,-optr->prec))) )
    {
    info = 1;
#ifdef LOGFILE
    fprintf(loga,"new root!\n");
    #ifdef LONLOG
      fprintf(loga,"the new root is:\n");
      for (i = 0;i<dim; i++)
	#ifdef LAPACKDOUBLE 
	fprintf(loga,"%lg,  ",bal1[i]); 
	#endif
      fprintf(loga,"\n"); 
    #endif
#endif
    #ifdef SOLVER_GSL
    /* GSL */
    gsl_permutation_free (p);
    gsl_vector_free (x);
    #endif
    free(newbal);
    /*******************************************
    We want to adecuate the stepsize of the continuation path in the following way:
    The idea is to get the root in each Newton iteraration process using always "iter_opt"
    number of iterations. We know that if the step is long we will need more iterations
    and with short steplength less iterations are needed. 
    with a given length, i.e. l, we have needed "i" iterations to get an accuracy of 10^acc
    so, we know that (K*l)^i = 10^acc ---> K = (10^(acc/i))/l

    but we want this other situation: (K * l_opt)^iter_opt = 10^acc
    sustituting K and reordering we obtain 
             l_opt = 10^(acc(1/iter_opt - 1/i)) * l
    so the factor for l to obtain l_opt is  10^(acc(1/iter_opt - 1/i))
    (take into account that we have the acc value as a positive value)
    ******************************************** */
    *factorptr = pow(10,optr->acc*(1.0/(double)iters - 1.0/(double)optr->iter_opt));
    return (1);
    }
  else 
    {  
    *factorptr = 0;
#ifdef LOGFILE
    fprintf(loga,"it has been considered there is no convergence! \n");
#endif
    }
//listfree(&listofvals);

#ifdef SOLVER_GSL
/* GSL */
gsl_permutation_free (p);
gsl_vector_free (x);
#endif
free(newbal);
if (nanfound)
    {
#ifdef LOGFILE
    fprintf(loga,"returning nanfound"); fflush(stdout);
#endif
    return(2);
    }
  else 
    {
    if (info) return(3);
      else return(0);
    }
}


int steptotaux(int dim, val_type *jac0, lapack_type *ipiv, time_type t0, time_type dt, 
                time_type t1, time_type *tnewptr, time_type *dtnewptr, 
                val_type *bal0, val_type *bal1, val_type *res, double *fnptr, double *ba, double *xa, double ta,
		double *preal, long numpreal, int *pint, long numpint,
                option_struct *optr)
{
/********************************************************
 *    Return values:
 *          = 1. Root has been found
 *          = 2. nanfound problem
 *          = 3  could not factorize jacobian
 *          = 4  out of memory
 *          = 0. maximal number of iterations reached.
 *               or (normdelta < 0.8*normdelta_prev)
 *
 *         AND
 *          bal1 points to new root
 *          res points to f(bal1,t)
 *          tnewptr points to the t for wich f(bal1,t) is considered to be 0
 *          dtnewptr points to the deltat proposed for the next step
 *
 *
 *********************************************************/
int i,j;
time_type dtused;
int rootfound, attempts;
double factor;
double dtreducefactor;
double dtfactor2reusejac;
double maxstepincr;

maxstepincr = 4.0;
dtreducefactor = 0.3;
if (!(optr->approxjacobian)) dtfactor2reusejac = 0.9;
  else dtfactor2reusejac = 0.6;
rootfound = -1;
if (dt >0) 
    {
    if ((t0 + dt) > t1) dt = t1-t0;
    }
  else
    {
    if ((t0 + dt) < t1) dt = t1-t0;
    }
for (attempts = 0, dtused = dt;
     (rootfound != 1) && (attempts < optr->maxstepattempts);
     attempts ++)
    {
    *tnewptr = t0 + dt;
#ifdef LOGFILE
    fprintf(loga,"    to findroot attempt number %d \n",attempts);
#endif
    rootfound = findrootbysnewton(dim, jac0,ipiv, t0+dt, bal0, bal1, res, fnptr, ba, xa, ta, preal, numpreal, pint, numpint, &factor, optr);
    if (rootfound == 0)  //lack of convergence  -> ((rootfound != 1)&&(rootfound != 2)&&(rootfound != 3))
        {
	if (optr->numstepswithsamejac == 0 )  // we reduce stepsize by dtreducefactor only if we are in the case with good jacobian and not convergence
		dt = dt * dtreducefactor;
	    else   // in other case  we reduce by 2*reducefactor
		{
		dt = dt * ((dtreducefactor > 0.5)?dtreducefactor:dtreducefactor*2.0);
		}
	optr->numstepswithsamejac= optr->maxnumstepswithsamejac;  // to force the evaluation of jacobian and LU factorization 
#ifdef LOGFILE
        fprintf(loga,"Newton iteration convergence failed! returned value = %d (we are going to try again with a shorter step: %lg : from %lg to %lg)\n",rootfound,dt,t0,t0+dt);
#endif
        }
      else  // rootfound== 1 or rootfound == 2...
        {
        if (rootfound != 1) // to force the evaluation of jacobian and LU factorization
            {
#ifdef LOGFILE
            fprintf(loga,"    ROOT NOT FOUND: NANFOUND or PROBLEM in FACTORIZATION... findroot returned %d value\n",rootfound);
#endif
            optr->numstepswithsamejac= optr->maxnumstepswithsamejac;  
            }
        }
    }
#ifdef LOGFILE
fprintf(loga,"findroot return value = %d (if 1 root found)\n", rootfound);
#endif
if (rootfound == 1)
    {
    if (factor > maxstepincr) factor = maxstepincr;
    *dtnewptr = dt *factor;
    if (optr->numstepswithsamejac < optr->maxnumstepswithsamejac) 
	// if in the next newton process we are going to reuse the jacobian we will use a shorter step
	*dtnewptr = (*dtnewptr) * dtfactor2reusejac;
    }
  else // we are not able to get the root at t0+dt so the last t for which we have a root is t0
    {
    *tnewptr = t0;
    }
return (rootfound);
}


int problem_conf_values(option_struct *optr)
/*int *maxNWTiter, int *maxINNERstep, int *precptr, int * accptr, 
int *iter_optptr, int *maxnumstepswithsamejacptr; int *maxdeltacondfailsptr, 
int approxjacobianptr, double *stepsizeptr, double *min_av_stepsizeptr, 
int *num_steps_for_stepsize_controlptr */
{
char word[81];
FILE *fitx;
struct stat stata;
int i,res;
char filename[30]; /* length of "./arcContinuation.conf" is 20 + end of string = 21 */

strcpy(filename,"./arcContinuation.conf");
optr->approxjacobian = 0;

if (stat(filename,&stata) == 0) 
    {
    fitx = fopen(filename ,"r");
    for (i = fscanf(fitx,"%80s",word);
     	 i != EOF;
     	 i = fscanf(fitx,"%80s",word))
         /* no more than 80 chars, it stops in a space or EOL */
    	{
    	if (i == 1)  /* it has read the word */
      	   {
           if ((word[0] == '\0') || (word[0] == '#'))
           		res=fscanf(fitx,"%*[^\n]\n"); /* comment line */
            else 
               {
               if (strcmp(word,"Newton")==0) 	//Newton iteration #
                 {
                 res=fscanf(fitx,"%d%*[^\n]\n",&(optr->maxnewtoniter));
                 }
                else
              	 {
	       	 if (strcmp(word,"Inner")==0)
		       	res=fscanf(fitx,"%d%*[^\n]\n",&(optr->maxinnersteps)); 
		  else
		   {
		   if (strcmp(word,"prec")==0)
		      res=fscanf(fitx,"%d%*[^\n]\n",&(optr->prec)); 
	            else
	              if (strcmp(word,"acc")==0)
		         res=fscanf(fitx,"%d%*[^\n]\n",&(optr->acc));
		      else
                        if (strcmp(word,"iter_opt")==0) 
			  res=fscanf(fitx,"%d%*[^\n]\n",&(optr->iter_opt));
		        else
                          if (strcmp(word,"maxnumstepswithsamejac")==0) 
			    res=fscanf(fitx,"%d%*[^\n]\n",&(optr->maxnumstepswithsamejac));
		           else
                            if (strcmp(word,"maxdeltacondfails")==0) 
			      res=fscanf(fitx,"%d%*[^\n]\n",&(optr->maxdeltacondfails));
			     else
                              if (strcmp(word,"approxjacobian")==0) 
				  optr->approxjacobian = 1;
				else
                                  if (strcmp(word,"stepsize")==0) 
					res=fscanf(fitx,"%lf%*[^\n]\n",&(optr->stepsize));
				    else
				      if (strcmp(word,"min_av_stepsize")==0)
					  res=fscanf(fitx,"%lg%*[^\n]\n",&(optr->min_av_stepsize));
				        else
					  if (strcmp(word,"num_steps_for_stepsize_control")==0)
					      res=fscanf(fitx,"%d%*[^\n]\n",&(optr->num_steps_for_stepsize_control));
					    else
					      res=fscanf(fitx,"%*[^\n]\n"); /* nothing of interest in the line, so read until EOL */
	            } // else Inner
           	 } // else Newton
               } //else  not comment
       	   } //else if 1==1
       	} /* for*/
    fclose(fitx);
    }
  else /* configuration file no present! */	
	{ 
	}
return(0);
}


void compute_average_stepsize(time_type *av_stepsizeptr, int nstep, time_type taux, time_type tnew, option_struct *optr)
{
time_type aux; // aux is the length of the last step
int num;       // num is the number of steps taken into account in the computation of the average

aux = tnew -taux;
if (aux < 0) aux = -aux;
if (nstep == 1) 
    *av_stepsizeptr = aux;
  else
    {
    if (nstep > optr->num_steps_for_stepsize_control) 
        {
        num = optr->num_steps_for_stepsize_control;
        }
      else  // we have to advance more to get the average of the first steps
        {   // this can be the case: we want the average of previous 10 steps but
            // we have advanced only 5 steps (so nstep is 5 and in the average_stepsize
	    // we have the average of the 4 previous steps
        num = nstep;
        }
    *av_stepsizeptr = ((*av_stepsizeptr * (num-1)) + aux) / ((double) num);
    }
}

// the minimum of the parabola formed by three points (the central point is near the minimum)
double minbyinterpolation(double a, double fa, double b, double fb, double c, double fc)
{
// x -> (c^2 (fa - fb) + a^2 (fb - fc) + b^2 (-fa + fc))/(2 (c (fa - fb) + a (fb - fc) + b (-fa + fc)))
return((c*c*(fa - fb) + a*a* (fb - fc) + b*b* (-fa + fc))/(2.0 * (c * (fa - fb) + a * (fb - fc) + b * (-fa + fc))));
}





/**************************************************************************************/
/*****************               MAIN FUNCTION                  ***********************/
/**************************************************************************************/

void arcContinuation(/* input */int dim, double *var0, double *b0, double init, double end, int action, double *preal, long numpreal, int *pint, long numpint,
		     /* output */ double *var1, double *b1, double *taptr, double *tbptr, int * infoptr)
{
int k;
int i,j,nstep;
int info;
time_type dt,dtnew,taux,tnew;
double t0,t1;
double *bal1;
double *res0;
double fn;   /* the las value returned by the user function rest is saved in this variable */
double * vara;
double ta;
double *var00;
double *b00;
double t00;
int is_closed;
double plane_eq;
double plane_eq_prev;
double b_diff;
double var_diff;
val_type *jac0;
lapack_type *ipiv;
int rootfound;
double norma3;
double auxfactor;
double average_stepsize;
option_struct opt_vals;
/*****************
int maxnewtoniter;
int maxinnersteps;
int maxstepattempts; // this value depends on maxinnersteps, but at first steps it is a higher value 
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
double stepsize;  // The user can set the initial stepsize. Starting from t0 the process goes to t1
		  // and the initial stepsize can be set in the .conf file
//
# the mecanism to control the minimun stepsize:
#    If the average_stepsize goes beyond min_av_stepsize the process will stop.
#    average_stepsize = \Sum_i^s stepsize_{-i}/s
#    where s is the number of steps taken into account for the estimation. (if s == 1 then 
#    the stepsize can not be beyond min_av_stepsize) s can be set in the configuration file:
#          num_steps_for_stepsize_control: is the number of steps taken into account to get the
#                average_stepsize
#          min_av_stepsize: The minimun average_stepsize acepted. Format #.#e#
*
double min_av_stepsize;
int num_steps_for_stepsize_control;
double stepsize;  // The user can set the initial stepsize. Starting from t0 the process goes to t1
		  // and the initial stepsize can be set in the .conf file
double iterationdecreasefactor;  // the minimum decrease of delta at each Newton iteration
int max_steps; 
*/
/************* variables used when action is 1 (look for the minimum)   */
int lastwasmin; 
double tmin,tminprev,tminpost,tmininterpolated;
double extraCond, newextraCond, lastextraCond, extraCondprev, extraCondpost;
double *varaux;
double *bmin;
/************* variables used when action is 2 (look for the absolute root)   */
int surechange; 
double prev_fn;
double troot;

#ifdef LOGFILE
loga = fopen("./arcContinuation.log","w");
fprintf(loga,"starting arcContinuation...\n");
#endif
t0 = init;
t1 = end;
opt_vals.maxnewtoniter = 20;   //???????????????
opt_vals.maxinnersteps = 5;
opt_vals.prec =8;
opt_vals.acc = 3;
opt_vals.iter_opt= 8;
opt_vals.maxdeltacondfails=2;
opt_vals.maxnumstepswithsamejac=0; // 0 -> the jacobian is computed for every newton method
opt_vals.stepsize=0;
opt_vals.min_av_stepsize = 1.0e-10;
opt_vals.num_steps_for_stepsize_control = 1;
opt_vals.numjac = 0;
opt_vals.iterationdecreasefactor = 0.8;
opt_vals.max_steps =400;
opt_vals.closure_tolerance = 1.0e-5;        //???????????????
problem_conf_values(&opt_vals);
if ((res0 = (val_type *)malloc (dim*sizeof(val_type)))==NULL)  // Vector of dimension dim 
    {
    t1 = t0;
#ifdef LOGFILE
    fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
    }
if ((ipiv = (lapack_type *)malloc (dim*sizeof(lapack_type)))==NULL)  // Vector of dimension dim 
    {
    t1 = t0;
#ifdef LOGFILE
    fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
    }
if ( ((t1>t0)&&(opt_vals.stepsize > 0)) || ((t1<t0)&&(opt_vals.stepsize < 0)) )
    dt = opt_vals.stepsize;
  else
    {
    dt = (t1-t0)/opt_vals.max_steps;
    }
if ((jac0 = (val_type *)malloc (dim*dim*sizeof(val_type)))==NULL) // n x n Matrix
    {
    t1 = t0;
#ifdef LOGFILE
    fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
    }
if ((vara = (val_type *)malloc ((dim)*sizeof(double)))==NULL)  // Vector of dimension dim 
    {
    t1 = t0;
#ifdef LOGFILE
    fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
    }

if ((var00 = (val_type *)malloc ((dim)*sizeof(double)))==NULL)  // Vector of dimension dim 
    {
    t1 = t0;
#ifdef LOGFILE
    fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
    }
if ((b00 = (val_type *)malloc ((dim)*sizeof(double)))==NULL)  // Vector of dimension dim 
    {
    t1 = t0;
#ifdef LOGFILE
    fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
    }
if (action == 1)  // in tha case of looking for minimum we need var min and b min to save both vectors.
    {
    if ((bmin = (val_type *)malloc ((dim)*sizeof(double)))==NULL)  // Vector of dimension dim 
        {
        t1 = t0;
#ifdef LOGFILE
        fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
        }
    }
if ((action == 1 || action == 2))
    {
    if ((varaux = (val_type *)malloc ((dim)*sizeof(double)))==NULL)  // Vector of dimension dim 
        {
        t1 = t0;
#ifdef LOGFILE
        fprintf(loga,"malloc problem! tfinal = %lf\n",t0);
#endif
        }
    }
/* first we have to obtain the orthogonal hiperplane for n-1 equations and as last equation b1 \dot (x-x_a) = t-t_a */
/* the return value is the vector b1 and the x_a vector and the t_a*/

hiperplaneortho(dim,var0, b0, &(b1[0]), t0, &(vara[0]), &ta, preal, numpreal, pint, numpint,&opt_vals);
/* Save the hiperplane equation, to control the closed paths*/
vectorcp(dim,b1,b00);
vectorcp(dim,vara,var00);
t00 = ta;
/* We evaluate the user function to get f_1(x), ...f_n-1(x) and f_n(x) 
 * But we have to move in the curve taking into account the first n-1 functions and the hiperplane
 * So we have to put the hiperplane equation as last element of the vector returned by rest. And we 
 * will save the last value returned by the user function. All this is made by the function called
 * evaluate_function(dim, var0, b0, b1, t0, &(res0[0]), &fn, preal, prelem, pint, pielem);
*/

 /* with this initialization the hiperplane equation value is 0 */
evaluate_function(dim, var0, vara, b1, t0, ta, &(res0[0]),&fn, preal, numpreal, pint, numpint);
/* rest(n, bal0, t0, &(res0[0]), preal, prelem, pint, pielem);    // res0 = F(bal0,t0)
/*******************************
 * MAIN FOR of the process
 *
 *   each iteration advances from taux to tnew. 
 *   if all goes OK in the iteration, at the end of the iteration we will have 
 *      a.- the new root at bal1,
 *      b.- the residual res0 = f(bal1,tnew)
 *      c.- the new stepsize proposed for the next step: dtnew
 *******************************/

switch (action)
    {
    case 0:  // do nothing!
	break;
    case 1:
	newextraCond= fn;    // extraCondNorm = || extra_conditions ||
	tmin = t0;
	extraCond = newextraCond;
	tminprev = tmin;
	extraCondprev = extraCond;
	lastwasmin = 1; // after the first step we have to save the time and the extraCondNorm for interpolation (3 points are needed)
	break;
    case 2:
	break;
    }

for (taux = t0, nstep = 0, opt_vals.numstepswithsamejac=opt_vals.maxnumstepswithsamejac,average_stepsize= ((dt>0)?dt:-dt),is_closed = 0,surechange=0;
     ((dt>0)?(taux-t1 < 0):(t1-taux < 0)) && (nstep < 10*opt_vals.max_steps) && (average_stepsize > opt_vals.min_av_stepsize)&&!is_closed && !surechange;
     nstep ++,dt=dtnew,compute_average_stepsize(&average_stepsize,nstep,taux,tnew, &opt_vals),taux = tnew
    )
    {
#ifdef LOGFILE
    fprintf(loga,"I' ll try to walk:  step= %d,  from t = %lg to t= %lg with dt = %lg \n",nstep, taux,taux+dt,dt);
#endif
    if (nstep <5) opt_vals.maxstepattempts = opt_vals.maxinnersteps*(10 - 2*nstep);
        else opt_vals.maxstepattempts = opt_vals.maxinnersteps;
    rootfound = steptotaux(dim, jac0, ipiv, taux, dt, t1, &tnew, &dtnew, var0, &(var1[0]),&(res0[0]), &fn, b1, vara, ta, preal,numpreal,pint,numpint, &opt_vals);
    if ((rootfound < 1) || (rootfound >=2))
        {
        #ifdef LOGFILE
	switch (rootfound)
	    {
	    case 0:
		fprintf(loga,"You have a problem: try to reduce the accuracy or increase working precision...\n");
	        break;
	    case 2:
                fprintf(loga,"Lapck reurned solution with NaN when t == %lg\n",taux);
	        break;
	    case 3:
                fprintf(loga,"The linear system could not been solved\n");
	        break;
	    case 4:
                fprintf(loga,"Out of memory\n");
	        break;
	    default:
		fprintf(loga,"some problem in the %d argument when solving the lineal system of a s-Newton step...\n", -rootfound);
		break;
	    }
        #endif
	t1 = taux;  // taux is the last t for which the root has been found. This root, bal0, will be the solution given.
        tnew = t1;  // this is a way to go out from the for
        }
      else  // rootfound
        {
#ifdef LOGFILE
        fprintf(loga,"arcContinuation: root found (else way) at %.12lg. Next step will advance %.12lg\n",tnew,dtnew);
#endif

	switch (action) // depending on the action we have to do different things
	    {
	    case 0:  // we want just move on the curve
		break;
	    case 1:  // we are moving on the curve but at the end we want the point for which f_n(x) has minimum value
		lastextraCond = newextraCond;  // we save the extraCond of previous root
		newextraCond= fn;     // newextraCond = f_n(x)
		if (lastwasmin) // we save the minimum, the previous and the next to get the minimum by polynomial interpolation
			{
			tminpost = tnew;
			extraCondpost = newextraCond;
			vectorcp(dim,b0,bmin);  // we have still the vector b of the hiperplane corresponding to the varaux into b0
			lastwasmin=0;
			}
		if (newextraCond < extraCond)     // this point has the minimun. So we have to save it.
			{
			tminprev = taux;    // taux is the time for the previous root
			extraCondprev = lastextraCond;   // we had saved the extraCond of the previous root
			tmin = tnew;
			lastwasmin = 1; 
			vectorcp(dim,var1,&(varaux[0]));
			extraCond = newextraCond;
			#ifdef LOGFILE
  			  fprintf(loga,"\n\n ALDAKETA!!! tmin (new) = %lf with norma = %lf\n",tmin,extraCond);
			#endif
			}
		break;
	    case 2:  // we are moving and we want to stop just when  f_n(x) is 0
		if ((fn*prev_fn)<0) // We have crossed the point that makes f_n(x) = 0
		    {
		    vectorcp(dim,var0,&(varaux[0])); // this has been the point where first time has changed (but with no accurate root)
		    vectorcp(dim,var1,&(var0[0])); 
		    opt_vals.acc=opt_vals.acc*4;
                    opt_vals.numstepswithsamejac= opt_vals.maxnumstepswithsamejac;
                    rootfound = findrootbysnewton(dim, jac0,ipiv, tnew, var0, &(var1[0]),&(res0[0]),&fn, b1, vara, ta, preal,numpreal,pint,numpint, &auxfactor,&opt_vals );
                    opt_vals.acc=opt_vals.acc/4.;
                    if (rootfound != 1)// tryin to get the more accurate root failed. so it is better to get the rest of the less accurate root
			{  
			vectorcp(dim,var0,&(var1[0])); //var1 must be a root
			evaluate_function(dim, var1, vara, b1, tnew, ta, &(res0[0]),&fn, preal, numpreal, pint, numpint);
			}
		      //else  We got a better root and we have the function evaluated at this root
			
		    // we have the root into var1 and the previous root in varaux
		    if ((fn*prev_fn)<0)  // after having got the more accurate root the change of the sign has been confirmed.
			{
			#ifdef LOGFILE
			    fprintf(loga," Sign changed! at t = %lf\n",tnew);
			#endif 
			surechange = 1;
                        /* we want to return the values tprev and tafter where:
			   tprev is the last value of t for which the sign of extra conditions were unchanged and 
			   tafter is the t for which some sign has changed. We want also return the root at tprev*/
			
			t0 = taux;  // tprev
			t1 = tnew;  // tafter
			vectorcp(dim,varaux,&(var0[0]));   // So we have the previous root in var0 and the more accurated root which changes the condition sign in var1
			}
		    else /* we don't have a "surechange" but probably next step there will be a sign change */
			{
			#ifdef LOGFILE
			fprintf(loga," False Sign change: at t = %lf, may be next step...\n",tnew);
			#endif 
			prev_fn = fn; // we save the f_n(x) (which has the same sign as prev_fn)
			}
		    }
		  else 
		    {
		    prev_fn = fn;   // we save the last f_n(x)
		    }
		
		break;
	    }  // end of switch(action)
	if (!surechange)  // there was not any surechange: we have to continue!!! 
			  // control posible closed curve and prepare for next step (recalculate the hiperplane orthogonal direction...)
	    {        
	    vectorcp(dim,var1,&(var0[0]));
	    vectorcp(dim,b1,&(b0[0])); /* we will compute a new vector which will be orthogonal to the curve at the new position. */
	    hiperplaneortho(dim,var0, b0, &(b1[0]), tnew, &(vara[0]), &ta, preal, numpreal, pint, numpint,&opt_vals);
	    // Control that we are not again in the initial point. We have the hiperplane equation saved.
	    for (i=0, plane_eq = 0; i<dim; i++) 
	        plane_eq += b00[i]*(var1[i]-var00[i]); 
	    plane_eq -= t00-tnew;
	    if (nstep > 1) 
	      {
	      if ((plane_eq_prev * plane_eq)< 0) //we have crossed the hiperplane
		{
		for (i=0, var_diff=0,b_diff=0; i<dim;i++) 
		    {
		    var_diff += pow(var1[i]-var00[i],2.0);
		    b_diff += pow(b1[i]-b00[i],2.0);
		    }
		var_diff = sqrt(var_diff);
		b_diff = sqrt(b_diff);
		if ((var_diff < opt_vals.closure_tolerance) && (b_diff < opt_vals.closure_tolerance))
		    is_closed = 1;
		}
	      }
	    plane_eq_prev = plane_eq;
	    /* We call a function which gives the posibility to apply some sort of projection of the solutions after each step.
            This is useful when working with homogenous variables, for performing projections onto the unit sphere */
	    //projection(n,&(bal0[0]),tnew,preal, prelem, pint, pielem);
	    if (dt > 0)
              {
              if ((t1 - tnew) < dtnew) dtnew = t1-tnew;
              }
            else
              {
              if ((t1 - tnew) > dtnew) dtnew = t1-tnew;
              }
	    } //end of surechange
        }  // end of if ( else way: rootfound)
#ifdef LOGFILE
    fprintf(loga," STEP %d .... rootfound = %d ...\n",nstep,rootfound);
#endif
    }  // end of the main FOR. The new (or may be the old) root is in bal0
       // end of the main FOR. The new (or may be the old) root is in bal0
       // end of the main FOR. The new (or may be the old) root is in bal0
       // end of the main FOR. The new (or may be the old) root is in bal0
if ((nstep >= 10*opt_vals.max_steps) || (average_stepsize <= opt_vals.min_av_stepsize))
	// in this case the proccess has stoped because the number of steps reaches the limit or 
	// because the stepsize required is too little.
	// The last t for wich the root has been obtained is taux and not t1, so we have to assign taux to t1
    {
    switch (action) // depending on the action we have to do different things
	    {
	    
	    case 1: 
	    t1 = taux;
		if (taux == tmin) *infoptr = -12;
		    else *infoptr = -11;
		break;
	    case 2:  // we are moving and we want to stop just when  f_n(x) is 0 
		break;
	    }
 
    }
switch (action) // depending on the action we have to do different things
    {
    case 0:  // we want just move on the curve	
	if (taux != init) // if we have been able to advance in t, then we will accurate the new root
	    {
	    if (average_stepsize <= opt_vals.min_av_stepsize) // the stepsize required is too little. we will not acurate the root.
		{
		*infoptr = -1;
		}
	      else  // we have moved and we havent stoped because we found a singularity. So it is posible to accurate the root.
	        {
	        opt_vals.acc=opt_vals.acc*4;
	        opt_vals.numstepswithsamejac= opt_vals.maxnumstepswithsamejac;
	        rootfound = findrootbysnewton(dim, jac0,ipiv, t1, var0, &(var1[0]),&(res0[0]),&fn,b1,var0,t1,preal,numpreal,pint,numpint, &auxfactor,&opt_vals);
	        opt_vals.acc = opt_vals.acc/4.;
	        /* We call a function which gives the posibility to apply some sort of projection of the solutions after each step.
	           This is useful when working with homogenous variables, for performing projections onto the unit sphere */
	       //projection(n,&(bal1[0]),t1,preal, numpreal, pint, numpint);  // <--- Anderrek
	        if (rootfound == 1) vectorcp(dim,var1,&(var0[0]));
	        #ifdef LOGFILE
	          else
	            {
	            fprintf(loga, "Warning: Last root with problems!!!!!\n");
	            }
	        fprintf(loga,"number of jacobian computations = %d\n",opt_vals.numjac);
	        #endif
		if (is_closed) *infoptr = 2; //We have returned to the init root.
		    else  if (nstep >= 10*opt_vals.max_steps) // We have stoped because max_steps was reached, not because t = tend. 
				{
				*infoptr = 3;
				t1 = taux;
				}
			    else *infoptr = 1; 
	        }
	    }
	  else    // we have not been able to move
	    {
	    vectorcp(dim, var0, var1);
	    *infoptr = -2;
	    }
	return_values(dim,var1,b1,t0,t1,preal,numpreal,pint,numpint);
	break;
    case 1: 
	if (lastwasmin) // we have to save the minimum, the previous and the next to get the minimum by polynomial interpolation
	    {	// the last root is the one with minimum norm for extra conditions, so we have not saved the point after it.
	    tminpost = tnew;  // the point after the minimum will be the same as the minimum (we repeat it --> linear interpolation!)
	    extraCondpost = newextraCond;
	    vectorcp(dim,b0,&(bmin[0]));     // because we have not saved the b vector of the last step (the min)
	    lastwasmin=0;
	    }
	// now we have the point with minimum value for extra conditions. and we have the previous point and the point after it with their extra conditions norm.
	// by interpolation we will obtain the point for wich the extra condition is the minimum:

	//     pol = InterpolatingPolynomial[{{a, fa}, {b, fb}, {c, fc}}, x]

	// fa + (-a +  x) ((-fa + fb)/(-a + b) + ((-((-fa + fb)/(-a + b)) + (-fb + fc)/(-b + c)) (-b + x))/(-a + c))

	//     eq = (Simplify[D[pol, x]] == 0)

	// (c^2 (fa - fb) + b^2 (-fa + fc) + a (fb - fc) (a - 2 x) + 2 c (-fa + fb) x + 2 b (fa - fc) x)/((a - b) (a - c) (b - c)) == 0

	//     Solve[eq, x] // Simplify

	// x -> (c^2 (fa - fb) + a^2 (fb - fc) + b^2 (-fa + fc))/(2 (c (fa - fb) + a (fb - fc) + b (-fa + fc)))

	tmininterpolated = minbyinterpolation(tminprev,extraCondprev, tmin, extraCond,tminpost,extraCondpost);
	//with this new tmin and the approximation of the root we recalculate the root at the new tmininterpolated with higher accuracy.
	if (taux != init) // if we have been able to advance in t, then we will accurate the new root
	    {
	    if ((tmininterpolated == tnew)&&(average_stepsize <= opt_vals.min_av_stepsize)) 
		{	// the minimum is reached at the end but the stepsize required is too little. we will not acurate the root.
		*infoptr = -12;
		*taptr = tmininterpolated;
		*tbptr = tnew;
		vectorcp(dim,bmin, &(b1[0]));
		// the return values are in var1
		}
	      else  // we have advanced and 
		    //    a.- if the minimum is the end point then  it is not a singularity
		    // or
		    //    b.- if the end point is a singularity then the minimum is in a previous point. 
		    // So, in both cases we can accurate the root.
	        {
		opt_vals.acc=opt_vals.acc*4;
		opt_vals.numstepswithsamejac = opt_vals.maxnumstepswithsamejac;   // we need a new jacobian!! 
		if ((tmin-tminprev)*(tmininterpolated -tmin) < 0)  // tmininterpolated is in the oposite direction of the direction given by bmin
			for (i=0; i< dim; i++) bmin[i] = -bmin[i];
		rootfound = findrootbysnewton(dim, jac0,ipiv, tmininterpolated, varaux, &(var1[0]),&(res0[0]),&fn,bmin,varaux,tmininterpolated,preal,numpreal,pint,numpint, &auxfactor,&opt_vals);
		if ((tmin-tminprev)*(tmininterpolated -tmin) < 0)  // to recover the direction we have used to advance in the curve
			for (i=0; i< dim; i++) bmin[i] = -bmin[i];	
		opt_vals.acc = opt_vals.acc/4.;
		if (rootfound == 1) 
			{
			// var1 has the vector to be returned
			*taptr = tmininterpolated;
			}
		      else
		        {
			*taptr = tmininterpolated;
			vectorcp(dim, varaux, &(var1[0]));  // we will return the non accurated value
			#ifdef LOGFILE
		        fprintf(loga, "Warning: root for point with minimum g(x) with problems!!!!!\n");
			#endif
		        }
		vectorcp(dim,bmin,&(b1[0]));  // the coefficients of the hiperplane for the point that minimizes the equation.
		*tbptr = taux;
		if (average_stepsize <= opt_vals.min_av_stepsize) *infoptr = -11; //We have stoped in a singularity but the minimum is a good point
		    else  // we have advanced in the curve and we have not stoped in a singularity
			{
			if (is_closed)
				{
				if ((tmininterpolated == t0)||(tmininterpolated == taux)) *infoptr = 20;
					else *infoptr = 21;
				}
			    else
				{
				if (tmininterpolated == t0) *infoptr = 10;
				    else if (tmininterpolated == taux) *infoptr = 12;
					    else *infoptr = 11;
				}
			}
		}
	    }
	  else    // we have not been able to move
	    {
	    *infoptr = -20;
	    *taptr = t0;
	    *tbptr = t0;
	    vectorcp(dim, var0,&(var1[0]));
	    vectorcp(dim, b0,&(b1[0]));
	    }
	break;
    case 2:  // we are moving and we want to stop just when  f_n(x) is 0
	/* if we got a shurechange we have the values tprev and tafter where:
	tprev is the last value of t for which the sign of extra conditions were unchanged and 
	tafter is the t for which some sign has changed. We want also return the root at tprev
	These values are saved into t0 and t1 and we have also the previous root
	in var0 (and varaux) and the more accurated root which changes the condition sign in var1.
	Now I am going to get a better approximation of the point where the sign changes, and after, 
	I will try to get the root where all conditions are satisfied: f_n(x) also! */ 
	
	if (surechange) 
	    {
	    //Lortu interpolazio lineala erabiliz t_root eta acc handiarekin beretzat lortu x_root
	    troot = t0 - prev_fn *(t1-t0)/(fn-prev_fn);
	    opt_vals.acc=opt_vals.acc*4;
            opt_vals.numstepswithsamejac= opt_vals.maxnumstepswithsamejac;
            rootfound = findrootbysnewton(dim, jac0,ipiv, troot, var0, &(var1[0]),&(res0[0]),&fn, b1, vara, ta, preal,numpreal,pint,numpint, &auxfactor,&opt_vals );
            opt_vals.acc=opt_vals.acc/4.;
	    if (rootfound) 
		{//ondoren saiatu f_n(x) ere betetzen duen findroot egiten eta lortzen badu hori itzuli. bestela aurrekoa.
		vectorcp(dim,var1, &(var0[0]));
		opt_vals.acc=opt_vals.acc*4;
        	opt_vals.numstepswithsamejac= opt_vals.maxnumstepswithsamejac;
		// we will try to get the root of the dim equations, without the hiperplane equation as last eqution.
		// to do it we call to findroot with parameters 0 as the pointers to the vector of previous values
        	rootfound = findrootbysnewton(dim, jac0,ipiv,troot, var0, &(var1[0]),&(res0[0]),&fn, (double *)0, (double *)0, ta, preal,numpreal,pint,numpint, &auxfactor,&opt_vals);
        	opt_vals.acc=opt_vals.acc/4.;
		if (rootfound)
		    {  // we have the root for all f_i(x) in the vector var1
			// the hiperplane coefficients are in vector b1
		    *infoptr = 111;
		    }
		  else
		    {  // we have the root for all f_i(x) but the jacobian may be singular we will return the vector obtained advancing in the curve.
			// the hiperplane coefficients are in vector b1
		    vectorcp(dim, var0, &(var1[0]));
		    *infoptr = 112;
		    }
		*taptr = troot;
		*tbptr = taux;

		}
	      else
		{ // there was a surechange but we could not get it at higher accuracy (why not? it should not happen this)
		  // the hiperplane coefficients are in vector b1
		vectorcp(dim, var0, &(var1[0]));
		*taptr = troot;
		*tbptr = taux;
		*infoptr = 112;
		}
	    }
	  else // we have not found any point in the curve for which the last equation is 0
  	    {
	    if (taux != init) // we could advance
		{    // we return the var1 and b1 vectors
		*taptr = init;
		*tbptr = taux;
		if (is_closed)
		    *infoptr = 20;
		  else
		    *infoptr = 10;
		}
	      else 
		{     // we have to return the initial vector
		vectorcp(dim, var0, &(var1[0]));
		vectorcp(dim, b00, &(b1[0]));
		*taptr = init;
		*tbptr = init;
		*infoptr = -20;
		}
	    }
	break;
    }
free(jac0);
free(ipiv);
free(res0);
free(vara);
free(var00);
free(b00);
if (action == 1) 
	{
	free(bmin);
	}
if ((action == 1) || (action == 2))
	{
	free(varaux);
	}
#ifdef LOGFILE
fprintf(loga,"end of arcContinuation!\n");
fclose(loga);
#endif
}






