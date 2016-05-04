#include <homotopia.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <complex.h>



void rest(int dim, val_type *var, time_type t, val_type *res, double *preal, long prelem, int *pint, long pielem)
{


<* If[Length[auxvars]>0, 
       SequenceForm["val_type aux[", Length[auxvars],"];"],
      " "] *>

<*auxaux//ColumnForm*>

<*resaux//ColumnForm*>
}


void projection(int dim,val_type *varvals, time_type t, double *preal, int numreal, int *pint, int numint) 
{
//hiperplaneortho(dim, varvals, t, preal, numreal, pint, numint);
 <* If[hvarsQ,
     {"int i;",
      "val_type norma;",
      "norma = Norm(dim,varvals);",
      "for (i=0; i<dim;i++) varvals[i]/=norma;"}//ColumnForm,
     " "]
   *>
}

int extra_len()
{
<* SequenceForm["return(", Length[resta]-Length[splicevars],");"] *>
}

void jacobian(int dim, val_type *var, time_type t, val_type *jac, double *preal, long prelem, int *pint, long pielem)
{

<* If[Length[auxvars]>0, 
       SequenceForm["val_type aux[", Length[auxvars],"];"],
      " "] *>

<* If[Length[dauxvars]>0, 
       SequenceForm["val_type daux[", Length[dauxvars],"];"],
      " "] *>

<*auxaux//ColumnForm*>

<*dauxaux//ColumnForm*>

<* jacaux//ColumnForm *>
}

