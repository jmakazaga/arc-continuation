void 
void mathArcContinuation P(( double *v, long n, double *b, long nb, double init, double end, int action, double *preal, long numpreal, int *pint, long numpint));

:Begin:
:Function:       mathArcContinuation
:Pattern:        mathArcContinuation[ v_List, b_List, init_Real, end_Real, action_Integer,preal_List, pint_List]
:Arguments:      {v, b, init, end, action, preal, pint}
:ArgumentTypes:  { RealList, RealList, Real, Real, Integer, RealList, IntegerList }
:ReturnType:     Manual
:End:

:Evaluate: mathArcContinuation::usage = "mathArcContinuation[{v1,v2,v3...vn}, {b0, b1, b2... bn}, t0, t1, action {preal0, preal1...prealm}, {pint1, pint2...pints}] starting with V0={v1,v2,v3...}, which is the solution of the n-1 equations f_i(V0)=0 (for i = 1...n-1), moves the vector V0 on the curve defined by the n-1 equations. This curve can be followed in two direcctions, the vector B implicitly defines the direction. At each step the process obtains the hiperplane defined by B_o(V- Vp) = t-tp which is orthogonal to the curve, so that the vector B_o indicates the direction for the movement. Using this direction the proces moves in the curve step by step. Action indicates the action to be done at each step:
0 means that we want just to move on the curve. 1 means we want to evaluate f_n(V) at each step and we will return the V value that minimizes it in the covered path. Action 2 means that we are looking for the V value that makes f_n(V) = 0  The return value is the list {V, B, t, tend, info}.  Depending on the action and the possibility to move on the curve they will be different values. In general if all was right they will be: if action is 0 the vector V is the point at the end of the covered curve and B represents the orthogonal direction to the hiperplane; info=1, t=t0 and tend=t1 (equal to input parameters). If action is 1 the vector V is the point with the minimum value for f_n(V), the vector B is the ortogonal direction for the hiperplane at this point, t is the parameter tp in the hiperplane equation, tend is equal to input value t1 and info = 11. If action is 2 the vector V is the point that makes f_n(V) = 0. The vector B and t are the coefficient vector and the tp value of the hiperplane at this point, tend is t1 (input parameter) and info=11. "

void realRest(double *v, long vlen, double t, double *preal, long prlen, int *pint, long pilen);

:Begin:
:Function:       realRest
:Pattern:        realRest[ v_List, t_Real, preal_List, pint_List]
:Arguments:      {v, t, preal, pint}
:ArgumentTypes:  { RealList, Real, RealList, IntegerList }
:ReturnType:     Manual
:End:

:Evaluate: realRest::usage = "realRest[{v1,v2,v3...vn}, t, {preal0, preal1...prealm}, {pint1, pint2...pints}] evaluates f(t,V): R^{n+1} -> R^n and returns the list {F, ||F(t,V)||}, where F= {F1, F2, F3,...Fn}. Be careful: t must be Real! not Integer."


void realJacobian(double *v, long vlen, double t, double *preal, long prlen, int *pint, long pilen);

:Begin:
:Function:       realJacobian
:Pattern:        realJacobian[ v_List, t_Real, preal_List, pint_List]
:Arguments:      {v, t, preal, pint}
:ArgumentTypes:  { RealList, Real, RealList, IntegerList }
:ReturnType:     Manual
:End:

:Evaluate: realJacobian::usage = "realJacobian[{v1,v2,v3...vn}, t, {preal0, preal1...prealm}, {pint1, pint2...pints}] evaluates d F(t,V)/d V and returns the Jacobian of F, which is a matrix of dimensions n x n."


void realFindroot(double *v, long vlen, double t, double *preal, long prlen, int *pint, long pilen, double acc, double prec, int approcjac);

:Begin:
:Function:       realFindroot
:Pattern:        realFindroot[ v_List, t_Real, preal_List, pint_List, acc_Real, prec_Real, approxjac_Integer]  
:Arguments:      {v, t, preal, pint, acc, prec, approxjac}
:ArgumentTypes:  { RealList, Real, RealList, IntegerList, Real, Real, Integer }
:ReturnType:     Manual
:End:

:Evaluate: realFindroot::usage = "realFindroot[{v1,v2,v3...vn}, t, {preal0, preal1...prealm}, {pint1, pint2...pints}, acc, prec, 0] gets a root of F(t,v) and returns the list {{r1, r2,...rn},retval} where r is the root and retvalue is 1 if the root has been found, otherwise has a different value. The input value acc is used to stop the Newton iterations of the process: when the root approximation has changed more than the first approximation's change multiplied by 10.0^-acc then the process stops and it returns the root approximation as the solution. It is possible that the last change does not reach the change requested by the accuracy, but, if the las iteration change is less than 10^-prec then the root is considered as a solution. The last input parameter guides the Newton iteration process in the sense of the use of the jacobien of F(t,v), that is, if the value is 0 then it uses the *jacobian* of F, but if this value is other number then the process approximates numerically the Jacobian."


