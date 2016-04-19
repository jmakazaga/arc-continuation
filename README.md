# arc-continuation
given n-1 equations f_1(x), ... f_{n-1}(x), where x \in R^n (yes, x has 1 more element than number of equations) and 
given x_0 for which f_i(x_0) = 0, there is a continuous curve for x values fao which f_i(x) = 0. 
We will move in this curve starting at x_0. This movement can be in two directions, we define an horthogonal 
hiperplane of the form
b(x-x_0) = s -s_1, 
where b,x,x_0 \in R^n and s,s_0 \in R. 
So that the vector b is hortoghonal to the curve.
With b we are able to move on the curve for which f_i(x) = 0 and we will get different values of x. 
While we are moving on the curve we can evaluate some f_n(x) and we can do several things with this value:
a.- we can return the minimum value for f_n(x) in the covered path.
b.- we can test the 0 value for f_n(x) (the sign cahnges in the step) and stop. So, we are near from the x value that
    is the root for F(x) (F: R^n -> R^n)
