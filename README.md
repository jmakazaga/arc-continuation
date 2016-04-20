# arc-continuation
Given n-1 equations f_1(x), ... f_{n-1}(x), where x \in R^n (yes, x has 1 more element than number of equations) and 
given x_0 for which f_i(x_0) = 0, there is a continuous curve for x values for which f_i(x) = 0. 
We will move in this curve starting at x_0. This movement can be in two directions, we define an horthogonal 
hiperplane of the form
b(x-x_0) = s -s_0, 
where b,x,x_0 \in R^n and s,s_0 \in R. 
So that the vector b is hortoghonal to the curve.
With b we are able to move on the curve for which f_i(x) = 0 and we will get different values of x. 
While we are moving on the curve we can evaluate some f_n(x) and we can do several things with this value:
a.- we can return the x that has minimum value for f_n(x) in the covered path.
b.- we can test the 0 value for f_n(x) (the sign cahnges in the step) and stop. So, we are near from the x value that
    is the root for F(x) (F: R^n -> R^n)

The function has as input: (X_0, b_0, s_0, s_1, action) where:
.- x_0 makes f_i(x_0) = 0 for i \in (1,n-1)
.- the vector b_0 is simply a vector that will help us stablishing the direction of the movement (it can be one component with value 1 and the rest with value 0). It will be used to get the orthogonal b_h that makes b_h(x-x_0) = s - s_0.
.- s_0 is the starting s value
.- s_1 is the final s value. These two values indicate the amount of movement on the curve.
.- action is the value taht tells the action to be done at each step: action = 0 if you just want to move from s_0 to s_1. action = 1 if you are looking for the x that makes f_n(x) be the minimun in the covered path. action = 2 if you are looking for the x value that makes f_n(x) =0. 

The output is: (x_e, b_e, s_a, s_e, info) where these values depend on the action defined as input:
If action was 0 (go from s_0 to s_1) then info can be:
    .- info = 1. So, all went right and  s_a = s_0, s_e = s_1, x_e and b_e are these that make f_i(x_e) = 0 for i \in (1,n-1) and we have been able to move until b_e(x-x_e) = s -s_e
    .- info = 2. So it has detected that the curve is a closed curve and we have arrived to the starting x_0. In this case we have x_e = x_0, s_a = s_0, s_e = s_p (where s_p is the length of the period) b_e is the vector orthoghonal to the curve at x_0.
    .- info = -2. So it was not posible to move on the curve, there was a singularity just at x_0. Return values are x_e = x_0, s_a = s_0, s_e = s_0 and b_e = b_h ( where b_h(x-x_0) = s -s_0).
    .- info = -1. The process has been able to move on the curve but it has stoped because it could not continue (a singularity in the patha s_s). The return values are the point until it has been able to move to: x_e and f_i(x_1) = 0, b_1 is the orthogonal direction at this point, s_a = s_0 and s_e = s_s a value between s_0 and s_1 (the point untill it was able to go to)
If action was 1 (go from s_0 to s_1 and give me the point that makes f_n(x) minimun)
    .- Info = 10, we arrived to s_1, but the minimum was the starting point, but we return the last point: x_e = x_1, b_e = b_1, s_a = s_0, s_e = s_1.
    .- info = 11, we arrived to s_1 and the minimum is a point in the midle of the path. return values: x_e = x_min for which f_n(x_min) is the minimum. s_a  = s_0, s_e = s_min, b_e = b_min
    .- info = 12, we arrived to s_1 and the minimum is at the last point. Return values are same as info = 10.
    ---
    .- info = 20. The curve is closed and the minimun is for the starting point. returns x_e = x_0, s_a = s_0, s_e = period of the curve, b_e = b_0 but ortoghonal.
    .- info = 21. The curve is closed and the minimun has been reached at s_min. return values x_e = x_min, s_a = s_min, s_b = period of the curve, b_e = b_min
    ---
    .- info = -20. It was not posible to move on the curve. return values are x_0, b_0, s_a=s_0, s_b = s_0
    ---
    .- info = -10. The process has been able to move on the curve but it has stoped because it could not continue (a singularity in the path at s_s) and the minimun value has been reached at the starting point. Reurn values correspond to the singularity point: x_e = x_s, b_e = b_s, s_a = s_0 s_e = s_s. 
    .- info = -11. The process has been able to move on the curve but it has stoped because it could not continue (a singularity in the path at s_s) and the minimun value has been reached at some point between s_0 and the singularity point s_s. Returns the min value: x_e = x_min, b_e = b_min, s_a = s_min, s_e = s_s. 
    .- info = -12. The process has been able to move on the curve but it has stoped because it could not continue (a singularity in the path at s_s) and the minimun value has been reached at the stoping point s_s. Returns the singularity point: x_e = x_s, b_e = b_s, s_a = s_0 s_e = s_s. 
if action was 2 ( go from s_0 to s_1 but stop when f_n(x) = 0)
    In this case info values are same as values for action 1. and return values are, when f_n(x) = 0 is reached, the point tha makes F(x) = 0. 
