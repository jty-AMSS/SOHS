# FSOS&SOHS

The new version of computation of FSOS and SOHS, which can compute function on every finite abelian group and easy to use.

and it is able to computing SOHS via FSOS on finite abelian groups.

usage:

Compute the SOHS of function on T^n, the main funtion is Fun_SOHS:

Example 1: verify 1-x^10-y^10-z^10>=0 for (x^2+y^2+z^2=1) with  relaxation order 50, Matlab code:
  
 et3=@(x,y,z)1-x^10-y^10-z^10;
 
 F=@(x)et3(1/2*(x(1)+conj(x(1))),(-1i)*1/2*(conj(x(1))-x(1))*(1/2)*(x(2)+x(2)'),(-1i)*1/2*(conj(x(1))-x(1))*(-i/2)*(x(2)-x(2)'));
 
 [Q,Index,SF,SOHS,err]=Fun_SOHS(F,[200,200])
 
 
Example 2: verify Motzkin polynomial is nonnegative on [-2,2]*[-2,2], Matlab code:

 M=@(x,y)x^4*y^2+x^2*y^4-3*x^2*y^2+1;
 
 F=@(x)M(x(1)+conj(x(1)),x(2)+conj(x(2)));
 
 [Q,Index,SF,SOHS,err]=Fun_SOHS(F,[8,8])
 


Other Examples:

MAX_SATExample.m  for example of MAX-3SAT

sphere_norm_p.m for computing the SOHS of 1-x^10-y^10 {(x,y):x^2+y^2=1} or 1-x^10-y^10-z^10 on {(x,y,z):x^2+y^2+z^1=1}

REFERENCES:
Computing sparse Fourier sum of squares on finite abelian groups in quasi-linear time." arXiv preprint arXiv:2201.03912 (2022)
