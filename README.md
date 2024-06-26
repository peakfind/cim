# CIM
 [`Julia`](http://julialang.org) code for the Contour Integral Method

# Reference
There are three references for the Contour Integral Method：
+ Wolf-J&uuml;rgen Beyn, [An integral method for solving nonlinear eigenvalue problems, **Linear Algebra Appl., 2012**](https://doi.org/10.1016/j.laa.2011.03.030)，
+ Wolf-J&uuml;rgen Beyn, Yuri Latushkin and Jens Rottmann-Matthes, [Finding eigenvalues of holomorphic Fredholm operator pencils using boundary value problems and contour integrals, **Integr. Equ. Oper. Theory, 2014**](https://doi.org/10.1007/s00020-013-2117-6)，
+ Reinhard Mennicken and Manfred M&ouml;ller, Non-self-adjoint Boundary Eigenvalue Problems, **Mathematics Studies 192, North-Holland, 2003**.

# Implementation
We plan to refer to some `Julia` codes which implement the Contour Integral Method:
+ the native `Julia` implementation in [**libBeyn**](https://github.com/HomerReid/libBeyn) writen by Homer Reid. 
+ [**NEP-PACK/NonlinearEigenproblems.jl**](https://github.com/nep-pack/NonlinearEigenproblems.jl).

A basic `MATLAB` implementation can be found in Figure 5.3 in Stefan G&uuml;ttel, Francoise Tisseur, [The nonlinear eigenvalue problem, **Acta Numerica, 2017**](https://doi.org/10.1017/S0962492917000034).