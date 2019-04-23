# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           steepest_descent

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c steepest_descent.f90```

and can be added to a program using

```$ gfortran program.f90 steepest_descent.o ``` 

**Description/Purpose:** This routine uses the steepest descent method to iteratively solve the system of equations given by:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" title="\mathbf{A}{\mathbf{x}} = \vv{\mathbf{b}}" /></a>

by choosing to minimize the residual vector, *r<sub>k</sub>* and varying the step size in the direction of steepest descent, i.e.:

<a href="https://www.codecogs.com/eqnedit.php?latex=x_{k&plus;1}&space;=&space;x_k&space;&plus;&space;\alpha_k{r_k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_{k&plus;1}&space;=&space;x_k&space;&plus;&space;\alpha_k{r_k}" title="x_{k+1} = x_k + \alpha_k{r_k}" /></a>

**Input:** 

*A* : REAL - the input coefficient matrix of size *n* x *n*

*n* : INTEGER - the rank of *A* and the length of *b* and *x*

*b* : REAL - the right hand size vector of length *n*

*tol* : REAL - the tolerance between the iterations

*maxiter* : INTEGER - a limit on the maximum number of iterations to be performed by the Jacobi iteration

*printit* : INTEGER - if equal to 1, the algorithm prints the final convergence error and iteration count

*x0* : REAL - the initial guess to the solution of Ax = bS

**Output:** 

*x0* : REAL - the solution to the system of equations given by *Ax* = *b*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, maxiter, i, printit
REAL*8, ALLOCATABLE :: A(:, :), b(:), x(:)
REAL*8 :: tol, error

n = 3
maxiter = 1000
tol = 10.D-16
printit = 1
ALLOCATE(A(n, n), b(n), x(n))
A = RESHAPE((/7.D0, 3.D0, 1.D0, &
& 3.D0, 10.D0, 2.D0, &
& 1.D0, 2.D0, 15.D0/), (/n, n/), ORDER=(/2, 1/))
b = (/28, 31, 22/)
x = 0.D0
CALL steepest_descent(A, n, b, tol, maxiter, printit, x)
WRITE(*,*) x
```

The outputs from the above code:

```fortran
   3.3287017925713278E-016          31
   2.9999999980826058        2.0000000016423951        1.0000000006619756  
```

This matches the answer given in Example 7.9 in Ascher (Ascher, U., and C. Greif. "A First Course in Numerical Methods. SIAM." *Society for* (2011).) of x = (3, 2, 1).

**Implementation/Code:** The code for steepest_descent can be seen below.

```fortran
SUBROUTINE steepest_descent(A, n, b, tol, maxiter, printit, x0)
	IMPLICIT NONE
	
	! Implements the steepest descent algorithm for solving a linear
	! system of equations.
	
	! Takes as inputs the coefficient matrix `A` of rank `n` as well
	! as the right-hand side vector `b` of length `n`. The argument
	! `tol` is given as the tolerance of convergence error and the
	! `maxiter` argument gives a limit to the number of iterations in
	! the algorithm. Finally, `printit` prints to the screen the
	! convergence error and iteration count if set to a value of 1.
	! The `x0` argument is passed in as the initial guess of the
	! solution to the system of equations and is passed out as the
	! final approximation at the end of the subroutine.
	INTEGER, INTENT(IN) :: n, maxiter, printit
	REAL*8, INTENT(IN) :: A(n, n), b(n), tol
	REAL*8, INTENT(INOUT) :: x0(n)
	
	! Sets up temporary variables that will be used in the algorithm.
	REAL*8 :: alpha, r0(n), r1(n), x1(n), d, s(n), temp
	INTEGER :: k
	
	! Initialize variables
	k = 0
	r0 = 0.D0
	r1 = 0.D0
	x1 = 0.D0
	d = 0.D0
	temp = 0.D0
	
	! Calculate the initial residual using the initial guess `x0`
	CALL mat_prod(A, x0, n, n, 1, r0)
	r0 = b - r0
	
	! Calculates the initial "error" using the residual
	CALL vec_dot_prod(r0, r0, n, d)
	
	! Iterates using the steepest descent algorithm until the "error" is
	! less than the given tolerance or the maximum number of iterations
	! is exceeded.
	DO WHILE (d > tol .AND. k < maxiter)
	
		! Calculate the product of `A` and the residual.
		CALL mat_prod(A, r0, n, n, 1, s)
		
		! Calculate the step size, `alpha`.
		CALL vec_dot_prod(r0, s, n, temp)
		alpha = d/temp
		
		! Increment `x0` and `r0` using the step size and the direction
		! `r0` that will drive the solution to a minimum residual.
		x1 = x0 + alpha*r0
		r1 = r0 - alpha*s
		
		! Calculate the "error" `d` (for delta).
		CALL vec_dot_prod(r1, r1, n, d)
		
		! Update variables and increment the counter.
		x0 = x1
		r0 = r1
		k = k + 1
	END DO
	
	! If the `printit` argument is provided as 1, then print the final
	! convergence error and number of iterations.
	IF (printit == 1) THEN
		WRITE(*,*) d, k
	END IF
	
END SUBROUTINE
```
