# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          K2_cond

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c K2_cond.f90```

and can be added to a program using

```$ gfortran program.f90 K2_cond.o ``` 

**Description/Purpose:** This subroutine uses the [power method](./power_method.md) and [inverse iteration](./inverse_iteration.md) method to calculate the highest and lowest eigenvalue of a ***symmetric positive definite matrix***. It then returns an approximation of the l2 condition number of the matrix by taking the absolute value of the highest over the lowest eigenvalue, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=||\kappa||_2&space;\approx&space;\frac{\lambda_1}{\lambda_n}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?||\kappa||_2&space;\approx&space;\frac{\lambda_1}{\lambda_n}" title="||\kappa||_2 \approx \frac{\lambda_1}{\lambda_n}" /></a>

**Input:** 

*A* : REAL - a symmetric positive definite coefficient matrix of size *n* x *n*

*n* : INTEGER - the rank of A and the length of the eigenvectors

*v0* : REAL - the initial guess of the eigenvetor

*tol* : REAL - the exit tolerance for the algorithm

*maxiter* : INTEGER - the maximum number of iterations before exit

*printit* : INTEGER - flag for printing final convergence information

**Output:** 

*cond* : REAL - the approximation of the l2 condition number

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, maxiter, i, printit
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam

n = 3
ALLOCATE(A(n, n), v0(n), v(n))
lam = 0.D0
tol = 10.D-16
maxiter = 10000
printit = 1
v0 = 1.D0
A = RESHAPE((/3.D0, 0.D0, 0.D0, &
			& 0.D0, 6.D0, 0.D0, &
			& 0.D0, 0.D0, 2.D0/), (/n, n/), ORDER=(/2, 1/))
CALL K2_cond(A, n, v0, tol, maxiter, printit, cond)
WRITE(*,*) cond
```

The outputs from the above code:

```fortran
   ! Power Method
   0.0000000000000000               28  ! Convergence error and number of iterations at exit
   6.0000000000000000     ! Eigenvalue found
   ! Inverse Iteration
   0.0000000000000000               19  ! Convergence error and number of iterations at exit
   2.0000000000000000     ! Eigenvalue found
   3.0000000000000000 ! Condition number
```

Additionally, using this routine on an 20 x 20 Hilbert matrix can be done with the following code:

```fortran
INTEGER :: n, maxiter, i, printit, j
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam

n = 20
ALLOCATE(A(n, n), v0(n), v(n))
lam = 0.D0
tol = 10.D-16
maxiter = 10000
printit = 1
v0 = 1.D0
DO i = 1, n
	DO j = 1, n
		A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
	END DO
END DO
CALL K2_cond(A, n, v0, tol, maxiter, printit, cond)
WRITE(*,*) cond
```

and the output is

```fortran
   0.0000000000000000               15
   1.1175345637892222E-017           1
   1.7065554679049347E+017           ! Awesome....
```



**Implementation/Code:** The code for K2_cond can be seen below.

```fortran
SUBROUTINE K2_cond(A, n, v0, tol, maxiter, printit, cond)
	IMPLICIT NONE
	
	! Finds an approximation of the l2 condition number using the `power
	! method` and `inverse iteration` subroutines.
	
	! Takes as an input the matrix `A` of rank `n` that contains the
	! system to be analyzed, an initial guess for the eigenvector, `v0`,
	! a tolerance, `tol`, for exiting the iterative solver as well as a
	! maximum number of iterations, `maxiter`. Finally, a flag to print
	! the final number of iterations and convergence error, `printit` is
	! an input. The output of the subroutine is an approximation of the
	! l2 condition number. ***Note that `A` needs to be a symmetric
	! positive definite matrix.
	INTEGER, INTENT(IN) :: n, maxiter, printit
	REAL*8, INTENT(IN) :: A(n, n), v0(n), tol
	REAL*8, INTENT(OUT) :: cond
	
	REAL*8 :: vhigh(n), vlow(n), lam_1, lam_n, alpha
	
	! Find the lowest eigenvalue using the inverse iteration method.
	alpha = 0.D0
	
	! Calculates the largest eigenvalue.
	CALL power_method(A, n, v0, tol, maxiter, printit, vhigh, lam_1)
	WRITE(*,*) lam_1
	
	! Calculates the smallest eigenvalue.
	CALL inverse_iteration(A, n, v0, alpha, tol, maxiter, printit, vlow, lam_n)
	WRITE(*,*) lam_n
	
	! Estimates the l2 condition number.
	cond = lam_1/lam_n
END SUBROUTINE
```