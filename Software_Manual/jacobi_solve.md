# Solve a Square Linear System Using the Jacobi Method

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           jacobi_solve

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c jacobi_solve.f90```

and can be added to a program using

```$ gfortran program.f90 jacobi_solve.o ``` 

**Description/Purpose:** This routine uses the Jacobi method (*simultaneous relaxation*) to iteratively solve the system of equations given by:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" title="\mathbf{A}{\mathbf{x}} = \vv{\mathbf{b}}" /></a>

by using the following algorithm

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;x_i^{(k&space;&plus;&space;1)}&space;=&space;\frac{1}{a_{ii}}\left[b_i&space;-&space;\sum\limits_{\substack{j=1&space;\\&space;j\neq&space;i}}^na_{ij}x_j^k&space;\right&space;],\qquad&space;i&space;=&space;1,...,n." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;x_i^{(k&space;&plus;&space;1)}&space;=&space;\frac{1}{a_{ii}}\left[b_i&space;-&space;\sum\limits_{\substack{j=1&space;\\&space;j\neq&space;i}}^na_{ij}x_j^k&space;\right&space;],\qquad&space;i&space;=&space;1,...,n." title="x_i^{(k + 1)} = \frac{1}{a_{ii}}\left[b_i - \sum\limits_{\substack{j=1 \\ j\neq i}}^na_{ij}x_j^k \right ],\qquad i = 1,...,n." /></a>

**Input:** 

*A* : REAL - the input coefficient matrix of size *n* x *n*

*n* : INTEGER - the size of *A* and the length of *b* and *x*

*b* : REAL - the right hand size vector of length *n*

*tol* : REAL - the tolerance between the iterations

*maxiter* : INTEGER - a limit on the maximum number of iterations to be performed by the Jacobi iteration

*printit* : INTEGER - if equal to 1, prints the convergence error and the iteration count

**Output:** 

*x* : REAL - the solution to the system of equations given by *Ax* = *b*

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, maxiter, i, printit
REAL*8, ALLOCATABLE :: A(:, :), b(:), x(:), Ab(:, :)
REAL*8 :: tol

n = 3
maxiter = 1000
tol = 1e-16
printit = 0
ALLOCATE(A(n, n), b(n), x(n), Ab(n, n + 1))
A = RESHAPE((/10.D0, 1.D0, 3.D0, &
& 1.D0, 10.D0, 0.D0, &
& 3.D0, 2.D0, 10.D0/), (/n, n/), ORDER=(/2, 1/))
b = (/2.D0, 4.D0, 1.D0/)
x = 0.D0
CALL jacobi_solve(A, n, b, tol, maxiter, x, printit)
WRITE(*,*) x
```

The outputs from the above code:

```fortran
  0.16997792494481237       0.38300220750551878       -2.7593818984547436E-002
```

To verify the routine, it can be compared to the [Gaussian Elimination](./direct_ge_bsin.md) subroutine as follows

```fortran
	IMPLICIT NONE

	INTEGER :: n, maxiter, i, printit
	REAL*8, ALLOCATABLE :: A(:, :), b(:), x(:), Ab(:, :)
	REAL*8 :: tol

	n = 3
	maxiter = 1000
	tol = 1e-16
	printit = 0
	ALLOCATE(A(n, n), b(n), x(n), Ab(n, n + 1))
	A = RESHAPE((/10.D0, 1.D0, 3.D0, &
				& 1.D0, 10.D0, 0.D0, &
				& 3.D0, 2.D0, 10.D0/), (/n, n/), ORDER=(/2, 1/))
	b = (/2.D0, 4.D0, 1.D0/)
	Ab(:n, :n) = A
	Ab(:, n + 1) = b
	x = 0.D0
	CALL jacobi_solve(A, n, b, tol, maxiter, x, printit)
	WRITE(*,*) "Jacobi"
	WRITE(*,*) x
	x = 0.D0
	CALL direct_ge_bsin(Ab, n, n + 1, x)
	WRITE(*,*) "Gauss-Elimination"
	WRITE(*,*) x
```

with the output

```fortran
 Jacobi
  0.16997792494481237       0.38300220750551878       -2.7593818984547436E-002
 Gauss-Elimination
  0.16997792494481237       0.38300220750551872       -2.7593818984547460E-002

```

In addition, with a larger, diagonally dominant matrix (n = 1000), the difference between the Jacobi iteration and the Gaussian-Elimination routine can be seen below for the first 20 elements (and it continues through the rest of the unknowns):

```fortran
2.3852447794681098E-018   0.0000000000000000        2.1684043449710089E-018   6.5052130349130266E-019   2.8189256484623115E-018   2.0599841277224584E-018   2.1684043449710089E-018   8.6736173798840355E-019   4.8789097761847700E-019   3.3881317890172014E-019   2.1684043449710089E-019   4.8789097761847700E-019   4.4045713257223618E-019   1.7347234759768071E-018   0.0000000000000000        1.0842021724855044E-019   2.3852447794681098E-018   4.3368086899420177E-019   1.5178830414797062E-018   5.9631119486702744E-019
```



**Implementation/Code:** The code for jacobi_solve can be seen below.

```fortran
SUBROUTINE jacobi_solve(A, n, b, tol, maxiter, x0, printit)
	IMPLICIT NONE
	
	! Takes as inputs the coefficient matrix, `A` of size `n` and the
	! right-hand side vector `b` of length `n`. Also as inputs are a
	! tolerance for exit of the algorithm, `tol` and a maximum number of
	! iterations `maxiter`. An initial guess, `x0` is input and refined
	! throughout the algorithm with each successive iteration. When the
	! algorithm exits, `x0` is output as the final approximation of x.
	! An input that tells the algorithm to print the final convergence
	! error and iteration count is also an input.
	INTEGER, INTENT(IN) :: n, maxiter, printit
	REAL*8, INTENT(IN) :: A(n, n), b(n), tol
	REAL*8, INTENT(INOUT) :: x0(n)
	
	INTEGER :: i, j, k
	REAL*8 :: error, x1(n), sum_ax
	
	x1 = 0.D0
	sum_ax = 0.D0
	error = 10.D0*tol
	k = 0
	
	! Iteration loop for the algorithm. Iterates until the error is
	! less than the given tolerance or the maximum number of iterations
	! is exceeded.
	DO WHILE (error > tol .AND. k < maxiter)
		
		! Implements the Jacobi algorithm of x_{k+1} = x_k + D^{-1}*r_k.
		DO i = 1, n
			sum_ax = b(i)
			DO j = 1, i - 1
				sum_ax = sum_ax - A(i, j)*x0(j)
			END DO
			DO j = i + 1, n
				sum_ax = sum_ax - A(i, j)*x0(j)
			END DO
			x1(i) = sum_ax/A(i, i)
		END DO
		
		! Increments the counter and calculates the absolute error
		! between x0 and x1 using the l2 norm. Then sets x0 = x1 for 
		! the next iteration.
		k = k + 1
		CALL abs_err_vecl2(x0, x1, n, error)
		x0 = x1
	END DO
	IF (printit == 1) THEN
		WRITE(*,*) error, k
	END IF
END SUBROUTINE
```
