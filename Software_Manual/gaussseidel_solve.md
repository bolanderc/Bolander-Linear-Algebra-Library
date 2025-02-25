# Solve a Square Linear System Using Gauss-Seidel

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           gaussseidel_solve

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c gaussseidel_solve.f90```

and can be added to a program using

```$ gfortran program.f90 gaussseidel_solve.o ``` 

**Description/Purpose:** This routine uses the Gauss-Seidel method to iteratively solve the system of equations given by:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}{\mathbf{x}}&space;=&space;\vv{\mathbf{b}}" title="\mathbf{A}{\mathbf{x}} = \vv{\mathbf{b}}" /></a>

by using the following algorithm

<a href="https://www.codecogs.com/eqnedit.php?latex=x_i^{(k&plus;1)}&space;=&space;\frac{1}{a_{ii}}\left[b_i&space;-&space;\sum\limits_{j<i}a_{ij}x_j^{(k&plus;1)}&space;-&space;\sum\limits_{j>i}a_{ij}x_j^{(k)}&space;\right&space;],&space;\qquad&space;i&space;=&space;1,...,n." target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_i^{(k&plus;1)}&space;=&space;\frac{1}{a_{ii}}\left[b_i&space;-&space;\sum\limits_{j<i}a_{ij}x_j^{(k&plus;1)}&space;-&space;\sum\limits_{j>i}a_{ij}x_j^{(k)}&space;\right&space;],&space;\qquad&space;i&space;=&space;1,...,n." title="x_i^{(k+1)} = \frac{1}{a_{ii}}\left[b_i - \sum\limits_{j<i}a_{ij}x_j^{(k+1)} - \sum\limits_{j>i}a_{ij}x_j^{(k)} \right ], \qquad i = 1,...,n." /></a>

**Input:** 

*A* : REAL - the input coefficient matrix of size *n* x *n*

*n* : INTEGER - the size of *A* and the length of *b* and *x*

*b* : REAL - the right hand size vector of length *n*

*tol* : REAL - the tolerance between the iterations

*maxiter* : INTEGER - a limit on the maximum number of iterations to be performed by the Jacobi iteration

*printit* : INTEGER - if equal to 1, the algorithm prints the final convergence error and iteration count

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
CALL gaussseidel_solve(A, n, b, tol, maxiter, x)
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
	CALL gaussseidel_solve(A, n, b, tol, maxiter, x, printit)
	WRITE(*,*) "Gauss-Seidel"
	WRITE(*,*) x
	x = 0.D0
	CALL direct_ge_bsin(Ab, n, n + 1, x)
	WRITE(*,*) "Gauss-Elimination"
	WRITE(*,*) x
```

with the output

```fortran
 Gauss-Seidel
  0.16997792494481237       0.38300220750551878       -2.7593818984547436E-002
 Gauss-Elimination
  0.16997792494481237       0.38300220750551872       -2.7593818984547460E-002

```

In addition, with a larger, diagonally dominant matrix (n = 1000), the difference between the Gauss-Seidel method and the Gaussian-Elimination routine can be seen below as the l2 norm of the absolute error between the two returned vector:

```fortran
4.1975798679441814E-017
```



**Implementation/Code:** The code for gaussseidel_solve can be seen below.

```fortran
SUBROUTINE gaussseidel_solve(A, n, b, tol, maxiter, x0, printit)
	IMPLICIT NONE
	
	! Takes as an input the square, coefficient matrix `A` of size `n`
	! as well as the right-hand side vector `b`. In addition, a
	! tolerance argument `tol` is used to determine sufficient
	! convergence and a `maxiter` argument specifies the maximum number
	! of iterations before the solver exits. An initial guess for the
	! solution, `x0`, is given as an input. It stores the iterations of
	! the Gauss-Seidel approximation to the solution x and is output as
	! soon as the solver exits. 
	! An input that tells the algorithm to print the final convergence
	! error and iteration count is also an input.
	INTEGER, INTENT(IN) :: n, maxiter, printit
	REAL*8, INTENT(IN) :: A(n, n), b(n), tol
	REAL*8, INTENT(INOUT) :: x0(n)
	
	INTEGER :: i, j, k
	REAL*8 :: ax_sum, x1(n), error
	
	x1 = x0
	k = 0
	error = 10.D0*tol
	
	! Starts the iterative Gauss-Seidel solver.
	DO WHILE (error > tol .AND. k < maxiter)
		
		! Finds the next iteration of the vector x.
		DO i = 1, n
			ax_sum = b(i)
			
			! Sets the Gauss-Seidel routine apart from the Jacobi in
			! that it uses any updated values previously calculated in 
			! the iterations to give a more 'accurate' approximation.
			DO j = 1, i - 1
				ax_sum = ax_sum - A(i, j)*x1(j)
			END DO
			DO j = i + 1, n
				ax_sum = ax_sum - A(i, j)*x0(j)
			END DO
		
			x1(i) = ax_sum/A(i, i)
			
			! The x vector approximation is updated, an error is
			! calculated, and the x0 vector is reset to the approximation
			! given in this iteration.
		END DO
		k = k + 1
		CALL abs_err_vecl2(x0, x1, n, error)
		x0 = x1
	END DO
	IF (printit == 1) THEN
		WRITE(*,*) error, k
	END IF
END SUBROUTINE
```
