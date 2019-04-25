# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          ls_solvejacobi

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c ls_solvejacobi.f90```

and can be added to a program using

```$ gfortran program.f90 ls_solvejacobi.o ``` 

**Description/Purpose:** This subroutine solves the least squares problem, given by

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{min}_\mathbf{\vv{x}}&space;||\mathbf{\vv{b}}&space;-&space;A\mathbf{\vv{x}}||" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{min}_\mathbf{\vv{x}}&space;||\mathbf{\vv{b}}&space;-&space;A\mathbf{\vv{x}}||" title="\mathrm{min}_\mathbf{\vv{x}} ||\mathbf{\vv{b}} - A\mathbf{\vv{x}}||" /></a>

using the the Jacobi iteration method. If *A* has full column rank, then there is a unique solution that satisfies the normal equations

<a href="https://www.codecogs.com/eqnedit.php?latex=\left(A^T&space;A&space;\right&space;)\mathbf{\vv{x}}&space;=&space;A^T&space;\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(A^T&space;A&space;\right&space;)\mathbf{\vv{x}}&space;=&space;A^T&space;\mathbf{\vv{b}}" title="\left(A^T A \right )\mathbf{\vv{x}} = A^T \mathbf{\vv{b}}" /></a>

The algorithm is as follows:

1. Form <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;B&space;=&space;A^T&space;A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;B&space;=&space;A^T&space;A" title="B = A^T A" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathbf{\vv{y}}&space;=&space;A^T\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathbf{\vv{y}}&space;=&space;A^T\mathbf{\vv{b}}" title="\mathbf{\vv{y}} = A^T\mathbf{\vv{b}}" /></a>.
2. Iteratively solve for **x** using Jacobi iteration on the system **Bx** = **y**

**Input:** 

*A* : REAL - a coefficient matrix of size *m* x *n*

*b* : REAL - the right hand side vector of size *m* x 1

*m* : INTEGER - number of rows in *A* and *b* (corresponds to the number of data points available)

*n* : INTEGER - number of columns in *A* (corresponds to the order of fit for the data)

*x0* : REAL - initial guess to be used in the [Jacobi iteration](./jacobi_solve.md)

*tol* : REAL - tolerance for convergence criteria

*maxiter* : INTEGER - maximum number of iterations until exit

*printit* : INTEGER - if equal to 1, then the final convergence error and number of iterations is printed

**Output:** 

*x0* : REAL - the final approximation to the solution to the least squares problem

**Usage/Example:**

This routine can be implemented in a program as follows (which follows example 6.1 in Ascher)

```fortran
INTEGER :: m, n, i, maxiter, printit
REAL*8, ALLOCATABLE :: A(:, :), b(:), x0(:), xstart(:)
REAL*8 :: tol, error

m = 5
n = 3
maxiter = 10000
tol = 10.D-16
printit = 1
ALLOCATE(A(1:m, 1:n), b(1:m), x0(1:n), xstart(n))
A = RESHAPE((/1.D0, 0.0D0, 1.0D0, &
& 2.D0, 3.0D0, 5.0D0, &
& 5.D0, 3.0D0, -2.0D0, &
& 3.D0, 5.0D0, 4.0D0, &
& -1.D0, 6.0D0, 3.0D0/), (/m, n/), ORDER=(/2, 1/))
! Ensures "diagonal dominance"
DO i = 1, n
	A(i, i) = A(i, i) + 50.D0
END DO
xstart = 1.D0
CALL mat_prod(A, xstart, m, n, 1, b)
x0 = 0.D0
CALL ls_solvejacobi(A, b, m, n, x0, tol, maxiter, printit)
WRITE(*,*) x0
CALL abs_err_vecl2(xstart, x0, n, error)
WRITE(*,*) error
```

The outputs from the above code:

```fortran
   3.8459253727671276E-016          26   ! Convergence error and number of iterations
   1.0000000000000000        1.0000000000000000        1.0000000000000000   !Final x0  
   0.0000000000000000  ! l2 norm of the absolute error between xstart and x0
```

**Implementation/Code:** The code for ls_solvejacobi can be seen below.

```fortran
SUBROUTINE ls_solvejacobi(A, b, m, n, x0, tol, maxiter, printit)
	IMPLICIT NONE
	
	! Takes as an input a coefficient matrix `A` of size `m` x `n`
	! (where m >> n) and a right-hand side matrix `b` of size `m` x 1.
	! Solves the least squares problem (min[x] ||b - Ax||) using the
	! normal equations.
	INTEGER, INTENT(IN) :: m, n, maxiter, printit
	REAL*8, INTENT(IN) :: A(1:m, 1:n), b(1:n), tol
	REAL*8, INTENT(OUT) :: x0(1:n)
	
	REAL*8 :: big_B(1:n, 1:n), y(1:n), z(1:n)
	INTEGER :: i
	
	! Form B = A^T * A
	CALL mat_prod(TRANSPOSE(A), A, n, m, n, big_B)
	! Form y = A^T * b
	CALL mat_prod(TRANSPOSE(A), b, n, m, 1, y)
	! Solve using the Jacobi Iteration method
	CALL jacobi_solve(big_B, n, y, tol, maxiter, x0, printit)
	
END SUBROUTINE
```