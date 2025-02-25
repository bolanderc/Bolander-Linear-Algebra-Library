# Solve the Least-Squares Problem Using Jacobi Iteration

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

As the Jacobi algorithm needs a diagonally dominant system, there we can create this system by adding a parameter to the diagonal of the *A<sup>T</sup>A* system such that

<a href="https://www.codecogs.com/eqnedit.php?latex=\left(A^TA&space;&plus;&space;\alpha{I}&space;\right&space;)\mathbf{\vv{x}}&space;=&space;A^T\mathbf{\vv{b}}&space;&plus;&space;\alpha{I}\mathbf{\vv{x}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(A^TA&space;&plus;&space;\alpha{I}&space;\right&space;)\mathbf{\vv{x}}&space;=&space;A^T\mathbf{\vv{b}}&space;&plus;&space;\alpha{I}\mathbf{\vv{x}}" title="\left(A^TA + \alpha{I} \right )\mathbf{\vv{x}} = A^T\mathbf{\vv{b}} + \alpha{I}\mathbf{\vv{x}}" /></a>

To set up the iterative process of the Jacobi method, we can set

<a href="https://www.codecogs.com/eqnedit.php?latex=\left(D_{A^TA}&space;&plus;&space;\alpha{I}&space;\right&space;)\mathbf{\vv{x}}_{k&plus;1}&space;=&space;A^T\mathbf{\vv{b}}&space;&plus;&space;\alpha{I}\mathbf{\vv{x}}_k&space;-&space;\left(L&space;&plus;&space;U&space;\right&space;)_{A^TA}\mathbf{\vv{x}}_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(D_{A^TA}&space;&plus;&space;\alpha{I}&space;\right&space;)\mathbf{\vv{x}}_{k&plus;1}&space;=&space;A^T\mathbf{\vv{b}}&space;&plus;&space;\alpha{I}\mathbf{\vv{x}}_k&space;-&space;\left(L&space;&plus;&space;U&space;\right&space;)_{A^TA}\mathbf{\vv{x}}_k" title="\left(D_{A^TA} + \alpha{I} \right )\mathbf{\vv{x}}_{k+1} = A^T\mathbf{\vv{b}} + \alpha{I}\mathbf{\vv{x}}_k - \left(L + U \right )_{A^TA}\mathbf{\vv{x}}_k" /></a>

Finally, rearranging we find

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vv{x}}_{k&plus;1}&space;=&space;\left[A^T\mathbf{\vv{b}}&space;&plus;&space;\alpha{I}\mathbf{\vv{x}}_k&space;-&space;\left(L&space;&plus;&space;U&space;\right&space;)_{A^TA}\mathbf{\vv{x}}_k\right]/\left(D_{A^TA}&space;&plus;&space;\alpha{I}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vv{x}}_{k&plus;1}&space;=&space;\left[A^T\mathbf{\vv{b}}&space;&plus;&space;\alpha{I}\mathbf{\vv{x}}_k&space;-&space;\left(L&space;&plus;&space;U&space;\right&space;)_{A^TA}\mathbf{\vv{x}}_k\right]/\left(D_{A^TA}&space;&plus;&space;\alpha{I}&space;\right&space;)" title="\mathbf{\vv{x}}_{k+1} = \left[A^T\mathbf{\vv{b}} + \alpha{I}\mathbf{\vv{x}}_k - \left(L + U \right )_{A^TA}\mathbf{\vv{x}}_k\right]/\left(D_{A^TA} + \alpha{I} \right )" /></a>

The algorithm is as follows:

1. Form <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;B&space;=&space;A^T&space;A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;B&space;=&space;A^T&space;A" title="B = A^T A" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathbf{\vv{y}}&space;=&space;A^T\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathbf{\vv{y}}&space;=&space;A^T\mathbf{\vv{b}}" title="\mathbf{\vv{y}} = A^T\mathbf{\vv{b}}" /></a>.
2. <a href="https://www.codecogs.com/eqnedit.php?latex=x_i^{(k&space;&plus;&space;1)}&space;=&space;\frac{1}{B_{ii}&space;&plus;&space;\alpha}\left[y_i&space;-&space;\sum\limits_{\substack{j=1&space;\\j\neq{i}}}^{n}B_{ij}x_j^{(k)}&space;&plus;&space;\alpha{x_i}^{(k)}\right&space;],\qquad{i&space;=&space;1,...,n.}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_i^{(k&space;&plus;&space;1)}&space;=&space;\frac{1}{B_{ii}&space;&plus;&space;\alpha}\left[y_i&space;-&space;\sum\limits_{\substack{j=1&space;\\j\neq{i}}}^{n}B_{ij}x_j^{(k)}&space;&plus;&space;\alpha{x_i}^{(k)}\right&space;],\qquad{i&space;=&space;1,...,n.}" title="x_i^{(k + 1)} = \frac{1}{B_{ii} + \alpha}\left[y_i - \sum\limits_{\substack{j=1 \\j\neq{i}}}^{n}B_{ij}x_j^{(k)} + \alpha{x_i}^{(k)}\right ],\qquad{i = 1,...,n.}" /></a>

**Input:** 

*A* : REAL - a coefficient matrix of size *m* x *n*

*b* : REAL - the right hand side vector of size *m* x 1

*m* : INTEGER - number of rows in *A* and *b* (corresponds to the number of data points available)

*n* : INTEGER - number of columns in *A* (corresponds to the order of fit for the data)

*x0* : REAL - initial guess to be used in the least squares Jacobi iteration

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
   8.0059320849734419E-016         169   ! Convergence error and number of iterations
   1.0000000000000013       0.99999999999999767        1.0000000000000024   !Final x0  
   3.6299369564111875E-015  ! l2 norm of the absolute error between xstart and x0
```

**Implementation/Code:** The code for ls_solvejacobi can be seen below.

```fortran
SUBROUTINE ls_solvejacobi(A, b, m, n, x0, tol, maxiter, printit)
	IMPLICIT NONE
	
	! Takes as an input a coefficient matrix `A` of size `m` x `n`
	! (where m >> n) and a right-hand side matrix `b` of size `m` x 1.
	! Solves the least squares problem (min[x] ||b - Ax||) using a
	! modified Jacobi iteration method.
	
	! *NOTE* This algorithm was developed in conjunction with Jackson
	! T. Reid.
	INTEGER, INTENT(IN) :: m, n, maxiter, printit
	REAL*8, INTENT(IN) :: A(m, n), b(n), tol
	REAL*8, INTENT(OUT) :: x0(n)
	
	REAL*8 :: big_B(n, n), y(n), error, x1(n), sum_ax, alpha(n)
	INTEGER :: i, j, k
	
	x1 = 0.D0
	sum_ax = 0.D0
	error = 10.D0*tol
	k = 0
	alpha = 0.D0
	
	! Determine the alpha for each row by summing all of the row values
	! in `A`. Alpha is the parameter that ensures diagonal dominance.
	DO i = 1, n
		DO j = 1, n
			alpha(i) = alpha(i) + 2.D0*A(i, j)
		END DO
	END DO
	
	! Form B = A^T * A
	CALL mat_prod(TRANSPOSE(A), A, n, m, n, big_B)
	! Form y = A^T * b
	CALL mat_prod(TRANSPOSE(A), b, n, m, 1, y)
	! Solve using the modified Jacobi Iteration method.
	! First, add Ialpha to the B matrix
	DO i = 1, n
		big_B(i, i) = big_B(i, i) + alpha(i)
	END DO
	
	! Iteration loop for the algorithm. Iterates until the error is
	! less than the given tolerance or the maximum number of iterations
	! is exceeded.
	DO WHILE (error > tol .AND. k < maxiter)
		
		! Implements the modified Jacobi algorithm.
		DO i = 1, n
			! Starts with y = A^Tb.
			sum_ax = y(i)
			! Subtracts L_ATA*x_k
			DO j = 1, i - 1
				sum_ax = sum_ax - big_B(i, j)*x0(j)
			END DO
			! Adds alpha*I*x_k
			sum_ax = sum_ax + alpha(i)*x0(i)
			! Subtracts U_ATA*x_k
			DO j = i + 1, n
				sum_ax = sum_ax - big_B(i, j)*x0(j)
			END DO
			! Divides by the diagonal of A^T*A + alpha*I
			x1(i) = sum_ax/big_B(i, i)
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