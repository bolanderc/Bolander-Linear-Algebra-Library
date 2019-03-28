# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          solve_normal_equations

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c solve_normal_equations.f90```

and can be added to a program using

```$ gfortran program.f90 solve_normal_equations.o ``` 

**Description/Purpose:** This subroutine solves the least squares problem, given by

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{min}_\mathbf{\vv{x}}&space;||\mathbf{\vv{b}}&space;-&space;A\mathbf{\vv{x}}||" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{min}_\mathbf{\vv{x}}&space;||\mathbf{\vv{b}}&space;-&space;A\mathbf{\vv{x}}||" title="\mathrm{min}_\mathbf{\vv{x}} ||\mathbf{\vv{b}} - A\mathbf{\vv{x}}||" /></a>

using the normal equations. If *A* has full column rank, then there is a unique solution that satisfies the normal equations

<a href="https://www.codecogs.com/eqnedit.php?latex=\left(A^T&space;A&space;\right&space;)\mathbf{\vv{x}}&space;=&space;A^T&space;\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(A^T&space;A&space;\right&space;)\mathbf{\vv{x}}&space;=&space;A^T&space;\mathbf{\vv{b}}" title="\left(A^T A \right )\mathbf{\vv{x}} = A^T \mathbf{\vv{b}}" /></a>

The algorithm is as follows (given in Ascher, U., and C. Greif. "A First Course in Numerical Methods. SIAM." *Society for* (2011).)

1. Form <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;B&space;=&space;A^T&space;A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;B&space;=&space;A^T&space;A" title="B = A^T A" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathbf{\vv{y}}&space;=&space;A^T\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathbf{\vv{y}}&space;=&space;A^T\mathbf{\vv{b}}" title="\mathbf{\vv{y}} = A^T\mathbf{\vv{b}}" /></a>.
2. Compute the Cholesky Factor satisfying <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;B&space;=&space;GG^T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;B&space;=&space;GG^T" title="B = GG^T" /></a>.
3. Solve <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;G\mathbf{\vv{z}}&space;=&space;\mathbf{\vv{y}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;G\mathbf{\vv{z}}&space;=&space;\mathbf{\vv{y}}" title="G\mathbf{\vv{z}} = \mathbf{\vv{y}}" /></a> for ***z.***
4. Solve <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;G^T\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{z}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;G^T\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{z}}" title="G^T\mathbf{\vv{x}} = \mathbf{\vv{z}}" /></a> for ***x***.

**Input:** 

*A* : REAL - a coefficient matrix of size *m* x *n*

*b* : REAL - the right hand side vector of size *m* x 1

*m* : INTEGER - number of rows in *A* and *b* (corresponds to the number of data points available)

*n* : INTEGER - number of columns in *A* (corresponds to the order of fit for the data)

**Output:** 

*x* : REAL - the solution to the least squares problem

**Usage/Example:**

This routine can be implemented in a program as follows (which follows example 6.1 in Ascher)

```fortran
INTEGER :: m, n, i
REAL*8, ALLOCATABLE :: A(:, :), b(:), x(:)

m = 5
n = 3
ALLOCATE(A(1:m, 1:n), b(1:m), x(1:n))
A = RESHAPE((/1.D0, 0.0D0, 1.0D0, &
			& 2.D0, 3.0D0, 5.0D0, &
			& 5.D0, 3.0D0, -2.0D0, &
			& 3.D0, 5.0D0, 4.0D0, &
			& -1.D0, 6.0D0, 3.0D0/), (/m, n/), ORDER=(/2, 1/))
b = (/4.D0, -2.D0, 5.D0, -2.D0, 1.D0/)
CALL solve_normal_equations(A, b, m, n, x)
DO i = 1, n
	WRITE(*,*) x(i)
END DO
```

The outputs from the above code:

```fortran
  0.34722617354196295     
  0.39900426742532002     
 -0.78591749644381215
```

**Implementation/Code:** The code for solve_normal_equations can be seen below.

```fortran
SUBROUTINE solve_normal_equations(A, b, m, n, x)
	IMPLICIT NONE
	
	! Takes as an input a coefficient matrix `A` of size `m` x `n`
	! (where m >> n) and a right-hand side matrix `b` of size `m` x 1.
	! Solves the least squares problem (min[x] ||b - Ax||) using the
	! normal equations.
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(IN) :: A(1:m, 1:n), b(1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	REAL*8 :: big_B(1:n, 1:n), y(1:n), z(1:n)
	INTEGER :: error, i
	
	! Form B = A^T * A
	CALL mat_prod(TRANSPOSE(A), A, n, m, n, big_B)
	! Form y = A^T * b
	CALL mat_prod(TRANSPOSE(A), b, n, m, 1, y)
	! Compute the Cholesky Factor, G, where B = G * G^T
	CALL cholesky_factor(big_B, n, error)
	! Solve G * z = y for z using forward substitution
	CALL forwardsub(n, big_B, y, z)
	! Solve G^T * x = z for x using backward substitution
	CALL backsub(n, TRANSPOSE(big_B), z, x)
	
END SUBROUTINE
```



