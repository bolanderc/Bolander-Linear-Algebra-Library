# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          ls_solveqrmod

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c ls_solveqrmod.f90```

and can be added to a program using

```$ gfortran program.f90 ls_solveqrmod.o ``` 

**Description/Purpose:** This subroutine solves a least-squares system of equations using the [modified Gram-Schmidt QR factorization](./qr_factor_modgs.md).

The algorithm is as follows (given in Ascher, U., and C. Greif. "A First Course in Numerical Methods. SIAM." *Society for* (2011).)

1. Decompose

   ​     (a) <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;A&space;=&space;QR" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;A&space;=&space;QR" title="A = QR" /></a>

2. Compute

   ​     (a) <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\vv{\mathbf{c}}&space;=&space;Q^T\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\vv{\mathbf{c}}&space;=&space;Q^T\vv{\mathbf{b}}" title="\vv{\mathbf{c}} = Q^T\vv{\mathbf{b}}" /></a>

3. Solve the upper triangular system <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;R\vv{\mathbf{x}}&space;=&space;\vv{\mathbf{c}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;R\vv{\mathbf{x}}&space;=&space;\vv{\mathbf{c}}" title="R\vv{\mathbf{x}} = \vv{\mathbf{c}}" /></a>

**Input:** 

*A* : REAL - a coefficient matrix of size *m* x *n*

*m* : INTEGER - number of rows in *A* and *Q* and the length of *rhs*

*n* : INTEGER - number of columns in *A*, *Q*, and the size of *R*

*rhs* : REAL - the right-hand side vector of size *m* x 1

*Q* : REAL (OPTIONAL) - the orthonormal column matrix of the *QR* factorization (size *m* x *n*)

*R* : REAL(OPTIONAL) - the upper triangular matrix of the *QR* factorization (size *n* x *n*)

*factor* : INTEGER - tells the subroutine whether to call the [qr_factor_modgs](./qr_factor_modgs.md) subroutine to decompose *A* into its *QR* factors

**Output:** 

*x* : REAL - the solution to the least-squares linear system of equations

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: m, n, i, j, factor, method
REAL*8, ALLOCATABLE :: A(:, :), Q(:, :), R(:, :), x(:), rhs(:)

m = 5
n = 3
factor = 1
method = 0
ALLOCATE(A(1:m, 1:n), Q(1:m, 1:n), R(1:n, 1:n), x(1:n), rhs(1:m))
A = RESHAPE((/1.D0, 0.0D0, 1.0D0, &
			& 2.D0, 3.0D0, 5.0D0, &
			& 5.D0, 3.0D0, -2.0D0, &
			& 3.D0, 5.0D0, 4.0D0, &
			& -1.D0, 6.0D0, 3.0D0/), (/m, n/), ORDER=(/2, 1/))
rhs = (/4.D0, -2.D0, 5.D0, -2.D0, 1.D0/)
CALL ls_solveqrmod(A, m, n, rhs, Q, R, x, factor)
WRITE(*,*) x
```

or, if Q and R have already been calculated,

```fortran
INTEGER :: m, n, i, j, factor, method
REAL*8, ALLOCATABLE :: A(:, :), Q(:, :), R(:, :), x(:), rhs(:)

m = 5
n = 3
factor = 0
method = 0
ALLOCATE(A(1:m, 1:n), Q(1:m, 1:n), R(1:n, 1:n), x(1:n), rhs(1:m))
A = RESHAPE((/1.D0, 0.0D0, 1.0D0, &
			& 2.D0, 3.0D0, 5.0D0, &
			& 5.D0, 3.0D0, -2.0D0, &
			& 3.D0, 5.0D0, 4.0D0, &
			& -1.D0, 6.0D0, 3.0D0/), (/m, n/), ORDER=(/2, 1/))
rhs = (/4.D0, -2.D0, 5.D0, -2.D0, 1.D0/)
CALL qr_factor_modgs(A, m, n, Q, R)
CALL ls_solveqrmod(A, m, n, rhs, Q, R, x, factor)
WRITE(*,*) x
```



The outputs from the above code:

```fortran
  0.34722617354196295       0.39900426742531991      -0.78591749644381215  
```

This matches the results of Example 6.1 in Ascher.

**Implementation/Code:** The code for ls_solveqrmod can be seen below.

```fortran
SUBROUTINE ls_solveqrmod(A, m, n, rhs, Q, R, x, factor)
	IMPLICIT NONE
	
	! Takes as an input a coefficient matrix `A` of size `m` x `n` and
	! its corresponding right-hand side vector `rhs` of length `m`. If
	! `factor` is set to a value of 1, then the `Q` and `R` factors of
	! `A` are calculated and then used to solve the least squares
	! problem.
	INTEGER, INTENT(IN) :: m, n, factor
	REAL*8, INTENT(IN) :: A(m, n), rhs(m)
	REAL*8, INTENT(INOUT) :: Q(m, n), R(n, n)
	REAL*8, INTENT(OUT) :: x(n)
	
	REAL*8 :: c(n)
	
	! If factor is set to 1 then the QR factorization is performed.
	! Otherwise it is assumed that the `Q` and `R` factors are passed
	! in as arguments to the function.
	IF (factor == 1) THEN
		! Calculates the QR factors of A using the modified Gram-Schmidt
		! orthogonalization method.
		CALL qr_factor_modgs(A, m, n, Q, R)
	END IF
	
	! Solves the system Q^T*b = c
	CALL mat_prod(TRANSPOSE(Q), rhs, n, m, 1, c)
			
	! Uses backward substitution to solve the system R*x = c
	CALL backsub(n, R, c, x)
	
END SUBROUTINE
```