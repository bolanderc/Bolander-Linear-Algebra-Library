# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          qr_sq_solve

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c qr_sq_solve.f90```

and can be added to a program using

```$ gfortran program.f90 qr_sq_solve.o ``` 

**Description/Purpose:** This subroutine solves a square, linear system of equations using the [classical Gram-Schmidt QR factorization](./qr_factor_gs.md).

The algorithm is as follows (given in Ascher, U., and C. Greif. "A First Course in Numerical Methods. SIAM." *Society for* (2011).)

1. Decompose

   ​     (a) <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;A&space;=&space;QR" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;A&space;=&space;QR" title="A = QR" /></a>

2. Compute

   ​     (a) <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\vv{\mathbf{c}}&space;=&space;Q^T\vv{\mathbf{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\vv{\mathbf{c}}&space;=&space;Q^T\vv{\mathbf{b}}" title="\vv{\mathbf{c}} = Q^T\vv{\mathbf{b}}" /></a>

3. Solve the upper triangular system <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;R\vv{\mathbf{x}}&space;=&space;\vv{\mathbf{c}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;R\vv{\mathbf{x}}&space;=&space;\vv{\mathbf{c}}" title="R\vv{\mathbf{x}} = \vv{\mathbf{c}}" /></a>

**Input:** 

*A* : REAL - a coefficient matrix of size *n* x *n*

*b* : REAL - the right hand side vector of size *n* x 1

*Q* : REAL (OPTIONAL) - the orthonormal column matrix of the *QR* factorization

*R* : REAL(OPTIONAL) - the upper triangular matrix of the *QR* factorization

*n* : INTEGER - number of rows and columns in *A*, *Q*, and *R*

*factor* : INTEGER - tells the subroutine whether to call the [qr_factor_gs](./qr_factor_gs.md) subroutine to decompose *A* into its *QR* factors

**Output:** 

*x* : REAL - the solution to the square, linear system of equations

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: m, n, i, j, factor
REAL*8, ALLOCATABLE :: A(:, :), Q(:, :), R(:, :), x(:), b(:)

n = 3
factor = 1
ALLOCATE(A(1:n, 1:n), Q(1:n, 1:n), R(1:n, 1:n), x(1:n), b(1:n))
A = RESHAPE((/1.D0, 2.0D0, 3.0D0, &
			& 4.D0, 8.0D0, 2.0D0, &
			& 6.D0, 4.0D0, 6.0D0/), (/n, n/), ORDER=(/2, 1/))
b = (/4.D0, 12.D0, 3.D0/)
CALL qr_sq_solve(A, n, b, Q, R, x, factor)
WRITE(*,*) x
```

or, if Q and R have already been calculated,

```fortran
INTEGER :: m, n, i, j, factor
REAL*8, ALLOCATABLE :: A(:, :), Q(:, :), R(:, :), x(:), b(:)

n = 3
factor = 0
ALLOCATE(A(1:n, 1:n), Q(1:n, 1:n), R(1:n, 1:n), x(1:n), b(1:n))
A = RESHAPE((/1.D0, 2.0D0, 3.0D0, &
			& 4.D0, 8.0D0, 2.0D0, &
			& 6.D0, 4.0D0, 6.0D0/), (/n, n/), ORDER=(/2, 1/))
b = (/4.D0, 12.D0, 3.D0/)
CALL qr_factor_gs(A, n, n, Q, R)
CALL qr_sq_solve(A, n, b, Q, R, x, factor)
WRITE(*,*) x
```



The outputs from the above code:

```fortran
  -1.2500000000000018        2.0250000000000004       0.40000000000000108 
```

This matches the results of the same matrices input into Numpy's [numpy.linalg.lstsq](<https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.linalg.lstsq.html#numpy.linalg.lstsq>) routine:

```python
A = np.array([[1, 2, 3],
              [4, 8, 2],
              [6, 4, 6]])

b = np.array([4, 12, 3])

x = np.linalg.lstsq(A, b)

print(x[0])
[-1.25   2.025  0.4  ]
```



**Implementation/Code:** The code for qr_sq_solve can be seen below.

```fortran
SUBROUTINE qr_sq_solve(A, n, b, Q, R, x, factor)
	IMPLICIT NONE
	
	! Takes as an input a square matrix `A` of size `n` and the 
	! corresponding right-hand side vector `b` of length `n`.
	! Optionally, the QR factors `Q` and `R` can be input if they are
	! already known (if not, the `factor` input can be set equal to
	! 1 to factor `A`). The solution to the system of equations, `x` is
	! output as a result of this subroutine.
	INTEGER, INTENT(IN) :: n, factor
	REAL*8, INTENT(IN) :: A(1:n, 1:n), b(1:n)
	REAL*8, INTENT(INOUT) :: Q(1:n, 1:n), R(1:n, 1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	! Create a temporary vector `c` that will be used in the algorithm.
	REAL*8 :: c(1:n)
	
	! Checks if the `Q` and `R` factors need to be found
	IF (factor == 1) THEN
		! If the factors needs to be found, performs a classical Gram-
		! Schmidt orthogonalization algorithm to find them.
		CALL qr_factor_gs(A, n, n, Q, R)
	ENDIF
	
	! Computes `c` = `Q`^T `b`.
	CALL mat_prod(TRANSPOSE(Q), b, n, n, 1, c)
	
	! Uses backward substitution to compute `R``x` = `c` and find `x`.
	CALL backsub(n, R, c, x)

END SUBROUTINE
```