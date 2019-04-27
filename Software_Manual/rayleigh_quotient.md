# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          rayleigh_quotient

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c rayleigh_quotient.f90```

and can be added to a program using

```$ gfortran program.f90 rayleigh_quotient.o ``` 

**Description/Purpose:** This subroutine uses the Rayleigh Quotient method to iteratively find any eigenvalue and the corresponding eigenvector given an appropriate initial eigenvector. This acts similarly to the [inverse_iteration](./inverse_iteration.md) subroutine, but it features a dynamically changing alpha to help the program converge to the eigenvalue associated most closely with the initial eigenvector given. This means that a good initial guess of the associated eigenvector is necessary.

**Input:** 

*A* : REAL - a coefficient matrix of size *n* x *n*

*n* : INTEGER - the rank of A and the length of the eigenvectors

*v0* : REAL - the initial guess of the eigenvetor

*tol* : REAL - the exit tolerance for the algorithm

*maxiter* : INTEGER - the maximum number of iterations before exit

*printit* : INTEGER - flag for printing final convergence information

**Output:** 

*v* : REAL - the final approximation of the eigenvector corresponding to the eigenvalue

*lam0* : REAL - the final approximation of the eigenvalue corresponding to the initial eigenvector

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, maxiter, i, printit, j, k
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam, alpha, cond
INTEGER :: s(1:4)

n = 3

tol = 10.D-15
maxiter = 10000
printit = 1
ALLOCATE(A(n, n), v0(n), v(n))
A = RESHAPE((/ 3.D0, 3.D0, 1.D0, &
			 & 3.D0, 4.D0, 4.D0, &
			 & 1.D0, 4.D0, 12.D0/), (/n, n/), ORDER=(/2, 1/))
v0 = 1.D0
CALL rayleigh_quotient(A, n, v0, tol, maxiter, printit, v, lam)
WRITE(*,*) lam, v
```

The outputs from the above code:

```fortran
   5.3290705182007514E-015           5 ! Final convergence error and number of iterations
   14.062770861175805  	! Eigenvalue most closely corresponding to initial eigenvector    
   -0.21562071895029772      -0.46178753555637664       -1.0000000000000000 ! Eigenvector
```

Additionally, using this routine on an 8 x 8 Hilbert matrix can be done with the following code:

```fortran
INTEGER :: n, maxiter, i, printit, j
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam
INTEGER :: s(1:4)

n = 8

tol = 10.D-15
maxiter = 10000
printit = 1
ALLOCATE(A(n, n), v0(n), v(n))
DO i = 1, n
	DO j = 1, n
		A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
	END DO
END DO
v0 = 1.D0
CALL rayleigh_quotient(A, n, v0, tol, maxiter, printit, v, lam)
WRITE(*,*) lam
WRITE(*,*)  v 
```

and the output is

```fortran
   5.5511151231257827E-017           8
  0.29812521131693076   ! Eigenvalue  
  -1.0000000000000000       0.19964107668627296       0.45500608109550816       0.52037881437800704       0.52756595376147575       0.51397652248247738       0.49288911462148283       0.46961742309873455
```

which matches precisely with the second eigenvalue presented in [this paper](<https://www.ams.org/journals/mcom/1967-21-099/S0025-5718-1967-0223075-0/S0025-5718-1967-0223075-0.pdf>).

**Implementation/Code:** The code for rayleigh_quotient can be seen below.

```fortran
SUBROUTINE rayleigh_quotient(A, n, v0, tol, maxiter, printit, v, lam0)
	IMPLICIT NONE
	
	! Implements the inverse iteration method for finding any eigenvalue
	! in the system (provided an accurate guess of `alpha` is given).
	
	! Takes as an input the matrix `A` of rank `n` that contains the
	! system to be analyzed, an initial guess for the eigenvector, `v0`,
	! the shift, `alpha` that is to be used to find the eigenvalue,
	! a tolerance, `tol`, for exiting the iterative solver as well as a
	! maximum number of iterations, `maxiter`. Finally, a flag to print
	! the final number of iterations and convergence error, `printit` is
	! an input. The output of the subroutine is the final eigenvector,
	! `v1` produced from the algorithm (scaled by the element with the
	! maximum absolute value) as well as the final eigenvalue, `lam0`.
	! If `alpha` is zero, the minimum eigenvalue will be found.
	INTEGER, INTENT(IN) :: n, maxiter, printit
	REAL*8, INTENT(IN) :: A(n, n), v0(n), tol
	REAL*8, INTENT(OUT) :: v(n), lam0
	
	REAL*8 :: v1(n), lam1, error, norm, ALU(n, n)
	INTEGER :: i, k
	
	! Initializes variables
	v = 0.D0
	lam1 = 0.D0
	error = 10.D0*tol
	norm = 0.D0
	k = 0
	CALL mat_prod(A, v0, n, n, 1, v)
	CALL vec_dot_prod(v0, v, n, lam0)
	v1 = v0
	
	! Iterate until the error or number of iterations reaches the given
	! limits
	DO WHILE (error > tol .AND. k < maxiter)
		ALU = A
		
		! Shift the `A` matrix by the eigenvalue of the previous iteration.
		DO i = 1, n
			ALU(i, i) = ALU(i, i) - lam0
		END DO
		
		! Perform an LU factorization on the shifted matrix.
		CALL lu_factor(ALU, n)
		! Solve the LU equation with the previous eigenvalue guess.
		CALL lu_solve(ALU, n, v1, v)
		! Normalize the newest guess.
		CALL l2_vec_norm(v, n, norm)
		v1 = v/norm
		
		! Calculate the new guess for the eigenvalue.
		CALL mat_prod(A, v1, n, n, 1, v)
		CALL vec_dot_prod(v1, v, n, lam1)
		
		! Calculate convergence error and increment values.
		error = DABS(lam1 - lam0)
		lam0 = lam1
		k = k + 1
	END DO
	
	! Normalizes the eigenvector according to the maximum absolute value
	! of the elements.
	v = v1/MAXVAL(DABS(v1))
	
	! Prints the error and number of iterations when exiting.
	IF (printit == 1) THEN
		WRITE(*,*) error, k
	END IF
END SUBROUTINE
```