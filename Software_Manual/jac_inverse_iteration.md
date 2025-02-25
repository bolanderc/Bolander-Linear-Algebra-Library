# Find Eigenvalues of a Square Linear System Using Inverse Iteration (Jacobi Iteration)

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          jac_inverse_iteration

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c jac_inverse_iteration.f90```

and can be added to a program using

```$ gfortran program.f90 jac_inverse_iteration.o ``` 

**Description/Purpose:** This subroutine uses the inverse iteration method with shifting to iteratively find any eigenvalue and the corresponding eigenvector given an appropriate shift (a shift of zero produces the smallest eigenvalue). The idea is that if the eigenvalues of *A* are <a href="https://www.codecogs.com/eqnedit.php?latex=\lambda_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda_j" title="\lambda_j" /></a>, then the eigenvalues of *A* - <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a>*I* are <a href="https://www.codecogs.com/eqnedit.php?latex=\lambda_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda_j" title="\lambda_j" /></a> - <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a> and the eigenvalues of *B* = (*A* - <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a>*I*)  are

<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_j&space;=&space;\frac{1}{\lambda_j&space;-&space;\alpha}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_j&space;=&space;\frac{1}{\lambda_j&space;-&space;\alpha}" title="\mu_j = \frac{1}{\lambda_j - \alpha}" /></a>

Meaning that a value can be found in the inverse problem that satisfies a given eigenvalue in the original problem. This routine differs from the [inverse_iteration](./inverse_iteration.md) routine because it uses Jacobi iteration to solve the inverse instead of LU decomposition. In general, this makes it less efficient than the other method, though if the initial guess for the vector is close to the actual vector it can be much more efficient.

**Input:** 

*A* : REAL - a coefficient matrix of size *n* x *n*

*n* : INTEGER - the rank of A and the length of the eigenvectors

*v0* : REAL - the initial guess of the eigenvetor

*alpha* : REAL - a shift given to the *A* matrix

*tol* : REAL - the exit tolerance for the algorithm

*maxiter* : INTEGER - the maximum number of iterations before exit

*printit* : INTEGER - flag for printing final convergence information

**Output:** 

*v* : REAL - the final approximation of the eigenvector corresponding to the largest eigenvalue

*lam0* : REAL - the final approximation of the largest eigenvalue

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
n = 3
ALLOCATE(A(n, n), v0(n), v(n))
lam = 0.D0
tol = 10.D-16
maxiter = 10000
printit = 1
v0 = 1.D0
alpha = 5.0D0
!~ 	CALL mat_dd(n, A)
A = RESHAPE((/3.D0, 0.D0, 0.D0, &
		 	& 0.D0, 6.D0, 0.D0, &
			& 0.D0, 0.D0, 2.D0/), (/n, n/), ORDER=(/2, 1/))
CALL inverse_iteration(A, n, v0, alpha, tol, maxiter, printit, v, lam)
WRITE(*,*) lam, v
v0 = 1.D0
CALL jac_inverse_iteration(A, n, v0, alpha, tol, maxiter, printit, v, lam)
WRITE(*,*) lam, v
```

The outputs from the above code:

```fortran
   0.0000000000000000               18
   6.0000000000000000        
   1.4551915228366852E-011   1.0000000000000000        2.5811747917131966E-009
   8.8817841970012523E-016          27
   6.0000000000000000       
   -7.4505805969238281E-009   1.0000000000000000       -1.3113726523970927E-013
```

Additionally, using this routine on a random 10 x 10 diagonally dominant matrix can be done with the following code:

```fortran
INTEGER :: n, maxiter, i, printit, j
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam, alpha

n = 20
ALLOCATE(A(n, n), v0(n), v(n))
lam = 0.D0
tol = 10.D-16
maxiter = 10000
printit = 1
v0 = 1.D0
alpha = 0.D0
DO i = 1, n
	DO j = 1, n
		A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
	END DO
END DO
CALL inverse_iteration(A, n, v0, tol, alpha, maxiter, printit, v, lam)
WRITE(*,*) lam, v
```

and the output is

```fortran
! Using inverse_iteration
   4.4408920985006262E-016          98
   3.1076841031614717      ! Eigenvalue    
   -1.1023089664893332E-003 -0.12139776077911797      -0.48109527394207219      -0.11306015810765833        3.4002745226522038E-002   4.4803443790737835E-002   5.4390219364371375E-003  -4.2721997152509295E-002   1.0000000000000000        4.4481088453212046E-002

! Using jac_inverse_iteration
   4.4408920985006262E-016          98
   3.1076828703134511       !Eigenvalue
   -1.8109729634054326E-003 -0.12141572470876785      -0.48074583030356316      -0.11301064663236836        3.3937948246448801E-002   4.4954366416069161E-002   5.5153767218869059E-003  -4.2705492596333501E-002   1.0000000000000000        4.4534131562955361E-002

```

Discrepancies between the two approaches stem from possible inaccuracies in the Jacobi iteration results.

**Implementation/Code:** The code for jac_inverse_iteration can be seen below.

```fortran
SUBROUTINE jac_inverse_iteration(A, n, v0, alpha, tol, maxiter, printit, v, lam0)
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
	REAL*8, INTENT(IN) :: A(n, n), v0(n), alpha, tol
	REAL*8, INTENT(OUT) :: v(n), lam0
	
	REAL*8 :: v1(n), lam1, error, norm, ALU(n, n), vguess(n)
	INTEGER :: i, k
	
	! Initializes variables
	v = 0.D0
	lam1 = 0.D0
	error = 10.D0*tol
	norm = 0.D0
	k = 0
	lam0 = 0.D0
	v1 = v0
	vguess = v1
	ALU = A
	
	! Shift the `A` matrix.
	DO i = 1, n
		ALU(i, i) = ALU(i, i) - alpha
	END DO
	
	
	! Iterate until the error or number of iterations reaches the given
	! limits
	DO WHILE (error > tol .AND. k < maxiter)
		! Solve for the new eigenvector using the Jacobi Iteration
		! The maximum iterations are fixed in an attempt to
		! limit the number of iterations that the jacobi uses so that
		! it is not as computationally inefficient as it could be.
		CALL jacobi_solve(ALU, n, v1, 10.D-16, 1000, vguess, 0)
		v = vguess
		! Normalize the newest guess
		CALL l2_vec_norm(v, n, norm)
		v1 = v/norm
		vguess = v1
		
		! Calculate the new guess for the eigenvalue
		CALL mat_prod(A, v1, n, n, 1, v)
		CALL vec_dot_prod(v1, v, n, lam1)
		
		! Calculate convergence error and increment values
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