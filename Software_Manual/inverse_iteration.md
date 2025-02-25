# Find Eigenvalues of a Square Linear System Using Inverse Iteration (LU Decomposition)

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          inverse_iteration

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c inverse_iteration.f90```

and can be added to a program using

```$ gfortran program.f90 inverse_iteration.o ``` 

**Description/Purpose:** This subroutine uses the inverse iteration method with shifting to iteratively find any eigenvalue and the corresponding eigenvector given an appropriate shift (a shift of zero produces the smallest eigenvalue). The idea is that if the eigenvalues of *A* are <a href="https://www.codecogs.com/eqnedit.php?latex=\lambda_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda_j" title="\lambda_j" /></a>, then the eigenvalues of *A* - <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a>*I* are <a href="https://www.codecogs.com/eqnedit.php?latex=\lambda_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda_j" title="\lambda_j" /></a> - <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a> and the eigenvalues of *B* = (*A* - <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a>*I*)  are

<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_j&space;=&space;\frac{1}{\lambda_j&space;-&space;\alpha}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_j&space;=&space;\frac{1}{\lambda_j&space;-&space;\alpha}" title="\mu_j = \frac{1}{\lambda_j - \alpha}" /></a>

Meaning that a value can be found in the inverse problem that satisfies a given eigenvalue in the original problem.

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
INTEGER :: n, maxiter, i, printit
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam, alpha

n = 3
ALLOCATE(A(n, n), v0(n), v(n))
lam = 0.D0
tol = 10.D-16
maxiter = 10000
printit = 1
v0 = 1.D0
alpha = 0.0D0
A = RESHAPE((/1.D0, 2.D0, 0.D0, &
			& 2.D0, 1.D0, 2.D0, &
            & 0.D0, 2.D0, 1.D0/), (/n, n/), ORDER=(/2, 1/))
CALL inverse_iteration(A, n, v0, alpha, tol, maxiter, printit, v, lam)
WRITE(*,*) lam, v
```

The outputs from the above code:

```fortran
   8.8817841970012523E-016          28  ! Convergence error at exit and the exit iterations
   -1.8284271247461898                   ! Approximation of largest eigenvalue
   0.70710678052474663   -1.0000000000000000   0.70710679885817773  ! Normalized eigenvector
```

Additionally, using this routine on an 20 x 20 Hilbert matrix can be done with the following code:

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
   8.4534686233368191E-020          31   ! Convergence error and number of iterations
   -1.8190197007123243E-018                    ! Eigenvalue
   ! Eigenvector
   7.1434857804492041E-009  -1.1494603201052072E-006   4.5632987553809310E-005  -7.7840548626109121E-004   7.0518472705676013E-003  -3.7471861751792644E-002  0.12163121197985835      -0.23992887217693432       0.27453697897244583      -0.20088876891016452       0.30726509218971626      -0.74604545697557068        1.0000000000000000      -0.62497068848836279        2.3580846900107453E-002  0.35574384095812017      -0.54925444346690810       0.52524084608927935      -0.27334781169413414        5.7591158369571113E-002 
```



**Implementation/Code:** The code for inverse_iteration can be seen below.

```fortran
SUBROUTINE inverse_iteration(A, n, v0, alpha, tol, maxiter, printit, v, lam0)
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
	
	REAL*8 :: v1(n), lam1, error, norm, ALU(n, n)
	INTEGER :: i, k
	
	! Initializes variables
	v = 0.D0
	lam1 = 0.D0
	error = 10.D0*tol
	norm = 0.D0
	k = 0
	lam0 = 0.D0
	v1 = v0
	ALU = A
	
	! Shift the `A` matrix.
	DO i = 1, n
		ALU(i, i) = ALU(i, i) - alpha
	END DO
	
	! Perform an LU factorization on the shifted matrix for efficiency.
	CALL lu_factor(ALU, n)
	
	! Iterate until the error or number of iterations reaches the given
	! limits
	DO WHILE (error > tol .AND. k < maxiter)
		! Solve the LU equation with the previous eigenvalue guess
		CALL lu_solve(ALU, n, v1, v)
		
		! Normalize the newest guess
		CALL l2_vec_norm(v, n, norm)
		v1 = v/norm
		
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