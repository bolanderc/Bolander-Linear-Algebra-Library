# Find the Largest Eigenvalue (and Eigenvector) Using the Power Method

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          power_method

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c power_method.f90```

and can be added to a program using

```$ gfortran program.f90 power_method.o ``` 

**Description/Purpose:** This subroutine uses the power method to iteratively find the largest eigenvalue and the corresponding eigenvector. The algorithm for finding the eigenvector is as follows:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vv{v}}_k&space;=&space;\gamma_k\lambda_1^k\sum_{j=1}^n\beta_j(\lambda_j/\lambda_1)^k\mathbf{\vv{x}}_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vv{v}}_k&space;=&space;\gamma_k\lambda_1^k\sum_{j=1}^n\beta_j(\lambda_j/\lambda_1)^k\mathbf{\vv{x}}_j" title="\mathbf{\vv{v}}_k = \gamma_k\lambda_1^k\sum_{j=1}^n\beta_j(\lambda_j/\lambda_1)^k\mathbf{\vv{x}}_j" /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma_k" title="\gamma_k" /></a> is some normalization coefficient that guarantees that the norm of the eigenvector is one. Essentially this shows that the eigenvector can be seen as a linear combination of the eigenvectors **x**<sub>j</sub> and that the eigenvalue term in the summation will go to zero for every j > 1. Knowing this, the eigenvalue can be found by

<a href="https://www.codecogs.com/eqnedit.php?latex=\lambda_k&space;=&space;\mathbf{\vv{v}}_k^TA\mathbf{\vv{v}}_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda_k&space;=&space;\mathbf{\vv{v}}_k^TA\mathbf{\vv{v}}_k" title="\lambda_k = \mathbf{\vv{v}}_k^TA\mathbf{\vv{v}}_k" /></a>

via a Rayleigh quotient.

**Input:** 

*A* : REAL - a coefficient matrix of size *n* x *n*

*n* : INTEGER - the rank of A and the length of the eigenvectors

*v0* : REAL - the initial guess of the eigenvetor

*tol* : REAL - the exit tolerance for the algorithm

*maxiter* : INTEGER - the maximum number of iterations before exit

*printit* : INTEGER - flag for printing final convergence information

**Output:** 

*v1* : REAL - the final approximation of the eigenvector corresponding to the largest eigenvalue

*lam0* : REAL - the final approximation of the largest eigenvalue

**Usage/Example:**

This routine can be implemented in a program as follows (which follows Example 4 [here](<http://ergodic.ugr.es/cphys/LECCIONES/FORTRAN/power_method.pdf>))

```fortran
INTEGER :: n, maxiter, i, printit
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam

n = 3
ALLOCATE(A(n, n), v0(n), v(n))
lam = 0.D0
tol = 10.D-16
maxiter = 10000
printit = 1
v0 = 1.D0
A = RESHAPE((/1.D0, 2.D0, 0.D0, &
			& -2.D0, 1.D0, 2.D0, &
			& 1.D0, 3.D0, 1.D0/), (/n, n/), ORDER=(/2, 1/))
CALL power_method(A, n, v0, tol, maxiter, printit, v, lam)
WRITE(*,*) lam, v
```

The outputs from the above code:

```fortran
   0.0000000000000000               34  ! Convergence error at exit and the exit iterations
   2.9999999999999996                   ! Approximation of largest eigenvalue
   0.49999999999999994  0.49999999999999994  0.99999999999999989  ! Normalized eigenvector
```

Additionally, using this routine on an 20 x 20 Hilbert matrix can be done with the following code:

```fortran
INTEGER :: n, maxiter, i, printit, j
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam

n = 20
ALLOCATE(A(n, n), v0(n), v(n))
lam = 0.D0
tol = 10.D-16
maxiter = 10000
printit = 1
v0 = 1.D0
DO i = 1, n
	DO j = 1, n
		A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
	END DO
END DO
CALL power_method(A, n, v0, tol, maxiter, printit, v, lam)
WRITE(*,*) lam, v
```

and the output is

```fortran
   0.0000000000000000               15   ! Convergence error and number of iterations
   1.9071347204072533                    ! Eigenvalue
   ! Eigenvector
   1.0000000000000000       0.63153893180834586       0.48170552470917699       0.39577939404225121       0.33864052058721722       0.29732839459896360       0.26579806044172338       0.24080108268105982       0.22041627505509773       0.20342569216373199       0.18901536311290157       0.17661823144515848       0.16582577118835690       0.15633539873115412       0.14791772253568922       0.14039535584329460       0.13362876033883500       0.12750652172661470       0.12193850695620556       0.11685094644504723 
```



**Implementation/Code:** The code for power_method can be seen below.

```fortran
SUBROUTINE power_method(A, n, v0, tol, maxiter, printit, v1, lam0)
	IMPLICIT NONE
	
	! Implements the power method for finding the largest eigenvalue in
	! a system.
	
	! Takes as an input the matrix `A` of rank `n` that contains the
	! system to be analyzed, an initial guess for the eigenvector, `v0`,
	! a tolerance, `tol`, for exiting the iterative solver as well as a
	! maximum number of iterations, `maxiter`. Finally, a flag to print
	! the final number of iterations and convergence error, `printit` is
	! an input. The output of the subroutine is the final eigenvector,
	! `v1` produced from the algorithm (scaled by the element with the
	! maximum absolute value) as well as the final eigenvalue, `lam0`. 
	INTEGER, INTENT(IN) :: n, maxiter, printit
	REAL*8, INTENT(IN) :: A(n, n), v0(n), tol
	REAL*8, INTENT(OUT) :: v1(n), lam0
	
	REAL*8 :: vt(n), lam1, error, norm
	INTEGER :: k
	
	! Initializes variables
	vt = 0.D0
	lam1 = 0.D0
	error = 10.D0*tol
	norm = 0.D0
	k = 0
	lam0 = 0.D0
	v1 = 0.D0
	
	! Find the first approximation of the eigenvector
	CALL mat_prod(A, v0, n, n, 1, vt)
	
	! Iterate until the error or number of iterations reaches the given
	! limits
	DO WHILE (error > tol .AND. k < maxiter)
		
		! Normalize the eigenvector approximation to prevent overflow.
		CALL l2_vec_norm(vt, n, norm)
		v1 = vt*(1.D0/norm)
		
		! Calculate the next approximation for the eigenvector
		CALL mat_prod(A, v1, n, n, 1, vt)
		
		! Find the corresponding eigenvalue
		CALL vec_dot_prod(v1, vt, n, lam1)
		
		! Calculate the error and increment values
		error = ABS(lam1 - lam0)
		lam0 = lam1
		k = k + 1
	END DO
	
	! Normalizes the eigenvector according to the maximum absolute value
	! of the elements.
	CALL l_inf_vec_norm(v1, n, norm)
	v1 = v1*(1.D0/norm)
	
	! Prints the error and number of iterations when exiting.
	IF (printit == 1) THEN
		WRITE(*,*) error, k
	END IF
END SUBROUTINE
```