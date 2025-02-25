# Calculate the L2 Condition Number of an SPD Matrix (Rayleigh Quotient)

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          rayleigh_cond

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c rayleigh_cond.f90```

and can be added to a program using

```$ gfortran program.f90 rayleigh_cond.o ``` 

**Description/Purpose:** This subroutine uses the Rayleigh Quotient method to find the l2 condition number of a symmetric positive definite matrix. Since the nature of finding an eigenvector that corresponds to the largest and smallest eigenvalue is rather difficult, a brute force approach is applied. A large eigenvector is first run through the Rayleigh Quotient algorithm to generate (hopefully) the largest eigenvalue. Then, random vectors are run through for a set number of iterations to try and find all possible eigenvalues in the system.

**Input:** 

*A* : REAL - a coefficient matrix of size *n* x *n*

*n* : INTEGER - the rank of A and the length of the eigenvectors

*tol* : REAL - the exit tolerance for the algorithm

*maxiter* : INTEGER - the maximum number of iterations before exit

**Output:** 

*cond* : REAL - the final approximation of the condition number of the system

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, maxiter, i,j
REAL*8, ALLOCATABLE :: A(:, :)
REAL*8 :: tol, cond

n = 3

tol = 10.D-16
maxiter = 10000

ALLOCATE(A(n, n))
A = RESHAPE((/ 3.D0, 3.D0, 1.D0, &
			& 3.D0, 4.D0, 4.D0, &
			& 1.D0, 4.D0, 12.D0/), (/n, n/), ORDER=(/2, 1/))
CALL rayleigh_cond(A, n, tol, maxiter, cond)
WRITE(*,*) cond
```

The outputs from the above code:

```fortran
 Maximum Eigenvalue:   14.062770861175807     
 Minimum Eigenvalue:  0.11804443416433479     
   119.13116413094421 
```

Additionally, using this routine on an 8 x 8 Hilbert matrix can be done with the following code:

```fortran
INTEGER :: n, maxiter, i, printit, j
REAL*8, ALLOCATABLE :: A(:, :)
REAL*8 :: tol, cond

n = 8

tol = 10.D-16
maxiter = 10000
ALLOCATE(A(n, n))
DO i = 1, n
	DO j = 1, n
		A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
	END DO
END DO
CALL rayleigh_cond(A, n, tol, maxiter, cond)
WRITE(*,*) cond
```

and the output is

```fortran
 Maximum Eigenvalue:   1.6959389969219500     
 Minimum Eigenvalue:   1.1115388618738692E-010
   15257577176.050144
```

which matches precisely with the second eigenvalue presented in [this paper](<https://www.ams.org/journals/mcom/1967-21-099/S0025-5718-1967-0223075-0/S0025-5718-1967-0223075-0.pdf>).

**Implementation/Code:** The code for rayleigh_quotient can be seen below.

```fortran
SUBROUTINE rayleigh_cond(A, n, tol, maxiter, cond)
	IMPLICIT NONE
	
	! Uses the Rayleigh Quotient iterative method to approximate the
	! condition number of a given symmetric positive definite matrix.
	
	! Takes as an input the matrix `A` of rank `n` that contains the
	! system to be analyzed, a tolerance, `tol`, for exiting the
	! iterative solver as well as a maximum number of iterations,
	! `maxiter`. The output of the subroutine is the approximation of 
	! the condition number.
	INTEGER, INTENT(IN) :: n, maxiter
	REAL*8, INTENT(IN) :: A(n, n), tol
	REAL*8, INTENT(OUT) :: cond
	
	INTEGER :: i
	REAL*8 :: lam_max, lam_min, lam_i, v0(n), v_i(n)
	
	! Runs through the Rayleigh Quotient algorithm once with a large
	! eigenvector to try and find the maximum off the bat.
	v0 = 10.D0
	CALL rayleigh_quotient(A, n, v0, tol, maxiter, 0, v_i, lam_max)
	
	! Superior to randomly assigning these values, as we start with
	! eigenvalues that we know exist in this system.
	lam_min = lam_max
	
	! Brute force approach that checks many random eigenvectors to see
	! if we can find the maximum and minimum.
	DO i = 1, INT(100*n)
		CALL rand_mat(n, 1, v0)
		
		! Useful to find extremely small eigenvalues
		IF (i < INT(50*n)) THEN
			v0 = v0/10.D0
		END IF
		
		! Stores maximum and minimum values found.
		CALL rayleigh_quotient(A, n, v0, tol, maxiter, 0, v_i, lam_i)
		IF (lam_i > lam_max) THEN
			lam_max = lam_i
		ELSE IF (lam_i < lam_min) THEN
			lam_min = lam_i
		END IF
	END DO
	
	! Finds the condition number and reports the eigenvalues found.
	WRITE(*,*) "Maximum Eigenvalue:", lam_max
	WRITE(*,*) "Minimum Eigenvalue:", lam_min
	cond = lam_max/lam_min
END SUBROUTINE
```