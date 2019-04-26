# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          eigen_search

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c eigen_search.f90```

and can be added to a program using

```$ gfortran program.f90 eigen_search.o ``` 

**Description/Purpose:** This subroutine performs searches between the largest and smallest eigenvalue (found using the [power method](./power_method.md) and [inverse_iteration](./inverse_iteration.md) respectively) to find additional eigenvalues. Once the interval is established, the inverse iteration method is used with variable alpha values to find additional eigenvalues between the minimum and maximum. The number of alpha locations given in each search are k = 1, 2, 4, 6, 8, ... etc.

**Input:** 

*A* : REAL - a coefficient matrix of size *n* x *n*

*n* : INTEGER - the rank of A and the length of the eigenvectors

*v0* : REAL - the initial guess of the eigenvetor

*tol* : REAL - the exit tolerance for the algorithm

*maxiter* : INTEGER - the maximum number of iterations before exit

*n_search* : INTEGER - the number of searches to perform

**Output:** 

There are no outputs from this code

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
INTEGER :: n, maxiter, i, printit, j, k
REAL*8, ALLOCATABLE :: A(:, :), v0(:), v(:)
REAL*8 :: tol, lam, alpha, cond

n = 6
tol = 10.D-18
maxiter = 10000
printit = 1
! Sets up a matrix with easily identifiable eigenvalues
ALLOCATE(A(n, n), v0(n))
lam = 5.D0
A(1, 1) = 1.D0
DO i = 2, n
	A(i, i) = lam*DBLE(i)
END DO
v0 = 1.D0
CALL eigen_search(A, n, v0, tol, maxiter, 5)
```

The outputs from the above code:

```fortran
 Minimum Eigenvalue:   1.0000000000000000     
 Maximum Eigenvalue:   29.999999999999993     
 Search           1
    Alpha                     Eigenvalue
   15.499999999999996        15.000000000000000     
 Search           2
    Alpha                     Eigenvalue
   10.666666666666664        10.000000000000000     
   20.333333333333329        20.000000000000000     
 Search           3
    Alpha                     Eigenvalue
   6.7999999999999989        10.000000000000000     
   12.599999999999998        14.999999999999973     
   18.399999999999999        20.000000000000000     
   24.199999999999996        25.000000000000000     
 Search           4
    Alpha                     Eigenvalue
   5.1428571428571415        10.000000000000000     
   9.2857142857142829        10.000000000000000     
   13.428571428571423        14.999999999999995     
   17.571428571428566        19.999999999999932     
   21.714285714285708        20.000000000000007     
   25.857142857142851        25.000000000000000     
 Search           5
    Alpha                     Eigenvalue
   4.2222222222222214        10.000000000000000     
   7.4444444444444429        10.000000000000000     
   10.666666666666664        10.000000000000000     
   13.888888888888886        15.000000000000000     
   17.111111111111107        15.000000000000000     
   20.333333333333329        20.000000000000000     
   23.555555555555550        25.000000000000000     
   26.777777777777771        25.000000000000000 
```

**Implementation/Code:** The code for eigen_search can be seen below.

```fortran
SUBROUTINE eigen_search(A, n, v0, tol, maxiter, n_search)
	IMPLICIT NONE
	
	! Searches for the maximum and minimum eigenvalues in a system and
	! then continues subdiving the space to search for additional
	! intermediate eigenvalues.
	
	! Takes as an input the matrix `A` of rank `n` that contains the
	! system to be analyzed, an initial guess for the eigenvector, `v0`,
	! a tolerance, `tol`, for exiting the iterative solver as well as a
	! maximum number of iterations, `maxiter`. Finally, the number of
	! searches to perform is given as `n_search`. There are no outputs
	! to the code, everything is written to the terminal.
	INTEGER, INTENT(IN) :: n, maxiter, n_search
	REAL*8, INTENT(IN) :: A(n, n), v0(n), tol
	
	REAL*8 :: lam_1, lam_n, lam_i, v(n), alpha, del_lam
	INTEGER :: i, j, k, delk(n_search)
	
	! Finds the largest eigenvalue using the power method.
	CALL power_method(A, n, v0, tol, maxiter, 0, v, lam_1)
	
	! Finds the minimum eigenvalues using the inverse iteration method.
	alpha = 0.D0
	CALL inverse_iteration(A, n, v0, alpha, tol, maxiter, 0, v, lam_n)
	
	WRITE(*,*) "Minimum Eigenvalue:", lam_n
	WRITE(*,*) "Maximum Eigenvalue:", lam_1
	
	! Sets up the manner of subdividing the interval from `lam_n` to
	! `lam_1`. The inverse iteration routine is used with a variable
	! shift to search in the interval at specific locations. The number
	! of locations (shifts) searched starts at 1, then 2 locations, and
	! then the number of locations increases by 2 for each additional
	! search. This avoids the repetition of having an alpha at the center
	! of the interval.
	k = 1
	delk = 1
	delk(2:) = 2
	
	! Executes each search group with 1 point, 2 points, 4 points, 6
	! points, etc.
	DO i = 1, n_search
		WRITE(*,*) "Search", i
		del_lam = (lam_1 - lam_n)/DBLE(k + 1)
		alpha = lam_n
		WRITE(*,*) "   Alpha                     Eigenvalue"
		
		! Does the inverse iteration for each alpha value in the current
		! search.
		DO j = 1, k
			alpha = alpha + del_lam
			CALL inverse_iteration(A, n, v0, alpha, tol, maxiter, 0, v, lam_i)
			WRITE(*,*) alpha, lam_i
		END DO
		
		! Increment k (the number of search points for the given search)
		k = k +  delk(i)
	END DO
END SUBROUTINE
```