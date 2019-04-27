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
