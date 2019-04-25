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
