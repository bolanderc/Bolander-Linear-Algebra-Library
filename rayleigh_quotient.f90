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
