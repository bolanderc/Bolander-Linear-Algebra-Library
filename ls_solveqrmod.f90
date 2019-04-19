SUBROUTINE ls_solveqrmod(A, m, n, rhs, Q, R, x, factor)
	IMPLICIT NONE
	
	! Takes as an input a coefficient matrix `A` of size `m` x `n` and
	! its corresponding right-hand side vector `rhs` of length `m`. If
	! `factor` is set to a value of 1, then the `Q` and `R` factors of
	! `A` are calculated and then used to solve the least squares
	! problem.
	INTEGER, INTENT(IN) :: m, n, factor
	REAL*8, INTENT(IN) :: A(m, n), rhs(m)
	REAL*8, INTENT(INOUT) :: Q(m, n), R(n, n)
	REAL*8, INTENT(OUT) :: x(n)
	
	REAL*8 :: c(n)
	
	! If factor is set to 1 then the QR factorization is performed.
	! Otherwise it is assumed that the `Q` and `R` factors are passed
	! in as arguments to the function.
	IF (factor == 1) THEN
		! Calculates the QR factors of A using the modified Gram-Schmidt
		! orthogonalization method.
		CALL qr_factor_modgs(A, m, n, Q, R)
	END IF
	
	! Solves the system Q^T*b = c
	CALL mat_prod(TRANSPOSE(Q), rhs, n, m, 1, c)
			
	! Uses backward substitution to solve the system R*x = c
	CALL backsub(n, R, c, x)
	
END SUBROUTINE
