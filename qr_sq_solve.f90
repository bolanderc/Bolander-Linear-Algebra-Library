SUBROUTINE qr_sq_solve(A, n, b, Q, R, x, factor)
	IMPLICIT NONE
	
	! Takes as an input a square matrix `A` of size `n` and the 
	! corresponding right-hand side vector `b` of length `n`.
	! Optionally, the QR factors `Q` and `R` can be input if they are
	! already known (if not, the `factor` input can be set equal to
	! 1 to factor `A`). The solution to the system of equations, `x` is
	! output as a result of this subroutine.
	INTEGER, INTENT(IN) :: n, factor
	REAL*8, INTENT(IN) :: A(1:n, 1:n), b(1:n)
	REAL*8, INTENT(INOUT) :: Q(1:n, 1:n), R(1:n, 1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	! Create a temporary vector `c` that will be used in the algorithm.
	REAL*8 :: c(1:n)
	
	! Checks if the `Q` and `R` factors need to be found
	IF (factor == 1) THEN
		! If the factors needs to be found, performs a classical Gram-
		! Schmidt orthogonalization algorithm to find them.
		CALL qr_factor_gs(A, n, n, Q, R)
	ENDIF
	
	! Computes `c` = `Q`^T `b`.
	CALL mat_prod(TRANSPOSE(Q), b, n, n, 1, c)
	
	! Uses backward substitution to compute `R``x` = `c` and find `x`.
	CALL backsub(n, R, c, x)

END SUBROUTINE
