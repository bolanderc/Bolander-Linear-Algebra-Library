SUBROUTINE ls_solveqr(A, m, n, rhs, Q, R, x, factor)
	IMPLICIT NONE
	
	! Takes as an input a coefficient matrix `A` of size `m` x `n` and
	! its corresponding right-hand side vector `rhs` of length `m`. If
	! `factor` is set to a value of 1, then the `Q` and `R` factors of
	! `A` are calculated and then used to solve the least squares
	! problem.
	INTEGER, INTENT(IN) :: m, n, factor
	REAL*8, INTENT(IN) :: A(m, n), rhs(m)
	REAL*8, INTENT(INOUT) :: Q(m, m), R(m, n)
	REAL*8, INTENT(OUT) :: x(n)
	
	REAL*8 :: c(n), c_d(m), temp(n, n), b(m)
	
	! If factor is set to 1 then the QR factorization is performed.
	! Otherwise it is assumed that the `Q` and `R` factors are passed
	! in as arguments to the function.
	IF (factor == 1) THEN
		! Finds the QR factors of A using the Householder Transformation
		! subroutine.
		b = rhs
		CALL qr_factor_hh(A, m, n, Q, R)
	ENDIF
	! Solves the system Q^T*b = c_d
	CALL mat_prod(TRANSPOSE(Q), rhs, m, m, 1, c_d)
	! Takes the first n values in c_d as c and takes the first n
	! rows of R.
	c = c_d(:n)
	temp = R(:n, :)
	
	! Uses backward substitution to solve the system R*x = c
	CALL backsub(n, temp, c, x)
	
END SUBROUTINE
