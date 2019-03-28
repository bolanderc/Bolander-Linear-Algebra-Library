SUBROUTINE solve_normal_equations(A, b, m, n, x)
	IMPLICIT NONE
	
	! Takes as an input a coefficient matrix `A` of size `m` x `n`
	! (where m >> n) and a right-hand side matrix `b` of size `m` x 1.
	! Solves the least squares problem (min[x] ||b - Ax||) using the
	! normal equations.
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(IN) :: A(1:m, 1:n), b(1:n)
	REAL*8, INTENT(INOUT) :: x(1:n)
	
	REAL*8 :: big_B(1:n, 1:n), y(1:n), z(1:n)
	INTEGER :: error, i
	
	! Form B = A^T * A
	CALL mat_prod(TRANSPOSE(A), A, n, m, n, big_B)
	! Form y = A^T * b
	CALL mat_prod(TRANSPOSE(A), b, n, m, 1, y)
	! Compute the Cholesky Factor, G, where B = G * G^T
	CALL cholesky_factor(big_B, n, error)
	! Solve G * z = y for z using forward substitution
	CALL forwardsub(n, big_B, y, z)
	! Solve G^T * x = z for x using backward substitution
	CALL backsub(n, TRANSPOSE(big_B), z, x)
	
END SUBROUTINE
	
