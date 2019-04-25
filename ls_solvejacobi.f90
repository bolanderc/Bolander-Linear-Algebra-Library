SUBROUTINE ls_solvejacobi(A, b, m, n, x0, tol, maxiter, printit)
	IMPLICIT NONE
	
	! Takes as an input a coefficient matrix `A` of size `m` x `n`
	! (where m >> n) and a right-hand side matrix `b` of size `m` x 1.
	! Solves the least squares problem (min[x] ||b - Ax||) using the
	! normal equations.
	INTEGER, INTENT(IN) :: m, n, maxiter, printit
	REAL*8, INTENT(IN) :: A(1:m, 1:n), b(1:n), tol
	REAL*8, INTENT(OUT) :: x0(1:n)
	
	REAL*8 :: big_B(1:n, 1:n), y(1:n), z(1:n)
	INTEGER :: i
	
	! Form B = A^T * A
	CALL mat_prod(TRANSPOSE(A), A, n, m, n, big_B)
	! Form y = A^T * b
	CALL mat_prod(TRANSPOSE(A), b, n, m, 1, y)
	! Solve using the Jacobi Iteration method
	CALL jacobi_solve(big_B, n, y, tol, maxiter, x0, printit)
	
END SUBROUTINE
	
