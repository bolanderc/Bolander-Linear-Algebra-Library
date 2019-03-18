SUBROUTINE lu_solve(LU, n, b, x)
	IMPLICIT NONE
	
	! Takes as inputs a square, LU-decomposed, coefficient matrix, `LU`
	! of size `n` x `n` and the right-hand side vector `b` to solve the
	! system of equations Ax=b for `x` using LU decomposition with
	! forward and backward substitution.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: b(1:n)
	REAL*8, INTENT(INOUT) :: LU(1:n, 1:n), x(1:n)
	
	INTEGER :: i, j, k
	REAL*8 :: fsum, backsum
	REAL*8 :: y(1:n)
	
	! Implement the forward substitution algorithm on Ly = b to find
	! `y`.
	
	y(1) = b(1)/LU(1, 1)
	
	DO k = 2, n
		fsum = 0.D0
		DO j = 1, k-1
			fsum = fsum + LU(k, j)*y(j)
		END DO
		y(k) = (b(k) - fsum)
	END DO
	
	! Implement the backward substitution algorithm on Ux = y to find
	! `x`.
	
	x(n) = y(n)/LU(n, n)
	
	DO k = n-1, 1, -1
		backsum = 0.D0
		DO j = k+1, n
			backsum = backsum + LU(k, j)*x(j)
		END DO
		x(k) = (y(k) - backsum)/LU(k, k)
	END DO
	
END SUBROUTINE
