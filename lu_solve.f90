SUBROUTINE lu_solve(A, n, b, x)
	IMPLICIT NONE
	
	! Takes as inputs a square coefficient matrix, `A` of size
	! `n` x `n` and the right-hand side vector `b` to solve the system
	! of equations Ax=b for `x` using LU decomposition with forward and
	! backward substitution.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: b(1:n)
	REAL*8, INTENT(INOUT) :: A(1:n, 1:n), x(1:n)
	
	INTEGER :: i, j, k
	REAL*8 :: factor, fsum, backsum
	REAL*8 :: y(1:n)
	
	! Find the LU decomposition of the coefficient matrix `A`.
	DO k = 1, n - 1
		
		DO i = k + 1, n
			
			factor = A(i, k)/A(k, k)
			
			DO j = k , n
				A(i, j) = A(i, j) - factor*A(k, j)
			END DO
			
			A(i, k) = factor
		END DO
	END DO
	
	! Implement the forward substitution algorithm on Ly = b to find
	! `y`.
	
	y(1) = b(1)/A(1, 1)
	
	DO k = 2, n
		fsum = 0.D0
		DO j = 1, k-1
			fsum = fsum + A(k, j)*y(j)
		END DO
		y(k) = (b(k) - fsum)
	END DO
	
	! Implement the backward substitution algorithm on Ux = y to find
	! `x`.
	
	x(n) = y(n)/A(n, n)
	
	DO k = n-1, 1, -1
		backsum = 0.D0
		DO j = k+1, n
			backsum = backsum + A(k, j)*x(j)
		END DO
		x(k) = (y(k) - backsum)/A(k, k)
	END DO
	
END SUBROUTINE
