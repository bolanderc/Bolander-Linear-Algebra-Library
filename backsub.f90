SUBROUTINE backsub(n, A, b, x)
	IMPLICIT NONE
	
	! Takes as an input the size of the square matrix `n`, the upper
	! triangular matrix `A`, and the right hand side vector `b`. Outputs
	! the solution vector `x`.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: A(1:n, 1:n), b(1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	! Initialize decrement variables and a variable to compute the sum
	! of previous solutions integrated into the algorithm.
	INTEGER :: k, j
	REAL*8 :: backsum
	
	! Calculate the last value in the solution vector `x`.
	x(n) = b(n)/A(n, n)
	
	! Loop through the remaining rows in `x` to calculate the solution
	! using the backward substitution algorithm.
	DO k = n-1, 1, -1
		backsum = 0.D0
		DO j = k+1, n
			backsum = backsum + A(k, j)*x(j)
		END DO
		x(k) = (b(k) - backsum)/A(k, k)
	END DO
END SUBROUTINE
