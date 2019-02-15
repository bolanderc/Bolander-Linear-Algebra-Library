SUBROUTINE forwardsub(n, A, b, x)
	IMPLICIT NONE
	
	! Takes as input the size of the square matrix, `n`, the lower
	! triangular matrix, `A`, and the right-hand-side vector `b`.
	! Outputs the solution vector `x`.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: A(1:n, 1:n), b(1:n)
	REAL*8, INTENT(OUT) :: x(1:n)
	
	! Initialize decrement variables and a variable to compute the sum
	! of previous solutions integrated into the algorithm.
	INTEGER k, j
	REAL*8 fsum
	
	! Calculate the first value in the solution vector `x`.
	x(1) = b(1)/A(1, 1)
	
	! Loop through the remaining rows in `x` to calculate the solution
	! using the forward substitution algorithm.
	DO k = 2, n
		fsum = 0.D0
		DO j = 1, k-1
			fsum = fsum + A(k, j)*x(j)
		END DO
		x(k) = (b(k) - fsum)/A(k, k)
	END DO
END SUBROUTINE
