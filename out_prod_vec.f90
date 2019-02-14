SUBROUTINE out_prod_vec(m, n, a, b, C)
	IMPLICIT NONE
	
	! Takes as inputs two vectors, `a` and `b`, of length `m` and `n`
	! respectively. Outputs the outer product of `a` and `b`, a matrix
	! `C` of size `m` x `n`.
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(IN) :: a(1:m), b(1:n)
	REAL*8, INTENT(OUT) :: C(1:m, 1:n)
	
	! Initialize two increment variables for looping.
	INTEGER :: i, j
	
	! Loop over all elements of `C` and multiply the appropriate
	! elements of `a` and `b` together to calculate the outer product.
	DO i = 1, m
		DO j = 1, n
			C(i, j) = a(i)*b(j)
		END DO
	END DO
END SUBROUTINE
