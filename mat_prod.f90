SUBROUTINE mat_prod(A, B, n, m, p, C)
	IMPLICIT NONE
	
	! Calculates the matrix product of A and B, which is C. A is of size
	! n x m and B is of size m x p, therefore C is of size n x p.
	INTEGER, INTENT(IN) :: n, m, p
	REAL*8, INTENT(IN) :: A(1:n, 1:m), B(1:m, 1:p)
	REAL*8, INTENT(OUT) :: C(1:n, 1:p)
	
	! Initialize i, j, and k as increments in do loops for the matrix
	! multiplication and also assign dotprodr_c to hold the value for
	! the dot product to be used in the calculation of the C element.
	INTEGER :: i, j, k
	REAL*8 :: dotprodr_c
	
	! Loop through the rows (i) and columns (j) of C and calculate the
	! value from the multiplication that belongs in that element.
	DO i = 1, n
		DO j = 1, p
			dotprodr_c = 0.D0
	
			! Calculate the dot product of the ith row of A and the jth
			! column of B and assign that value to C_i, j
			DO k = 1, m
				dotprodr_c = dotprodr_c + A(i, k)*B(k, j)
			END DO
			C(i, j) = dotprodr_c
		END DO
	END DO
END SUBROUTINE
