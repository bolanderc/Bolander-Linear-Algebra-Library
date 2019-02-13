SUBROUTINE s_mult_mat(s, n, m, A)
	IMPLICIT NONE
	
	! Takes as inputs a matrix A of size n x m and a scalar s. Outputs
	! the product of s and A as a modified version of A.
	INTEGER, INTENT(IN) :: n, m
	REAL*8, INTENT(IN) :: s
	REAL*8, INTENT(INOUT) :: A(1:n, 1:m)
	
	! Initializes two increment variables to loop over A.
	INTEGER :: i, j
	
	! Loops over every value in A and multiplies that value by s.
	DO i = 1, n
		DO j = 1, m
			A(i, j) = s*A(i, j)
		END DO
	END DO
END SUBROUTINE
