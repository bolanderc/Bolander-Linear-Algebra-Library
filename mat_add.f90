SUBROUTINE mat_add(n, m, A, B, C)
	IMPLICIT NONE
	
	! Takes as inputs: two arrays, A and B of size n x m to be added
	! together. Outputs an array C that is the sum of A and B and has
	! size n x m. n represents the number of rows and m the number of
	! columns.
	INTEGER, INTENT(IN) :: n, m
	REAL*8, INTENT(IN) :: A(1:n, 1:m), B(1:n, 1:m)
	REAL*8, INTENT(OUT) :: C(1:n, 1:m)
	
	! Initializes increment variables i and j to loop over the matrices.
	INTEGER :: i, j
	
	! Loops over the matrix C and fills its elements with the sum of the
	! corresponding elements in A and B.
	DO i = 1, n
		DO j = 1, m
			C(i, j) = A(i, j) + B(i, j)
		END DO
	END DO
END SUBROUTINE
