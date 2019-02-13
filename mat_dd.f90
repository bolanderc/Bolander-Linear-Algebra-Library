SUBROUTINE mat_dd(n, mat)
	IMPLICIT NONE
	
	! Takes in a value for the number of rows and columns in the
	! diagonally dominant matrix and outputs the matrix.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(OUT) :: mat(1:n, 1:n)
	
	! Intialize do loop incrementers and a row summation variable to
	! help enforce diagonal dominance.
	INTEGER :: i, j
	REAL*8 :: row_sum
	
	! Fills the matrix with random numbers in all elements from 0 to 1
	! non-inclusive.
	CALL RANDOM_NUMBER(mat)
	
	! Loops through all rows and calculates a value for the summation
	! of all elements in that row.
	DO i = 1, n
		row_sum = 0.D0
		DO j = 1, n
			row_sum = row_sum + mat(i, j)
		END DO
		! Adds the value of the row summation to the diagonal element.
		mat(i, i) = row_sum
	END DO
END SUBROUTINE
