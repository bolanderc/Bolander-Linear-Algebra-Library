SUBROUTINE sym_dd_mat_gen(n, mat)
	IMPLICIT NONE
	
	! Takes as input the size of the square, symmetric matrix, `n` and
	! the matrix `mat` to be created.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(OUT) :: mat(n, n)
	
	! Initializes increment variables `i` and `j` as well as a `row_sum`
	! variable to add to the diagonal elements to make them dominant
	INTEGER :: i, j
	REAL*8 :: row_sum
	
	!~ Calculates a symmetric, square matrix by looping over the matrix
	DO i = 1, n
		!~ goes along all rows, i, and columns where j >= i
		DO j = i, n
			!~ A random number is assigned to the entries along the diagonal and to 
			!~ the right on each row (for example, a 3x3 matrix would assign random
			!~ values to all elements with an x in the following diagram
			!~ |x x x|
			!~ |- x x|
			!~ |- - x|
			CALL RANDOM_NUMBER(mat(i, j))
			!~ The random number assigned above is then copied across the diagonal.
			!~ The following diagram gives an example of how the values would match
			!~ |x 1 2|
			!~ |1 x 3|
			!~ |2 3 x|
			!~ The diagonal values are unique.
			mat(j, i) = mat(i, j)
		END DO
		row_sum = 0.D0
		! Loops through the columns and adds up all of the elements in
		! the ith row.
		DO j = 1, n
			row_sum = row_sum + mat(i, j)
		END DO
		! Adds the value of the row summation to the diagonal element.
		mat(i, i) = row_sum
	END DO
END SUBROUTINE
