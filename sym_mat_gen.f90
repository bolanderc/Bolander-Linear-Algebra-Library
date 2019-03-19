SUBROUTINE sym_mat_gen(n, mat)
	IMPLICIT NONE
	
	! Takes as an input the size of the matrix `n` and outputs a
	! symmetric, square matrix of that size.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(OUT) :: mat(n, n)
	INTEGER :: i, j

	!~ Fills `mat` with random values
	CALL RANDOM_NUMBER(mat)

	!~ goes along all rows, i, and columns where j >= i
	!~ The random number assigned to the upper triangular elements is
	!~ then copied across the diagonal.
	!~ The following diagram gives an example of how the values would match
	!~ |x 1 2|
	!~ |1 x 3|
	!~ |2 3 x|
	!~ The diagonal values are unique.
	DO i = 1, n
		DO j = i, n
			mat(j, i) = mat(i, j)
		END DO
	END DO

END SUBROUTINE
