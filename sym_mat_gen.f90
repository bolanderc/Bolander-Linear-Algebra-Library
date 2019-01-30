SUBROUTINE sym_mat_gen(n, mat)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(OUT) :: mat(n, n)
INTEGER :: i, j

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
END DO

END SUBROUTINE
