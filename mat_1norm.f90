SUBROUTINE mat_1norm(A, n, norm)
	IMPLICIT NONE
	
	! I/O variables
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: A(1:n, 1:n)
	REAL*8, INTENT(OUT) :: norm
	
	! Initialize subroutine variables
	INTEGER :: i
	REAL*8 :: col_1norm
	
	norm = 0.D0
	col_1norm = 0.D0
	
	! Loop through each column in the matrix A and find the l1 vector
	! norm for that column.
	DO i = 1, n
		CALL l1_vec_norm(A(:, i), n, col_1norm)
		! If a new maximum value for the column l1 norm has been found,
		! save it
		IF (col_1norm > norm) THEN
			norm = col_1norm
		ENDIF
	END DO
	
END SUBROUTINE
	
