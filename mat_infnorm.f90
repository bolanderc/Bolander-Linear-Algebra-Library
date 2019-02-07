SUBROUTINE mat_infnorm(A, n, norm)
	IMPLICIT NONE
	
	! I/O variables
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(IN) :: A(1:n, 1:n)
	REAL*8, INTENT(OUT) :: norm
	
	! Initialize subroutine variables
	INTEGER :: i
	REAL*8 :: row_1norm
	
	norm = 0.D0
	row_1norm = 0.D0
	
	! Loop through each row of the matrix A and calculate the l1 norm
	! treating that row as a vector.
	DO i = 1, n
		CALL l1_vec_norm(A(i, :), n, row_1norm)
		! If a new maximum value is found from the row l1 norm, then
		! store it.
		IF (row_1norm > norm) THEN
			norm = row_1norm
		ENDIF
	END DO
END SUBROUTINE
