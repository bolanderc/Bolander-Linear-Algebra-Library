SUBROUTINE rand_mat(r, c, mat)
	IMPLICIT NONE
	INTEGER, INTENT(in) :: r, c
	REAL*8, INTENT(out) ::  mat(r, c)
	CALL RANDOM_NUMBER(mat)
END SUBROUTINE
