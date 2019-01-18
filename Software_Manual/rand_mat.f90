SUBROUTINE rand_mat(r, c, mat)
	IMPLICIT NONE
	INTEGER, INTENT(in) :: r, c
	REAL*8, INTENT(out) ::  mat(r, c)
	! Fill a matrix with `r` rows and `c` columns with a random
	! number between 0 and 1 (not including 1).
	CALL RANDOM_NUMBER(mat)
END SUBROUTINE
