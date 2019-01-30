SUBROUTINE l_inf_vec_norm(a, n, norm)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n)
REAL*8, INTENT(OUT) :: norm
INTEGER :: i

norm = 0.0D0

! The element in the vector with the maximum absolute value is found
! and represents the l_infinity norm.
DO i = 1, n
	IF (ABS(a(i)) > norm) THEN
		norm = ABS(a(i))
	ENDIF
END DO
END SUBROUTINE
