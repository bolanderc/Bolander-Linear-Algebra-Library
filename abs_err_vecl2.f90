SUBROUTINE abs_err_vec(a, approx, n, error)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n), approx(1:n)
REAL*8, INTENT(OUT) :: error

REAL*8 :: diff(1:n)
INTEGER :: i
REAL*8 :: norm

DO i = 1, n
	diff(i) = approx(i) - a(i)
END DO

CALL l2_vec_norm(diff, n, error)

END SUBROUTINE
