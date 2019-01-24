SUBROUTINE abs_err_n(x, y, e_abs)
IMPLICIT NONE
REAL*8, INTENT(IN) :: x, y
REAL*8, INTENT(OUT) :: e_abs

! Calculates the absolute error between an approximation, y, and
! the actual value, x. These are single, double precision values.
e_abs = ABS(y - x)

END SUBROUTINE

