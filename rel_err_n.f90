SUBROUTINE rel_err_n(x, y, e_rel)
IMPLICIT NONE
REAL*8, INTENT(IN) :: x, y
REAL*8, INTENT(OUT) :: e_rel

! Calculates the relative error between two numbers, x and y, where
! y is the approximation to x.
e_rel = ABS(x - y)/ABS(x)

END SUBROUTINE
