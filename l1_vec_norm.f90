SUBROUTINE l1_vec_norm(a, n, norm)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n)
REAL*8, INTENT(OUT) :: norm
INTEGER :: i

norm = 0.0D0

! Sum up the absolute value of each value in the vector a to find
! the l1-norm of the vector.
DO i = 1, n
	norm = norm + ABS(a(i))
END DO
END SUBROUTINE
