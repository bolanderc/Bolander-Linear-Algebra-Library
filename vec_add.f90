SUBROUTINE vec_add(a, b, c, n)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n), b(1:n)
REAL*8, INTENT(OUT) :: c(1:n)
INTEGER :: i

! Add two vectors together element-wise.
DO i = 1, n
	c(i) = a(i) + b(i)
END DO
END SUBROUTINE

