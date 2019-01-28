SUBROUTINE s_mult_vec(s, a, c, n)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n)
REAL*8, INTENT(IN) :: s
REAL*8, INTENT(OUT) :: c(1:n)
INTEGER :: i

! Multiply a vector, a, by a scalar, s, element-wise.
DO i = 1, n
	c(i) = s*a(i)
END DO
END SUBROUTINE

