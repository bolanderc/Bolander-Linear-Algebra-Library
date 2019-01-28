SUBROUTINE l2_vec_norm(a, n, norm)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n)
REAL*8, INTENT(OUT) :: norm
REAL*8 :: summation
INTEGER :: i

summation = 0.0D0

! Sum up the squares of each value in the vector a.
DO i = 1, n
	summation = summation + (a(i)*a(i))
END DO

! Find the l2 norm of the vector by taking the square root of the sum
! of the squares
norm = SQRT(summation)
END SUBROUTINE
