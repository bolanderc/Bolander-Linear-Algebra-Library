SUBROUTINE abs_err_vecl1(a, approx, n, error)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(1:n), approx(1:n)
REAL*8, INTENT(OUT) :: error
INTEGER :: i
REAL*8 :: diff(1:n)

!~ Find the difference between *a* and its approximation for each
!~ element.
DO i = 1, n
	diff(i) = approx(i) - a(i)
END DO

!~ Calculate the l1 vector norm using the difference between the vectors
!~ to find the error.
CALL l1_vec_norm(diff, n, error)

END SUBROUTINE
