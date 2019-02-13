SUBROUTINE s_mult_vec(s, a, n)
	IMPLICIT NONE
	
	! Takes a vector a, its length, n, and a scalar, s as inputs. s is
	! multiplied into a and then a is returned.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(INOUT) :: a(1:n)
	REAL*8, INTENT(IN) :: s
	
	! An increment variable to loop over the vector a
	INTEGER :: i

	! Multiply a vector, a, by a scalar, s, element-wise.
	DO i = 1, n
		a(i) = s*a(i)
	END DO
END SUBROUTINE

