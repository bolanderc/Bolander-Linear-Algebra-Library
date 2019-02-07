SUBROUTINE vec_cross_prod3(a, b, crossprod)
	IMPLICIT NONE
	
	! I/O variable initiation
	REAL*8, INTENT(IN) :: a(1:3), b(1:3)
	REAL*8, INTENT(OUT) :: crossprod(1:3)
	
	! Calculate the cross product for two vectors of length 3.
	crossprod(1) = a(2)*b(3) - a(3)*b(2)
	crossprod(2) = -(a(1)*b(3) - a(3)*b(1))
	crossprod(3) = a(1)*b(2) - a(2)*b(1)
END SUBROUTINE
	
