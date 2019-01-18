PROGRAM main
IMPLICIT NONE
!~ INTEGER :: r, c
!~ REAL*8, ALLOCATABLE :: mat(:, :)

!~ r = 5
!~ c = 5

!~ ALLOCATE(mat(r, c))
!~ CALL rand_mat(r, c, mat)
!~ WRITE(*, *) mat
INTEGER :: i
REAL*8 :: mach_eps
CALL dmaceps(mach_eps, i)
WRITE(*,*) mach_eps, i
END PROGRAM
