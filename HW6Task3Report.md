# Homework 6 Task 3: QR factorization Using Householder Transformations on Hilbert Matrices

*Task: Create a third version of the QR-factorization algorithm using Householder Transformations. As usual, document you code in your software manual. Use the third incarnation of the code on the same matrices as the previous two QR-factorization and compare/explain your results.*

This process can be executed in Fortran with the following code.

```fortran
INTEGER :: m, n, i, j, k
REAL*8, ALLOCATABLE :: A(:, :), Q(:, :), R(:, :), QTQ(:, :)
INTEGER :: s(1:4)

s = (/4, 6, 8, 10/)
DO k = 1, 4
	IF(ALLOCATED(A)) THEN
		DEALLOCATE(A, Q, R, QTQ)
	END IF
	WRITE(*,*) "Hilbert", s(k), "x", s(k)
	m = s(k)
	n = s(k)
	ALLOCATE(A(1:m, 1:n), Q(1:m, 1:n), R(1:n, 1:n), QTQ(1:n, 1:n))
	DO i = 1, m
		DO j = 1, n
			A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
		END DO
	END DO
	CALL qr_factor_hh(A, m, n, Q, R)
	WRITE(*,*) "A"
	DO i = 1, m
		WRITE(*,*) A(i, :), "/"
	END DO
	WRITE(*,*) "Q"
	DO i = 1, m
		WRITE(*,*) Q(i, :), "/"
	END DO
	WRITE(*,*) "R"
	DO i = 1, n
		WRITE(*,*) R(i, :), "/"
	END DO
	CALL mat_prod(TRANSPOSE(Q), Q, n, m, n, QTQ)
	WRITE(*,*) "QTQ"
	DO i = 1, n
		WRITE(*,*) QTQ(i, :), "/"
	END DO
	CALL mat_prod(Q, R, m, n, n, A)
	WRITE(*,*) "A = QR"
	DO i = 1, m
		WRITE(*,*) A(i, :), "/"
	END DO
END DO
```



The output of running the Hilbert matrices can be seen [here](./HilbertHouseholderQR.txt). As can be seen, the QTQ section of each of the successive matrices maintain themselves very nearly equal to the identity matrix. When compared to the modified Gram-Schmidt orthogonalization method, as executed in the following code,

```fortran
INTEGER :: m, n, i, j, k
REAL*8, ALLOCATABLE :: A(:, :), Q(:, :), R(:, :), QTQ(:, :)
REAL*8, ALLOCATABLE :: Am(:, :), Qm(:, :), Rm(:, :), QTQm(:, :)
INTEGER :: s(1:4)

s = (/4, 6, 8, 10/)
DO k = 1, SIZE(s)
	IF(ALLOCATED(A)) THEN
		DEALLOCATE(A, Q, R, QTQ, Qm, Rm, QTQm)
	END IF
	WRITE(*,*) "Hilbert", s(k), "x", s(k)
	m = s(k)
	n = s(k)
	ALLOCATE(A(1:m, 1:n), Q(1:m, 1:n), R(1:n, 1:n), QTQ(1:n, 1:n), Qm(1:m, 1:n), Rm(1:n, 1:n), QTQm(1:n, 1:n))
	DO i = 1, m
		DO j = 1, n
			A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
		END DO
	END DO
	CALL qr_factor_modgs(A, m, n, Q, R)
	CALL qr_factor_hh(A, m, n, Qm, Rm)
	CALL mat_prod(TRANSPOSE(Q), Q, n, m, n, QTQ)
	CALL mat_prod(TRANSPOSE(Qm), Qm, n, m, n, QTQm)
	WRITE(*,*) "QTQ Mod_GS - QTQ Householder"
	DO i = 1, n
		WRITE(*,*) QTQ(i, :) - QTQm(i, :), "/"
	END DO
END DO
```

the modified algorithm performs much better, as seen [here](./HilbertModHouseDiff.txt). This indicates a much more stable algorithm, as there is a large difference between the QTQ matrices of the two methods. The software manual entry for the Householder transformation routine can be seen [here](./Software_Manual/qr_factor_hh.md).