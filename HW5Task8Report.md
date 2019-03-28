# Homework 5 Task 8: QR factorization Using Gram-Schmidt on Hilbert Matrices

*Task: Try out your QR-factorization method from the previous task on the Hilbert matrices of sizes n=4,6,8,10. Test to see if the orthogonal matrix is really orthogonal by multiplying QtQ and comparing the result to the identity matrix. Explain the results you obtain.*

This process can be executed in Fotran with the following code.

```fortran
INTEGER :: m, n, i, j
REAL*8, ALLOCATABLE :: A(:, :), Q(:, :), R(:, :), QTQ(:, :)

m = 8
n = 8
ALLOCATE(A(1:m, 1:n), Q(1:m, 1:n), R(1:n, 1:n), QTQ(1:n, 1:n))
DO i = 1, m
	DO j = 1, n
		A(i, j) = 1.0D0/(REAL(i) + REAL(j) - 1.0D0)
	END DO
END DO
CALL qr_factor_modgs(A, m, n, Q, R)
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
```



The output of running the Hilbert matrices can be seen [here](./HilbertQR.txt). As can be seen, the QTQ section of each of the successive matrices becomes less and less similar to the identity matrix as the size of the matrix increases. This essentially means that Q is becoming less and less orthogonal due to round off error. If this same process were to be implemented using the classical Gram-Schmidt implementation, this would likely be even worse. This is where the Gram-Schmidt algorithm becomes less viable, as it is less stable for these matrices with small values in them where round off error is more likely to be significant.