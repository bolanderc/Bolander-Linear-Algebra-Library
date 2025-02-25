# Find the Orthogonal QR Factorization of a Matrix (Householder)

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          qr_factor_hh

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c qr_factor_hh.f90```

and can be added to a program using

```$ gfortran program.f90 qr_factor_hh.o ``` 

**Description/Purpose:** This subroutine finds the orthogonal QR factorization of a matrix using the Householder transformations algorithm to be used in solving least squares problems. That is, it decomposes *A* such that

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\begin{pmatrix}&space;a_{11}&space;&a_{12}&space;\\&space;a_{21}&space;&a_{22}&space;\\&space;a_{31}&space;&a_{32}&space;\end{pmatrix}&space;=&space;\begin{pmatrix}&space;q_{11}&space;&q_{12}&space;\\&space;q_{21}&space;&q_{22}&space;\\&space;q_{31}&space;&q_{32}&space;\end{pmatrix}\begin{pmatrix}&space;r_{11}&space;&r_{12}&space;\\&space;0&space;&r_{22}&space;\end{pmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\begin{pmatrix}&space;a_{11}&space;&a_{12}&space;\\&space;a_{21}&space;&a_{22}&space;\\&space;a_{31}&space;&a_{32}&space;\end{pmatrix}&space;=&space;\begin{pmatrix}&space;q_{11}&space;&q_{12}&space;\\&space;q_{21}&space;&q_{22}&space;\\&space;q_{31}&space;&q_{32}&space;\end{pmatrix}\begin{pmatrix}&space;r_{11}&space;&r_{12}&space;\\&space;0&space;&r_{22}&space;\end{pmatrix}" title="\begin{pmatrix} a_{11} &a_{12} \\ a_{21} &a_{22} \\ a_{31} &a_{32} \end{pmatrix} = \begin{pmatrix} q_{11} &q_{12} \\ q_{21} &q_{22} \\ q_{31} &q_{32} \end{pmatrix}\begin{pmatrix} r_{11} &r_{12} \\ 0 &r_{22} \end{pmatrix}" /></a>

and returns the *Q* and *R* matrices. Since *Q* is orthogonal, <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Q^TQ" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;Q^TQ" title="Q^TQ" /></a> should produce the identity matrix and *Q* x *R* should reproduce *A*.

**Input:** 

*A* : REAL - a coefficient matrix of size *m* x *n*

*m* : INTEGER - the number of rows in *A* and *Q*

*n* : INTEGER - the number of columns in *A* and the size of the square matrix *R*

**Output:** 

*Q* : REAL - the orthogonal matrix 'Q' of size *m* x *n* of the QR factorization

*R* : REAL - the upper triangular 'R' factor of size *n* x *n* of the QR factorization

**Usage/Example:**

This routine can be implemented in a program as follows (which follows example 6.1 in Ascher)

```fortran
INTEGER :: m, n, i
REAL*8, ALLOCATABLE :: A(:, :), Q(:, :), R(:, :), QTQ(:, :)

m = 5
n = 3
ALLOCATE(A(1:m, 1:n), Q(1:m, 1:n), R(1:n, 1:n), QTQ(1:n, 1:n))
A = RESHAPE((/1.D0, 0.0D0, 1.0D0, &
			& 2.D0, 3.0D0, 5.0D0, &
			& 5.D0, 3.0D0, -2.0D0, &
			& 3.D0, 5.0D0, 4.0D0, &
			& -1.D0, 6.0D0, 3.0D0/), (/m, n/), ORDER=(/2, 1/))
CALL qr_factor_modgs(A, m, n, Q, R)
WRITE(*,*) "A"
DO i = 1, m
	WRITE(*,*) A(i, :)
END DO
WRITE(*,*) "Q"
DO i = 1, m
	WRITE(*,*) Q(i, :)
END DO
WRITE(*,*) "R"
DO i = 1, n
	WRITE(*,*) R(i, :)
END DO
CALL mat_prod(TRANSPOSE(Q), Q, n, m, n, QTQ)
WRITE(*,*) "QTQ"
DO i = 1, n
	WRITE(*,*) QTQ(i, :)
END DO
CALL mat_prod(Q, R, m, n, n, A)
WRITE(*,*) "A = QR"
DO i = 1, m
	WRITE(*,*) A(i, :)
END DO
```

The outputs from the above code:

```fortran
 A
   1.0000000000000000        0.0000000000000000        1.0000000000000000     
   2.0000000000000000        3.0000000000000000        5.0000000000000000     
   5.0000000000000000        3.0000000000000000       -2.0000000000000000     
   3.0000000000000000        5.0000000000000000        4.0000000000000000     
  -1.0000000000000000        6.0000000000000000        3.0000000000000000     
 Q
 -0.15811388300841900        9.9778515785660965E-002  0.25545570859468680      -0.30866541581688367       0.89694609080623544 /    
 -0.31622776601683794      -0.19955703157132196       0.69185921077727630      -0.47549548008128917      -0.39422312468366333 /    
 -0.79056941504209477        9.9778515785660896E-002 -0.54639137671641325      -0.24538502820116251       -7.9289968925886023E-002 /
 -0.47434164902525688      -0.36585455788075655       0.26609969645279874       0.74270133159821761       0.13687996956399021  /   
  0.15811388300841897      -0.89800664207094805      -0.29448366407443055      -0.25847752219062209       0.12268990550144958 /    
 R
  -6.3245553203367573       -4.7434164902525691       -1.5811388300841895     
  -1.2176676335913021E-016  -7.5166481891864541       -5.2550018313781415     
  -4.1468336586232195E-016  -6.1895392470044310E-016   4.9884823095017978     
   1.5775447872901429E-016  -2.3830241699585675E-016   1.3877787807814457E-017
  -6.8497333484492075E-016   6.3108851118458503E-016   0.0000000000000000     
 QTQ
  0.99999999999999989        2.7755575615628914E-017  -3.4694469519536142E-017   1.3877787807814457E-017   6.9388939039072284E-018 /
   2.7755575615628914E-017   1.0000000000000000       -5.5511151231257827E-017   1.1102230246251565E-016   5.5511151231257827E-017 /
  -3.4694469519536142E-017  -5.5511151231257827E-017   1.0000000000000004        0.0000000000000000        1.8735013540549517E-016 /
   1.3877787807814457E-017   1.1102230246251565E-016   0.0000000000000000       0.99999999999999978        3.4694469519536142E-017 /
   6.9388939039072284E-018   5.5511151231257827E-017   1.8735013540549517E-016   3.4694469519536142E-017   1.0000000000000004 /    
 A = QR
  0.99999999999999922        3.7403564420058577E-017   1.0000000000000004     
   1.9999999999999996        3.0000000000000004        5.0000000000000009     
   4.9999999999999982        3.0000000000000000       -2.0000000000000009     
   2.9999999999999991        5.0000000000000000        4.0000000000000000     
 -0.99999999999999967        6.0000000000000000        3.0000000000000004 
```

**Implementation/Code:** The code for qr_factor_hh can be seen below.

```fortran
SUBROUTINE qr_factor_hh(A, m, n, Q, R)
	IMPLICIT NONE
	
	! Takes as an input the matrix `A` of size `m` x `n` and returns the
	! QR factorization of that matrix, where `Q` is size `m` x `m` and
	! `R` is size `m` x `n`.
	! SOURCES
	! ------------
	! Algorithm taken from:
	! https://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(IN) :: A(m, n)
	REAL*8, INTENT(OUT) :: Q(m, m), R(m, n)
	
	! Creates subroutine-level variables for use in the algorithm.
	INTEGER :: i, j, sub_len
	REAL*8 :: normx, tau, u1, s, w(m), H(m, m), eye(m, m), temp(m, m)
	
	! Initialize `Q` as the identity matrix and `R` as a copy of `A`.
	! These will be modified throughout the algorithm from this initial
	! condition.
	DO i = 1, m
		Q(i, i) = 1.D0
	END DO
	eye = Q
	R = A
	
	! Loops through the lower-triangular part of the `R` matrix
	DO j = 1, n
		! Finds the number of rows to be modified from the diagonal.
		sub_len = SIZE(w(j:))
		
		! Find H = I - tau*w*w^T to put zeros below the diagonal of R.
		CALL l2_vec_norm(R(j:, j), sub_len, normx)
		s = -SIGN(1.D0, R(j, j))
		u1 = R(j, j) - s*normx
		w(j:) = R(j:, j)/u1
		w(j) = 1.D0
		tau = -s*u1/normx
		H = eye
		CALL out_prod_vec(sub_len, sub_len, tau*w(j:), w(j:), H(j:, j:))
		H(j:, j:) = eye(j:, j:) - H(j:, j:)
		
		! Finds R = H*R and Q = Q*H
		CALL mat_prod(H(j:, j:), R(j:, :), sub_len, sub_len, n, temp(j:, :))
		R(j:, :) = temp(j:, :)
		CALL mat_prod(Q(:, j:), H(j:, j:), m, sub_len, sub_len, temp(:, j:))
		Q(:, j:) = temp(:, j:)
	END DO

END SUBROUTINE
```



