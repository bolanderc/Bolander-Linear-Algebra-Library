# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**          qr_factor_gs

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c qr_factor_gs.f90```

and can be added to a program using

```$ gfortran program.f90 qr_factor_gs.o ``` 

**Description/Purpose:** This subroutine finds the orthogonal QR factorization of a matrix to be used in solving least squares problems using the classical Gram Schmidt orthogonalization algorithm. That is, it decomposes *A* such that

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
CALL qr_factor_gs(A, m, n, Q, R)
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
   1.0000000000000000        0.0000000000000000        1.0000000000000000      /
   2.0000000000000000        3.0000000000000000        5.0000000000000000      /
   5.0000000000000000        3.0000000000000000       -2.0000000000000000      /
   3.0000000000000000        5.0000000000000000        4.0000000000000000      /
  -1.0000000000000000        6.0000000000000000        3.0000000000000000      /
 Q
  0.15811388300841897       -9.9778515785660896E-002  0.25545570859468664      /
  0.31622776601683794       0.19955703157132179       0.69185921077727630      /
  0.79056941504209477       -9.9778515785660840E-002 -0.54639137671641314      /
  0.47434164902525688       0.36585455788075660       0.26609969645279863      /
 -0.15811388300841897       0.89800664207094805      -0.29448366407443044      /
 R
   6.3245553203367590        4.7434164902525691        1.5811388300841895      /
   0.0000000000000000        7.5166481891864541        5.2550018313781406      /
   0.0000000000000000        0.0000000000000000        4.9884823095017978      /
 QTQ
  0.99999999999999989        2.7755575615628914E-017   2.0816681711721685E-017 /
   2.7755575615628914E-017  0.99999999999999989       -5.5511151231257827E-017 /
   2.0816681711721685E-017  -5.5511151231257827E-017  0.99999999999999989      /
 A = QR
   1.0000000000000000        0.0000000000000000        1.0000000000000000      /
   2.0000000000000000        3.0000000000000000        5.0000000000000000      /
   5.0000000000000000        3.0000000000000000       -2.0000000000000000      /
   3.0000000000000000        5.0000000000000000        4.0000000000000000      /
  -1.0000000000000000        6.0000000000000000        3.0000000000000000      /
```

**Implementation/Code:** The code for qr_factor_gs can be seen below.

```fortran
SUBROUTINE qr_factor_gs(A, m, n, Q, R)
	IMPLICIT NONE
	
	! Takes as a input a coefficient matrix `A` of size `m` x `n` and
	! finds the QR factorization of that matrix using the classical
	! Gram-Schmidt orthogonalization algorithm. Outputs `Q` and `R`,
	! which are of size `m` x `n` and `n` x `n` respectively. Uses the
	! `vec_dot_prod` and `l2_vec_norm` subroutines.
	INTEGER, INTENT(IN) :: m, n
	REAL*8, INTENT(IN) :: A(1:m, 1:n)
	REAL*8, INTENT(OUT) :: Q(1:m, 1:n), R(1:n, 1:n)
	
	INTEGER i, j
	
	! Iterates over each column to perform the QR factorization
	DO j = 1, n
		Q(:, j) = A(:, j)
		
		! Calculates the dot product of the columns to the left of the
		! current column to find R(i, j).
		DO i = 1, j - 1
			! This line is what makes this the "classical" method.
			CALL vec_dot_prod(A(:, j), Q(:, i), m, R(i, j))
			
			! Updates Q(j) using the previous columns and R(i, j)
			Q(:, j) = Q(:, j) - R(i, j)*Q(:, i)
		END DO
		
		! Computes the new Q column using R.
		CALL l2_vec_norm(Q(:, j), m, R(j, j))
		Q(:, j) = Q(:, j)/R(j, j)
	END DO
	
END SUBROUTINE
```



