# Cholesky Factorization of a Symmetric, Positive Definite Matrix

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           cholesky_factor

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by

```$ gfortran -c cholesky_factor.f90```

and can be added to a program using

```$ gfortran program.f90 cholesky_factor.o ``` 

**Description/Purpose:** This routine takes a symmetric, positive definite matrix and finds its Cholesky factorization, ***G***, such that

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vv{A}}&space;=&space;\mathbf{\vv{G}}\mathbf{\vv{G}}^T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vv{A}}&space;=&space;\mathbf{\vv{G}}\mathbf{\vv{G}}^T" title="\mathbf{\vv{A}} = \mathbf{\vv{G}}\mathbf{\vv{G}}^T" /></a>

where ***G*** is a lower triangular matrix. The output of this method is that ***G*** is stored in the lower triangular part of the input ***A***, like

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vv{A}}_{output}&space;=&space;\begin{bmatrix}&space;g_{11}&space;&a_{12}&space;&a_{13}&space;\\&space;g_{21}&space;&g_{22}&space;&a_{23}&space;\\&space;g_{31}&space;&g_{32}&space;&g_{33}&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vv{A}}_{output}&space;=&space;\begin{bmatrix}&space;g_{11}&space;&a_{12}&space;&a_{13}&space;\\&space;g_{21}&space;&g_{22}&space;&a_{23}&space;\\&space;g_{31}&space;&g_{32}&space;&g_{33}&space;\end{bmatrix}" title="\mathbf{\vv{A}}_{output} = \begin{bmatrix} g_{11} &a_{12} &a_{13} \\ g_{21} &g_{22} &a_{23} \\ g_{31} &g_{32} &g_{33} \end{bmatrix}" /></a>

The Cholesky factorization is useful in least squares problems as well as being a valuable method to determine if a matrix is symmetric and positive definite.

**Input:** 

*n* : INTEGER - number of rows and columns in the matrix

*A* : REAL - a symmetric, positive definite array of size (n, n)

**Output:** 

*A* : REAL - the Cholesky factorization of *A* as the lower triangular and diagonal with the upper triangular (without the diagonal) filled with the original elements of *A*.

**Usage/Example:**

This routine can be implemented in a program as follows

```fortran
IMPLICIT NONE

INTEGER :: n, i, error
REAL*8, ALLOCATABLE :: mat(:, :)

n = 3

ALLOCATE(mat(n, n))

CALL spd_mat_gen(mat, n)

WRITE(*, *) "SPD MATRIX"
DO i = 1, n
	WRITE(*, *) mat(i, :)
END DO

CALL cholesky_factor(mat, n, error)

WRITE(*,*) "CHOLESKY FACTORIZATION"
DO i = 1, n
	WRITE(*, *) mat(i, :)
END DO
```

The outputs from the above code for several iterations:

```fortran
n = 3
 SPD MATRIX
   1.4562893868127005       0.41766112582738740        1.1072287345244631     
  0.41766112582738740       0.18262626225367931       0.29113699430590895     
   1.1072287345244631       0.29113699430590895       0.86936702617046580     
 CHOLESKY FACTORIZATION
   1.2067681578549794       0.41766112582738740        1.1072287345244631     
  0.34609889489442336       0.25068270224835659       0.29113699430590895     
  0.91751570284432538      -0.10536896347408566       0.12817699770608676 

n = 4
 SPD MATRIX
  0.51158588923659165  0.83823970663441538  0.18292133423321394  0.50232653457365195     
  0.83823970663441538  1.4799031965642322   0.30653084541539288  0.83881474976684034     
0.18292133423321394 0.30653084541539288  8.4358977183772194E-002 0.13737083047203774     
  0.50232653457365195  0.83881474976684034  0.13737083047203774  0.60147574581917396     
 CHOLESKY FACTORIZATION
  0.71525232557230578  0.83823970663441538  0.18292133423321394  0.50232653457365195     
   1.1719496416368893  0.32624719773723548  0.30653084541539288  0.83881474976684034     
0.25574378117100183 2.0879911822754612E-002 0.13608131708065188  0.13737083047203774     
0.70230674772251556 4.8265881097165980E-002 -0.31780659536203587 7.0073893775801283E-002

```

An additional, large example that is verified with another solver can be seen [here](../HW5Task5Report.md).

**Implementation/Code:** The code for cholesky_factor can be seen below.

```fortran
SUBROUTINE cholesky_factor(A, n)
	IMPLICIT NONE
	
	! Takes as inputs a symmetric, positive definite matrix `A` of size
	! `n` and outputs a modified version of `A` with the Cholesky
	! factorization stored in the lower triangular and diagonal parts of
	! the matrix and the original elements in the upper triangular part.
	INTEGER, INTENT(IN) :: n
	REAL*8, INTENT(INOUT) :: A(1:n, 1:n)
	
	INTEGER :: i, j, k
	REAL*8 :: factor
	
	! The Cholesky factorization method first loops through all of the
	! diagonal elements and takes their square root. If this process fails
	DO k = 1, n - 1
		
		A(k, k) = SQRT(A(k, k))
		
		! Then loops through the lower triangular components to compute
		! the Cholesky decomposition.
		DO i = k + 1, n ! rows
			A(i, k) = A(i, k)/A(k, k)
		END DO
		DO j = k + 1, n ! rows
			DO i = j, n ! columns
				A(i, j) = A(i, j) - A(i, k)*A(j, k)
			END DO
		END DO
	END DO
	! The last diagonal value is simply factored to its square root.
	A(n, n) = SQRT(A(n, n))
	
END SUBROUTINE
```
