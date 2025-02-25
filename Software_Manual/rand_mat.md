# Generate a Random Matrix

This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           rand_mat

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c rand_mat.f90```

and can be added to a program using

```$ gfortran program.f90 rand_mat.o ``` 

**Description/Purpose:** This routine takes a matrix with a given number of rows and columns and fills each element of the matrix with a random number from 0 to 1 not including one. It can be used to quickly generate a matrix to be tested with other linear algebra routines.

**Input:** 

*r* : INTEGER - number of rows in the matrix

*c* : INTEGER - number of columns in the matrix

**Output:** 

*mat* : REAL - an array of size (r, c) containing random numbers

**Usage/Example:**

This routine has two inputs, `r` and `c`, and can be implemented in a program as follows

```fortran
r = 4
c = 4
CALL rand_mat(r, c, mat)
WRITE(*,*) mat
```

The outputs from the above code:

```fortran
  0.19553478045561368       0.64591609807036332       0.28778363522487882       0.77345893123498888     
  0.51572313511517420       0.53781067742119892       0.58312487473801156       0.48948996514044674     
  0.87174012219643315       0.21245055839531113       0.64895831046172725       0.82264052755663353     
  0.14500803756514291        8.3044258666368442E-002  0.42828671870518820       0.28045865572420892     
```

**Implementation/Code:** The code for rand_mat can be seen below.

```fortran
SUBROUTINE rand_mat(r, c, mat)
	IMPLICIT NONE
	INTEGER, INTENT(in) :: r, c
	REAL*8, INTENT(out) ::  mat(r, c)
	! Fill a matrix with `r` rows and `c` columns with a random
	! number between 0 and 1 (not including 1).
	CALL RANDOM_NUMBER(mat)
END SUBROUTINE
```



