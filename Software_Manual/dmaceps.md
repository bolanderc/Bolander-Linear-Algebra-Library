# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual
This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           dmaceps

**Author:** Christian Bolander

**Language:** Fortran. This code can be compiled using the GNU Fortran compiler by
```$ gfortran -c dmaceps.f90```

and can be added to a program using

```$ gfortran program.f90 dmaceps.o ``` 

**Description/Purpose:** This routine will compute the double precision value for the machine epsilon. This is a routine for analyzing the behavior of any computer. This
usually will need to be run one time for each computer.

**Input:** There are no inputs needed.

**Output:** 

*mach_eps* : REAL - double precision machine epsilon value

*i* : INTEGER - the number of binary digits that represent machine epsilon

**Usage/Example:**

This routine has no inputs and can be implemented in a program as follows

```fortran
CALL dmaceps(mach_eps, i)
WRITE(*,*) mach_eps, i
```

The outputs from the above code:

```fortran
   1.1102230246251565E-016          53
```

The number of decimal digits that can be represented is roughly sixteen (E-16 on the
end of the first value) and the number of binary digits that represent machine epsilon are 53.

**Implementation/Code:** The code for dmaceps can be seen [here](dmaceps.f90).

