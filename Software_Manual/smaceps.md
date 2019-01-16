# Math 5610 Computational Linear Algebra and Solution of Systems of Equations Software Manual
This is a part of the student software manual project for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. 

**Routine Name:**           smaceps

**Author:** Christian Bolander

**Language:** Python. This code was written in Python 3.

**Description/Purpose:** This routine will compute the single precision value for the machine epsilon. This is a routine for analyzing the behavior of any computer. This
usually will need to be run one time for each computer.

**Input:** There are no inputs needed.

**Output:** This routine returns a single precision value for the number of decimal digits that can be represented on the
computer being used to run the program.

**Usage/Example:**

This routine has no inputs and can be implemented as follows

      import mach_prec as mp
      

      mp.smaceps()

If the value returned by mp.smaceps() is printed, the output is:

      5.960464477539063e-08

This represents the decimal value of machine epsilon for single precision. The number of decimal digits that can be represented is roughly eight (E-08 on the
end of the second value).

**Implementation/Code:** The following is the code for dmaceps()

      def smaceps()
          mach_eps = 2**-(24 - 1)/2
          return mach_eps

