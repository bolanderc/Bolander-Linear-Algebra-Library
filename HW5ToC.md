# Homework 5 Table of Contents

- [x] [Task 1](./Software_Manual/direct_ge_bsin.md)
- [x] [Task 2](./Software_Manual/lu_factor.md)
- [x] [Task 3](./Software_Manual/lu_solve.md)
- [x] [Task 4](./HW5Task4Report.md)
- [x] [Task 5](./HW5Task5Report.md)
- [x] [Task 6](./Software_Manual/solve_normal_equations.md) 
- [x] [Task 7](./Software_Manual/qr_factor_gs.md)
- [x] [Task 8](./HW5Task8Report.md)
- [x] [Task 9](./Software_Manual/mat_dd.md)
- [x] [Task 10](./HW5Task10Report.md)

## Description of Tasks

### Task 1
Implement a method that will return the approximate solution of a square linear system of equations where previous methods are not used. That is, inline the row reduction operations and the backsubstitution methods. Test the speed of the code you generated in this problem and the code that references your previous methods. Try this for increasing sizes of the linear system. You will likely need to use large systems of linear equations - possibly 10,000 by 10,000 to see any kind of time. Use cpu_timing methods in the language you have chosen to do your coding. Add a manual page to document the inline version of the solution process. Report any differences you see in the time it takes to solve the linear systems in the two approaches.

------

### Task 2
Implement a method that returns the LU-factorization of a square matrix. Add an entry to your software manual to document the method you have created. Hint: You can actually modify the Gaussian elimination code in two lines to come up with the new method.

------

### Task 3
Use the LU factorization method you created in the previous step, along with the forward substitution method (for lower triangular square systems) and the back substitution method (for upper triangular square systems) to create a method that will solve a square linear system of equations. Document the method in your software manual.

------

### Task 4
Write a code that will generate a symmetric, positive definite matrix for a given integer, n. Make sure that you add an entry to your software manual with a couple of examples. Your examples should be relatively small for your examples, but you should include a large example in the task solution write-up.

------

### Task 5
Implement the Cholesky factorization method for square matrices. Do not include any pivoting in the algorithm. Document the algorithm in your software manual. Test the code on at least 2 or 3 matrices of different sizes. At least one example should involve a matrix that is bigger than 100×100 in size. Use output from the method you created in the previous task.

------

### Task 6
Write a routine/method that will return an approximate solution of the least squares problem using the normal equation approach. Create an entry in your software manual for the method. Also, make sure you use the Cholesky algorithm that you created in a previous task.

------

### Task 7
Implement the QR factorization of a square matrix. Use the Gram-Schmidt process to create the orthogonal vectors for the orthogonal matrix. Document the method in your software manual. Include examples showing the orthogonal matrix and the other factor which should be upper triangular.

------

### Task 8
Try out your QR-factorization method from the previous task on the Hilbert matrices of sizes *n* = 4, 6, 8, 10. Test to see if the orthogonal matrix is really orthogonal by multiplying *Q*<sup>T</sup>*Q* and comparing the result to the identity matrix. Explain the results you obtain.

------

### Task 9
Implement a method that will return a square diagonally dominant matrix. Document this method in your software manual.

------

### Task 10
Search the internet for sites that discuss the use of direct methods for the approximate solution of linear systems of systems of equations. Note that direct methods include factorization methods and the standard Gaussian elimination with back substitution. Find at least a couple of sites where limitations of direct methods are discussed. As usual, cite your sites.
