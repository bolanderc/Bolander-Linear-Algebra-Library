# Homework 6 Table of Contents

- [x] [Task 1](./Software_Manual/qr_sq_solve.md)
- [x] [Task 2](./HW6Task2Report.md)
- [x] [Task 3](./HW6Task3Report.md)
- [x] [Task 4](./Software_Manual/ls_solveqrmod.md)
- [x] [Task 5](./Software_Manual/ls_solveqr.md)
- [x] [Task 6](./Software_Manual/jacobi_solve.md) 
- [x] [Task 7](./Software_Manual/gaussseidel_solve.md)
- [x] [Task 8](./HW6Task8Report.md)
- [ ] [Task 9](https://bolanderc.github.io/math5610)
- [ ] [Task 10](https://bolanderc.github.io/math5610)

## Description of Tasks

### Task 1
Implement a method that will compute the solution of a square linear system of equations using the QR-factorization of the matrix. Give examples and document the code in your software manual.

------

### Task 2
Create another version of the QR-factorization algorithm using the Modified Gram-Schmidt process. Document the code in your software manual. For examples, use the same matrices used when testing the modified version. Compare the results to the first version of the QR-factorization.

------

### Task 3
Create a third version of the QR-factorization algorithm using Householder Transformations. As usual, document you code in your software manual. Use the third incarnation of the code on the same matrices as the previous two QR-factorization and compare/explain your results.

------

### Task 4
Build a code that will solve the least squares problem using QR factorization. Document the code in your software manual. Use the modified Gram-Schmidt algorithm to compute the QR factorization.

------

### Task 5

Build a code that will solve the least squares problem using QR factorization. Document the code in your software manual. (I used the Householder transformation QR factorization here).

------

### Task 6

Implement the Jacobi Iteration algorithm for computing a sequence of approximate solutions for the linear system equations Ax=b. Include a software manual entry for the code you write. Include at least one example that solves a system of equations with 1000 equations in 1000 unknowns. You can use the code you developed to create a diagonally dominant system.

------

### Task 7

Repeat the previous task using the Gauss-Seidel algorithm.

------

### Task 8

Compare the Jacobi and Gauss-Seidel in terms of the number of iterations needed to converge to a given tolerance. For example, compute the number of iterations needed to produce a solution to within four digits of accuracy. Tabulate and/or plot the number of iterations needed for the two methods as the size of the system changes. Do this for large systems of equations - greater than 500 by 500.

------

### Task 9

Do an internet search for pages that discuss the difference between the solution of the least squares problem using the normal equations and the solution via QR factorization of the matrix. Make sure that you cite the sites you use and include a couple of paragraphs in your own words.

------

### Task 10

Complete an internet search for sites that discuss the stability of various algorithms used in computing the QR factorization of both rectangular and square matrices. Give a brief description of what you found and include citations for the pages you find.