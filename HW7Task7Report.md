# Homework 7 Task 7: Iterative Methods

*Task: Do an internet search on the use of iterative methods for the solution of linear systems of equations. Provide a list of at least three methods that do not use preconditioning of the system.*



## Iterative Methods

In contrast to direct methods, iterative methods aim to save time and space by iterating on a an approximate solution to a system of equations and reducing some form of error<sup>[1]</sup>. The tradeoff is that these methods are far less predictable than direct methods and may not converge to the correct solution. To improve the convergence performance of a specific system, preconditioning can be used by taking the system and pre-multiplying both sides by a matrix. The hope is that this new system (which has the same solution as the original system) now has a coefficient matrix on the left hand side that has more favorable spectral properties<sup>[2]</sup>. A method that we have not covered in this class is the generalized minimal residual method (GMRES) that is used to solve nonsymmetric systems of linear equations. This method finds orthonormal vectors to form a basis for the Krylov subspace and then solves a linear least squares problem<sup>[3]</sup>.

There are some methods that do not use preconditioning of the system in their most general form. These methods include:

* Jacobi Iteration<sup>[1]</sup>
* Gauss-Seidel<sup>[1]</sup>
* Richardson Iteration<sup>[4]</sup>

Each of these methods can, however, work as a good preconditioner for more advanced iterative solvers.


* [Source 1](<http://www.netlib.org/utk/people/JackDongarra/WEB-PAGES/SPRING-2003/lect08.pdf>)
* [Source 2](<http://matrixeditions.com/Woz.1-2.pdf>)
* [Source 3](<https://en.wikipedia.org/wiki/Generalized_minimal_residual_method>)
* [Source 4](<https://en.wikipedia.org/wiki/Modified_Richardson_iteration>)