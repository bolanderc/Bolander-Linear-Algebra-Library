# Homework 5 Task 10: Direct Methods

*Task: Search the internet for sites that discuss the use of direct methods for the approximate solution of linear systems of systems of equations. Note that direct methods include factorization methods and the standard Gaussian elimination with back substitution. Find at least a couple of sites where limitations of direct methods are discussed. As usual, cite your sites.*



## Direct Methods for Approximating Solutions to Linear Systems

Direct methods for approximating the solution of a linear system of equations are useful in many circumstances. A direct method has a predetermined number of steps that it must complete to get a solution and that solution is guaranteed to be within machine precision of the actual solution for a stable algorithm and well-conditioned problem. Direct methods also allow for a cheaper solution for multiple right-hand side vectors with the same coefficient matrix. Direct methods are generally very useful, but have a few disadvantages as well when compared to an iterative method. Direct methods may require extra storage when compared to an indirect method and do not have the ability to stop when the approximation is "good enough". In addition, direct methods also are limited when dealing with very large matrices, as assembling the A matrix in full form requires memory that may not be available.


* [Source 1](http://www.cs.nthu.edu.tw/~cchen/CS6331/ch2.pdf)
* [Source 2 (and the paper it references)](https://math.stackexchange.com/questions/1842630/direct-vs-iterative-solvers-choice)
* [Source 3](https://www.simscale.com/blog/2016/08/how-to-choose-solvers-for-fem/)
* [Source 4](<https://www.researchgate.net/post/what_is_the_best_method_for_solving_system_of_linear_equations>)