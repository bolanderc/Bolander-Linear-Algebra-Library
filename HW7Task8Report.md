# Homework 7 Task 8: Preconditioning Strategies

*Task: Look for internet sites that include descriptions of preconditioning of systems of equations. Document at least 3 different preconditioning strategies.*



## Preconditioning

Iterative methods, and the convergence of such methods, are heavily reliant on well-conditioned systems. If a system is poorly conditioned, then the system may not converge to a solution at all (or it will converge to an incorrect system). Thus, a preconditioning matrix can be applied to the original system of equations i.e.:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vv{B}}\mathbf{\vv{A}}\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{B}}\mathbf{\vv{b}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vv{B}}\mathbf{\vv{A}}\mathbf{\vv{x}}&space;=&space;\mathbf{\vv{B}}\mathbf{\vv{b}}" title="\mathbf{\vv{B}}\mathbf{\vv{A}}\mathbf{\vv{x}} = \mathbf{\vv{B}}\mathbf{\vv{b}}" /></a>

The condition number for **BA** may then be smaller than the condition number of **A** alone<sup>[1]</sup>. This also presents the requirement for preconditioners that they must<sup>[2]</sup>:

* Have better convergence than the original system
* Operations with **B** should be easy to perform

Different preconditioning strategies include:

* Jacobi preconditioner<sup>[2]</sup>
* Gauss-Seidel, SOR, and SSOR<sup>[2]</sup>
* Incomplete LU factorization<sup>[2]</sup>
* Multigrid preconditioning<sup>[3]</sup>


* [Source 1](<<https://en.wikipedia.org/wiki/Preconditioner>)
* [Source 2](<http://users.tem.uoc.gr/~spapadem/NLA_2016/Preconditioning.pdf>)
* [Source 3](<http://www.hpcs.cs.tsukuba.ac.jp/~tatebe/research/paper/CM93-tatebe.pdf>)