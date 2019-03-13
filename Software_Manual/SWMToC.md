# Software Manual Table of Contents
This is the Software Manual created for Math 5610: Computational Linear Algebra and Solution of Systems of Equations. Contained within are a list of documented subroutines created for solving systems of equations computationally.

## Subroutines
The subroutines contained in this manual will be split according to their general purpose in solving computational linear algebra problems

### Basic Linear Algebra
- [abs_err_n](./abs_err_n.md) : Calculates the absolute error for a value and its approximation.
- [abs_err_vecl1](./abs_err_vecl1.md) : Calculates the absolute error between two vectors using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_1" title="\ell_1" /></a>-norm.
- [abs_err_vecl2](./abs_err_vecl2.md) : Calculates the absolute error between two vectors using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_2" title="\ell_2" /></a>-norm.
- [abs_err_vecl_inf](./abs_err_vecl_inf.md) : Calculates the absolute error between two vectors using the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm.
- [dmaceps](dmaceps.md) : Finds the machine epsilon for double precision.
- [l1_vec_norm](./l1_vec_norm.md) : Calculates the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_1" title="\ell_1" /></a>-norm of a vector.
- [l2_vec_norm](./l2_vec_norm.md) : Calculates the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_2" title="\ell_2" /></a>-norm of a vector.
- [l_inf_vec_norm](./l_inf_vec_norm.md) : Calculates the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm of a vector.
- [mat_1norm](./mat_1norm.md) : Calculates the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_1" title="\ell_1" /></a>-norm of a square matrix.
- [mat_add](./mat_add.md) : Computes the sum of two matrices of equal size.
- [mat_dd](./mat_dd.md) : Generates a random, square, diagonally dominant matrix with n rows and columns.
- [mat_infnorm](./mat_infnorm.md) : Calculates the <a href="https://www.codecogs.com/eqnedit.php?latex=\ell_\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell_\infty" title="\ell_\infty" /></a>-norm of a square matrix.
- [mat_prod](./mat_prod.md) : Calculates the product of two matrices with equal inner dimension.
- [mat_row_ech](./mat_row_ech.md) : Transforms an arbitrary matrix of size (m, n) into row echelon form.
- [out_prod_vec](./out_prod_vec.md) : Calculates the outer product of two vectors.
- [rand_mat](rand_mat.md) : Generates a random matrix of size (r, c).
- [rel_err_n](./rel_err_n.md) : Calculates the relative error for a value and its approximation.
- [smaceps](smaceps.md) : Finds the machine epsilon for single precision.
- [s_mult_mat](./s_mult_mat.md) : Multiplies a matrix by a scalar value.
- [s_mult_vec](./s_mult_vec.md) : Multiplies a vector by a scalar value.
- [sym_dd_mat_gen](./sym_dd_mat_gen.md) : Generates a symmetric, diagonally dominant matrix.
- [sym_mat_gen](./sym_mat_gen.md) : Generates a symmetric matrix.
- [vec_add](./vec_add.md) : Adds two vectors together.
- [vec_cross_prod3](./vec_cross_prod3.md) : Calculates the cross product of two vectors of length 3.
- [vec_dot_prod](./vec_dot_prod.md) : Calculates the dot product of two vectors of the same size.

### Direct Methods
- [backsub](./backsub.md) : Computes the solution of a square, linear system of equations where the coefficient matrix is an upper-triangular matrix using the backward substitution algorithm.
- [direct_ge_bs](./direct_ge_bs.md) : Solves a square linear system of equations using Gaussian elimination and backward substitution.
- [direct_ge_bsin](./direct_ge_bsin.md) : Solves a square linear system of equations using Gaussian elimination and backward substitution in-line. More computationally efficient than `direct_ge_bs`.
- [forwardsub](./forwardsub.md) : Computes the solution of a square, linear system of equations where the coefficient matrix is a lower-triangular matrix using the forward substitution algorithm.
- [lss_diag](./lss_diag.md) : Computes the solution of a square, linear system of equations where the coefficient matrix is a diagonal matrix.
- [lu_factor](./lu_factor.md) : Uses LU decomposition on a square coefficient matrix to decompose the matrix into its L and U components.
- [lu_solve](./lu_solve.md) : Solves a square linear system of equations using LU decomposition, forward substitution, and backward substitution.
