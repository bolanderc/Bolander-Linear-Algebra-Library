# Homework 1 Task 6: Creation of the mylib shared library

*Task: To complete this problem, you will need to go to the Engineering Computer Lab on the third floor of the Engineering Building. Login to one of the computers and open up a Cygwin window. A Linux shell window will pop up for you to use. Complete the following steps:*

1. *Log onto a computer (Engineering Lab) and open a command terminal to work in.*

2. *Upload/copy the routines that you created in the first problem.*

3. *Compile the routines into object modules (.o files). For example, put the files you have uploaded into a folder, say hw1_prob3, using the command*

   ```
               % mkdir hw1_prob3
   ```

   *and in a Cygwin/Linux/Unix operating system. Note that the "%" is the command prompt that may appear in the command terminal. Then use*

   ```
               % mv *.f hw1_prob3
               % cd hw1_prob3
   ```

   *to move all files ending with a .f suffix to the folder just created and change the working folder to the folder just created. Finally, compile the files using*

   ```
   % gfortran -c *.f
   ```

   *or*

   ```
               % gcc -c *.c
   ```

   *using the C-compiler in Cygwin. The result will be a bunch of object files with a suffix of ".o".*

4. *The last step in this problem is to create a shared library from the routines you have created.*

   ```
               % ar rcv mylib *.o
   ```

   *or*

   ```
               % ar rcv mylib *.o
               % ranlib mylib.a
   ```

----------

By following the steps outlined in this task, I have created the mylib shared library, with the contents outlined here

```
christian@Joberta:~/Python Projects/math5610$ ar t mylib
abs_err_n.o
abs_err_vecl1.o
abs_err_vecl2.o
abs_err_vecl_inf.o
dmaceps.o
l1_vec_norm.o
l2_vec_norm.o
l_inf_vec_norm.o
mat_1norm.o
mat_infnorm.o
mat_prod.o
rand_mat.o
rel_err_n.o
smaceps.o
s_mult_vec.o
sym_mat_gen.o
vec_add.o
vec_cross_prod3.o
vec_dot_prod.o
mat_dd.o
s_mult_mat.o
mat_add.o
out_prod_vec.o
lss_diag.o
backsub.o
forwardsub.o
mat_row_ech.o
direct_ge_bs.o
sym_dd_mat_gen.o
```

