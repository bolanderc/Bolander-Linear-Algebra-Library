# Homework 1 Table of Contents

- [x] [Task 1](Software_Manual/SWMToC.md)
- [x] [Task 2](README.md)
- [x] [Task 3](https://bolanderc.github.io/math5610)
- [x] [Task 4](Software_Manual/SWMToC.md)
- [x] [Task 5 - HW](HomeworkToC.md) [Task 5 - Software Manual](Software_Manual/SWMToC.md)
- [x] [Task 6](mylib) 
- [x] [Task 7](./HW1Task7Report.md)
- [x] [Task 8](./HW1Task8Report.md)
- [x] [Task 9](./rand_mat.f90)
- [x] [Task 10](./HW1Task10Report.md)

## Description of Tasks

1. Task: Write a code that will return machine precision for your computer (or any other computer for that matter) in single precision arithmetic. Give the method a name that is descriptive. For example, smaceps(), or something like that. The routine should return the default machine precision in terms of a decimal number. Create a second routine that will return the machine precision in double precision computations. Give the routine a unique name, say dmaceps(). Make sure that your code is fully internally documented with you as the author and so on. An example of a Fortran code is included in this repository. You can translate this into a Python, C, or C++ code to do the work. You can also modify the Fortran and use the given code directly.

------

2. Task: Create a repository on [Github](https://www.github.com) and create a repositry for this homework. Name the homework repository for your work as follows:

   ```
               username/math5610
            
   ```

   Note that your username will be set when you create your Github account. It is important for appropriate grading to name the repository as above.

------

3. Task: Create a repository for Github Pages. You can do this by using the documentation at the following link:

   ```
              https://pages.github.com/
            
   ```

   The information on the web page given above will show you how to build web pages from Jeykll css files. This repository allows access to projects as you create them in your repository. I will use the link created in this repository to grade your Software Manual.

------

4. Task: Create a folder in your "math5610" repository for your for your software manual. Download a copy of the following markdown file for your software manual.

   Put a copy of the file in the folder for the software manual for the codes you write. At this point in the homework, you have two routines that need to be documented in the software manual. These are the single precision and double precision machine epsilon codes from the tasks above.

   > Note: The software manual must use the template to be graded!. The theme can be changed to one you like more, but the format should be consistent.

------

5. Task: Create a table of contents for the homework tasks using markdown. We will go through the basics of markdown and how to set up a table of contents for your work. Also build a markdown file that will serve as a table of contents for the entries in the software manual.

   ------

6. Task: To complete this problem, you will need to go to the Engineering Computer Lab on the third floor of the Engineering Building. Login to one of the computers and open up a Cygwin window. A Linux shell window will pop up for you to use. Complete the following steps:

    1. Log onto a computer (Engineering Lab) and open a command terminal to work in.

    2. Upload/copy the routines that you created in the first problem.

    3. Compile the routines into object modules (.o files). For example, put the files you have uploaded into a folder, say hw1_prob3, using the command

       ```
                   % mkdir hw1_prob3
                 
       ```

       and in a Cygwin/Linux/Unix operating system. Note that the "%" is the command prompt that may appear in the command terminal. Then use

       ```
                   % mv *.f hw1_prob3
                   % cd hw1_prob3
                 
       ```

       to move all files ending with a .f suffix to the folder just created and change the working folder to the folder just created. Finally, compile the files using

       ```
                   % gfortran -c *.f
                 
       ```

       or

       ```
                   % gcc -c *.c
                 
       ```

       using the C-compiler in Cygwin. The result will be a bunch of object files with a suffix of ".o".

    4. The last step in this problem is to create a shared library from the routines you have created.

       ```
                   % ar rcv mylib *.o
                 
       ```

       or

       ```
                   % ar rcv mylib *.o
                   % ranlib mylib.a
                 
       ```

------

7. Task: In this problem you will learn a bit more about your computer and how many processes you can run on the cores included in your computer. You will use OpenMP to do this work.

    1. To start, write code like the following.

       ```
                   program main
                   print *, "hello world!"
                   stop
                   end
                 
       ```

       to implement the good old Hello World introductory program. Compile and run the code you have created. The code should produce the string in the print statement. That is,

       ```
                   hello world!
                 
       ```

       Now, let's do a bit more with this example. We will need to add some code to the example to have each of the cores on your computer write out the same statement. In addition, the way that the code is compiled will also change. That we will do on the other side of the following code.

       ```
                   |      program main
                   |      integer id, nthrds
                   |      integer omp_get_thread_num, omp_get_num_threads
                   |C$OMP PARALLEL PRIVATE(id)
                   |      id = omp_get_thread_num()
                   |      print *, 'hello world from thread', id
                   |C$OMP BARRIER
                   |      if(id .eq. 0) then
                   |        nthrds = omp_get_num_threads()
                   |        print *, 'There are', nthrds, ' threads!'
                   |      end if
                   |C$OMP END PARALLEL
                   |      stop
                   |      end
                 
       ```

       First, there are some things that need to be described in the code:

       1. The line on the left indicates how the code needs to be written in the column. In fortran, the comment line must start as the first character of the line typed in. So, the pipe, "|", is like the edge of the window and is not actually a character on your screen, just the border of the document.

       2. The string "C$OMP" tells the compiler that the line is a line associated with OpenMP. The first letter must appear in the first column and indicates that when compiled without OpenMP extensions, will be ignored. If the OpenMP library is available, the compiler will include the parameters contained in the comments. These strings are the start of a "directive" to the compiler to do something.

       3. For the present time, (1) save the code above into a file named "hello.f", (2) compile the code above, using

          ```
                          % gfortran -fopenmp hello.f -o hello
                        
          ```

          and then run the compiled version of the code using the command

          ```
                          % ./hello.exe
                        
          ```

          Report the results of running the code.

------

8. Task: Read the three disaster articles at the site

    ```
               http://www-users.math.umn.edu/~arnold//disasters/
             
    ```

    Write a brief paragraph on each of the disasters describing the particular problem as described.

------

9. Task: Write a routine that will generate a random matrix of a given size. That is, input the number of rows and columns and output the matrix created by setting each entry in the matrix to a random value between zero and one.

------

10. Task: Search the internet for sites that discuss linear algebra packages for solving linear algebra problems. Find a couple of sites that most closely line up with what you think we will be doing in the class. Reference the sites in a brief discussion.