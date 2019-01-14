      program main
      integer id, nthrds
      integer omp_get_thread_num, omp_get_num_threads
C$OMP PARALLEL PRIVATE(id)
      id = omp_get_thread_num()
      print *, 'hello world from thread', id
C$OMP BARRIER
      if(id .eq. 0) then
        nthrds = omp_get_num_threads()
        print *, 'There are', nthrds, ' threads!'
      end if
C$OMP END PARALLEL
      stop
      end
