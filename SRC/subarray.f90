program main

  implicit none
  include 'mpif.h' 

  double precision, dimension(:,:), allocatable :: subarray,subarray2
  integer :: filetype, rank, ierror, thefile, nb_tasks
  integer :: sizes(2), subsizes(2), starts(2) 
  
  call MPI_init(ierror)
  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_tasks, ierror)
 call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror) 
  sizes(1)=100 
  sizes(2)=100 
  subsizes(1)=100 
  subsizes(2)=100/nb_tasks 
  allocate(subarray(100,subsizes(2)))
  allocate(subarray2(100,subsizes(2)))

  starts(1)=0 
  starts(2)=rank*subsizes(2) 
  
  call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, & 
       MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,       & 
       filetype, ierror)



  call MPI_FILE_OPEN(MPI_COMM_WORLD, 'testfile', & 
       MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
       MPI_INFO_NULL, thefile, ierror) 


  call MPI_FILE_SET_VIEW(thefile, 0, MPI_DOUBLE_PRECISION, & 
       filetype, 'native',   MPI_INFO_NULL, ierror) 
  call MPI_FILE_WRITE_ALL(thefile, subarray, 100*subsizes(2), MPI_DOUBLE_PRECISION, & 
       MPI_STATUS_IGNORE, ierror) 

  call MPI_FILE_CLOSE(thefile, ierror) 


  call MPI_FILE_OPEN(MPI_COMM_WORLD, 'testfile10', & 
       MPI_MODE_RDONLY , & 
       MPI_INFO_NULL, thefile, ierror) 


  call MPI_FILE_READ_ALL(thefile, subarray2, 100*subsizes(2), MPI_DOUBLE_PRECISION, & 
       MPI_STATUS_IGNORE, ierror) 

  call MPI_FILE_CLOSE(thefile, ierror) 


  
!  call MPI_TYPE_FREE(filetype,ierror)

  call MPI_FINALIZE(ierror) 
  

  deallocate(subarray)
  deallocate(subarray2)

end  program main
