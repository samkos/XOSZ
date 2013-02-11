program main

  implicit none
  include 'mpif.h' 

  double precision, dimension(:,:), allocatable :: subarray,subarray2
  integer :: filetype, rank, ierror, thefile, nb_tasks, nb_i_blocks,nb_k_blocks, icolumn, iline, i,k,k4
  integer :: sizes(2), subsizes(2), starts(2) 
  character(len=10) :: filename
  integer(kind=MPI_OFFSET_KIND) :: disp

  filename = "test2D"

  call MPI_init(ierror)
  
 call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_tasks, ierror)
 call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror) 


 nb_k_blocks=1
 nb_i_blocks=nb_tasks
 
 nb_k_blocks=nb_tasks
 nb_i_blocks=1
 
 do while (nb_k_blocks.gt.nb_i_blocks)
    nb_i_blocks=nb_i_blocks*2
    nb_k_blocks=nb_k_blocks/2
 enddo


 icolumn=mod(rank,nb_i_blocks)
 iline=rank/nb_i_blocks



  sizes(1)=10 
  sizes(2)=10 
  subsizes(1)=10/nb_i_blocks  
  subsizes(2)=10/nb_k_blocks  
  allocate(subarray(subsizes(1),subsizes(2)))
  allocate(subarray2(subsizes(1),subsizes(2)))

  starts(1)=icolumn*subsizes(1)
  starts(2)=iline*subsizes(2) 

  print '(13(I4,X))',rank,iline,icolumn,sizes(1),sizes(2),subsizes(1),subsizes(2),starts(1),starts(2)
  
  call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, & 
       MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,       & 
       filetype, ierror)

  call MPI_TYPE_COMMIT(filetype,ierror)


  if (nb_tasks.eq.1) then 

     call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, & 
          MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
          MPI_INFO_NULL, thefile, ierror) 

     disp = 0_MPI_OFFSET_KIND
     call MPI_FILE_SET_VIEW(thefile, disp, MPI_DOUBLE_PRECISION, & 
          filetype, 'native',   MPI_INFO_NULL, ierror) 

     
     print *,"filling the array : ",subsizes(1),subsizes(2)
     do k=1,subsizes(2)
        do i=1,subsizes(1)
           subarray(i,k) = 10000+i+k*100
        end  do
     end do
     
     
     call MPI_FILE_WRITE_ALL(thefile, subarray, subsizes(1)*subsizes(2), MPI_DOUBLE_PRECISION, & 
          MPI_STATUS_IGNORE, ierror) 

     call MPI_FILE_CLOSE(thefile, ierror) 
  end if

  call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, & 
       MPI_MODE_RDONLY , & 
       MPI_INFO_NULL, thefile, ierror) 
  print *,"ierror open :" ,ierror

  disp = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_VIEW(thefile, disp, MPI_DOUBLE_PRECISION, & 
       filetype, 'native',   MPI_INFO_NULL, ierror) 
  print *,"ierror set_view :" ,ierror,subsizes(1)*subsizes(2)
  
  call MPI_FILE_READ_ALL(thefile, subarray2, subsizes(1)*subsizes(2), MPI_DOUBLE_PRECISION, & 
       MPI_STATUS_IGNORE, ierror) 
  print *,"ierror read_all :" ,ierror

  call MPI_FILE_CLOSE(thefile, ierror) 
  
!!$  do k=subsizes(2),1,-1
!!$     print '(I5,"b",I2,"->",100(F6.0,X))',rank,k,(subarray(i,k),i=1,subsizes(1))
!!$  end do
!!$
  do k=subsizes(2),1,-1
     print '(I5,"o",I2,"->",100(F6.0,X))',rank,k,(subarray2(i,k),i=1,subsizes(1))
     call flush(6)
  end do

  call MPI_TYPE_FREE(filetype,ierror)
     
  !  call MPI_TYPE_FREE(filetype,ierror)
  
  call MPI_FINALIZE(ierror) 
  

  
  deallocate(subarray)
  deallocate(subarray2)
  
end  program main
