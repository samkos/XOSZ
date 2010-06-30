module selectprec
  integer, parameter :: prec=kind(1.d0)
end module selectprec

module uncol
  
  use selectprec

  integer, parameter :: max_tasks=128, max_msg_size=200 &
       , user_timers=6, system_timers=3                 &
       , all_timers=user_timers + system_timers
  real(8), save    ::  tchrono(0:all_timers), uchrono(0:all_timers)
  logical, save ::  is_on(0:all_timers)=.false.

  integer, save :: nb_tasks,my_task, master_task        &
       , my_tid, host_tid

  integer, save :: tid( 0:max_tasks)

  logical, save :: power_of_two_tasks

  character*76, save ::  host_name, node_name
  
  integer,save :: nb_msg=0


  interface snd_msg
     module procedure snd_real_msg0d, snd_real_msg1d, snd_real_msg2d &
                     ,snd_int_msg0d
  end interface


  interface rcv_msg
     module procedure rcv_real_msg0d, rcv_real_msg1d, rcv_real_msg2d &
                     ,rcv_int_msg0d
  end interface

  interface bcast
     module procedure real_bcast, int_bcast
  end interface

  interface global_add
     module procedure global_int_add, global_real_add&
          & , global_int_madd1d, global_real_madd1d
  end interface

  interface global_min
     module procedure global_int_min, global_real_min
  end interface

  interface global_max
     module procedure global_int_max, global_real_max
  end interface
  
  include "mpif.h"

contains

  !***********************************************************************
  !  begin uncol_routines.f
  !***********************************************************************
  !       
  !  UNiversal COmmunication Library  (compact version)
  !    Samuel Kortas & Bertrand Melty
  !    version 2.3.0     xx Feb 1994
  !
  !***********************************************************************
  !  New in this version :
  !    - Initialization only for structured mesh and 2^p processors
  !***********************************************************************

  !***********************************************************************
  !***********************************************************************
  !                                                                      *
  !     INITIALIZATION ROUTINES                                          *
  !                                                                      *
  !***********************************************************************
  !***********************************************************************

  subroutine initialize_processors

    !  This routine initializes the tasks.
    implicit none

    integer i
    integer info, error
    character*1 arch

    master_task=0

    call MPI_init(error)

    call MPI_COMM_RANK(MPI_COMM_WORLD, my_task, error)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_tasks, error)

    if (my_task .eq. 0) then
       print *, ' Initialization successfully completed for ',nb_tasks, ' tasks'
    endif
    
    return

  end subroutine initialize_processors

  !***********************************************************************

  subroutine release_processors

    !  This routine terminates all tasks but the master_task, in case of
    !  PVM, or finishes the tracing, in case of PICL.

    implicit none

    integer error

!!    call MPI_finalize(error)

    return

  end subroutine release_processors


  !***********************************************************************
  !   reduction operations
  !***********************************************************************


  function global_int_add(i) result (j)
    implicit none

    integer          :: i,j,error
    real(kind=prec) :: r,s


    call timer_start(2)
    r=real(i)
    call MPI_Allreduce(r,s,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
    j=int(s)
    call timer_stop(2)

  end function global_int_add

  !***********************************************************************

  function global_int_max(i) result (j)
    implicit none

    integer          :: i,j,error
    real(kind=prec) :: r,s

    call timer_start(2)
    r=real(i)
    call MPI_Allreduce(r,s,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)
    j=int(s)
    call timer_stop(2)
    
  end function global_int_max

  !***********************************************************************

  function global_int_min(i) result (j)
    implicit none

    integer          :: i,j,error
    real(kind=prec) :: r,s

    call timer_start(2)
    r = i
    j = global_real_min(r)
    call timer_stop(2)
    
  end function global_int_min

  !***********************************************************************

  function global_real_add(r) result (s)
    implicit none

    real(kind=prec), intent(in)  :: r
    real(kind=prec) :: s
    integer :: error
    
    call timer_start(2)
    call MPI_Allreduce(r,s,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
    call timer_stop(2)

  end function global_real_add


  !***********************************************************************

  function global_real_max(r) result(s)
    implicit none

    real(kind=prec) :: r,s
    integer :: error

    call timer_start(2)
    call MPI_Allreduce(r,s,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)
    call timer_stop(2)


  end function global_real_max

  !***********************************************************************

  function global_real_min (r) result (s)
    implicit none

    real(kind=prec), intent(in)   :: r
    real(kind=prec) :: s
    integer :: error

    call timer_start(2)
    print *,"avant",r
    call MPI_Allreduce(r,s,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,error)
    s=r
    print *,"apres",r
    call timer_stop(2)

    return
  end function global_real_min

  !***********************************************************************
  
  function global_real_madd1d (r) result(s)

    implicit none

    real(kind=prec), dimension(:),      intent(in)  :: r
    real(kind=prec), dimension(size(r))             :: s

    integer i, nb, error

    call timer_start(2)
    nb = size(r)
    call MPI_Allreduce(r(1),s(1),nb,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)
    call timer_stop(2)

  end function global_real_madd1d

  !***********************************************************************

  function global_int_madd1d (a) result(b)

    implicit none

    integer, dimension(:),      intent(in)  :: a
    integer, dimension(size(a))             :: b
    real(kind=prec), dimension(size(a))     :: ra,rb

    call timer_start(2)
    ra=a
    rb=global_real_madd1d(ra)
    b=rb
    call timer_stop(2)

    return
  end function global_int_madd1d

  !***********************************************************************

  !***********************************************************************
  !   send and receives...
  !***********************************************************************

  !***********************************************************************

  subroutine snd_real_msg0d( task, tag, buf, nb)
    implicit none

    integer          :: task, tag
    real(kind=prec) :: buf
    integer          :: error,nb

    call timer_start(1)
    call MPI_send(buf,nb,MPI_DOUBLE_PRECISION,task,tag,MPI_COMM_WORLD,error)
    call timer_stop(1)

  end subroutine snd_real_msg0d

  !***********************************************************************

  subroutine rcv_real_msg0d( task, tag, buf,nb)
    implicit none

    integer          :: task, tag
    real(kind=prec) :: buf
    integer          :: error,nb,status
    
    call timer_start(1)
    call MPI_recv(buf,nb,MPI_DOUBLE_PRECISION,task,tag,MPI_COMM_WORLD,status,error)
!!$    nb_msg = nb_msg+1
!!$    if (nb_msg.ge.50) buf=0.
    call timer_stop(1)

    return
  end subroutine rcv_real_msg0d

  !***********************************************************************

  subroutine snd_real_msg1d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:) :: buf
    integer                        :: error,nb

    call timer_start(1)
    nb=size(buf)
    call MPI_send(buf(1),nb,MPI_DOUBLE_PRECISION,task,tag,MPI_COMM_WORLD,error)
    call timer_stop(1)
    
    return
  end subroutine snd_real_msg1d

  !***********************************************************************

  subroutine rcv_real_msg1d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:) :: buf
    integer                        :: error,nb,status
    
    call timer_start(1)
    nb=size(buf)
    call MPI_recv(buf(1),nb,MPI_DOUBLE_PRECISION,task,tag,MPI_COMM_WORLD,status,error)
!!$    nb_msg = nb_msg+1
!!$    if (nb_msg.ge.50) buf=0.
    call timer_stop(1)

    return
  end subroutine rcv_real_msg1d

  !***********************************************************************

  subroutine snd_real_msg2d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:,:) :: buf
    real(kind=prec), dimension(size(buf)) :: buf1d
    integer                        :: error,nb

    call timer_start(1)
    nb=size(buf)
    buf1d = reshape(buf,(/size(buf)/))
    !!print *,my_task,'sends ',buf1d
    call MPI_send(buf1d(1),nb,MPI_DOUBLE_PRECISION,task,tag,MPI_COMM_WORLD,error)
    call timer_stop(1)

    return
  end subroutine snd_real_msg2d

  !***********************************************************************

  subroutine rcv_real_msg2d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:,:) :: buf
    real(kind=prec), dimension(size(buf)) :: buf1d
    integer                        :: error,nb,status
    
    call timer_start(1)
    nb=size(buf)
    call MPI_recv(buf1d(1),nb,MPI_DOUBLE_PRECISION,task,tag,MPI_COMM_WORLD,status,error)
    buf = reshape(buf1d,(/size(buf,1),size(buf,2)/))
!!$    nb_msg = nb_msg+1
!!$    if (nb_msg.ge.50) buf=0.
    !!print *,my_task,'receives ',buf
    call timer_stop(1)

    return
  end subroutine rcv_real_msg2d

  !***********************************************************************

  subroutine snd_int_msg0d( task, tag, buf)
    implicit none

    integer          :: task, tag
    integer          :: buf
    real(kind=prec) :: bufr

    call timer_start(1)
    bufr=buf
    call snd_real_msg0d(task,tag,bufr,1)
    call timer_stop(1)

    return
  end subroutine snd_int_msg0d

  !***********************************************************************

  subroutine rcv_int_msg0d( task, tag, buf)
    implicit none

    integer          :: task, tag
    integer          :: buf
    real(kind=prec) :: bufr
    
    call timer_start(1)
    call rcv_real_msg0d(task,tag,bufr,1)
    buf=bufr
    call timer_stop(1)
    
    return
  end subroutine rcv_int_msg0d

  !***********************************************************************

  subroutine real_bcast(task,value)
    implicit none

    integer task
    real(kind=prec) value

    call timer_start(2)
    if (my_task.ne.task) value=0._prec
    value=global_real_add(value)
    call timer_stop(2)

    return
  end subroutine real_bcast



  subroutine int_bcast(task,value)
    implicit none

    integer task,value
    real(kind=prec) rvalue

    call timer_start(2)
    rvalue=value
    if (my_task.ne.task) rvalue=0._prec
    value=global_real_add(rvalue)
    call timer_stop(2)

    return
  end subroutine int_bcast


  !***********************************************************************
  !***********************************************************************
  !                                                                      *
  !     TIMERS                                                           *
  !                                                                      *
  !***********************************************************************
  !***********************************************************************

  subroutine timer_clear(n)
    implicit none

    integer n

    tchrono(n) = 0.
    is_on(n)   = .false.

  end subroutine timer_clear

  !***********************************************************************

  subroutine timer_start(n)
    implicit none

    integer n
    integer ic
    real(8) s
    real(8) :: second

    is_on(n)=.true.
!!    call get_time(s)
!!    call system_clock(ic)
    s=second()
    uchrono(n) = s
    is_on(n)   = .true.

  end subroutine timer_start


  !***********************************************************************

  subroutine timer_stop(n)
    implicit none

    integer n,ic
    real(8) s
    real(8) :: second

    if (is_on(n)) then
       s = second()
!!       call get_time(s)
       tchrono(n) = tchrono(n) + s - uchrono(n)
       is_on(n)=.false.
    endif

  end subroutine timer_stop

  !***********************************************************************

  subroutine timer_print(n)
    implicit none

    integer, intent(in) :: n

    if (my_task.eq.0) &
         & print '(A,I3,A,I3,A,F10.2,A,E10.5,A)&
         &',' task ',my_task&
         & ,' Timer ',n,' Cpu Time : ',tchrono(n),' seconds'
    

  end subroutine timer_print

  !***********************************************************************
  !***********************************************************************
  !                                                                      *
  !     MISCELLANEOUS                                                    *
  !                                                                      *
  !***********************************************************************
  !***********************************************************************

  subroutine synchronize_processors( )
    implicit none
    integer error

    call MPI_Barrier(MPI_COMM_WORLD)

  end subroutine synchronize_processors



  !***********************************************************************


end module uncol
