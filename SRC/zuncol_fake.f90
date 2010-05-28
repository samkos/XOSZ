module selectprec
  integer, parameter :: prec=kind(1.d0)
end module selectprec

module uncol
  
  use selectprec

  integer, parameter :: max_tasks=128, max_msg_size=200 &
       , user_timers=6, system_timers=3                 &
       , all_timers=user_timers + system_timers
  real, save    ::  tchrono(0:all_timers), uchrono(0:all_timers)
  logical, save ::  is_on(0:all_timers)=.false.

  integer, save :: nb_tasks,my_task, master_task        &
       , my_tid, host_tid

  integer, save :: tid( 0:max_tasks)

  logical, save :: power_of_two_tasks

  character*76, save ::  host_name, node_name

  integer, parameter :: mask_msg_type=256, add_tag=1*mask_msg_type  &
       , init_tag=2*mask_msg_type,  rel_tag=3*mask_msg_type         &
       , shift_tag=4*mask_msg_type, max_tag=5*mask_msg_type         &
       , size_real=8, size_integer=8                                &
       , PVMDEFAULT = 0, PVMRAW = 1, BYTE1 = 1, REAL8 = 6           &
       , PVMSHOWTIDS = 14                                           &
       , pvm_mode = PVMDEFAULT, pvm_data = PVMRAW                   &
       , pvm_real = REAL8, pvm_int = REAL8                          



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

!!$    call pvmfcatchout(1,i)
!!$    call pvmfsetopt(PvmShowTids, 0 ,I )

    node_name = 'node_name_uncol'
    arch = '*'
    host_name = node_name

    nb_tasks=1
    my_task=0
    master_task=0

    if (my_task .eq. 0) then
       print *, ' Initialization successfully completed for ',nb_tasks, ' tasks'
    endif

    !  is nb_tasks a power of 2
    power_of_two_tasks = .false.
    do i=1,20
       if (2**i .eq. nb_tasks) power_of_two_tasks = .true.
    enddo



  end subroutine initialize_processors

  !***********************************************************************

  subroutine release_processors

    !  This routine terminates all tasks but the master_task, in case of
    !  PVM, or finishes the tracing, in case of PICL.

    implicit none

    integer i
    integer error

  end subroutine release_processors


  !***********************************************************************
  !   reduction operations
  !***********************************************************************


  function global_int_add(i) result (j)
    implicit none

    integer          :: i,j
    j=i

    return
  end function global_int_add

  !***********************************************************************

  function global_int_max(i) result (j)
    implicit none

    integer          :: i,j

    j=i
    
    return
  end function global_int_max

  !***********************************************************************

  function global_int_min(i) result (j)
    implicit none

    integer          :: i,j
    real(kind=prec) :: r,s

    j=i
    
    return
  end function global_int_min

  !***********************************************************************

  function global_real_add(r) result (s)
    implicit none

    real(kind=prec), intent(in)  :: r
    real(kind=prec) :: s
    
    s=r

    return
  end function global_real_add


  !***********************************************************************

  function global_real_max(r) result(s)
    implicit none

    real(kind=prec) :: r,s,a(1),b(1)
    
    s=r

    return
  end function global_real_max

  !***********************************************************************

  function global_real_min (r) result (a)
    implicit none

    real(kind=prec), intent(in)   :: r
    real(kind=prec) :: a

    a = r

    return
  end function global_real_min



  !***********************************************************************
  
  function global_real_madd1d (a) result(b)

    implicit none

    real(kind=prec), dimension(:),      intent(in)  :: a
    real(kind=prec), dimension(size(a))             :: b,c

    integer i, stride, other_task, tag

    
    b=a

  end function global_real_madd1d

  !***********************************************************************

  function global_int_madd1d (a) result(b)

    implicit none

    integer, dimension(:),      intent(in)  :: a
    integer, dimension(size(a))             :: b
    real(kind=prec), dimension(size(a))     :: ra,rb

    b=a

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

  end subroutine snd_real_msg0d

  !***********************************************************************

  subroutine rcv_real_msg0d( task, tag, buf,nb)
    implicit none

    integer          :: task, tag
    real(kind=prec) :: buf
    integer          :: error,nb
    
  end subroutine rcv_real_msg0d

  !***********************************************************************

  subroutine snd_real_msg1d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:) :: buf
    integer                        :: error,nb


  end subroutine snd_real_msg1d

  !***********************************************************************

  subroutine rcv_real_msg1d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:) :: buf
    integer                        :: error,nb
    

  end subroutine rcv_real_msg1d

  !***********************************************************************

  subroutine snd_real_msg2d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:,:) :: buf
    integer                        :: error,nb


  end subroutine snd_real_msg2d

  !***********************************************************************

  subroutine rcv_real_msg2d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:,:) :: buf
    integer                        :: error,nb
    

  end subroutine rcv_real_msg2d

  !***********************************************************************

  subroutine snd_int_msg0d( task, tag, buf)
    implicit none

    integer          :: task, tag
    integer          :: buf
    real(kind=prec) :: bufr


  end subroutine snd_int_msg0d

  !***********************************************************************

  subroutine rcv_int_msg0d( task, tag, buf)
    implicit none

    integer          :: task, tag
    integer          :: buf
    real(kind=prec) :: bufr
    

  end subroutine rcv_int_msg0d

  !***********************************************************************

  subroutine real_bcast(task,value)
    implicit none

    integer task
    real(kind=prec) value


    return
  end subroutine real_bcast



  subroutine int_bcast(task,value)
    implicit none

    integer task
    integer value

    value = value

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
    real(kind=prec) s

    real(8) :: second

    is_on(n)=.true.
    !    call get_time(s)
    !!    call system_clock(ic)
    s=second()
    uchrono(n) = s
    is_on(n)   = .true.

  end subroutine timer_start


  !***********************************************************************

  subroutine timer_stop(n)
    implicit none

    integer n,ic
    real(kind=prec) s
    
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
         & print '(A,I3,A,F10.2,A)',' task ',my_task&
         & ,' Cpu Time : ',tchrono(n),' seconds'
    

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



  end subroutine synchronize_processors

  !***********************************************************************
  !***********************************************************************

  subroutine display_error_msg( comment, error_code)
    implicit none
    character(len=*), intent(in) :: comment
    integer,          intent(in):: error_code


  end subroutine display_error_msg

  !***********************************************************************

!SP2      subroutine flush(n)
!SP2      return
!SP2      end subroutine flush

!SP2      subroutine get_time(s)
!SP2      implicit none
!SP2      real(kind=prec), intent(out) :: s
!SP2      integer ::  mclock
!SP2      external mclock 
!SP2      
!SP2      s=0.01*mclock()
!SP2      
!SP2      end subroutine get_time


!SGI      subroutine get_time(s)
!SGI      implicit none
!SGI      real(kind=prec), intent(out) :: s
!SGI      integer ::  mclock
!SGI      external mclock 
!SGI      
!SGI      s=0.01*mclock()
!SGI      
!SGI      end subroutine get_time


  !***********************************************************************


end module uncol
