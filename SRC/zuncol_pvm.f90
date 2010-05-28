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

    call pvmfmytid( my_tid)
    if ( my_tid .lt. 0 ) then
       print *, 'cannot start main program ; error code :', my_tid
       call display_error_msg( ' pvmfmytid ', my_tid)
    endif

    call pvmfjoingroup( 'global', error)
    if (error .lt. 0) then
       print *, ' Error in join_group ; error code :', error
       call display_error_msg( ' pvmfjoingroup ', error)
    endif

    !   Code for a general purpose (non-T3D) PVM
    call pvmfparent( info)
    if ( info .lt. 0 ) then
       !**     Code for the master node        ***
       my_task = 0
       master_task = my_task
       host_tid = my_tid
       tid( 0) = my_tid

       print *, ' number of tasks :'
       read  *, nb_tasks
       print *, '   ', nb_tasks

       !       Initiate nb_tasks instances of node program 
       !       Remember that the numbering starts at 0
       do i=1,nb_tasks-1
          !          If arch is set to "*" then ANY type of machine is good
          !          otherwise arch should be set to architecture type you use.
          call pvmfspawn( node_name, pvm_mode, arch,1,tid(i),error)
          if ( error .le. 0) then
             print *, 'cannot start node ', i   &
                              ,', error code  :', error
             call display_error_msg( ' pvmfspawn ', error)
          endif
       enddo

       !   Back to common 1  section (T3D or not)

       do i=1,nb_tasks-1
          call pvmfinitsend( pvm_data, error )
          call pvmfpack( pvm_int, master_task, 1, 1, error )
          call pvmfpack( pvm_int, nb_tasks, 1, 1, error )
          call pvmfpack( pvm_int, host_tid, 1, 1, error )
          call pvmfpack( pvm_int, tid(i), 1, 1, error )
          call pvmfpack( pvm_int, tid, nb_tasks, 1, error )
          call pvmfsend(tid(i), init_tag, error )
       enddo
       
    else
       !**     Code for a slave node            ***

       call pvmfrecv( -1, init_tag, error)
       call pvmfunpack( pvm_int, master_task, 1, 1, error )
       call pvmfunpack( pvm_int, nb_tasks, 1, 1, error )
       call pvmfunpack( pvm_int, host_tid, 1, 1, error )
       call pvmfunpack( pvm_int, my_tid, 1, 1, error )
       call pvmfunpack( pvm_int, tid, nb_tasks, 1, error )

       !       What is my task number
       do i=1, nb_tasks-1
          if ( tid(i) .eq. my_tid ) my_task = i
       enddo

    endif

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

    if (my_task .eq. master_task) then
       !       receive notice of end from the slave nodes
       do  i=1,nb_tasks-1
          call pvmfrecv( tid(i), rel_tag, error)
       enddo
    else
       !       send notice of end to the master node
       call pvmfinitsend( pvm_data, error )
       call pvmfsend( tid(0), rel_tag, error )
    endif
    !    program finished, leave PVM before exiting 

    call pvmflvgroup( 'global', error)

    call pvmfexit( error)
  end subroutine release_processors


  !***********************************************************************
  !   reduction operations
  !***********************************************************************


  function global_int_add(i) result (j)
    implicit none

    integer          :: i,j
    real(kind=prec) :: r,s

    r=real(i)
    s=global_real_add(r)
    j=int(s)
    
    return
  end function global_int_add

  !***********************************************************************

  function global_int_max(i) result (j)
    implicit none

    integer          :: i,j
    real(kind=prec) :: r,s

    r=real(i)
    s=global_real_max(r)
    j=int(s)
    
    return
  end function global_int_max

  !***********************************************************************

  function global_int_min(i) result (j)
    implicit none

    integer          :: i,j
    real(kind=prec) :: r,s

    r=real(i)
    s=global_real_min(r)
    j=int(s)
    
    return
  end function global_int_min

  !***********************************************************************

  function global_real_add(r) result (s)
    implicit none

    real(kind=prec), intent(in)  :: r
    real(kind=prec) :: s
    
    real(kind=prec) :: a(1),b(1)
    integer i, stride, other_task, tag

    a(1) = r
    if ( power_of_two_tasks ) then
       i = 0
       do while (2**i<nb_tasks)
          stride = 2**i
          if (mod(my_task, stride*2) .ge. stride) stride = -stride
          other_task = my_task + stride
          tag = add_tag + i
          call snd_msg( other_task, tag, a)
          call rcv_msg( other_task, tag, b)
          a(1) = b(1) + a(1)
          i = i + 1
       enddo
    else
       if (nb_tasks.ne.1) stop 'global operation not implemented for nb_tasks<>2^p'
    endif

    s=a(1)

    return
  end function global_real_add


  !***********************************************************************

  function global_real_max(r) result(s)
    implicit none

    real(kind=prec) :: r,s,a(1),b(1)
    integer i, stride, other_task, tag

    a(1) = r
    if ( power_of_two_tasks ) then

       i = 0
       do while (2**i<nb_tasks) 
          stride = 2**i
          if (mod(my_task, stride*2) .ge. stride) stride = -stride
          other_task = my_task + stride
          tag = add_tag + i
          call snd_msg( other_task, tag, a)
          call rcv_msg( other_task, tag, b)
          a(1) = max(b(1),a(1))
          i    = i + 1
       enddo
    else
       if (nb_tasks.ne.1) stop 'global operation not implemented for nb_tasks<>2^p'
    endif

    s=a(1)

    return
  end function global_real_max

  !***********************************************************************

  function global_real_min (r) result (a)
    implicit none

    real(kind=prec), intent(in)   :: r
    real(kind=prec) :: a

    a = - global_real_max(-r)

    return
  end function global_real_min



  !***********************************************************************
  
  function global_real_madd1d (a) result(b)

    implicit none

    real(kind=prec), dimension(:),      intent(in)  :: a
    real(kind=prec), dimension(size(a))             :: b,c

    integer i, stride, other_task, tag

    
    b=a
    if ( power_of_two_tasks ) then

       i = 0
       do while (2**i<nb_tasks)
          stride = 2**i
          if (mod(my_task, stride*2) .ge. stride) stride = -stride
          other_task = my_task + stride
          tag = add_tag + i
          call snd_msg( other_task, tag, b)
          call rcv_msg( other_task, tag, c)
          b = b+c
          i = i + 1
       enddo
    else
       if (nb_tasks.ne.1) stop 'global operation not implemented for nb_tasks<>2^p'
    endif

  end function global_real_madd1d

  !***********************************************************************

  function global_int_madd1d (a) result(b)

    implicit none

    integer, dimension(:),      intent(in)  :: a
    integer, dimension(size(a))             :: b
    real(kind=prec), dimension(size(a))     :: ra,rb

    ra=a
    rb=global_real_madd1d(ra)
    b=rb

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

    call pvmfinitsend( pvm_data, error)
    call pvmfpack( pvm_real, buf, nb,1, error)
    call pvmfsend( tid(task), tag, error)

  end subroutine snd_real_msg0d

  !***********************************************************************

  subroutine rcv_real_msg0d( task, tag, buf,nb)
    implicit none

    integer          :: task, tag
    real(kind=prec) :: buf
    integer          :: error,nb
    
    call pvmfrecv ( tid(task), tag, error)
    call pvmfunpack( pvm_real, buf, nb, 1, error)

  end subroutine rcv_real_msg0d

  !***********************************************************************

  subroutine snd_real_msg1d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:) :: buf
    integer                        :: error,nb

    nb=size(buf)
    call pvmfinitsend( pvm_data, error)
    call pvmfpack( pvm_real, buf, nb,1, error)
    call pvmfsend( tid(task), tag, error)

  end subroutine snd_real_msg1d

  !***********************************************************************

  subroutine rcv_real_msg1d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:) :: buf
    integer                        :: error,nb
    
    nb=size(buf)
    call pvmfrecv ( tid(task), tag, error)
    call pvmfunpack( pvm_real, buf, nb, 1, error)

  end subroutine rcv_real_msg1d

  !***********************************************************************

  subroutine snd_real_msg2d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:,:) :: buf
    integer                        :: error,nb

    nb=size(buf)
    call pvmfinitsend( pvm_data, error)
    call pvmfpack( pvm_real, buf, nb,1, error)
    call pvmfsend( tid(task), tag, error)

  end subroutine snd_real_msg2d

  !***********************************************************************

  subroutine rcv_real_msg2d( task, tag, buf)
    implicit none

    integer                        :: task, tag
    real(kind=prec), dimension(:,:) :: buf
    integer                        :: error,nb
    
    nb=size(buf)
    call pvmfrecv ( tid(task), tag, error)
    call pvmfunpack( pvm_real, buf, nb, 1, error)

  end subroutine rcv_real_msg2d

  !***********************************************************************

  subroutine snd_int_msg0d( task, tag, buf)
    implicit none

    integer          :: task, tag
    integer          :: buf
    real(kind=prec) :: bufr

    bufr=buf
    call snd_real_msg0d(task,tag,bufr,1)

  end subroutine snd_int_msg0d

  !***********************************************************************

  subroutine rcv_int_msg0d( task, tag, buf)
    implicit none

    integer          :: task, tag
    integer          :: buf
    real(kind=prec) :: bufr
    
    call rcv_real_msg0d(task,tag,bufr,1)
    buf=bufr

  end subroutine rcv_int_msg0d

  !***********************************************************************

  subroutine real_bcast(task,value)
    implicit none

    integer task
    real(kind=prec) value

    if (my_task.ne.task) value=0._prec
    value=global_real_add(value)

    return
  end subroutine real_bcast



  subroutine int_bcast(task,value)
    implicit none

    integer task,value
    real(kind=prec) rvalue

    rvalue=value
    if (my_task.ne.task) rvalue=0._prec
    value=global_real_add(rvalue)

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
    real(kind=prec) s

    is_on(n)=.true.
    call get_time(s)
    uchrono(n) = s
    is_on(n)   = .true.

  end subroutine timer_start


  !***********************************************************************

  subroutine timer_stop(n)
    implicit none

    integer n
    real(kind=prec) s

    if (is_on(n)) then
       call get_time(s)
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

    call pvmfbarrier( 'global', nb_tasks, error)

  end subroutine synchronize_processors

  !***********************************************************************
  !***********************************************************************

  subroutine display_error_msg( comment, error_code)
    implicit none
    character(len=*), intent(in) :: comment
    integer,          intent(in):: error_code

    if ( error_code .lt. 0) then
       call pvmfperror( comment, error_code)
    endif

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
