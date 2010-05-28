module selectprec
  integer, parameter :: prec=kind(1.e0)
end module selectprec



module uncol

  use selectprec

  integer, parameter :: max_tasks=128, max_msg_size=200 &
       , user_timers=6, system_timers=3                 &
       , all_timers=user_timers + system_timers
  real, save    ::  tchrono(0:all_timers), uchrono(0:all_timers)

  integer, save :: nb_tasks,my_task, master_task        &
       , my_tid, host_tid, tid( 0:max_tasks)

  logical, save :: power_of_two_tasks

  character*31, save ::  host_name, node_name

  integer, parameter :: mask_msg_type=256, add_tag=1*mask_msg_type  &
       , init_tag=2*mask_msg_type,  rel_tag=3*mask_msg_type         &
       , shift_tag=4*mask_msg_type, max_tag=5*mask_msg_type         &
       , size_real=8, size_integer=8                                &
       , PVMDEFAULT = 0, PVMRAW = 1, BYTE1 = 1, REAL8 = 6           &
       , PVMSHOWTIDS = 14                                           &
       , pvm_mode = PVMDEFAULT, pvm_data = PVMRAW                   &
       , pvm_real = REAL8, pvm_int = REAL8                          

  character*8 PVMALL
  common /PVMALL/ PVMALL


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
     module procedure global_int_add, global_real_add, global_real_madd1d
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
  !    Samuel Kortas & Bertrand Meltz
  !    version 2.3.0     xx Feb 1994
  !
  !***********************************************************************
  !  New in this version :
  !    - Initialization only for structured mesh and 2^p processors
  !***********************************************************************

!!$  subroutine start_network
!!$
!!$    !  main routine. It does nothing, except calling other routines.
!!$
!!$    implicit none
!!$
!!$    call initialize_processors
!!$
!!$
!!$    !*********************************
!!$    !     call your_program here     *
!!$    !*********************************
!!$    call programme
!!$    !*********************************
!!$
!!$    call release_processors
!!$
!!$  end subroutine start_network

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
    integer other_task, other_tid
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
!   Code for the 1
!    Make sure the task number is the same as the PE number
!     my_task = info
      call pvmfgetpe( my_tid, my_task)
  
!    Get the number of processors (specified as "-npes xx" in the
!    mppexec command line)
      call pvmfgsize( PVMALL, nb_tasks)
!    Make sure every task has finished starting
      call pvmfbarrier( PVMALL, nb_tasks, error)
      if (my_task .eq. 0) then
!**     Code for the master node        ***
         master_task = my_task
         host_tid  = my_tid
         tid( 0) = my_tid
!       Get the tids of all other tasks
         do i=1,nb_tasks-1
            call pvmfgettid( PVMALL, i, other_tid)
            call pvmfgetpe( other_tid, other_task)
            tid( other_task) = other_tid
         enddo
      endif
!   Back to common 1 section (1 or not)
      if (my_task .eq. 0 .and. nb_tasks .gt. 1) then
         call pvmfinitsend( pvm_data, error )
         call pvmfpack( pvm_int, master_task, 1, 1, error )
         call pvmfpack( pvm_int, nb_tasks, 1, 1, error )
         call pvmfpack( pvm_int, host_tid, 1, 1, error )
         call pvmfpack( pvm_int, tid(0), nb_tasks, 1, error )
         call pvmfmcast( nb_tasks - 1, tid(1), init_tag, error )
      else if (nb_tasks .gt. 1) then
!**     Code for a slave node            ***
         call pvmfrecv( -1, init_tag, error)
         call pvmfunpack( pvm_int, master_task, 1, 1, error )
         call pvmfunpack( pvm_int, nb_tasks, 1, 1, error )
         call pvmfunpack( pvm_int, host_tid, 1, 1, error )
         call pvmfunpack( pvm_int, tid(0), nb_tasks, 1, error )
      endif
      if (my_task .eq. master_task) then
         print *, ' Initialization successfully completed for ',nb_tasks, ' tasks'
         print *
      endif

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

!    call pvmflvgroup( 'global', error)

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

  end subroutine timer_clear

  !***********************************************************************

  subroutine timer_start(n)
    implicit none

    integer n
    integer irtc
    external irtc

    uchrono(n) = real(irtc())
      
  end subroutine timer_start


  !***********************************************************************

  subroutine timer_stop(n)
    implicit none

    integer n
    real(kind=prec) s
    integer irtc
    external irtc

    s = real(irtc())
!    call get_time(s)
    tchrono(n) = tchrono(n) + s - uchrono(n)

  end subroutine timer_stop

  !***********************************************************************

  subroutine timer_print(n)
    implicit none

    integer, intent(in) :: n
    real*8,  parameter  :: temps_cycle=3.333e-9

    if (my_task==0) &
         print '(A,I3,A,F15.5,A)',' task ',my_task,' Cpu Time : ',tchrono(n)*temps_cycle,' seconds'

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

    call pvmfbarrier( PVMALL, nb_tasks, error)

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

end module uncol
