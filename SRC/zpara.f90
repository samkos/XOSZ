
module para

  use uncol


  logical, save :: is_north, is_south, is_west, is_east, is_mpp

  integer, save ::  &
        p_north,  p_south,  p_west,  p_east &
       ,nb_i_blocks, nb_k_blocks,iline,icolumn&
       ,mpi_comm_line ,mpi_comm_column


  integer, parameter ::                                  &
       no_proc=-99,add_tag_local=78965                   &
       ,stn_n=8,stn_e=4,stn_w=2,stn_s=1                  &  
       ,stn_news=stn_n+stn_e+stn_w+stn_s                 &
       ,stn_ne=stn_n+stn_e,stn_ws=stn_w+stn_s            &
       ,stn_ns=stn_n+stn_s,stn_ew=stn_e+stn_w            & 
       ,stn_black=128,stn_red=64,stn_carre=32,stn_mos=16 &
       ,from_n=990,from_s=991,from_w=992,from_e=993

  
  logical, save  :: is_output
  integer, save  :: tag_outdon_nobord=0, tag_outdon_bord=0 &
       , tag_save_nobord=0, tag_output=0, nunit
  character(137), save  :: nom_file

  integer, parameter :: p_save=1, p_retrieve=2, p_save_txt=3

  interface outdon
     module procedure outdon_bord, outdon_nobord 
  end interface


  interface rfr_stn
     module procedure nonblocking_rfr_stn
  end interface

  interface rfr_ddm_stn
     module procedure nonblocking_rfr_ddm_stn
  end interface

contains

  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine init_para
  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
#ifndef SEQ
    use mpi
#else
#endif
    implicit none
    integer error

      nb_k_blocks=1
      nb_i_blocks=nb_tasks

      nb_k_blocks=nb_tasks
      nb_i_blocks=1

      do while (nb_k_blocks.gt.nb_i_blocks)
         nb_i_blocks=nb_i_blocks*2
         nb_k_blocks=nb_k_blocks/2
      enddo
      
!SK===================================================================
!SK  processeurs voisins   : no_proc = aucun
!SK===================================================================
!SK
      
      p_east = mod(my_task+1,nb_i_blocks)  &
           + int(my_task/nb_i_blocks)*nb_i_blocks
      p_west = mod(my_task+nb_i_blocks-1,nb_i_blocks) &
           + int(my_task/nb_i_blocks)*nb_i_blocks
      p_south = mod(my_task+nb_tasks-nb_i_blocks,nb_tasks)
      p_north = mod(my_task+nb_i_blocks,nb_tasks)
      
      if (mod(my_task,nb_i_blocks).eq.0) p_west=no_proc
      if (mod(my_task,nb_i_blocks).eq.nb_i_blocks-1) p_east=no_proc
      if (my_task.lt.nb_i_blocks) p_south=no_proc
      if (my_task.ge.nb_tasks-nb_i_blocks) p_north=no_proc

      is_west=p_west.eq.no_proc
      is_east=p_east.eq.no_proc
      is_north=p_north.eq.no_proc
      is_south=p_south.eq.no_proc


      icolumn=mod(my_task,nb_i_blocks)
      iline=my_task/nb_i_blocks

      is_mpp = (nb_tasks>1)
 
#ifndef SEQ
     call MPI_Comm_split( MPI_COMM_WORLD, iline, my_task, mpi_comm_line, error)
     call MPI_Comm_split( MPI_COMM_WORLD, icolumn, my_task, mpi_comm_column, error)
#else
#endif



    return
  end subroutine init_para


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine blocking_rfr_stn(INP,direction,nx0,ny0)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,0:)  :: INP
    integer          :: direction,nx,ny,nx0,ny0
    integer :: lm,nm

    
    lm=size(INP,1)-2*nx0; nm=size(INP,2)-2*ny0

    nx=nx0-1
    ny=ny0-1

    !!print *,'nx,ny',nx,ny

    if (iand(direction,stn_e)/=0.and.(.not.is_east))  &
         call snd_msg(p_east,from_w,INP(lm:nx+lm,:))

    if (iand(direction,stn_w)/=0.and.(.not.is_west))  &
         call snd_msg(p_west,from_e,INP(nx+1:2*nx+1,:))

    if (iand(direction,stn_w)/=0.and.(.not.is_west))  &
         call rcv_msg(p_west,from_w,INP(0:nx,:))

    if (iand(direction,stn_e)/=0.and.(.not.is_east))  &
         call rcv_msg(p_east,from_e,INP(nx+lm+1:2*nx+lm+1,:))

    if (iand(direction,stn_n)/=0.and.(.not.is_north)) &
         call snd_msg(p_north,from_s,INP(:,nm:ny+nm))

    if (iand(direction,stn_s)/=0.and.(.not.is_south)) &
         call snd_msg(p_south,from_n,INP(:,ny+1:2*ny+1))

    if (iand(direction,stn_s)/=0.and.(.not.is_south)) &
         call rcv_msg(p_south,from_s,INP(:,0:ny))

    if (iand(direction,stn_n)/=0.and.(.not.is_north)) &
         call rcv_msg(p_north,from_n,INP(:,ny+nm+1:2*ny+nm+1))


    return
  end subroutine blocking_rfr_stn


  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine nonblocking_rfr_stn(INP,direction,nx0,ny0)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
#ifndef SEQ
    use mpi
#endif
    implicit none
    real(kind=prec), dimension(0:,0:)  :: INP
    integer          :: direction,nx,ny,nx0,ny0
#ifndef SEQ
    integer :: lm,nm
    real(kind=prec), dimension(:), allocatable ::&
         & buf_to_north,buf_to_east,buf_to_west,buf_to_south,&
         & buf_from_north,buf_from_east,buf_from_west,buf_from_south
    integer, dimension(2) :: req_north,req_east,req_west,req_south
    integer, dimension(MPI_STATUS_SIZE,2) :: statuses
    integer error, status, nb
    integer, save :: tagplus=0,tag0
    tag0 = tagplus*1000
    tagplus = tagplus+1


    lm=size(INP,1)-2*nx0; nm=size(INP,2)-2*ny0

    nx=nx0-1
    ny=ny0-1

    !!print *,'nx,ny',nx,ny


    if (iand(direction,stn_n)/=0.and.(.not.is_north)) then
       nb = size(INP(:,ny+nm+1:2*ny+nm+1))
       allocate(buf_from_north (nb))
       allocate(buf_to_north (nb))
       call MPI_IRECV (buf_from_north(1),nb,MPI_DOUBLE_PRECISION,p_north,&
            & tag0+from_n,MPI_COMM_WORLD,req_north(1),error)
       !call rcv_msg(p_north,from_n,INP(:,ny+nm+1:2*ny+nm+1))

       buf_to_north  = reshape(INP(:,nm:ny+nm),(/nb/))
       call MPI_ISEND (buf_to_north(1),nb,MPI_DOUBLE_PRECISION,p_north,&
            & tag0+from_s,MPI_COMM_WORLD,req_north(2),error)
       !call snd_msg(p_north,from_s,INP(:,nm:ny+nm))
    end if

    if (iand(direction,stn_e)/=0.and.(.not.is_east))  then
       nb = size(INP(nx+lm+1:2*nx+lm+1,:))
       allocate(buf_from_east (nb))
       allocate(buf_to_east (nb))
       call MPI_IRECV (buf_from_east(1),nb,MPI_DOUBLE_PRECISION,p_east,&
            & tag0+from_e,MPI_COMM_WORLD,req_east(1),error)
       !call rcv_msg(p_east,from_e,INP(nx+lm+1:2*nx+lm+1,:))
       buf_to_east  = reshape(INP(lm:nx+lm,:),(/nb/))
       call MPI_ISEND (buf_to_east(1),nb,MPI_DOUBLE_PRECISION,p_east,&
            & tag0+from_w,MPI_COMM_WORLD,req_east(2),error)
       !call snd_msg(p_east,from_w,INP(lm:nx+lm,:))
    end if

    if (iand(direction,stn_w)/=0.and.(.not.is_west))  then
       nb = size(INP(0:nx,:))
       allocate(buf_from_west (nb))
       allocate(buf_to_west (nb))
       call MPI_IRECV (buf_from_west(1),nb,MPI_DOUBLE_PRECISION,p_west,&
            & tag0+from_w,MPI_COMM_WORLD,req_west(1),error)
       !call rcv_msg(p_west,from_w,INP(0:nx,:))
       buf_to_west  = reshape(INP(nx+1:2*nx+1,:),(/nb/))
       call MPI_ISEND (buf_to_west(1),nb,MPI_DOUBLE_PRECISION,p_west,&
            & tag0+from_e,MPI_COMM_WORLD,req_west(2),error)
       !call snd_msg(p_west,from_e,INP(nx+1:2*nx+1,:))
    end if

    if (iand(direction,stn_s)/=0.and.(.not.is_south)) then
       nb = size(INP(:,0:ny))
       allocate(buf_from_south (nb))
       allocate(buf_to_south (nb))
       call MPI_IRECV (buf_from_south(1),nb,MPI_DOUBLE_PRECISION,p_south,&
            & tag0+from_s,MPI_COMM_WORLD,req_south(1),error)
       !call rcv_msg(p_south,from_s,INP(:,0:ny))
       buf_to_south  = reshape(INP(:,ny+1:2*ny+1),(/nb/))
       call MPI_ISEND (buf_to_south(1),nb,MPI_DOUBLE_PRECISION,p_south,&
            & tag0+from_n,MPI_COMM_WORLD,req_south(2),error)
       !call snd_msg(p_south,from_n,INP(:,ny+1:2*ny+1))
    end if


    if (iand(direction,stn_n)/=0.and.(.not.is_north)) then
       call MPI_WAITALL(2, req_north, statuses, error)
       INP(:,ny+nm+1:2*ny+nm+1) = &
            &reshape(buf_from_north,(/size(INP(:,ny+nm+1:2*ny+nm+1),1),size(INP(:,ny+nm+1:2*ny+nm+1),2)/))
       deallocate(buf_to_north)
       deallocate(buf_from_north)       
    end if

    if (iand(direction,stn_e)/=0.and.(.not.is_east))  then
       call MPI_WAITALL(2, req_east, statuses, error)
       INP(nx+lm+1:2*nx+lm+1,:) = &
            &reshape(buf_from_east,(/size(INP(nx+lm+1:2*nx+lm+1,:),1),size(INP(nx+lm+1:2*nx+lm+1,:),2)/))
       deallocate(buf_to_east)
       deallocate(buf_from_east)       
    end if

    if (iand(direction,stn_w)/=0.and.(.not.is_west))  then
       call MPI_WAITALL(2, req_west, statuses, error)
       INP(0:nx,:) = &
            &reshape(buf_from_west,(/size(INP(0:nx,:),1),size(INP(0:nx,:),2)/))
       deallocate(buf_to_west)
       deallocate(buf_from_west)       
    end if

    if (iand(direction,stn_s)/=0.and.(.not.is_south)) then
       call MPI_WAITALL(2, req_south, statuses, error)
       INP(:,0:ny) = reshape(buf_from_south,(/size(INP(:,0:ny),1),size(INP(:,0:ny),2)/))
       deallocate(buf_to_south)
       deallocate(buf_from_south)       
    end if


#endif
    return
  end subroutine nonblocking_rfr_stn


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine nonblocking_rfr_ddm_stn(INP,ideb,ifin,kdeb,kfin,rec0x,rec0y)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK 
    !SK Avril 97 : nouvelle version a 2 directions grace a Robert Cimrman
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
#ifndef SEQ
    use mpi
#endif
    use debug
    implicit none
    real(kind=prec), dimension(0:,0:)  :: INP

    integer          :: ideb,ifin,kdeb,kfin,rec0x,rec0y,lm0,lm1,nm0,nm1 &
                        ,recx,recy

    real(kind=prec), dimension(:), allocatable ::&
         & buf_to_north,buf_to_east,buf_to_west,buf_to_south,&
         & buf_from_north,buf_from_east,buf_from_west,buf_from_south
    integer, dimension(2) :: req_north,req_east,req_west,req_south
    integer, dimension(MPI_STATUS_SIZE,2) :: statuses
    integer error, status, nb, tag0


    !debug_stn_ddm = .true.

    tag0 = 0
    if (debug_stn_ddm) print *,my_task,'in nonblocking_ddm_stn'

    lm0=ideb+1; lm1=ifin-1
    nm0=kdeb+1; nm1=kfin-1
    recx=rec0x-1
    recy=rec0y-1

    if (.not.is_east)  then
       nb = size(INP(lm1-recx:lm1,:))
       allocate(buf_from_east (nb))
       allocate(buf_to_east (nb))       
       call MPI_IRECV (buf_from_east(1),nb,MPI_DOUBLE_PRECISION,p_east,&
            & tag0+from_w,MPI_COMM_WORLD,req_east(1),error)
       buf_to_east  = reshape(INP(lm1-recx:lm1,:),(/nb/))
       call MPI_ISEND (buf_to_east(1),nb,MPI_DOUBLE_PRECISION,p_east,&
            & tag0+from_e,MPI_COMM_WORLD,req_east(2),error)
       !call snd_msg(p_east ,from_w,INP(lm1-recx:lm1,:))
    endif

    if (.not.is_west)  then
       nb = size(INP(lm0:lm0+recx,:))
       allocate(buf_from_west (nb))
       allocate(buf_to_west (nb))       
       call MPI_IRECV (buf_from_west(1),nb,MPI_DOUBLE_PRECISION,p_west,&
            & tag0+from_e,MPI_COMM_WORLD,req_west(1),error)
       buf_to_west  = reshape(INP(lm0:lm0+recx,:),(/nb/))
       call MPI_ISEND (buf_to_west(1),nb,MPI_DOUBLE_PRECISION,p_west,&
            & tag0+from_w,MPI_COMM_WORLD,req_west(2),error)
       !call snd_msg(p_west ,from_e,INP(lm0:lm0+recx,:))
    end if




    if (debug_stn_ddm) print *,my_task,'waiting nonblocking_ddm_stn'

    if (.not.is_west)  then
       !call rcv_msg(p_west ,from_w,INP(0:recx,:))
       call MPI_WAITALL(2, req_west, statuses, error)
       INP(0:recx,:)= &
            &reshape(buf_from_west,(/size(INP(0:recx,:),1),size(INP(0:recx,:),2)/))
       deallocate(buf_to_west)
       deallocate(buf_from_west)       
    end if

    if (.not.is_east) then
       !call rcv_msg(p_east ,from_e,INP(lm1+1:lm1+1+recx,:))
       call MPI_WAITALL(2, req_east, statuses, error)
       INP(lm1+1:lm1+1+recx,:)= &
            &reshape(buf_from_east,(/size(INP(lm1+1:lm1+1+recx,:),1),size(INP(lm1+1:lm1+1+recx,:),2)/))
       deallocate(buf_to_east)
       deallocate(buf_from_east)       
     end if

    if (debug_stn_ddm) print *,my_task,'excanging north/south nonblocking_ddm_stn'


    if (.not.is_north) then
       nb = size(INP(:,nm1-recy:nm1))
       allocate(buf_from_north (nb))
       allocate(buf_to_north (nb))       
       call MPI_IRECV (buf_from_north(1),nb,MPI_DOUBLE_PRECISION,p_north,&
            & tag0+from_s,MPI_COMM_WORLD,req_north(1),error)
       buf_to_north  = reshape(INP(:,nm1-recy:nm1),(/nb/))
       call MPI_ISEND (buf_to_north(1),nb,MPI_DOUBLE_PRECISION,p_north,&
            & tag0+from_n,MPI_COMM_WORLD,req_north(2),error)
       !call snd_msg(p_north,from_s,INP(:,nm1-recy:nm1))
    end if

    if (.not.is_south) then
       nb = size(INP(:,nm0:nm0+recy))
       allocate(buf_from_south (nb))
       allocate(buf_to_south (nb))       
       call MPI_IRECV (buf_from_south(1),nb,MPI_DOUBLE_PRECISION,p_south,&
            & tag0+from_n,MPI_COMM_WORLD,req_south(1),error)
       buf_to_south  = reshape(INP(:,nm0:nm0+recy),(/nb/))
       call MPI_ISEND (buf_to_south(1),nb,MPI_DOUBLE_PRECISION,p_south,&
            & tag0+from_s,MPI_COMM_WORLD,req_south(2),error)
       !call snd_msg(p_south,from_n,INP(:,nm0:nm0+recy))
    end if


    if (debug_stn_ddm) print *,my_task,'waiting nonblocking_ddm_stn'

    if (.not.is_south) then
       !call rcv_msg(p_south,from_s,INP(:,0:recy))
       call MPI_WAITALL(2, req_south, statuses, error)
       INP(:,0:recy)= &
            &reshape(buf_from_south,(/size(INP(:,0:recy),1),size(INP(:,0:recy),2)/))
       deallocate(buf_to_south)
       deallocate(buf_from_south)       
    end if

    if (.not.is_north) then
       !call rcv_msg(p_north,from_n,INP(:,nm1+1:nm1+1+recy))
       call MPI_WAITALL(2, req_north, statuses, error)
       INP(:,nm1+1:nm1+1+recy) = &
            &reshape(buf_from_north,(/size(INP(:,nm1+1:nm1+1+recy),1),size(INP(:,nm1+1:nm1+1+recy),2)/))
       deallocate(buf_to_north)
       deallocate(buf_from_north)       
    end if


    if (debug_stn_ddm) print *,my_task,'out nonblocking_ddm_stn'

   return
  end subroutine nonblocking_rfr_ddm_stn

  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine blocking_rfr_ddm_stn(INP,ideb,ifin,kdeb,kfin,rec0x,rec0y)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK 
    !SK Avril 97 : nouvelle version a 2 directions grace a Robert Cimrman
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    implicit none
    real(kind=prec), dimension(0:,0:)  :: INP
    integer          :: ideb,ifin,kdeb,kfin,rec0x,rec0y,lm0,lm1,nm0,nm1 &
                        ,recx,recy

    lm0=ideb+1; lm1=ifin-1
    nm0=kdeb+1; nm1=kfin-1
    recx=rec0x-1
    recy=rec0y-1

    if (.not.is_east)  call snd_msg(p_east ,from_w,INP(lm1-recx:lm1,:))

    if (.not.is_west)  call snd_msg(p_west ,from_e,INP(lm0:lm0+recx,:))

    if (.not.is_west)  call rcv_msg(p_west ,from_w,INP(0:recx,:))

    if (.not.is_east)  call rcv_msg(p_east ,from_e,INP(lm1+1:lm1+1+recx,:))

    if (.not.is_north) call snd_msg(p_north,from_s,INP(:,nm1-recy:nm1))

    if (.not.is_south) call snd_msg(p_south,from_n,INP(:,nm0:nm0+recy))

    if (.not.is_south) call rcv_msg(p_south,from_s,INP(:,0:recy))

    if (.not.is_north) call rcv_msg(p_north,from_n,INP(:,nm1+1:nm1+1+recy))


    return
  end subroutine blocking_rfr_ddm_stn


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine rfr_ddm_stn_old(INP,ideb,ifin,kdeb,kfin,rec0)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,0:)  :: INP
    integer          :: ideb,ifin,kdeb,kfin,rec0,lm0,lm1,nm0,nm1,rec

    lm0=ideb+1; lm1=ifin-1
    nm0=kdeb+1; nm1=kfin-1
    rec=rec0-1

    if (.not.is_east)  call snd_msg(p_east ,from_w,INP(lm1-rec:lm1,:))

    if (.not.is_west)  call snd_msg(p_west ,from_e,INP(lm0:lm0+rec,:))

    if (.not.is_west)  call rcv_msg(p_west ,from_w,INP(0:rec,:))

    if (.not.is_east)  call rcv_msg(p_east ,from_e,INP(lm1+1:lm1+1+rec,:))

    if (.not.is_north) call snd_msg(p_north,from_s,INP(:,nm1-rec:nm1))

    if (.not.is_south) call snd_msg(p_south,from_n,INP(:,nm0:nm0+rec))

    if (.not.is_south) call rcv_msg(p_south,from_s,INP(:,0:rec))

    if (.not.is_north) call rcv_msg(p_north,from_n,INP(:,nm1+1:nm1+1+rec))


    return
  end subroutine rfr_ddm_stn_old

!!$
  !************************************************************************

  function global_maxval(TAB) result(r)
    implicit none
    real(kind=prec), dimension(:,:) :: TAB
    real(kind=prec) :: r

    r=global_max(maxval(TAB))
    
    return
  end function global_maxval

  !************************************************************************

  function global_sum(TAB) result(r)
    implicit none
    real(kind=prec), dimension(:,:) :: TAB
    real(kind=prec) :: r

    r=global_add(sum(TAB))
    
    return
  end function global_sum


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine outdon_nobord(solution,chaine,nx,ny,is_int)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only : ntype_solver,ntype_precond,is_print,nom_fic_output
    implicit none
    real(kind=prec), dimension(0:,0:)            :: solution
    real(kind=prec), dimension(:,:), allocatable :: buffer,buffer2
    character(len=*)  :: chaine
    integer           :: nx,ny
    logical, optional :: is_int
    integer           :: i,k,task,i0,k0,i1,k1,lm,nm,lm_all,nm_all,ok
    integer, dimension(0:nb_tasks-1) :: lmtask,nmtask

!SP2
     return

    lm=size(solution,1)-2
    nm=size(solution,2)-2
    
    if (my_task==0.and.tag_outdon_nobord==0) then
       open (file=trim(nom_fic_output)//'outdon_save0' &
            ,unit=76,form='unformatted',iostat=ok)
       if (ok/=0) stop 'pb ouverture fichier save'  
    endif

    if (.not.is_mpp) then
       if (is_output) then
          print *,chaine,' sur ',nb_tasks,' proc'
          do k=nm+1,0,-1
          if (present(is_int)) then
             print '(100I2)',int(solution(:,k))
          else
             print '(100E11.3)',solution(:,k)
          endif
          enddo
       endif
       write(76) lm+2,nm+2
       write(76) solution
       ! print *,'enr : lm_global,nm_global ',lm+2,nm+2
       print *,'   ',chaine,' sauvegarde...'
    else
       i0=0; k0=0; i1=lm+1; k1=nm+1
       if (.not.is_west)  i0=nx
       if (.not.is_east)  i1=lm+1-nx
       if (.not.is_south) k0=ny
       if (.not.is_north) k1=nm+1-ny
       if (my_task>0) then
          call snd_msg(0,tag_outdon_nobord+1000*my_task,i1-i0+1)
          call snd_msg(0,tag_outdon_nobord+3000*my_task,k1-k0+1)
          call snd_msg(0,tag_outdon_nobord+100*my_task,solution(i0:i1,k0:k1))
       else
          do k=0,nb_k_blocks-1
             do i=0,nb_i_blocks-1
                if ((i+k)==0) then
                   lmtask(0)=i1-i0+1; nmtask(0)=k1-k0+1
                   cycle
                endif
                task=k*nb_i_blocks+i
                call rcv_msg(task,tag_outdon_nobord+1000*task,lmtask(task))
                call rcv_msg(task,tag_outdon_nobord+3000*task,nmtask(task))
             enddo
          enddo
          lm_all=sum(lmtask(0:nb_i_blocks-1)); 
          nm_all=sum(nmtask(0:nb_tasks-1:nb_i_blocks)); 

          allocate(buffer (0:lm_all-1,0:nm_all-1),stat=ok); 
          if (ok/=0) stop 'buffer : error alloc'
          allocate(buffer2(0:lm_all-1,0:nm_all-1),stat=ok); 
          if (ok/=0) stop 'buffer2 : error alloc'

          buffer(0:lm+1,0:nm+1)=solution
          k0=0
          do k=0,nb_k_blocks-1
             i0=0
             do i=0,nb_i_blocks-1
                task=k*nb_i_blocks+i
                if ((i+k)>0) call rcv_msg(task,tag_outdon_nobord+100*task &
                     ,buffer(i0:i0+lmtask(task)-1,k0:k0+nmtask(task)-1))
                i0=i0+lmtask(task)                
             enddo
             k0=k0+nmtask(task)
          enddo

          read(76) i,k
          if (i/=size(buffer,1).or.k/=size(buffer,2)) then
             print *,'pour ',chaine&
                  & ,' : tableau non conformant avec tableau sauve dans outdon'
             print *,'ici : lm,nm,lm_global,nm_global '&
                  & ,lm,nm,size(buffer,1),size(buffer,2)
             print *,'enr : lm_global,nm_global ',i,k
             stop
          endif
          read(76) buffer2
          buffer2= abs(buffer2-buffer)

          if (is_output) then
             print *,chaine,' sur ',nb_tasks,' procs'
             do k=k0-1,0,-1
             if (present(is_int)) then
                print '(100I2)',int(buffer(:,k))
             else
                print '(100E11.3)',buffer(:,k)
             endif
             enddo
             print *,'Differences ...'
             do k=k0-1,0,-1
                if (present(is_int)) then
                   print '(100I2)',int(buffer2(:,k))
                else
                   print '(100E11.3)',buffer2(:,k)
                endif
             enddo
          else
             buffer2=-log(max(buffer2,9.*(0.1_prec)**10))/log(10._prec)
             print *,chaine,'difference de ',nb_tasks,' procs  a 1 proc'
             do k=k0-1,0,-1
                print '(3X,100I1)',(int(buffer2(:,k)))
             enddo
          endif

          deallocate(buffer2,stat=ok); if (ok/=0) stop 'buffer2 : error desalloc'
          deallocate(buffer,stat=ok);  if (ok/=0) stop 'buffer  : error desalloc'
       endif
    endif


    tag_outdon_nobord=tag_outdon_nobord+1

    return
  end subroutine outdon_nobord

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine outdon_bord(solution,chaine,is_int)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,0:)            :: solution
    real(kind=prec), dimension(:,:), allocatable :: buffer
    character(len=*)                   :: chaine
    integer                            :: i,k,task,i0,k0,lm,nm,ok,lm_all,nm_all
    integer, dimension(0:nb_tasks-1)   :: lmtask,nmtask
    logical, optional                  :: is_int
    
    
!SP2  
   return

    lm=size(solution,1)-2
    nm=size(solution,2)-2

    if (.not.is_mpp) then
       print *,chaine,' sur ',nb_tasks,' proc'
       do k=nm+1,0,-1
          if (present(is_int)) then
             print '(100I2)',int(solution(:,k))
          else
             print '(100E11.3)',solution(:,k)
          endif
       enddo
    else
       if (my_task>0) then
          call snd_msg(0,tag_outdon_bord+1000*my_task,lm+2)
          call snd_msg(0,tag_outdon_bord+3000*my_task,nm+2)
          call snd_msg(0,tag_outdon_bord+100*my_task,solution)
       else
          do k=0,nb_k_blocks-1
             do i=0,nb_i_blocks-1
                if ((i+k)==0) then
                   lmtask(0)=lm+2; nmtask(0)=nm+2;
                   cycle
                endif
                task=k*nb_i_blocks+i
                call rcv_msg(task,tag_outdon_bord+1000*task,lmtask(task))
                call rcv_msg(task,tag_outdon_bord+3000*task,nmtask(task))
             enddo
          enddo
          lm_all=sum(lmtask(0:nb_i_blocks-1)); 
          nm_all=sum(nmtask(0:nb_tasks-1:nb_i_blocks)); 

          allocate(buffer(0:lm_all-1,0:nm_all-1),stat=ok); 
          if (ok/=0) stop 'error alloc'
          buffer(0:lm+1,0:nm+1)=solution

          k0=0
          do k=0,nb_k_blocks-1
             i0=0
             do i=0,nb_i_blocks-1
                task=k*nb_i_blocks+i
                if ((i+k)>0) call rcv_msg(task,tag_outdon_bord+100*task&
                     ,buffer(i0:i0+lmtask(task)-1,k0:k0+nmtask(task)-1))
                i0=i0+lmtask(task)
             enddo
             k0=k0+nmtask(task)
          enddo

          print *,chaine,' sur ',nb_tasks,' proc'
          do k=(size(solution,2))*nb_k_blocks-1,0,-1
             if (present(is_int)) then
                print '(100I2)',int(buffer(:,k))
             else
                print '(100E11.3)',buffer(:,k)
             endif
          enddo

          deallocate(buffer,stat=ok); if (ok/=0) stop 'error alloc'
       endif
       
       tag_outdon_bord=tag_outdon_bord+1
    
    endif

    return
  end subroutine outdon_bord


  !************************************************************************
  
  subroutine output(s,i1,i2,i3,i4,i5,r1,r2,r3,r4,l1,l2,l3,l4)
    use data,  only : nom_fic_output
    implicit none
    character(len=*), optional :: s
    integer, optional :: i1,i2,i3,i4,i5
    real(kind=prec), optional :: r1,r2,r3,r4
    logical, optional :: l1,l2,l3,l4
    character(6) :: sformat
    integer :: ok

!SP2
    return
!T3D    return

    if (tag_output==0) then
       sformat='(A,I1)'
       if (my_task>9) sformat='(A,I2)'
       write (nom_file,sformat) trim(nom_fic_output)//'output0',my_task
       nunit=10+my_task
       open(file=nom_file,unit=nunit,iostat=ok)
       if (ok/=0) then
          print *,'pb ouverture fichier',nom_file
          stop
       endif
       write (nunit,'(A,i2,A)') 'proc ',my_task &
            ,' ==========================================================='
    else
       open (file=nom_file,unit=nunit,position='append',iostat=ok)
       if (ok/=0) stop 'pb reouverture fichier'  
    endif
    
    write (nunit,'(a,i2,A,I4,A)',advance='no')&
         & 'proc ',my_task,' output ',tag_output,' ---> '
    if (present(s))  write (nunit,'(1X,A)',advance='no') s
    if (present(i1)) write (nunit,'(1X,I5)'   ,advance='no') i1
    if (present(i2)) write (nunit,'(1X,I5)'   ,advance='no') i2
    if (present(i3)) write (nunit,'(1X,I5)'   ,advance='no') i3
    if (present(i4)) write (nunit,'(1X,I5)'   ,advance='no') i4
    if (present(i5)) write (nunit,'(1X,I5)'   ,advance='no') i5
    if (present(r1)) write (nunit,'(1X,E11.3)',advance='no') r1
    if (present(r2)) write (nunit,'(1X,E11.3)',advance='no') r2
    if (present(r3)) write (nunit,'(1X,E11.3)',advance='no') r3
    if (present(r4)) write (nunit,'(1X,E11.3)',advance='no') r4
    if (present(l1)) write (nunit,'(1X,L1)'    ,advance='no') l1
    if (present(l2)) write (nunit,'(1X,L1)'    ,advance='no') l2
    if (present(l3)) write (nunit,'(1X,L1)'    ,advance='no') l3
    if (present(l4)) write (nunit,'(1X,L1)'    ,advance='no') l4
    write (nunit,'(1X)')
    call flush(nunit)
    close(nunit)
    
       tag_output=tag_output+1

    return
  end subroutine output

  !************************************************************************  ! start_out_light

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine save_or_retrieve(nom_solution,solution,nx,ny,unit,flag)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data
    use debug
    implicit none
#ifndef SEQ
    include 'mpif.h' 
    real(kind=prec), dimension(0:,0:)            :: solution
    real(kind=prec), dimension(:,:), allocatable :: solution_buf
    character (len=*) :: nom_solution
    integer           :: nx,ny,unit,nsize
    integer           :: flag
    logical           :: is_save, is_save_txt

    integer           :: i,k,task,i0,k0,i1,k1,k2,k3,k4,lm,nm,lm_all,nm_all,ok
    integer :: i0_task,i1_task,k0_task,k1_task,lm0,nm0,i4
    logical :: is_task_north,is_task_east,is_task_west,is_task_south
    integer :: p_task_north,p_task_east,p_task_west,p_task_south
    integer, dimension(0:nb_tasks-1) :: lmtask,nmtask
    integer, save :: filetype
    integer :: ierror
    integer :: sizes(2), subsizes(2), starts(2) 
    integer(kind=MPI_OFFSET_KIND) :: disp

    debug_save=.false.
    check_save = .false.

    if (debug_save) then
       print *,my_task,' in save_or_retrieve ',nom_solution,size(solution,1),size(solution,2)
       call flush(6)
    endif

    is_save     = (flag==p_save).or.(flag==p_save_txt)
    is_save_txt = (flag==p_save_txt)

    !print *,my_task, is_save,is_mpp
    call flush(6)

    if (.not.is_save) solution=0._prec


    !print *,'global',lm_global,nm_global,nb_i_blocks,nb_k_blocks
    sizes(1) = lm_global+2
    sizes(2) = nm_global+2
    lm0=(lm_global)/nb_i_blocks
    nm0=(nm_global)/nb_k_blocks
    


    lm=size(solution,1)-2
    nm=size(solution,2)-2

    i0=0; k0=0; i1=lm+1; k1=nm+1
    if (.not.is_west)  i0=nx
    if (.not.is_east)  i1=lm+1-nx
    if (.not.is_south) k0=ny
    if (.not.is_north) k1=nm+1-ny


    allocate(solution_buf(i0:i1,k0:k1))
    subsizes(1) = i1-i0+1
    subsizes(2) = k1-k0+1
        


   do k=0,nb_k_blocks-1
      do i=0,nb_i_blocks-1
         task=k*nb_i_blocks+i

         if (task.eq.my_task) then
            ! getting rid of the ghost cells that need not to be saved
            starts(1) = (lm0+1)*i
            starts(2) = (nm0+1)*k
            subsizes(1) = i1-i0+1
            subsizes(2) = k1-k0+1
            if (.not.is_west)  then
               starts(1) = starts(1)+1
            end if
            if (.not.is_south) then
               starts(2) = starts(2)+1
            end if
            if (check_save) then
               print '(13(I4,X),4(L,X))',my_task,sizes(1),sizes(2),subsizes(1),subsizes(2),starts(1),starts(2),&
                    & i0,i1,k0,k1,size(solution_buf,1),size(solution_buf,2),is_north,is_east,is_west,is_south
               call flush(6)
            end if
         end if
      enddo
   enddo

   
   call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, & 
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,       & 
        filetype, ierror)
   call MPI_TYPE_COMMIT(filetype,ierror)

   if (unit.gt.90) then
      unit = unit-100
      disp = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_VIEW(unit, disp, MPI_DOUBLE_PRECISION, & 
           filetype, 'native',   MPI_INFO_NULL, ierror) 
      if (check_save) print *,"herrrrrrrre"
   end if

   if (is_save) then
      if (check_save) then
         do k4=0,nm_global+1
            do i=0,lm_global+1
               solution(i,k4) = i+k4*100
            end  do
         end do
         print *,lm_global,nm_global
         do k4=nm_global+1,0,-1
            print '(I5,"b",I2,"->",100(F6.0,X))',my_task,k4,(solution(i,k4),i=0,lm_global+1)
         end do
      end if
          
      solution_buf(i0:i1,k0:k1) = solution(i0:i1,k0:k1)
      call MPI_FILE_WRITE_ALL(unit, solution_buf, subsizes(1)*subsizes(2), MPI_DOUBLE_PRECISION, & 
           MPI_STATUS_IGNORE, ierror) 

      
   else
      call MPI_FILE_READ_ALL(unit, solution_buf(i0:i1,k0:k1), subsizes(1)*subsizes(2), MPI_DOUBLE_PRECISION, & 
               MPI_STATUS_IGNORE, ierror) 
      if (check_save) then
         print *,"xxxxxxx",my_task,i0,i1,k0,k1
         do k=k1,k0,-1
            print '(I5,"x",I2,"->",100(F6.0,X))',my_task,k,(solution_buf(i,k),i=i0,i1)
         end do
      end if
       solution(i0:i1,k0:k1) = solution_buf(i0:i1,k0:k1)
   end if
       

   deallocate(solution_buf)

   if (debug_save) then
      print *,my_task,' out of  save_or_retrieve ',nom_solution
      call flush(6)
   endif

   if (check_save) then
      do k=size(solution,2)-1,0,-1
          print '(I5,"o",I2,"->",100(F6.0,X))',my_task,k,(solution(i,k),i=0,size(solution,1)-1)
       end do
       call parallel_stop
    end if

   call MPI_TYPE_FREE(filetype,ierror)
#else
  write (*,*) "save_or_retrieved not implemented in sequentiel"
#endif

    return
  end subroutine save_or_retrieve                                        ! end_out_light


  !************************************************************************

  function global_minval(INP) result(a)
    implicit none
    real(kind=prec), dimension(0:,0:)  :: INP
    real(kind=prec) :: a
    integer :: lm,nm,i0,i1,k0,k1
    
    lm=size(INP,1)-2; nm=size(INP,2)-2
    i0=1; k0=1; i1=lm; k1=nm
    
    if (is_west)  i0=0
    if (is_east)  i1=lm+1
    if (is_south) k0=0
    if (is_north) k1=nm+1

    a=global_min(minval(INP(i0:i1,k0:k1)))

    return
  end function global_minval

  !************************************************************************

  subroutine parallel_stop(s)
    implicit none
    character(len=*), optional :: s

    if (present(s)) print *,s
    print *
    print *,'*************** PARALLEL STOP CALLED ***************'
    print *
    call flush(6)

    print *,'releasing network...'
    call release_processors
    print *,'network free'

    stop
    return
  end subroutine parallel_stop

  !************************************************************************

end module para

