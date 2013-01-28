module dac

  use constante
  use para
  use uncol

  integer, parameter ::              &
        iu_a=0,iu_b=1,iu_c=2,iu_d=0  &
       ,id_a=3,id_b=4,id_c=5,id_d=1

  type compact_part_matrix
     real(kind=prec), dimension (:,:), pointer :: P,Pi,Q,Qi
  end type compact_part_matrix

  type (compact_part_matrix), dimension(:), allocatable :: &
       p_df4x,p_df4y,p_vf4x,p_vf4y     ! coefficients partitiones

  integer, save :: tag_dacx=0,tag_dacy=1000

contains

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine initialise_df4_part_coefx(lmg,lm)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc,     only : dx
    implicit none
    integer :: lm,lmg,ok

    allocate(p_df4x(lmg)%P (0:lm+1,-1:3)     ,stat=ok);   if (ok/=0) stop 'p_df4x%P  : error alloc'
    allocate(p_df4x(lmg)%Pi(-1:nb_i_blocks-1,4),stat=ok); if (ok/=0) stop 'p_df4x%Pi : error alloc'
    allocate(p_df4x(lmg)%Q (0:lm+1,-1:3),stat=ok);        if (ok/=0) stop 'p_df4x%Q  : error alloc'
    allocate(p_df4x(lmg)%Qi(-1:nb_i_blocks-1,4),stat=ok); if (ok/=0) stop 'p_df4x%Qi : error alloc'

    call coef_p_df4_seconde (p_df4x(lmg)%P(:,:));          p_df4x(lmg)%P(:,:)=p_df4x(lmg)%P(:,:)*dx*dx
    call coef_p_df4_premiere(p_df4x(lmg)%Q(:,:));          p_df4x(lmg)%Q(:,:)=p_df4x(lmg)%Q(:,:)*dx
    
    if (.not.is_west) p_df4x(lmg)%P(0,:)=0._prec 
    if (.not.is_east) p_df4x(lmg)%P(lm+1,:)=0._prec 
    if (.not.is_west) p_df4x(lmg)%Q(0,:)=0._prec 
    if (.not.is_east) p_df4x(lmg)%Q(lm+1,:)=0._prec 
    
    call facto_dacx(p_df4x(lmg)%P(:,:),p_df4x(lmg)%Pi(:,:)); 
    call facto_dacx(p_df4x(lmg)%Q(:,:),p_df4x(lmg)%Qi(:,:))

    return
  end subroutine initialise_df4_part_coefx


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine release_df4_part_coefx(lmg)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc,     only : dx
    implicit none
    integer :: lm,lmg,ok

    deallocate(p_df4x(lmg)%P ,stat=ok); if (ok/=0) stop 'p_df4x%P  : error dealloc'
    deallocate(p_df4x(lmg)%Pi,stat=ok); if (ok/=0) stop 'p_df4x%Pi : error dealloc'
    deallocate(p_df4x(lmg)%Q ,stat=ok); if (ok/=0) stop 'p_df4x%Q  : error dealloc'
    deallocate(p_df4x(lmg)%Qi,stat=ok); if (ok/=0) stop 'p_df4x%Qi : error dealloc'
 
    return
  end subroutine release_df4_part_coefx


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine initialise_df4_part_coefy(nmg,nm)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc,     only : dy
    implicit none
    integer :: nmg,nm,ok

    allocate(p_df4y(nmg)%P (0:nm+1,-1:3)     ,stat=ok);   if (ok/=0) stop 'p_df4y%P  :error alloc'
    allocate(p_df4y(nmg)%Pi(-1:nb_k_blocks-1,4),stat=ok); if (ok/=0) stop 'p_df4y%Pi :error alloc'
    allocate(p_df4y(nmg)%Q (0:nm+1,-1:3),stat=ok);        if (ok/=0) stop 'p_df4y%Q  :error alloc'
    allocate(p_df4y(nmg)%Qi(-1:nb_k_blocks-1,4),stat=ok); if (ok/=0) stop 'p_df4y%Qi :error alloc'
    
    call coef_p_df4_seconde (p_df4y(nmg)%P(:,:));         p_df4y(nmg)%P(:,:)=p_df4y(nmg)%P(:,:)*dy*dy
    call coef_p_df4_premiere(p_df4y(nmg)%Q(:,:));         p_df4y(nmg)%Q(:,:)=p_df4y(nmg)%Q(:,:)*dy

    if (.not.is_south) p_df4y(nmg)%P(0,:)=0._prec 
    if (.not.is_north) p_df4y(nmg)%P(nm+1,:)=0._prec 
    if (.not.is_south) p_df4y(nmg)%Q(0,:)=0._prec 
    if (.not.is_north) p_df4y(nmg)%Q(nm+1,:)=0._prec 
    
    call facto_dacy(p_df4y(nmg)%P(:,:),p_df4y(nmg)%Pi(:,:)); 
    call facto_dacy(p_df4y(nmg)%Q(:,:),p_df4y(nmg)%Qi(:,:))

    return
  end subroutine initialise_df4_part_coefy

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine release_df4_part_coefy(lmg)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc,     only : dy
    implicit none
    integer :: lm,lmg,ok

    deallocate(p_df4y(lmg)%P ,stat=ok); if (ok/=0) stop 'p_df4y%P  : error dealloc'
    deallocate(p_df4y(lmg)%Pi,stat=ok); if (ok/=0) stop 'p_df4y%Pi : error dealloc'
    deallocate(p_df4y(lmg)%Q ,stat=ok); if (ok/=0) stop 'p_df4y%Q  : error dealloc'
    deallocate(p_df4y(lmg)%Qi,stat=ok); if (ok/=0) stop 'p_df4y%Qi : error dealloc'

    return
  end subroutine release_df4_part_coefy

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine initialise_vf4_part_coefx(lmg,lm)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc,     only : dx
    implicit none
    integer :: lmg,lm2,lm,ok

    lm2=lm+2; if (.not.is_east) lm2=lm+1
    allocate(p_vf4x(lmg)%P (0:lm2,-1:3)     ,stat=ok);    if (ok/=0) stop 'p_vf4x%P  : error alloc'
    allocate(p_vf4x(lmg)%Pi(-1:nb_i_blocks-1,4),stat=ok); if (ok/=0) stop 'p_vf4x%Pi : error alloc'
    allocate(p_vf4x(lmg)%Q (0:lm2,-1:3),stat=ok);         if (ok/=0) stop 'p_vf4x%Q  : error alloc'
    allocate(p_vf4x(lmg)%Qi(-1:nb_i_blocks-1,4),stat=ok); if (ok/=0) stop 'p_vf4x%Qi : error alloc'

    call coef_p_vf4_premiere (p_vf4x(lmg)%P(:,:),is_west,is_east); p_vf4x(lmg)%P(:,:)=p_vf4x(lmg)%P(:,:)*dx
    call coef_p_vf4_milieu   (p_vf4x(lmg)%Q(:,:),is_west,is_east); p_vf4x(lmg)%Q(:,:)=p_vf4x(lmg)%Q(:,:)

    if (.not.is_west) p_vf4x(lmg)%P(0,:)=0._prec 
    if (.not.is_east) p_vf4x(lmg)%P(lm2,:)=0._prec 
    if (.not.is_west) p_vf4x(lmg)%Q(0,:)=0._prec 
    if (.not.is_east) p_vf4x(lmg)%Q(lm2,:)=0._prec 

    call facto_dacx(p_vf4x(lmg)%P(:,:),p_vf4x(lmg)%Pi(:,:)); 
    call facto_dacx(p_vf4x(lmg)%Q(:,:),p_vf4x(lmg)%Qi(:,:))

    return
  end subroutine initialise_vf4_part_coefx

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine release_vf4_part_coefx(lmg)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc,     only : dx
    implicit none
    integer :: lmg,lm2,lm,ok

    lm2=lm+2; if (.not.is_east) lm2=lm+1
    deallocate(p_vf4x(lmg)%P ,stat=ok); if (ok/=0) stop 'p_vf4x%P  : error dealloc'
    deallocate(p_vf4x(lmg)%Pi,stat=ok); if (ok/=0) stop 'p_vf4x%Pi : error dealloc'
    deallocate(p_vf4x(lmg)%Q ,stat=ok); if (ok/=0) stop 'p_vf4x%Q  : error dealloc'
    deallocate(p_vf4x(lmg)%Qi,stat=ok); if (ok/=0) stop 'p_vf4x%Qi : error dealloc'

    return
  end subroutine release_vf4_part_coefx

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine initialise_vf4_part_coefy(nmg,nm)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc,     only : dy
    implicit none
    integer :: nmg,nm,nm2,ok

    nm2=nm+2; if (.not.is_north) nm2=nm+1
    allocate(p_vf4y(nmg)%P (0:nm2,-1:3)     ,stat=ok);    if (ok/=0) stop 'p_vf4y%P  : error alloc'
    allocate(p_vf4y(nmg)%Pi(-1:nb_k_blocks-1,4),stat=ok); if (ok/=0) stop 'p_vf4y%Pi : error alloc'
    allocate(p_vf4y(nmg)%Q (0:nm2,-1:3),stat=ok);         if (ok/=0) stop 'p_vf4y%Q  : error alloc'
    allocate(p_vf4y(nmg)%Qi(-1:nb_k_blocks-1,4),stat=ok); if (ok/=0) stop 'p_vf4y%Qi : error alloc'

    call coef_p_vf4_premiere(p_vf4y(nmg)%P(:,:),is_south,is_north); p_vf4y(nmg)%P(:,:)=p_vf4y(nmg)%P(:,:)*dy
    call coef_p_vf4_milieu  (p_vf4y(nmg)%Q(:,:),is_south,is_north); p_vf4y(nmg)%Q(:,:)=p_vf4y(nmg)%Q(:,:)

    if (.not.is_south) p_vf4y(nmg)%P(0,:)=0._prec 
    if (.not.is_north) p_vf4y(nmg)%P(nm2,:)=0._prec 
    if (.not.is_south) p_vf4y(nmg)%Q(0,:)=0._prec 
    if (.not.is_north) p_vf4y(nmg)%Q(nm2,:)=0._prec 
    
    call facto_dacy(p_vf4y(nmg)%P(:,:),p_vf4y(nmg)%Pi(:,:)); 
    call facto_dacy(p_vf4y(nmg)%Q(:,:),p_vf4y(nmg)%Qi(:,:))

    return
  end subroutine initialise_vf4_part_coefy


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine release_vf4_part_coefy(lmg)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc,     only : dy
    implicit none
    integer :: lmg,lm2,lm,ok

    lm2=lm+2; if (.not.is_east) lm2=lm+1
    deallocate(p_vf4y(lmg)%P ,stat=ok); if (ok/=0) stop 'p_vf4y%P  : error dealloc'
    deallocate(p_vf4y(lmg)%Pi,stat=ok); if (ok/=0) stop 'p_vf4y%Pi : error dealloc'
    deallocate(p_vf4y(lmg)%Q ,stat=ok); if (ok/=0) stop 'p_vf4y%Q  : error dealloc'
    deallocate(p_vf4y(lmg)%Qi,stat=ok); if (ok/=0) stop 'p_vf4y%Qi : error dealloc'

    return
  end subroutine release_vf4_part_coefy
  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coef_p_df4_seconde(P)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none 
    real(kind=prec), dimension(0:,-1:) :: P
    integer :: lm

    lm = size(P,1)-2

    P=0._prec
    P(0,0)       = 1._prec
    P(1:lm,-1:1) = spread( (/   1._prec,   10._prec, 1._prec /) / 12._prec, 1, lm )
    P(lm+1,0)    = P(0,0)

    return
  end subroutine coef_p_df4_seconde

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coef_p_df4_premiere(Q)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none 
    real(kind=prec), dimension(0:,-1:) :: Q
    integer :: lm

    lm = size(Q,1)-2

    Q=0._prec
    Q(0,0)       = 1._prec
    Q(1:lm,-1:1) = spread( (/  1._prec,  4._prec, 1._prec /) /  6._prec, 1, lm)
    Q(lm+1,0)    = Q(0,0)

    return
  end subroutine coef_p_df4_premiere

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coef_p_vf4_premiere(P,is_deb,is_fin)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none 
    real(kind=prec), dimension(0:,-1:) :: P
    logical :: is_deb,is_fin
    integer :: lm,lm0,lm1

    lm = size(P,1)-3

    P=0._prec
    if (is_deb) P(0,0)       = 1._prec
    if (is_deb) P(1,-1:1)    = (/ 9._prec, 54._prec, 42._prec /) /27._prec
    if (is_fin) P(lm+1,-1:1) = (/42._prec, 54._prec, 9._prec /) /27._prec
    if (is_fin) P(lm+2,0)    = 1._prec
    lm0=0; lm1=lm+2; 
    if (is_deb) lm0=2
    if (is_fin) lm1=lm
    P(lm0:lm1,-1:1) = spread( (/  27._prec,   186._prec,  27._prec /) / 186._prec, 1, lm1-lm0+1 )

    return
  end subroutine coef_p_vf4_premiere

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coef_p_vf4_milieu(Q,is_deb,is_fin)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

    implicit none 
    real(kind=prec), dimension(0:,-1:) :: Q
    logical :: is_deb,is_fin
    integer :: lm,lm0,lm1

    lm = size(Q,1)-3

    Q=0._prec
    if (is_deb) Q(0,0)       = 1._prec
    if (is_deb) Q(1,0:1)     = (/ 1._prec, 5._prec/3._prec /)
    if (is_fin) Q(lm+1,-1:0) = (/ 5._prec/3._prec, 1._prec /)
    if (is_fin) Q(lm+2,0)    = 1._prec
    lm0=0; lm1=lm+2; 
    if (is_deb) lm0=2
    if (is_fin) lm1=lm
    Q(lm0:lm1,-1:1) = spread( (/  0.3_prec,  1._prec, 0.3_prec /) , 1, lm1-lm0+1)

    return
  end subroutine coef_p_vf4_milieu

  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine tridiacx(Q,Qi,H,INPX,OUTX)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(-6:6,-2:2)                   :: H
    real(kind=prec), dimension(0:,-1:),target               :: Q
    real(kind=prec), dimension(-1:,:)  ,target    :: Qi
    real(kind=prec), dimension(0:,0:)                       :: INPX, OUTX
    real(kind=prec), dimension(-2:size(INPX,1)+1,0:size(INPX,2)-1)  :: INPX2
    real(kind=prec), dimension(0:size(OUTX,1)-1,0:size(OUTX,2)-1)   :: DLX

    real(kind=prec), dimension(:), pointer :: &
         alx,blx,clx,flx,glx,aix,bix,cix,fix

    real(kind=prec), dimension(0:1,0:size(OUTX,2)-1,0:nb_i_blocks-1) :: buffer
    real(kind=prec), dimension(0:1,0:size(OUTX,2)-1) :: buffer2d
    real(kind=prec), dimension(2*size(OUTX,2)) :: buffer1d
    real(kind=prec), dimension(-1:nb_i_blocks-1,0:size(OUTX,2)-1)    :: vec,dix

    integer :: lm,lmu,nm,i,j,k,i0,i1,k0,k1

    lm =size(INPX,1)-2
    nm =size(INPX,2)-2
    lmu=size(OUTX,1)-2

    INPX2=0._prec
    INPX2(0:lm+1,:)=INPX
    OUTX=0._prec


    alx => Q(:,-1);    aix => Qi(:,1)
    blx => Q(:,0) ;    bix => Qi(:,2)
    clx => Q(:,1) ;    cix => Qi(:,3)
    flx => Q(:,2) ;    fix => Qi(:,4)
    glx => Q(:,3)      ! on doit rajouter 1 a tous les indices 
                       ! car un pointeur commence forcement a 1!!!!!!!!!!

    call rfr_stn(INPX2,stn_ew,3,1)
    i0=1; i1=lm; k0=1; k1=nm; if (is_south) k0=0;  if (is_north) k1=nm+1

    if (is_west) then
       do k=k0,k1
          DLX(0,k)   =dot_product(H(0:6,-2),INPX2(0:6,k))
          DLX(1,k)   =dot_product(H(-1:6,-1),INPX2(0:7,k))
       enddo
       i0=2
    endif

    if (is_east) then
       do k=k0,k1
          DLX(lmu,k)  =dot_product(H(-6:1,1),INPX2(lmu-6:lmu+1,k))
          DLX(lmu+1,k)=dot_product(H(-6:0,2),INPX2(lmu-5:lmu+1,k))
       enddo
       i1=lmu-1
    endif

    do k=k0,k1
       do i=i0,i1
          DLX(i,k)=dot_product(H(-2:2,0),INPX2(i-2:i+2,k))
       enddo
    enddo

    i0=1; i1=lm; if (is_west) i0=0; if (is_east) i1=lmu+1

    if (ncheck==df4e.or.ncheck==vf4e) then
       OUTX(i0:i1,k0:k1)=DLX(i0:i1,k0:k1)
    else
       do k=k0,k1
          do i=i0+1,i1
             DLX(i,k)=DLX(i,k)-flx(1+i)*DLX(i-1,k)
          enddo
       enddo
       do k=k0,k1
          do i=i1-2,i0,-1
             DLX(i,k)=DLX(i,k)-glx(1+i)*DLX(i+1,k)
          enddo
       enddo

       !SK     
       !SK     collecte information systeme interface
       !SK     
       buffer(iu_d,:,icolumn)=DLX(i0,:)   ! iu_d=0
       buffer(id_d,:,icolumn)=DLX(i1,:)   ! id_d=1

       if (nb_i_blocks.ne.1) then
          j=iline*nb_i_blocks
          do i=j,j+nb_i_blocks-1
             if (i.ne.my_task) then
                buffer2d = buffer(:,:,icolumn)
                buffer1d = reshape(buffer(:,:,icolumn),(/size(buffer,1)*size(buffer,2)/))
                call snd_msg(i,tag_dacx+my_task,buffer1d)
             end if
          enddo
          do i=j,j+nb_i_blocks-1
             if (i.ne.my_task) then
                call rcv_msg(i,tag_dacx+i,buffer1d)
                !buffer(:,:,i-j) = buffer2d
                buffer(:,:,i-j) = reshape(buffer1d,(/size(buffer,1),size(buffer,2)/))
             end if
          enddo
       endif

       !SK     
       !SK     formation second membre systeme interface
       !SK     


       do k=k0,k1
          do i=0,nb_i_blocks-2
             dix(i,k)=buffer(id_d,k,i)-fix(2+i)*buffer(iu_d,k,i+1)
          enddo

          if (icolumn.ne.nb_i_blocks-1) then
             DLX(i1,k)=dix(icolumn,k)
          endif

          dix(-1,k)=buffer(iu_d,k,0)

          dix(nb_i_blocks-1,k)=buffer(id_d,k,nb_i_blocks-1)

          !SK     
          !SK     resolution systeme interface, facto, desc, montee
          !SK     

          vec(-1,k)=dix(-1,k)
          do i=0,nb_i_blocks-1
             vec(i,k)=dix(i,k)-aix(2+i)*vec(i-1,k)
          enddo

          vec(nb_i_blocks-1,k)=vec(nb_i_blocks-1,k)/bix(2+nb_i_blocks-1)
          do i=nb_i_blocks-2,-1,-1
             vec(i,k)=(vec(i,k)-cix(2+i)*vec(i+1,k))/bix(2+i)
          enddo

          !SK     
          !SK     resolution systeme local
          !SK     

          OUTX(i0,k)=vec(icolumn-1,k)
          OUTX(i1,k)=vec(icolumn,k)
          do i=i1-1,i0,-1
             OUTX(i,k)=(DLX(i,k)-alx(1+i)*OUTX(i0,k)-clx(1+i)*OUTX(i1,k))/blx(1+i)
          enddo
       enddo
    endif

    tag_dacx=tag_dacx+1

    return
  end subroutine tridiacx

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine tridiacy(Q,Qi,H,INPY,OUTY)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(-6:6,-2:2) :: H
    real(kind=prec), dimension(0:,-1:), target    :: Q
    real(kind=prec), dimension(-1:,:), target    :: Qi
    real(kind=prec), dimension(0:,0:)     :: INPY, OUTY
    real(kind=prec), dimension(0:size(INPY,1)-1,-2:size(INPY,2)+1)  :: INPY2
    real(kind=prec), dimension(0:size(OUTY,1)-1,0:size(OUTY,2)-1)   :: DLY

    real(kind=prec), dimension(:), pointer :: &
         aly,bly,cly,fly,gly,aiy,biy,ciy,fiy

    real(kind=prec), dimension(0:1,0:size(OUTY,1)-1,0:nb_k_blocks-1) :: buffer
    real(kind=prec), dimension(0:size(OUTY,1)-1,-1:nb_k_blocks-1)   :: vec,diy

    integer :: lm,nm,nmv,i,j,k,i0,i1,k0,k1

    lm =size(INPY,1)-2
    nm =size(INPY,2)-2
    nmv=size(OUTY,2)-2

    INPY2=0._prec
    INPY2(:,0:nm+1)=INPY
    OUTY=0._prec

    aly => Q(:,-1);    aiy => Qi(:,1)
    bly => Q(:,0) ;    biy => Qi(:,2)
    cly => Q(:,1) ;    ciy => Qi(:,3)
    fly => Q(:,2) ;    fiy => Qi(:,4)
    gly => Q(:,3)      ! on doit rajouter 1 a tous les indices 
                       ! car un pointeur commence forcement a 1!!!!!!!!!!

    call rfr_stn(INPY2,stn_ns,1,3)
    i0=1; i1=lm; k0=1; k1=nm; if (is_west) i0=0; if (is_east) i1=lm+1

    if (is_south) then
       do i=i0,i1
          DLY(i,0)   =dot_product(H(0:6,-2),INPY2(i, 0:+6))
          DLY(i,1)   =dot_product(H(-1:6,-1),INPY2(i,0:+7))
       enddo
       k0=2
    endif

    if (is_north) then
       do i=i0,i1
          DLY(i,nmv)  =dot_product(H(-6:1, 1),INPY2(i,nmv-6:nmv+1))
          DLY(i,nmv+1)=dot_product(H(-6:0, 2),INPY2(i,nmv-5:nmv+1))
       enddo
       k1=nmv-1
    endif

    do k=k0,k1
       do i=i0,i1
          DLY(i,k)=dot_product(H(-2:2,0),INPY2(i,k-2:k+2))
       enddo
    enddo

    k0=1; k1=nm; if (is_south) k0=0; if (is_north) k1=nmv+1

    if (ncheck==df4e.or.ncheck==vf4e) then
       OUTY(i0:i1,k0:k1)=DLY(i0:i1,k0:k1)
    else
       !SK    
       !SK    elimination diagonale inferieure
       !SK    
       do k=k0+1,k1
          do i=i0,i1
             DLY(i,k)=DLY(i,k)-fly(1+k)*DLY(i,k-1)
          enddo
       enddo
       !SK    
       !SK    elimination diagonale superieure
       !SK    
       do k=k1-2,k0,-1
          do i=i0,i1
             DLY(i,k)=DLY(i,k)-gly(1+k)*DLY(i,k+1)
          enddo
       enddo

       !SK    
       !SK    collecte information systeme interface
       !SK    
       buffer(iu_d,:,iline)=DLY(:,k0)  ! iu_d=0
       buffer(id_d,:,iline)=DLY(:,k1)  ! id_d=1

       if (nb_k_blocks.ne.1) then
          j=icolumn
          do i=j,nb_tasks-1,nb_i_blocks
             if (i.ne.my_task) call snd_msg(i,tag_dacy+my_task,buffer(:,:,iline))
          enddo
          k=0
          do i=j,nb_tasks-1,nb_i_blocks
             if (i.ne.my_task) call rcv_msg(i,tag_dacy+i,buffer(:,:,k))
             k=k+1
          enddo
       endif

       !SK    
       !SK    formation second membre systeme interface
       !SK    
       do i=i0,i1
          do k=0,nb_k_blocks-2
             diy(i,k)=buffer(id_d,i,k)-fiy(2+k)*buffer(iu_d,i,k+1)
          enddo

          if (iline.ne.nb_k_blocks-1) DLY(i,k1)=diy(i,iline)
          diy(i,-1)=buffer(iu_d,i,0)
          diy(i,nb_k_blocks-1)=buffer(id_d,i,nb_k_blocks-1)

          !SK    
          !SK    resolution a l'interface
          !SK    

          vec(i,-1)=diy(i,-1)
          do k=0,nb_k_blocks-1
             vec(i,k)=diy(i,k)-aiy(2+k)*vec(i,k-1)
          enddo

          vec(i,nb_k_blocks-1)=vec(i,nb_k_blocks-1)/biy(2+nb_k_blocks-1)
          do k=nb_k_blocks-2,-1,-1
             vec(i,k)=(vec(i,k)-ciy(2+k)*vec(i,k+1))/biy(2+k)
          enddo

          !SK    
          !SK    resolution systeme local
          !SK    

          OUTY(i,k0)=vec(i,iline-1)
          OUTY(i,k1)=vec(i,iline)
          do k=k1-1,k0,-1
             OUTY(i,k)=(DLY(i,k)-aly(1+k)*OUTY(i,k0)-cly(1+k)*OUTY(i,k1))/bly(1+k)
          enddo
       enddo
    endif

    tag_dacy=tag_dacy+1

    return
  end subroutine tridiacy

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine facto_dacx(Q,Qi)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,-1:), target :: Q
    real(kind=prec), dimension(:,:) , target   :: Qi

    real(kind=prec), dimension(:), pointer :: &
         alx,blx,clx,flx,glx,aix,bix,cix,fix
    real(kind=prec), dimension(0:5,0:nb_i_blocks-1) :: buffer
    integer          :: lm,i0,i1,i,j

    lm=size(Q,1)-2
    i0=1; i1=lm; if (is_west) i0=0;  if (is_east) i1=lm+1

    alx => Q(:,-1);    aix => Qi(:,1)
    blx => Q(:,0) ;    bix => Qi(:,2)
    clx => Q(:,1) ;    cix => Qi(:,3)
    flx => Q(:,2) ;    fix => Qi(:,4)
    glx => Q(:,3)

    Qi=0._prec;  flx=0._prec; glx=0._prec

    !SK     
    !SK     elimination diagonale inferieure
    !SK     
    do i=i0+1,i1
       flx(1+i)=alx(1+i)/blx(1+i-1)
       alx(1+i)=-flx(1+i)*alx(1+i-1)
       blx(1+i)=blx(1+i)-flx(1+i)*clx(1+i-1)
    enddo
    !SK     
    !SK     elimination diagonale superieure
    !SK     
    do i=i1-2,i0,-1
       glx(1+i)=clx(1+i)/blx(1+i+1)
       alx(1+i)=alx(1+i)-glx(1+i)*alx(1+i+1)
       clx(1+i)=-glx(1+i)*clx(1+i+1)
    enddo

    !SK     
    !SK     collecte information systeme interface
    !SK     
    buffer(iu_a,icolumn)=alx(1+i0)  ! iu_a = 0
    buffer(iu_b,icolumn)=blx(1+i0)  ! iu_b = 1
    buffer(iu_c,icolumn)=clx(1+i0)  ! iu_c = 2
    buffer(id_a,icolumn)=alx(1+i1)  ! id_a = 3
    buffer(id_b,icolumn)=blx(1+i1)  ! id_b = 4
    buffer(id_c,icolumn)=clx(1+i1)  ! id_c = 5

    if (nb_i_blocks.ne.1) then
       j=iline*nb_i_blocks
       do i=j,j+nb_i_blocks-1
          if (i/=my_task) call snd_msg(i,tag_dacx+my_task,buffer(:,icolumn))
       enddo
       do i=j,j+nb_i_blocks-1
          if (i/=my_task) call rcv_msg(i,tag_dacx+i,buffer(:,i-j))
       enddo
    endif

    !SK     
    !SK     formation systeme interface
    !SK     
    do i=0,nb_i_blocks-2
       aix(2+i)=buffer(id_a,i)
       fix(2+i)=buffer(id_c,i)/buffer(iu_b,i+1)
       bix(2+i)=buffer(id_b,i)-fix(2+i)*buffer(iu_a,i+1)
       cix(2+i)=-fix(2+i)*buffer(iu_c,i+1)
    enddo

    if (icolumn.ne.nb_i_blocks-1) then
       blx(1+i1)=bix(2+icolumn)
       clx(1+i1)=cix(2+icolumn)
    endif


    aix(2-1)=buffer(iu_a,0)
    bix(2-1)=buffer(iu_b,0)
    cix(2-1)=buffer(iu_c,0)

    aix(2+nb_i_blocks-1)=buffer(id_a,nb_i_blocks-1)
    bix(2+nb_i_blocks-1)=buffer(id_b,nb_i_blocks-1)
    cix(2+nb_i_blocks-1)=buffer(id_c,nb_i_blocks-1)

    do i=0,nb_i_blocks-1
       aix(2+i)=aix(2+i)/bix(2+i-1)
       bix(2+i)=bix(2+i)-cix(2+i-1)*aix(2+i)
    enddo

    tag_dacx=tag_dacx+1

    return 
  end subroutine facto_dacx

  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine facto_dacy(Q,Qi)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,-1:), target   :: Q
    real(kind=prec), dimension(:,:),    target   :: Qi

    real(kind=prec), dimension(:), pointer :: &
         aly,bly,cly,fly,gly,aiy,biy,ciy,fiy
    real(kind=prec), dimension(0:5,0:nb_k_blocks-1) :: buffer

    integer :: nm, i, j, k, k0, k1

    nm=size(Q,1)-2
    k0=1; k1=nm; if (is_south) k0=0;  if (is_north) k1=nm+1

    aly => Q(:,-1);    aiy => Qi(:,1)
    bly => Q(:,0) ;    biy => Qi(:,2)
    cly => Q(:,1) ;    ciy => Qi(:,3)
    fly => Q(:,2) ;    fiy => Qi(:,4)
    gly => Q(:,3)

    !SK    
    !SK    elimination diagonale inferieure
    !SK    
    do k=k0+1,k1
       fly(1+k)= aly(1+k)/bly(1+k-1)
       aly(1+k)=-fly(1+k)*aly(1+k-1)
       bly(1+k)=bly(1+k)-fly(1+k)*cly(1+k-1)
    enddo
    !SK    
    !SK    elimination diagonale superieure
    !SK    
    do k=k1-2,k0,-1
       gly(1+k)=cly(1+k)/bly(1+k+1)
       aly(1+k)=aly(1+k)-gly(1+k)*aly(1+k+1)
       cly(1+k)=-gly(1+k)*cly(1+k+1)
    enddo

    !SK    
    !SK    collecte information systeme interface
    !SK    
    buffer(iu_a,iline)=aly(1+k0)  ! iu_a = 0
    buffer(iu_b,iline)=bly(1+k0)  ! iu_b = 1
    buffer(iu_c,iline)=cly(1+k0)  ! iu_c = 2
    buffer(id_a,iline)=aly(1+k1)  ! id_a = 3
    buffer(id_b,iline)=bly(1+k1)  ! id_b = 4
    buffer(id_c,iline)=cly(1+k1)  ! id_c = 5

    if (nb_k_blocks/=1) then
       j=icolumn
       do i=j,nb_tasks-1,nb_i_blocks
          if (i/=my_task) call snd_msg(i,tag_dacy+my_task,buffer(:,iline))
       enddo
       k=0
       do i=j,nb_tasks-1,nb_i_blocks
          if (i/=my_task) call rcv_msg(i,tag_dacy+i,buffer(:,k))
          k=k+1
       enddo
    endif

    !SK    
    !SK    formation systeme interface
    !SK    

    do k=0,nb_k_blocks-2
       aiy(2+k)=buffer(id_a,k)
       fiy(2+k)=buffer(id_c,k)/buffer(iu_b,k+1)
       biy(2+k)=buffer(id_b,k)-fiy(2+k)*buffer(iu_a,k+1)
       ciy(2+k)=-fiy(2+k)*buffer(iu_c,k+1)
    enddo

    if (iline.ne.nb_k_blocks-1) then
       bly(1+k1)=biy(2+iline)
       cly(1+k1)=ciy(2+iline)
    endif

    aiy(2-1)=buffer(iu_a,0)
    biy(2-1)=buffer(iu_b,0)
    ciy(2-1)=buffer(iu_c,0)

    aiy(2+nb_k_blocks-1)=buffer(id_a,nb_k_blocks-1)
    biy(2+nb_k_blocks-1)=buffer(id_b,nb_k_blocks-1)
    ciy(2+nb_k_blocks-1)=buffer(id_c,nb_k_blocks-1)



    !SK    
    !SK    resolution systeme interface, facto, desc, montee
    !SK    

    do k=0,nb_k_blocks-1
       aiy(2+k)=aiy(2+k)/biy(2+k-1)
       biy(2+k)=biy(2+k)-ciy(2+k-1)*aiy(2+k)
    enddo

    tag_dacy=tag_dacy+1

    return
  end subroutine facto_dacy

  !************************************************************************

end module dac







