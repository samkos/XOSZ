module navier

  use selectprec
  use para

contains

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine ajout_compress(PRE,VTU,VTV,rho_ca,maxdiv)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use blas, only : mavecxy
    use drapeaux
    implicit none
    real(kind=prec), dimension(:,:)                     :: PRE,VTU,VTV
    real(kind=prec), dimension(size(PRE,1),size(PRE,2)) :: DIV
    real(kind=prec), dimension(size(VTU,1),size(VTU,2)) :: TMPX
    real(kind=prec), dimension(size(VTV,1),size(VTV,2)) :: TMPY
    real(kind=prec)                                     :: rho_ca,maxdiv
    integer :: old_switch

    old_switch=current_switch
    call set_switch(div_code)
    call mavecxy(VTU,VTV,PRE,TMPX,TMPY,DIV)
    maxdiv=maxval(abs(div))
    PRE = PRE - rho_ca*div
    call set_switch(old_switch)

    return
  end subroutine ajout_compress

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine ajout_grad(PRE,SMU,SMV,SMU0,SMV0)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use blas, only : mavecxy
    use drapeaux
    implicit none
    real(kind=prec), dimension(:,:) :: PRE,SMU,SMV
    real(kind=prec), dimension(size(SMU,1),size(SMU,2)) :: SMU0
    real(kind=prec), dimension(size(SMV,1),size(SMV,2)) :: SMV0
    real(kind=prec), dimension(size(PRE,1),size(PRE,2)) :: TMPP
    integer :: old_switch

    old_switch=current_switch
    call set_switch(grdp_code*u_code*v_code)
    call mavecxy(SMU,SMV,PRE,SMU0,SMV0,TMPP)
    SMU0=SMU-SMU0
    SMV0=SMV-SMV0
    call set_switch(old_switch)

    return
  end subroutine ajout_grad

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine sortie_fichier(VTU,VTV,PRE)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use drapeaux
    use plot_flags
    use io,    only : sortie_champ,read_plot_flags
    use para,  only : is_mpp
    use data,  only : nom_fic_output
    use uncol, only : timer_stop,timer_start
    implicit none

    real(kind=prec), dimension(:,:)              :: VTU,VTV,PRE
    real(kind=prec), dimension(:,:), allocatable :: PSI
    integer :: ok

    if (.not.is_ns) return

    call timer_stop(0)


    call read_plot_flags

    if (is_ns) call sortie_champ(PRE,'p')

    if (is_plot_psi.or.is_plot_w) then
       allocate(PSI(size(VTV,1),size(VTU,2)),stat=ok)
       if (ok/=0) stop 'PSI : allocation error!' 

       call calcul_psi(VTU,VTV,PSI)
       call sortie_champ(PSI,'v')

        deallocate(PSI,stat=ok); if (ok/=0) stop 'PSI : error desalloc'
    endif


    call timer_start(0)

    return
  end subroutine sortie_fichier


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine calcul_psi(VTU,VTV,PSI)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use plot_flags
    use calcul
    use data, only : is_print,nu,nom_fic_output
    use drapeaux
    use disc
    use io, only : sortie_champ
    implicit none
      
    integer lm,nm,i,k,niv_solve_old
    
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PSI
    real(kind=prec), dimension(0:size(PSI,1)-1,0:size(PSI,2)-1) :: SMPSI,FAKE
    real(kind=prec) :: nu_old,cx,cy
    integer :: old_switch
    logical :: old_is_precond

    lm=size(PSI,1)-2
    nm=size(PSI,2)-2

    if (is_decal) then
       do k=1,nm
          do i=1,lm
             SMPSI(i,k)=0.25_prec*(VTV(i+1,k)+VTV(i+1,k+1)-VTV(i-1,k)-VTV(i-1,k+1))/dx &
                       -0.25_prec*(VTU(i,k+1)+VTU(i+1,k+1)-VTU(i,k-1)-VTU(i+1,k-1))/dy
          enddo
       enddo
    else
       do k=1,nm
          do i=1,lm
             SMPSI(i,k)=(VTV(i+1,k)-VTV(i-1,k))/2._prec/dx-(VTU(i,k+1)-VTU(i,k-1))/2._prec/dy
          enddo
       enddo
     endif


     if (is_plot_w) call sortie_champ(SMPSI(1:lm,1:nm),'w')

     if (.not.is_plot_psi) return

      PSI=0._prec

     if (is_decal) then
        cx=0.25_prec*dx; cy=0.25_prec*dy; 
        do i=1,lm+1; PSI(i,0)   = PSI(i-1,0)-cx*(VTV(i,0)+VTV(i,1)+VTV(i-1,0)+VTV(i-1,1));            enddo
        do k=1,nm+1; PSI(0,k)   = PSI(0,k-1)+cy*(VTU(0,k)+VTU(0,k-1)+VTU(1,k)+VTU(1,k-1));            enddo;
        do i=1,lm+1; PSI(i,nm+1)= PSI(i-1,nm+1)-cx*(VTV(i,nm+1)+VTV(i-1,nm+1)+VTV(i,nm+2)+VTV(i-1,nm+2)); enddo
        do k=1,nm;   PSI(lm+1,k)= PSI(lm+1,k-1)+cy*(VTU(lm+1,k)+VTU(lm+1,k-1)+VTU(lm+2,k)+VTU(lm+2,k-1)); enddo
     else
        cx=0.5_prec*dx; cy=0.5_prec*dy; 
        do i=1,lm+1; PSI(i,0)   = PSI(i-1,0)-cx*(VTV(i,0)+VTV(i-1,0));          enddo
        do k=1,nm+1; PSI(0,k)   = PSI(0,k-1)+cy*(VTU(0,k)+VTU(0,k-1));          enddo;
        do i=1,lm+1; PSI(i,nm+1)= PSI(i-1,nm+1)-cx*(VTV(i,nm+1)+VTV(i-1,nm+1)); enddo
        do k=1,nm;   PSI(lm+1,k)= PSI(lm+1,k-1)+cy*(VTU(lm+1,k)+VTU(lm+1,k-1)); enddo
     endif

     if (is_print>1) &
        print *,'difference d''integration des cond. aux limites' &
          ,PSI(lm+1,nm)-dy*VTU(lm+1,nm)-PSI(lm+1,nm+1)

      SMPSI(0,0:nm+1)=PSI(0,0:nm+1)
      SMPSI(lm+1,0:nm+1)=PSI(lm+1,0:nm+1)
      SMPSI(1:lm,0)=PSI(1:lm,0)
      SMPSI(1:lm,nm+1)=PSI(1:lm,nm+1)

!!$     call outdon(VTU,'VTU')
!!$     call outdon(VTV,'VTV')
!!$     call outdon(SMPSI,'SMPSI')
!!$
!SK
!SK   conditions de  dirichlet
!SK
      nu_old=nu
      nu=1._prec
      old_switch=current_switch
      call set_switch(u_code*diff_code)
      old_is_precond=is_precond
      is_precond=.false.
      niv_solve_old=niv_solve
      niv_solve=10
      call solve_bicgstab(SMPSI,FAKE,FAKE,PSI,FAKE,FAKE,200)
      niv_solve=niv_solve_old
      nu=nu_old
      is_precond=old_is_precond
      call set_switch(old_switch)

      return
      end subroutine calcul_psi

end module navier
