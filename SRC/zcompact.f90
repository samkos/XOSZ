module compact

  use constante
  use para
  use dac

  type compact_matrix
     real(kind=prec), dimension (:,:), pointer :: P,Q
  end type compact_matrix

  type (compact_matrix), dimension(:), allocatable :: &
       df4x,df4y,vf4x,vf4y             ! coefficients monodomaines

  real(kind=prec), dimension(-6:6,-2:2) :: Gdf4, Hdf4, Gvf4,Hvf4  &
                                          ,Gcdf4x,Gcdf4y,Hcdf4x,Hcdf4y &
                                          ,Gcvf4x,Gcvf4y,Hcvf4x,Hcvf4y 


  real(kind=prec), dimension(-3:3,-2:2) :: Mvf4
  real(kind=prec), dimension(-5:5,-2:2) :: Mavf4
  real(kind=prec), dimension(-5:5,-2:2) :: Mdvf4

  integer :: tag_derm0x=0, tag_derm1x=500, tag_derm0y=0, tag_derm1y=500

contains


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine initialise_general_coef
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use champs
    use disc,     only : lmu_global,nmv_global,dx,dy
    use data,     only : ncheck_precond
    use drapeaux, only : is_decal, is_ns,is_precond
    implicit none
    integer :: lm0,lm1,nm0,nm1,ok

    Gdf4=0._prec
    Gdf4(0:5,-2)  = (/ 225._prec, -770._prec, 1070._prec, -780._prec, 305._prec, -50._prec /) / 60._prec
    Gdf4(-1:1,-1) = (/   1._prec,   -2._prec, 1._prec /)
    Gdf4(-1:1,0)  = (/   1._prec,   -2._prec, 1._prec /)
    Gdf4(-1:1,1)  = (/   1._prec,   -2._prec, 1._prec /)
    Gdf4(-5:0,2)  = Gdf4(5:0:-1,-2)

    Gcdf4x=0._prec                            ! coefficients explicites --> P=Id
    Gcdf4x(0:5,-2)  = Gdf4(0:5,-2)  
    Gcdf4x(-1:3,-1) = (/ 11._prec,-20._prec,   6._prec,  4._prec, -1._prec /) / 12._prec
    Gcdf4x(-2:2,0)  = (/ -1._prec, 16._prec, -30._prec, 16._prec, -1._prec /) / 12._prec
    Gcdf4x(-3:1,1)  = Gcdf4x(3:-1:-1,-1)
    Gcdf4x(-5:0,2)  = Gcdf4x(5:0:-1,-2)
    Gcdf4y=Gcdf4x/dy/dy
    Gcdf4x=Gcdf4x/dx/dx

    Hdf4=0._prec
    Hdf4(0:3,-2)  = (/ -11._prec, 18._prec, -9._prec, 2._prec /) / 6._prec
    Hdf4(-1:1,-1) = (/ -1._prec,   0._prec,  1._prec /) / 2._prec
    Hdf4(-1:1,0)  = (/ -1._prec,   0._prec,  1._prec /) / 2._prec
    Hdf4(-1:1,1)  = (/ -1._prec,   0._prec,  1._prec /) / 2._prec
    Hdf4(-3:0,2)  = - Hdf4(3:0:-1,-2)

    Hcdf4x=0._prec                            ! coefficients explicites --> P=Id
    Hcdf4x(0:3,-2)  = Hdf4(0:3,-2)  
    Hcdf4x(-1:3,-1) = (/ -0.75_prec, -2.5_prec,  4.5_prec,  -1.5_prec, 0.25_prec /) / 3._prec
    Hcdf4x(-2:2,0)  = (/ 1._prec, -8._prec, 0._prec, 8._prec, -1._prec  /) / 12._prec
    Hcdf4x(-3:1,1)  = -Hcdf4x(3:-1:-1,-1)
    Hcdf4x(-3:0,2)  = - Hcdf4x(3:0:-1,-2)
    Hcdf4y=Hcdf4x/dy
    Hcdf4x=Hcdf4x/dx

    Gvf4=0._prec
    Gvf4(0:5,-2)  = (/ -6.85_prec/3._prec, 5._prec, -5._prec, 10._prec/3._prec, -1.25_prec, 0.20_prec /)
    Gvf4(-1:2,-1) = (/ -66.5_prec, 27._prec, 40.5_prec, -1._prec /) /27._prec
    Gvf4(-2:1,0)  =  (/ -17._prec, -189._prec, 189._prec, 17._prec /)/186._prec
    Gvf4(-3:0,1)  = -Gvf4(2:-1:-1,-1)
    Gvf4(-6:-1,2) = -Gvf4(5:0:-1,-2)

    Gcvf4x=0._prec
    Gcvf4x(0:5,-2)  = Gvf4(0:5,-2)
    Gcvf4x(-1:3,-1) = (/ -22._prec, 17._prec, 9._prec, -5._prec, 1._prec /) / 24._prec
    Gcvf4x(-2:1,0)  =  (/ 1._prec, -27._prec, 27._prec, -1._prec /)/24._prec
    Gcvf4x(-4:0,1)  = -Gcvf4x(3:-1:-1,-1)
    Gcvf4x(-6:-1,2) = -Gcvf4x(5:0:-1,-2)
    Gcvf4y  = Gcvf4x/dy
    Gcvf4x  = Gcvf4x/dx
    

    Hvf4=0._prec
    Hvf4(0,-2)    = 1._prec
    Hvf4(-1:2,-1) = (/ 0.625_prec, 5.625_prec, 1.875_prec, -0.125_prec /) / 3._prec
    Hvf4(-2:1,0)  = (/ 0.05_prec,   0.75_prec,  0.75_prec,   0.05_prec /)
    Hvf4(-3:0,1)  = Hvf4(2:-1:-1,-1)
    Hvf4(-1,2)    = Hvf4(0,-2)


    Hcvf4x=0._prec
    Hcvf4x(0,-2)    = 1._prec
    Hcvf4x(-1:3,-1) = (/ 0.273437499995_prec, 1.09375_prec, -0.546874999993_prec,&
         & 0.218749999997_prec, -0.0390624999994_prec /)
    Hcvf4x(-2:1,0)  = (/ -.25_prec, 1.25_prec, 1.25_prec, -.25_prec /)/2._prec
    Hcvf4x(-4:0,1)  = Hcvf4x(3:-1:-1,-1)
    Hcvf4x(-1,2)    = Hcvf4x(0,-2)
    Hcvf4y          = Hcvf4x

    Mvf4=0._prec
    Mvf4(0:3,-2)  = (/ 29.75_prec, 26.75_prec, -10.75_prec, 2.25_prec /) /48._prec
    Mvf4(-1:1,-1) = (/ 1._prec/24._prec, 11._prec/12._prec, 1._prec/24._prec /)
    Mvf4(-2:2,0)  = (/ 1._prec, -1._prec, 72._prec, -1._prec, 1._prec /) /72._prec
    Mvf4(-1:1,1)  = Mvf4(1:-1:-1,-1)
    Mvf4(-3:0,2)  = Mvf4(3:0:-1,-2)

    Mavf4=0._prec
    Mavf4(0:5,-2)  = &
(/ 166.725_prec, -78.05_prec, 108.95_prec, -79.80_prec, 31.325_prec, -5.15_prec /) / 144._prec
!!$(/ 1.0081747298_prec,   -0.0106752494115_prec, -7.86109888581_prec*0.001_prec, &
!!$                        0.0163041617756_prec, -7.1976782564*0.001_prec, 1.25513497638_prec*0.001_prec/)
!!$(/ 29.75_prec, 26.75_prec, -10.75_prec, 2.25_prec, 0._prec, 0._prec /) /24._prec
    Mavf4(-1:1,-1) = (/ 1._prec/24._prec, 11._prec/12._prec, 1._prec/24._prec /)
    Mavf4(-2:2,0)  = (/ 1._prec, -1._prec, 72._prec, -1._prec, 1._prec /) /72._prec
    Mavf4(-1:1,1)  = Mavf4(1:-1:-1,-1)
    Mavf4(-5:0,2)  = Mavf4(5:0:-1,-2)

    Mdvf4=0._prec
    Mdvf4(0:5,-2)  = (/ 3800._prec,   11416._prec,  -6384._prec, 3856._prec, -1384._prec, 216._prec /) /256._prec/45._prec
    Mdvf4(-1:4,-1) = (/ -432._prec, 10192._prec, 16352._prec, -4128._prec, 1232._prec, -176._prec /)/256._prec/90._prec
    Mdvf4(-2:3,0)  = (/ 22._prec, -186._prec, 1604._prec, 1604._prec, -186._prec, 22._prec /) / 2880._prec
    Mdvf4(-3:2,1)  = Mdvf4(4:-1:-1,-1)
    Mdvf4(-4:1,2)  = Mdvf4(5:0:-1,-2)


    lm1=size(VTUS,1)-2
    nm1=size(VTVS,2)-2
    lm0=lm1; if (is_decal.or.is_ns) lm0=lm1-1
    nm0=nm1; if (is_decal.or.is_ns) nm0=nm1-1
    
    if (is_precond.and.(ncheck_precond==df4.or.ncheck_precond==df4e.or.ncheck_precond==vf4)) then
       lm1=size(VTUp,1)-2
       nm1=size(VTVp,2)-2
    endif

    allocate(df4x(lm0:lm1),stat=ok); if (ok/=0) stop 'df4x : error alloc'
    allocate(df4y(nm0:nm1),stat=ok); if (ok/=0) stop 'df4y : error alloc'
    allocate(vf4x(lm0:lm1),stat=ok); if (ok/=0) stop 'vf4x : error alloc'
    allocate(vf4y(nm0:nm1),stat=ok); if (ok/=0) stop 'vf4y : error alloc'

    if (is_mpp) then
       lm1=lmu_global
       nm1=nmv_global
       lm0=lm1; if (is_decal.or.is_ns) lm0=lm1-1
       nm0=nm1; if (is_decal.or.is_ns) nm0=nm1-1

       allocate(p_df4x(lm0:lm1),stat=ok); if (ok/=0) stop 'p_df4x : error alloc'
       allocate(p_df4y(nm0:nm1),stat=ok); if (ok/=0) stop 'p_df4y : error alloc'
       allocate(p_vf4x(lm0:lm1),stat=ok); if (ok/=0) stop 'p_vf4x : error alloc'
       allocate(p_vf4y(nm0:nm1),stat=ok); if (ok/=0) stop 'p_vf4y : error alloc'
    endif
    
    return
  end subroutine initialise_general_coef

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine release_general_coef
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use champs
    use disc,     only : lmu_global,nmv_global,dx,dy
    use data,     only : ncheck_precond
    use drapeaux, only : is_decal, is_ns,is_precond
    implicit none
    integer ok

    deallocate(df4x,stat=ok); if (ok/=0) stop 'df4x : error dealloc'
    deallocate(df4y,stat=ok); if (ok/=0) stop 'df4y : error dealloc'
    deallocate(vf4x,stat=ok); if (ok/=0) stop 'vf4x : error dealloc'
    deallocate(vf4y,stat=ok); if (ok/=0) stop 'vf4y : error dealloc'

    if (is_mpp) then
       deallocate(p_df4x,stat=ok); if (ok/=0) stop 'p_df4x : error dealloc'
       deallocate(p_df4y,stat=ok); if (ok/=0) stop 'p_df4y : error dealloc'
       deallocate(p_vf4x,stat=ok); if (ok/=0) stop 'p_vf4x : error dealloc'
       deallocate(p_vf4y,stat=ok); if (ok/=0) stop 'p_vf4y : error dealloc'
    endif

    return
  end subroutine release_general_coef

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine initialise_df4_coef(VTUS,VTVS)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use drapeaux, only : is_decal, is_ns
    use disc,     only : dx,dy
    implicit none
    real(kind=prec), dimension(:,:) :: VTUS,VTVS
    integer :: lm0,nm0,lm1,nm1,lm,nm,ok

    lm1=size(VTUS,1)-2
    nm1=size(VTVS,2)-2
    lm0=lm1; if (is_decal.or.is_ns) lm0=lm1-1
    nm0=nm1; if (is_decal.or.is_ns) nm0=nm1-1

    do lm=lm0,lm1
!       if (.not.allocated(df4x(lm)%P)) then
          allocate(df4x(lm)%P(-1:1,0:lm+1),stat=ok); if (ok/=0) stop 'df4x%P : error alloc'
          allocate(df4x(lm)%Q(-1:1,0:lm+1),stat=ok); if (ok/=0) stop 'df4x%Q : error alloc'
          call coef_df4_seconde (df4x(lm)%P);    df4x(lm)%P=df4x(lm)%P*dx*dx
          call coef_df4_premiere(df4x(lm)%Q);    df4x(lm)%Q=df4x(lm)%Q*dx
          call factorise(df4x(lm)%P); call factorise(df4x(lm)%Q)
!       endif
    enddo

    do nm=nm0,nm1
!       if (.not.allocated(df4y(lm)%P(:,:))) then
          allocate(df4y(nm)%P(-1:1,0:nm+1),stat=ok); if (ok/=0) stop 'df4y%P : error alloc'
          allocate(df4y(nm)%Q(-1:1,0:nm+1),stat=ok); if (ok/=0) stop 'df4y%Q : error alloc'
          call coef_df4_seconde (df4y(nm)%P);    df4y(nm)%P=df4y(nm)%P*dy*dy
          call coef_df4_premiere(df4y(nm)%Q);    df4y(nm)%Q=df4y(nm)%Q*dy
          call factorise(df4y(nm)%P); call factorise(df4y(nm)%Q)
!       endif
    enddo

    return
  end subroutine initialise_df4_coef

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine release_df4_coef(VTUS,VTVS)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use drapeaux, only : is_decal, is_ns
    use disc,     only : dx,dy
    implicit none
    real(kind=prec), dimension(:,:) :: VTUS,VTVS
    integer :: lm0,nm0,lm1,nm1,lm,nm,ok

    lm1=size(VTUS,1)-2
    nm1=size(VTVS,2)-2
    lm0=lm1; if (is_decal.or.is_ns) lm0=lm1-1
    nm0=nm1; if (is_decal.or.is_ns) nm0=nm1-1

    do lm=lm0,lm1
          deallocate(df4x(lm)%P,stat=ok); if (ok/=0) stop 'df4x%P : error dealloc'
          deallocate(df4x(lm)%Q,stat=ok); if (ok/=0) stop 'df4x%Q : error dealloc'
    enddo

    do nm=nm0,nm1
          deallocate(df4y(nm)%P,stat=ok); if (ok/=0) stop 'df4y%P : error dealloc'
          deallocate(df4y(nm)%Q,stat=ok); if (ok/=0) stop 'df4y%Q : error dealloc'
    enddo

    return
  end subroutine release_df4_coef

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine initialise_vf4_coef(VTUS,VTVS)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use drapeaux, only : is_decal, is_ns
    use disc,    only : dx,dy
    implicit none
    real(kind=prec), dimension(:,:) :: VTUS,VTVS
    integer :: lm,nm,lm0,nm0,lm1,nm1,ok

    lm1=size(VTUS,1)-2
    nm1=size(VTVS,2)-2
    lm0=lm1; if (is_decal.or.is_ns) lm0=lm1-1
    nm0=nm1; if (is_decal.or.is_ns) nm0=nm1-1

    do lm=lm0,lm1
!       if (.not.allocated(vf4x(lm)%P(:,:))) then
          allocate(vf4x(lm)%P(-1:1,0:lm+2),stat=ok); if (ok/=0) stop 'vf4x%P : error alloc'
          allocate(vf4x(lm)%Q(-1:1,0:lm+2),stat=ok); if (ok/=0) stop 'vf4x%Q : error alloc'
          call coef_vf4_premiere(vf4x(lm)%P);    vf4x(lm)%P=vf4x(lm)%P*dx
          call coef_vf4_milieu  (vf4x(lm)%Q);   
          ! xsamxsamx
          ! if (.not.is_west) vf4x(lm)%P(:,1)    = vf4x(lm)%P(:,2)
          ! if (.not.is_east) vf4x(lm)%P(:,lm+1) = vf4x(lm)%P(:,2)
          ! if (.not.is_west) vf4x(lm)%Q(:,1)    = vf4x(lm)%Q(:,2)
          ! if (.not.is_east) vf4x(lm)%Q(:,lm+1) = vf4x(lm)%Q(:,2)
          ! xsamxsamx
          call factorise(vf4x(lm)%P); call factorise(vf4x(lm)%Q)
!       endif
    enddo

    do nm=nm0,nm1
!       if (.not.allocated(vf4y(lm)%P(:,:))) then
          allocate(vf4y(nm)%P(-1:1,0:nm+2),stat=ok); if (ok/=0) stop 'vf4y%P : error alloc'
          allocate(vf4y(nm)%Q(-1:1,0:nm+2),stat=ok); if (ok/=0) stop 'vf4y%Q : error alloc'
          call coef_vf4_premiere(vf4y(nm)%P);    vf4y(nm)%P=vf4y(nm)%P*dy
          call coef_vf4_milieu  (vf4y(nm)%Q);  
          ! xsamxsamx
          ! if (.not.is_south) vf4y(nm)%P(:,1)    = vf4y(nm)%P(:,2)
          ! if (.not.is_north) vf4y(nm)%P(:,nm+1) = vf4y(nm)%P(:,2)
          ! if (.not.is_south) vf4y(nm)%Q(:,1)    = vf4y(nm)%Q(:,2)
          ! if (.not.is_north) vf4y(nm)%Q(:,nm+1) = vf4y(nm)%Q(:,2)
          ! xsamxsamx
          call factorise(vf4y(nm)%P); call factorise(vf4y(nm)%Q)
!       endif
    enddo

    return
  end subroutine initialise_vf4_coef




  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine release_vf4_coef(VTUS,VTVS)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use drapeaux, only : is_decal, is_ns
    use disc,    only : dx,dy
    implicit none
    real(kind=prec), dimension(:,:) :: VTUS,VTVS
    integer :: lm,nm,lm0,nm0,lm1,nm1,ok

    lm1=size(VTUS,1)-2
    nm1=size(VTVS,2)-2
    lm0=lm1; if (is_decal.or.is_ns) lm0=lm1-1
    nm0=nm1; if (is_decal.or.is_ns) nm0=nm1-1

    do lm=lm0,lm1
       deallocate(vf4x(lm)%P,stat=ok); if (ok/=0) stop 'vf4x%P : error dealloc'
       deallocate(vf4x(lm)%Q,stat=ok); if (ok/=0) stop 'vf4x%Q : error dealloc'
    enddo

    do nm=nm0,nm1
          deallocate(vf4y(nm)%P,stat=ok); if (ok/=0) stop 'vf4y%P : error dealloc'
          deallocate(vf4y(nm)%Q,stat=ok); if (ok/=0) stop 'vf4y%Q : error dealloc'
    enddo

    return
  end subroutine release_vf4_coef




  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coef_df4_seconde(P)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none 
    real(kind=prec), dimension(-1:,0:) :: P
    integer :: lm

    lm = size(P,2)-2

    P=0._prec
    P(0,0)       = 1._prec
    P(-1:1,1:lm) = spread( (/   1._prec,   10._prec, 1._prec /) / 12._prec, 2, lm )
    P(0,lm+1)    = P(0,0)

    return
  end subroutine coef_df4_seconde

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coef_df4_premiere(Q)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none 
    real(kind=prec), dimension(-1:,0:) :: Q
    integer :: lm

    lm = size(Q,2)-2

    Q=0._prec
    Q(0,0)       = 1._prec
    Q(-1:1,1:lm) = spread( (/  1._prec,  4._prec, 1._prec /) /  6._prec, 2, lm)
    Q(0,lm+1)    = Q(0,0)

    return
  end subroutine coef_df4_premiere

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coef_vf4_premiere(P)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none 
    real(kind=prec), dimension(-1:,0:) :: P
    integer :: lm

    lm = size(P,2)-3

    P=0._prec
    P(0,0)       = 1._prec
    P(-1:1,1)    = (/ 9._prec, 54._prec, 42._prec /) /27._prec
    P(-1:1,2:lm) = spread( (/  27._prec,   186._prec,  27._prec /) / 186._prec, 2, lm-1 )
    P(-1:1,lm+1) = P(1:-1:-1,1)
    P(0,lm+2)    = P(0,0)

    return
  end subroutine coef_vf4_premiere

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coef_vf4_milieu(Q)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

    implicit none 
    real(kind=prec), dimension(-1:,0:) :: Q
    integer :: lm

    lm = size(Q,2)-3

    Q=0._prec
    Q(0,0)       = 1._prec
    Q(0:1,1)     = (/ 1._prec, 5._prec/3._prec /)
    Q(-1:1,2:lm) = spread( (/  0.3_prec,  1._prec, 0.3_prec /) , 2, lm-1)
    Q(-1:0,lm+1) = Q(1:0:-1,1)
    Q(0,lm+2)    = Q(0,0)

    return
  end subroutine coef_vf4_milieu

  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine factorise(P)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

    implicit none
    real(kind=prec), dimension(-1:,0:) :: P
    integer :: lm,i

    lm=size(P,2)-2

    do i=1,lm+1
       P(-1,i)=P(-1,i)/P(0,i-1)
       P(0,i)=P(0,i)-P(1,i-1)*P(-1,i)
    enddo

    do i=lm,0,-1
       P(1,i)=P(1,i)/P(0,i+1)
    enddo

    return
  end subroutine factorise

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine tridiagx(Q,H,INPX,OUTX)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only : ncheck
    implicit none
    real(kind=prec), dimension(-6:6,-2:2) :: H
    real(kind=prec), dimension(-1:,0:)    :: Q
    real(kind=prec), dimension(0:,0:)     :: INPX, OUTX
    real(kind=prec), dimension(-2:size(INPX,1)+1,0:size(INPX,2)-1)  :: INPX2
    real(kind=prec), dimension(0:size(OUTX,1)-1,0:size(OUTX,2)-1)   :: DLX
    integer :: lm,lmu,nm,i,k,i0,i1
    real(kind=prec) :: a

    lm =size(OUTX,1)-2
    nm =size(OUTX,2)-2
    lmu=size(INPX,1)-2

    INPX2=0._prec
    INPX2(0:lmu+1,:)=INPX
    OUTX=0._prec

    do k=0,nm+1
       DLX(0,k)   =dot_product(H(0:6,-2),INPX2(0:6,k))
       DLX(1,k)   =dot_product(H(-1:6,-1),INPX2(0:7,k))
       DLX(lm,k)  =dot_product(H(-6:1,1),INPX2(lm-6:lm+1,k))
       DLX(lm+1,k)=dot_product(H(-6:0,2),INPX2(lm-5:lm+1,k))
    enddo
    i0=2; i1=lm-1

    do k=0,nm+1
       do i=i0,i1
          DLX(i,k)=dot_product(H(-2:2,0),INPX2(i-2:i+2,k))
       enddo
    enddo

    if (ncheck==df4e.or.ncheck==vf4e) then
       OUTX=DLX
    else
       !
       !     elimination diagonale inferieure
       !     
       do k=0,nm+1
          do i=1,lm+1
             DLX(i,k)=DLX(i,k)-Q(-1,i)*DLX(i-1,k)
          enddo
       enddo
       !     
       !     elimination diagonale superieure
       !     
       do k=0,nm+1
          OUTX(lm+1,k)=DLX(lm+1,k)/Q(0,lm+1)
       enddo
       do i=lm,0,-1
          a=Q(1,i)
          do k=0,nm+1
             DLX(i,k)=DLX(i,k)-a*DLX(i+1,k)
          enddo
       enddo
       do i=lm,0,-1
          a=1._prec/Q(0,i)
          do k=0,nm+1
             OUTX(i,k)=a*DLX(i,k)
          enddo
       enddo
    endif

    return
  end subroutine tridiagx

  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine tridiagy(Q,H,INPY,OUTY)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only : ncheck
    implicit none
    real(kind=prec), dimension(-6:6,-2:2) :: H
    real(kind=prec), dimension(-1:,0:)    :: Q
    real(kind=prec), dimension(0:,0:)     :: INPY, OUTY
    real(kind=prec), dimension(0:size(INPY,1)-1,-2:size(INPY,2)+1)  :: INPY2
    real(kind=prec), dimension(0:size(OUTY,1)-1,0:size(OUTY,2)-1)   :: DLY
    integer :: lm,nm,nmv,i,k,k0,k1
    real(kind=prec) :: a

    lm =size(OUTY,1)-2
    nm =size(OUTY,2)-2
    nmv=size(INPY,2)-2

    INPY2=0._prec
    INPY2(:,0:nmv+1)=INPY
    OUTY=0._prec

    do i=0,lm+1
       DLY(i,0)   =dot_product(H(0:6,-2),INPY2(i, 0:+6))
       DLY(i,1)   =dot_product(H(-1:6,-1),INPY2(i,0:+7))
       DLY(i,nm)  =dot_product(H(-6:1, 1),INPY2(i,nm-6:nm+1))
       DLY(i,nm+1)=dot_product(H(-6:0, 2),INPY2(i,nm-5:nm+1))
    enddo
    k0=2; k1=nm-1

    do k=k0,k1
       do i=0,lm+1
          DLY(i,k)=dot_product(H(-2:2,0),INPY2(i,k-2:k+2))
       enddo
    enddo

    if (ncheck==df4e.or.ncheck==vf4e) then
       OUTY=DLY
    else
       !     
       !     elimination diagonale inferieure
       !     
       do k=1,nm+1
          do i=0,lm+1
             DLY(i,k)=DLY(i,k)-Q(-1,k)*DLY(i,k-1)
          enddo
       enddo
       !     
       !     elimination diagonale superieure
       !     
       do i=0,lm+1
          OUTY(i,nm+1)=DLY(i,nm+1)/Q(0,nm+1)
       enddo
       do k=nm,0,-1
          a=Q(1,k)
          do i=0,lm+1
             DLY(i,k)=DLY(i,k)-a*DLY(i,k+1)
          enddo
       enddo
       do k=nm,0,-1
          a=1._prec/Q(0,k)
          do i=0,lm+1
             OUTY(i,k)=a*DLY(i,k)
          enddo
       enddo
    endif

    return
  end subroutine tridiagy

  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function integ_vf(INP) result(OUT)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(:,:) :: INP
    real(kind=prec), dimension(size(INP,1),size(INP,2)) :: OUT,TMP         

    TMP=integ_x(INP)
    OUT=integ_y(TMP)

    return
  end function integ_vf

  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function integ_x(INP) result(OUT)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,0:) :: INP
    real(kind=prec), dimension(-2:size(INP,1)+1,0:size(INP,2)-1) :: INP2
    real(kind=prec), dimension(0:size(INP,1)-1 ,0:size(INP,2)-1) :: OUT
    integer :: lm,nm,i,k,i0,i1,k0,k1

    lm=size(INP,1)-2
    nm=size(OUT,2)-2

    INP2=0._prec
    INP2(0:lm+1,:)=INP
    OUT=0._prec

    if (is_mpp) call rfr_stn(INP2,stn_ew,3,0)
    i0=1; i1=lm; k0=1; k1=nm
    if (is_south.or..not.is_mpp) k0=0;   
    if (is_north.or..not.is_mpp) k1=nm+1

    if (is_west.or..not.is_mpp) then
       do k=k0,k1
          OUT(0,k)   = dot_product(Mvf4(0:3,-2) ,INP2(0:3,k))
          OUT(1,k)   = dot_product(Mvf4(-1:3,-1) ,INP2(0:4,k))
       enddo
       i0=2
    endif

    if (.not.is_mpp.or.is_east) then
       do k=k0,k1
          OUT(lm,k)  = dot_product(Mvf4(-3:1,1) ,INP2(lm-3:lm+1,k))
          OUT(lm+1,k)= dot_product(Mvf4(-3:0,2) ,INP2(lm-2:lm+1,k))
       enddo
       i1=lm-1
    endif

    do k=k0,k1
       do i=i0,i1
          OUT(i,k)=dot_product(Mvf4(-2:2,0),INP2(i-2:i+2,k))
       enddo
    enddo

    return
  end function integ_x

  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function integ_y(INP) result(OUT)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,0:) :: INP
    real(kind=prec), dimension(0:size(INP,1)-1,-2:size(INP,2)+1) :: INP2
    real(kind=prec), dimension(0:size(INP,1)-1,0:size(INP,2)-1 ) :: OUT
    integer :: lm,nm,i,k,i0,i1,k0,k1

    lm=size(INP,1)-2
    nm=size(OUT,2)-2

    INP2=0._prec
    INP2(:,0:nm+1)=INP
    OUT=0._prec

    if (is_mpp) call rfr_stn(INP2,stn_ns,0,3)
    i0=1; i1=lm; k0=1; k1=nm
    if (is_west.or..not.is_mpp) i0=0; 
    if (is_east.or..not.is_mpp) i1=lm+1


    if (.not.is_mpp.or.is_south) then
       do i=i0,i1
          OUT(i,0)   =dot_product(Mvf4(0:3,-2),INP2(i,0:+3))
          OUT(i,1)   =dot_product(Mvf4(-1:3,-1),INP2(i,0:4))
       enddo
       k0=2
    endif

    if (.not.is_mpp.or.is_north) then
       do i=i0,i1
          OUT(i,nm)  =dot_product(Mvf4(-3:1, 1),INP2(i,nm-3:nm+1))
          OUT(i,nm+1)=dot_product(Mvf4(-3:0, 2),INP2(i,nm-2:nm+1))
       enddo
       k1=nm-1
    endif

    do k=k0,k1
       do i=i0,i1
          OUT(i,k)=dot_product(Mvf4(-2:2,0),INP2(i,k-2:k+2))
       enddo
    enddo

    return
  end function integ_y

  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function intega_x(INP) result(OUT)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,0:) :: INP
    real(kind=prec), dimension(-2:size(INP,1)+1,0:size(INP,2)-1) :: INP2
    real(kind=prec), dimension(0:size(INP,1)-1 ,0:size(INP,2)-1) :: OUT
    integer :: lm,nm,i,k,i0,i1,k0,k1

    lm=size(INP,1)-2
    nm=size(OUT,2)-2

    INP2=0._prec
    INP2(0:lm+1,:)=INP
    OUT=0._prec

    if (is_mpp) call rfr_stn(INP2,stn_ew,3,0)
    i0=0; i1=lm+1; k0=1; k1=nm
    if (is_south.or..not.is_mpp) k0=0;   
    if (is_north.or..not.is_mpp) k1=nm+1
    k0=0; k1=nm+1

    if (is_west.or..not.is_mpp) then
       do k=k0,k1
          OUT(0,k)   = dot_product(Mavf4(0:3,-2) ,INP2(0:3,k))
          OUT(1,k)   = dot_product(Mavf4(-1:3,-1) ,INP2(0:4,k))
       enddo
       i0=2
    endif

    if (.not.is_mpp.or.is_east) then
       do k=k0,k1
          OUT(lm,k)  = dot_product(Mavf4(-3:1,1) ,INP2(lm-3:lm+1,k))
          OUT(lm+1,k)= dot_product(Mavf4(-3:0,2) ,INP2(lm-2:lm+1,k))
       enddo
       i1=lm-1
    endif

    do k=k0,k1
       do i=i0,i1
          OUT(i,k)=dot_product(Mavf4(-2:2,0),INP2(i-2:i+2,k))
       enddo
    enddo

    return
  end function intega_x

  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function intega_y(INP) result(OUT)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,0:) :: INP
    real(kind=prec), dimension(0:size(INP,1)-1,-2:size(INP,2)+1) :: INP2
    real(kind=prec), dimension(0:size(INP,1)-1,0:size(INP,2)-1 ) :: OUT
    integer :: lm,nm,i,k,i0,i1,k0,k1

    lm=size(INP,1)-2
    nm=size(OUT,2)-2

    INP2=0._prec
    INP2(:,0:nm+1)=INP
    OUT=0._prec

    if (is_mpp) call rfr_stn(INP2,stn_ns,0,3)
    i0=1; i1=lm; k0=1; k1=nm
    if (is_west.or..not.is_mpp) i0=0; 
    if (is_east.or..not.is_mpp) i1=lm+1
    i0=0; i1=lm+1

    if (.not.is_mpp.or.is_south) then
       do i=i0,i1
          OUT(i,0)   =dot_product(Mavf4(0:3,-2),INP2(i,0:+3))
          OUT(i,1)   =dot_product(Mavf4(-1:3,-1),INP2(i,0:4))
       enddo
       k0=2
    endif

    if (.not.is_mpp.or.is_north) then
       do i=i0,i1
          OUT(i,nm)  =dot_product(Mavf4(-3:1, 1),INP2(i,nm-3:nm+1))
          OUT(i,nm+1)=dot_product(Mavf4(-3:0, 2),INP2(i,nm-2:nm+1))
       enddo
       k1=nm-1
    endif

    do k=k0,k1
       do i=i0,i1
          OUT(i,k)=dot_product(Mavf4(-2:2,0),INP2(i,k-2:k+2))
       enddo
    enddo

    return
  end function intega_y

  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function integd_x(INP) result(OUT)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,0:) :: INP
    real(kind=prec), dimension(-3:size(INP,1)+2,0:size(INP,2)-1) :: INP2
    real(kind=prec), dimension(0:size(INP,1)-2 ,0:size(INP,2)-1) :: OUT
    integer :: lm,nm,i,k,i0,i1,k0,k1

    lm=size(INP,1)-2
    nm=size(INP,2)-2

    INP2=0._prec
    INP2(0:lm+1,:)=INP
    OUT=0._prec

    if (is_mpp) call rfr_stn(INP2,stn_ew,4,0)
    i0=1; i1=lm; k0=1; k1=nm
    if (is_south.or..not.is_mpp) k0=0;   
    if (is_north.or..not.is_mpp) k1=nm+1

    if (is_west.or..not.is_mpp) then
       do k=k0,k1
          OUT(0,k)   = dot_product(Mdvf4(0:5,-2)  ,INP2(0:5,k))
          OUT(1,k)   = dot_product(Mdvf4(-1:4,-1) ,INP2(0:5,k))
       enddo
       i0=2
    endif

    if (.not.is_mpp.or.is_east) then
       do k=k0,k1
          OUT(lm-1,k)  = dot_product(Mdvf4(-3:2,1) ,INP2(lm-4:lm+1,k))
          OUT(lm,k)= dot_product(Mdvf4(-4:1,2) ,INP2(lm-4:lm+1,k))
       enddo
       i1=lm-2
    endif

    do k=k0,k1
       do i=i0,i1
          OUT(i,k)=dot_product(Mdvf4(-2:3,0),INP2(i-2:i+3,k))
       enddo
    enddo

    return
  end function integd_x

  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function integd_y(INP) result(OUT)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(0:,0:) :: INP
    real(kind=prec), dimension(0:size(INP,1)-1,-3:size(INP,2)+2) :: INP2
    real(kind=prec), dimension(0:size(INP,1)-1,0:size(INP,2)-2 ) :: OUT
    integer :: lm,nm,i,k,i0,i1,k0,k1

    lm=size(INP,1)-2
    nm=size(INP,2)-2

    INP2=0._prec
    INP2(:,0:nm+1)=INP
    OUT=0._prec

    if (is_mpp) call rfr_stn(INP2,stn_ns,0,4)
    i0=1; i1=lm; k0=1; k1=nm
    if (is_west.or..not.is_mpp) i0=0; 
    if (is_east.or..not.is_mpp) i1=lm+1


    if (.not.is_mpp.or.is_south) then
       do i=i0,i1
          OUT(i,0)   =dot_product(Mdvf4( 0:5,-2),INP2(i,0:5))
          OUT(i,1)   =dot_product(Mdvf4(-1:4,-1),INP2(i,0:5))
       enddo
       k0=2
    endif

    if (.not.is_mpp.or.is_north) then
       do i=i0,i1
          OUT(i,nm-1)=dot_product(Mdvf4(-3:2, 1),INP2(i,nm-4:nm+1))
          OUT(i,nm)  =dot_product(Mdvf4(-4:1, 2),INP2(i,nm-4:nm+1))
       enddo
       k1=nm-2
    endif

    do k=k0,k1
       do i=i0,i1
          OUT(i,k)=dot_product(Mdvf4(-2:3,0),INP2(i,k-2:k+3))
       enddo
    enddo

    return
  end function integd_y

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function der1x(INPX,lmu) result(OUTX)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(:,:) :: INPX
    real(kind=prec), dimension(size(INPX,1),size(INPX,2)) :: OUTX
    integer :: lm,lmu


    if (is_mpp) then
       if (ncheck==df4e) then
          call tridiacx(p_df4x(lmu)%Q(:,:),p_df4x(lmu)%Qi(:,:),Hcdf4x,INPX,OUTX)
       else
          call tridiacx(p_df4x(lmu)%Q(:,:),p_df4x(lmu)%Qi(:,:),Hdf4,INPX,OUTX)
       endif
    else
       lm=size(INPX,1)-2
       if (ncheck==df4e) then
          call tridiagx(df4x(lm)%Q(:,:),Hcdf4x,INPX,OUTX)
       else
          call tridiagx(df4x(lm)%Q(:,:),Hdf4,INPX,OUTX)
       endif
    endif
    
    return
  end function der1x

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function der2x(INPX,lmu) result(OUTX)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(:,:) :: INPX
    real(kind=prec), dimension(size(INPX,1),size(INPX,2)) :: OUTX
    integer :: lm,lmu


    if (is_mpp) then
       if (ncheck==df4e) then
          call tridiacx(p_df4x(lmu)%P(:,:),p_df4x(lmu)%Pi(:,:),Gcdf4x,INPX,OUTX)
       else
          call tridiacx(p_df4x(lmu)%P(:,:),p_df4x(lmu)%Pi(:,:),Gdf4,INPX,OUTX)
       endif
    else
       lm=size(INPX,1)-2
       if (ncheck==df4e) then
          call tridiagx(df4x(lm)%P(:,:),Gcdf4x,INPX,OUTX)
       else
          call tridiagx(df4x(lm)%P(:,:),Gdf4,INPX,OUTX)
       endif
    endif

    return
  end function der2x

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function der1y(INPY,nmu) result(OUTY)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(:,:) :: INPY
    real(kind=prec), dimension(size(INPY,1),size(INPY,2)) :: OUTY
    integer :: nm,nmu

    if (is_mpp) then
       if (ncheck==df4e) then
          call tridiacy(p_df4y(nmu)%Q(:,:),p_df4y(nmu)%Qi(:,:),Hcdf4y,INPY,OUTY)
       else
          call tridiacy(p_df4y(nmu)%Q(:,:),p_df4y(nmu)%Qi(:,:),Hdf4,INPY,OUTY)
       endif
    else
       nm=size(INPY,2)-2
       if (ncheck==df4e) then
          call tridiagy(df4y(nm)%Q(:,:),Hcdf4y,INPY,OUTY)
       else
          call tridiagy(df4y(nm)%Q(:,:),Hdf4,INPY,OUTY)
       endif
    endif

    return
  end function der1y

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function der2y(INPY,nmu) result(OUTY)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(:,:) :: INPY
    real(kind=prec), dimension(size(INPY,1),size(INPY,2)) :: OUTY
    integer :: nm,nmu

    if (is_mpp) then
       if (ncheck==df4e) then
          call tridiacy(p_df4y(nmu)%P(:,:),p_df4y(nmu)%Pi(:,:),Gcdf4y,INPY,OUTY)
       else
          call tridiacy(p_df4y(nmu)%P(:,:),p_df4y(nmu)%Pi(:,:),Gdf4,INPY,OUTY)
       endif
    else
       nm=size(INPY,2)-2
       if (ncheck==df4e) then
          call tridiagy(df4y(nm)%P(:,:),Gcdf4y,INPY,OUTY)
       else
          call tridiagy(df4y(nm)%P(:,:),Gdf4,INPY,OUTY)
       endif
    endif


    return
  end function der2y

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine der2(INP,lmu,nmu,OUTX,OUTY)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !     OUTX = der2x(INP),   OUTY = der2y(INP)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(:,:) :: INP,OUTX,OUTY
    integer :: lm,lmu,nm,nmu


    if (is_mpp) then
       if (ncheck==df4e) then
          call tridiacx(p_df4x(lmu)%P(:,:),p_df4x(lmu)%Pi(:,:),Gcdf4x,INP,OUTX)
       else
          call tridiacx(p_df4x(lmu)%P(:,:),p_df4x(lmu)%Pi(:,:),Gdf4,INP,OUTX)
       endif
    else
       lm=size(INP,1)-2
       if (ncheck==df4e) then
          call tridiagx(df4x(lm)%P(:,:),Gcdf4x,INP,OUTX)
       else
          call tridiagx(df4x(lm)%P(:,:),Gdf4,INP,OUTX)
       endif
    endif


    if (is_mpp) then
       if (ncheck==df4e) then
          call tridiacy(p_df4y(nmu)%P(:,:),p_df4y(nmu)%Pi(:,:),Gcdf4y,INP,OUTY)
       else
          call tridiacy(p_df4y(nmu)%P(:,:),p_df4y(nmu)%Pi(:,:),Gdf4,INP,OUTY)
       endif
    else
       nm=size(INP,2)-2
       if (ncheck==df4e) then
          call tridiagy(df4y(nm)%P(:,:),Gcdf4y,INP,OUTY)
       else
          call tridiagy(df4y(nm)%P(:,:),Gdf4,INP,OUTY)
       endif
    endif



    return
  end subroutine der2



  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function derm0x(INPX,lmu) result(OUTX)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(:,:) :: INPX
    real(kind=prec), dimension(0:size(INPX,1),size(INPX,2)) :: OUTX
    integer :: lm,lmu

    lm=size(INPX,1)-2
    if (is_mpp) then
       if (ncheck==vf4e) then
          call tridiacx(p_vf4x(lmu)%Q(:,:),p_vf4x(lmu)%Qi(:,:),Hcvf4x,INPX,OUTX)
       else
          call tridiacx(p_vf4x(lmu)%Q(:,:),p_vf4x(lmu)%Qi(:,:),Hvf4,INPX,OUTX)
       endif
       if (.not.is_west) call snd_msg(p_west,tag_derm0x,OUTX(1,:))
       if (.not.is_east) call rcv_msg(p_east,tag_derm0x,OUTX(lm+1,:))
       tag_derm0x=tag_derm0x+1
    else
       if (ncheck==vf4e) then
          call tridiagx(vf4x(lm)%Q(:,:),Hcvf4x,INPX,OUTX)
       else
          call tridiagx(vf4x(lm)%Q(:,:),Hvf4,INPX,OUTX)
       endif
    endif
    
    return
  end function derm0x

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function derm1x(INPX,lmu) result(OUTX)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(:,:) :: INPX
    real(kind=prec), dimension(0:size(INPX,1),size(INPX,2)) :: OUTX
    integer :: lm,lmu

    lm=size(INPX,1)-2
    if (is_mpp) then
       if (ncheck==vf4e) then
          call tridiacx(p_vf4x(lmu)%P(:,:),p_vf4x(lmu)%Pi(:,:),Gcvf4x,INPX,OUTX)
       else
          call tridiacx(p_vf4x(lmu)%P(:,:),p_vf4x(lmu)%Pi(:,:),Gvf4,INPX,OUTX)
       endif
       if (.not.is_west) call snd_msg(p_west,tag_derm1x,OUTX(1,:))
       if (.not.is_east) call rcv_msg(p_east,tag_derm1x,OUTX(lm+1,:))
       tag_derm1x=tag_derm1x+1
    else
       if (ncheck==vf4e) then
          call tridiagx(vf4x(lm)%P(:,:),Gcvf4x,INPX,OUTX)
       else
          call tridiagx(vf4x(lm)%P(:,:),Gvf4,INPX,OUTX)
       endif
    endif

    return
  end function derm1x

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function derm0y(INPY,nmu) result(OUTY)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(:,:) :: INPY
    real(kind=prec), dimension(size(INPY,1),0:size(INPY,2)) :: OUTY
    integer :: nm,nmu

    nm=size(INPY,2)-2
    if (is_mpp) then
       if (ncheck==vf4e) then
          call tridiacy(p_vf4y(nmu)%Q(:,:),p_vf4y(nmu)%Qi(:,:),Hcvf4y,INPY,OUTY)
       else
          call tridiacy(p_vf4y(nmu)%Q(:,:),p_vf4y(nmu)%Qi(:,:),Hvf4,INPY,OUTY)
       endif
       if (.not.is_south) call snd_msg(p_south,tag_derm0y,OUTY(:,1))
       if (.not.is_north) call rcv_msg(p_north,tag_derm0y,OUTY(:,nm+1))
       tag_derm0y=tag_derm0y+1
    else
       if (ncheck==vf4e) then
          call tridiagy(vf4y(nm)%Q(:,:),Hcvf4y,INPY,OUTY)
       else
          call tridiagy(vf4y(nm)%Q(:,:),Hvf4,INPY,OUTY)
       endif
    endif

    return
  end function derm0y

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function derm1y(INPY,nmu) result(OUTY)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only  : ncheck
    implicit none
    real(kind=prec), dimension(:,:) :: INPY
    real(kind=prec), dimension(size(INPY,1),0:size(INPY,2)) :: OUTY
    integer :: nm,nmu

    nm=size(INPY,2)-2
    if (is_mpp) then
       if (ncheck==vf4e) then
          call tridiacy(p_vf4y(nmu)%P(:,:),p_vf4y(nmu)%Pi(:,:),Gcvf4y,INPY,OUTY)
       else
          call tridiacy(p_vf4y(nmu)%P(:,:),p_vf4y(nmu)%Pi(:,:),Gvf4,INPY,OUTY)
       endif
       if (.not.is_south) call snd_msg(p_south,tag_derm1y,OUTY(:,1))
       if (.not.is_north) call rcv_msg(p_north,tag_derm1y,OUTY(:,nm+1))
       tag_derm1y=tag_derm1y+1
    else
       if (ncheck==vf4e) then
          call tridiagy(vf4y(nm)%P(:,:),Gcvf4y,INPY,OUTY)
       else
          call tridiagy(vf4y(nm)%P(:,:),Gvf4,INPY,OUTY)
       endif
    endif

    return
  end function derm1y

  !************************************************************************


end module compact

