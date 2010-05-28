module test

  use blas
  use compact

contains

  !************************************************************************

 subroutine test_df4
   use disc
   use para
   use drapeaux
   use data, only : ncheck
   implicit none
   real(kind=prec), dimension(size(xu,1),size(xu,2)) :: inpx
   real(kind=prec), dimension(size(yu,1),size(yu,2)) :: inpy
   real(kind=prec), dimension(size(xu,1),size(xu,2)) :: outx
   real(kind=prec), dimension(size(yu,1),size(yu,2)) :: outy


   real(kind=prec), dimension(size(xu,1),size(xu,2)) :: inpu
   real(kind=prec), dimension(size(xv,1),size(xv,2)) :: inpv
   real(kind=prec), dimension(size(xp,1),size(xp,2)) :: inpp
   real(kind=prec), dimension(size(xu,1),size(xu,2)) :: outu
   real(kind=prec), dimension(size(xv,1),size(xv,2)) :: outv
   real(kind=prec), dimension(size(xp,1),size(xp,2)) :: outp


   if (ncheck/=df4e.and.ncheck/=df4 ) stop 'Verification  DF4 necessite ncheck=df4,df4e'

    print *,'verification derivee premiere DF4'
    inpu=xu*xu/2._prec
    outu=der1x(inpu,lmu_global)-xu
    call outdon(outu,'der 1 en u',1,1)
    call outdon(outu,'der 1 en u')
    inpu=yu*yu/2._prec
    outu=der1y(inpu,nmu_global)-yu
    call outdon(outu,'der 1 en u',1,1)
    call outdon(outu,'der 1 en u')

    print *,'verification dérivée seconde DF4'
    inpu=xu*xu/2._prec
    outu=der2x(inpu,lmu_global)-1.
    call outdon(outu,'der 2 en u',1,1)
    call outdon(outu,'der 2 en u')
    inpu=yu*yu/2._prec
    outu=der2y(inpu,nmu_global)-1.
    call outdon(outu,'der 2 en u',1,1)
    call outdon(outu,'der 2 en u')

    print *,'verification derivee premiere DF4'
    inpv=xv*xv/2._prec
    outv=der1x(inpv,lmv_global)-xv
    call outdon(outv,'der 1 en v',1,1)
    call outdon(outv,'der 1 en v')
    inpv=yv*yv/2._prec
    outv=der1y(inpv,nmv_global)-yv
    call outdon(outv,'der 1 en v',1,1)
    call outdon(outv,'der 1 en v')

    print *,'verification dérivée seconde DF4'
    inpv=xv*xv/2._prec
    outv=der2x(inpv,lmv_global)-1.
    call outdon(outv,'der 2 en v',1,1)
    call outdon(outv,'der 2 en v')
    inpv=yv*yv/2._prec
    outv=der2y(inpv,nmv_global)-1.
    call outdon(outv,'der 2 en v',1,1)
    call outdon(outv,'der 2 en v')



!!$
!!$    print *,'verification Couplage pression/vitesse pour DF4'
!!$    inpp=(xp-.5*dx)**2/2._prec
!!$    inpu=0._prec; inpv=0._prec
!!$    call set_switch(u_code*v_code*div_code*grdp_code)
!!$    call mavecxy(inpu,inpv,inpp,outu,outv,outp)
!!$    call outdon(xp-.5*dx,'inpp')
!!$    call outdon(outv,'outv')
!!$    call outdon(outu,'outu')
!!$
!!$    inpp=(yp-.5*dy)**2/2._prec
!!$    inpu=0._prec; inpv=0._prec
!!$    call set_switch(u_code*v_code*div_code*grdp_code)
!!$    call mavecxy(inpu,inpv,inpp,outu,outv,outp)
!!$    call outdon(yp-.5*dy,'inpp')
!!$    call outdon(outv,'outv')
!!$    call outdon(outu,'outu')


!!$    inpu=(xu)**2/2._prec
!!$    inpp=0._prec; inpv=0._prec
!!$    call set_switch(u_code*v_code*div_code*grdp_code)
!!$    call mavecxy(inpu,inpv,inpp,outu,outv,outp)
!!$    call outdon(inpu,'inpu')
!!$    call outdon(outp,'outp')
!!$
!!$
!!$    inpv=(yv)**2/2._prec
!!$    inpp=0._prec; inpu=0._prec
!!$    call set_switch(u_code*v_code*div_code*grdp_code)
!!$    call mavecxy(inpu,inpv,inpp,outu,outv,outp)
!!$    call outdon(inpv,'inpv')
!!$    call outdon(outp,'outp')
!!$
!!$



    call parallel_stop

    return
  end subroutine test_df4


  !************************************************************************

!!$  subroutine test_vf4
!!$    use disc
!!$    use para
!!$    use data, only : ncheck
!!$    implicit none
!!$    real(kind=prec), dimension(size(xu,1),size(xu,2))   :: inpx,outix
!!$    real(kind=prec), dimension(size(yu,1),size(yu,2))   :: inpy
!!$    real(kind=prec), dimension(size(xu,1)+1,size(xu,2)) :: outx
!!$    real(kind=prec), dimension(size(yu,1),size(yu,2)+1) :: outy
!!$
!!$    if (ncheck/=vf4.and.ncheck/=vf4e) stop 'Verification  VF4 necessite ncheck=vf4'
!!$
!!$    print *,'verification milieu VF4'
!!$    inpx=xu
!!$    outx=derm0x(inpx,lmu_global)/dx
!!$    call outdon(outx,'outx')
!!$    inpy=yu
!!$    outy=derm0y(inpy,nmu_global)/dy
!!$    call outdon(outy,'outy')
!!$
!!$    print *,'verification dérivée première VF4'
!!$    inpx=xu
!!$    outx=derm1x(inpx,lmu_global)
!!$    call outdon(outx,'outx')
!!$    inpy=yu
!!$    outy=derm1y(inpy,nmu_global)
!!$    call outdon(outy,'outy')
!!$
!!$!LIN
!!$!LIN    print *,'verification intégration VF4'
!!$!LIN    inpx=XU/dx+YU/dy
!!$!LIN    outix=integ_vf(inpx)
!!$!LIN    call outdon(inpx,'inpx',1,1)
!!$!LIN    call outdon(outix,'outx',1,1)
!!$!LIN    call outdon(outix,'outx')
!!$!LIN
!!$    call parallel_stop
!!$
!!$    return
!!$  end subroutine test_vf4

  !************************************************************************

end module test
