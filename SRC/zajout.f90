module ajout

  use para, only : outdon, output, parallel_stop

  use calcul
  
  use data,     only : is_print,epsv,epsva, ninterne,ndirection
  use drapeaux, only : is_u, is_v, is_div, is_precond
  use blas,     only : mavecxy, global_ddot, fixe_cl
  use disc,     only : it

contains


  !SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !SK  Routines calculs
  !SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !************************************************************************

!!$  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!!$  recursive &
!!$       subroutine solve_qmr(SMU,SMV,SMP,VTU,VTV,PRE,ncgmax)
!!$    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!!$    !SK
!!$    !SK   Méthode QMR
!!$    !SK
!!$    !SK================================================================
!!$    !SK
!!$    implicit none
!!$
!!$    integer ncgmax,l
!!$    double precision, dimension(:,:) :: VTU,VTV,PRE,SMU,SMV,SMP
!!$
!!$    double precision, dimension(size(VTU,1),size(VTV,2)) :: pu,qu,ru,su,uu,vu,yu,wu
!!$    double precision, dimension(size(VTV,1),size(VTV,2)) :: pv,qv,rv,sv,uv,vv,yv,wv
!!$    double precision, dimension(size(PRE,1),size(PRE,2)) :: pp,qp,rp,sp,up,vp,yp,wp
!!$
!!$
!!$    double precision :: err0,err
!!$    integer :: ncg
!!$
!!$    call mavecxy(VTU,VTV,PRE,ru,rv,rp)
!!$
!!$    if (is_u)   then; ru=SMU-ru;   vu=ru; wu=ru; endif
!!$    if (is_v)   then; rv=SMV-rv;   vv=rv; wv=rv; endif
!!$    if (is_div) then; rp=SMP-rp;   vp=rp; wp=rp; endif
!!$
!!$    if (is_precond) then
!!$       call precond(ru,rv,rp,yu,yv,yp)
!!$    else
!!$       if (is_u) yu=ru;   if (is_v) yv=rv;   if (is_div) yp=rp
!!$    endif
!!$
!!$    rho=global_ddot(yu,yv,yp,yu,yv,yp)
!!$    phi=global_ddot(wu,wv,wp,wu,wv,wp)
!!$
!!$    gamma=1_prec;  nu=-1_prec
!!$
!!$    err=sqrt(phi)
!!$    err0=err
!!$    ncg=1
!!$    
!!$    do while (err>epsv*err0.and.err>epsva.and.ncg.lt.abs(ncgmax))
!!$
!!$       if (abs(rho)<1.d-30.or.abs(phi)<1.d-30) stop 'QMR Method fails'
!!$       alpha=1_prec/rho
!!$       if (is_u)   then;  vu=alpha*vu;  yu=alpha*yu;  wu=alpha*wu;  endif
!!$       if (is_v)   then;  vv=alpha*vv;  yv=alpha*yv;  wv=alpha*wv;  endif
!!$       if (is_div) then;  vp=alpha*vp;  yp=alpha*yp;  wp=alpha*wp;  endif 
!!$
!!$       delta=global_ddot(wu,wv,wp,yu,yv,yp)
!!$       if (abs(delta)<1.d-30) stop 'QMR Method fails'
!!$
!!$       if (is_precond) then
!!$          call precond(wu,wv,wp,zu,zv,zp)
!!$       else
!!$          if (is_u) zu=wu;   if (is_v) zv=wv;   if (is_div) zp=wp
!!$       endif
!!$
!!$       if (ncg==1) then
!!$          if (is_u) pu=yu;   if (is_v) pv=yv;   if (is_div) pp=yp
!!$          if (is_u) qu=zu;   if (is_v) qv=zv;   if (is_div) qp=zp
!!$       else
!!$          alpha=phi*delta/epsilon
!!$          if (is_u) pu=yu+alpha*pu;  if (is_v) pv=yv+alpha*pv;  if (is_div) pp=yp+alpha*pp
!!$          alpha=rho*delta/epsilon
!!$          if (is_u) qu=zu+alpha*qu;  if (is_v) qv=zv+alpha*qv;  if (is_div) qp=zp+alpha*qp
!!$       endif
!!$
!!$       call mavecxy(pu,pv,pp,zu,zv,zp)
!!$       epsilon=global_ddot(qu,qv,qp,zu,zv,zp)
!!$       if (abs(epsilon)<1.d-30) stop 'QMR Method fails'
!!$       beta=epsilon/delta
!!$       if (is_u) vu=zu-beta*vu;  if (is_v) vv=zv-beta*vv;  if (is_div) vp=zp-beta*vp
!!$
!!$       if (is_precond) then
!!$          call precond(vu,vv,vp,yu,yv,yp)
!!$       else
!!$          if (is_u) yu=vu;   if (is_v) yv=vv;   if (is_div) yp=vp
!!$       endif
!!$       rho=global_ddot(yu,yv,yp,yu,yv,yp)
!!$
!!$   BESOIN DE A'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$       err=sqrt(global_ddot(ru,rv,rp,ru,rv,rp))
!!$
!!$       if (is_print>=3*niv_solve-1) print *,'niveau ',niv_solve,'erreur QMR : ',ncg,err
!!$       if (is_cv_file.and.niv_solve==1) write (78,*) 'niveau ',niv_solve,'erreur QMR : ',ncg,err
!!$
!!$       ncg=ncg+1
!!$    enddo
!!$
!!$    if (is_u)   VTU=xu(:,:,0); 
!!$    if (is_v)   VTV=xv(:,:,0); 
!!$    if (is_div) PRE=xp(:,:,0); 
!!$    
!!$    if (is_print>=3*niv_solve-2) then
!!$       if (ncg.eq.ncgmax) then
!!$          print *,'Sorry! QMR n''a pas converge : ncg =',ncg,' err=',err
!!$       else
!!$          print *,'Bingo! QMR a converge : ncg =',ncg,' err=',err
!!$       endif
!!$    endif
!!$    
!!$    return
!!$
!!$  end subroutine solve_qmr
!!$
  !************************************************************************

!!$  !************************************************************************
!!$
!!$  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!!$  subroutine mavecxy_vf2(INPU,INPV,INPP,OUTU,OUTV,OUTP)
!!$    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!!$    !SK
!!$    !SK   Produit Matrice/Vecteur
!!$    !SK     l'equation de Poisson 2D -/\u=f  est discretisee
!!$    !SK     au 2e ordre en differences finies, soit
!!$    !SK             
!!$    !SK             -in    + 2 in   - in       -in    + 2 in   - in 
!!$    !SK               i,j-1     i,j    i,j+1     i-1,j     i,j    i+1,j
!!$    !SK    out    = ------------------------ + ------------------------
!!$    !SK      i,j                dx2                     dy2
!!$    !SK
!!$    !SK================================================================
!!$    !SK
!!$    use drapeaux
!!$    use data,   only : nu
!!$    use disc,   only : dx,dy,invdx,invdy,invdx2,invdy2,p5invdx,p5invdy &
!!$         ,coeft,dx,dy,U_EN_U, U_EN_V, V_EN_V, V_EN_U
!!$    implicit none
!!$    double precision, dimension(0:,0:), target  :: OUTU,OUTV,OUTP,INPU,INPV,INPP
!!$    double precision, dimension(0:size(INPU,1)-1,0:size(INPU,2)-1),target :: TMPU
!!$    double precision, dimension(0:size(INPV,1)-1,0:size(INPV,2)-1),target :: TMPV
!!$    double precision, dimension(:,:), pointer :: &
!!$         INPU_n,INPU_e,INPU_w,INPU_s,TMPU_u,TMPU_v,TMPU_n,TMPU_e,TMPU_w,TMPU_s &
!!$        ,INPV_n,INPV_e,INPV_w,INPV_s,TMPV_u,TMPV_v,TMPV_n,TMPV_e,TMPV_w,TMPV_s &
!!$        ,OUTU_c,OUTV_c,INPU_c,INPV_c
!!$
!!$    integer :: lmu,nmu,lmv,nmv
!!$    double precision :: coef
!!$
!!$    coef=2_prec/dx/dx+2_prec/dy/dy
!!$
!!$    lmu=size(INPU,1)-2;  nmu=size(INPU,2)-2
!!$    lmv=size(INPV,1)-2;  nmv=size(INPV,2)-2
!!$
!!$    OUTU_c   => OUTU(1:lmu,1:nmu);        OUTV_c   => OUTV(1:lmv,1:nmv)
!!$
!!$    INPU_c   => INPU(1:lmu  ,1:nmu);      INPV_c   => INPV(1:lmv  ,1:nmv);  
!!$    INPU_n   => INPU(1:lmu,2:nmu+1);      INPV_n   => INPV(1:lmv  ,2:nmv+1);
!!$    INPU_e   => INPU(2:lmu+1,1:nmu);      INPV_e   => INPV(2:lmv+1,1:nmv);  
!!$    INPU_w   => INPU(0:lmu-1,1:nmu);      INPV_w   => INPV(0:lmv-1,1:nmv);  
!!$    INPU_s   => INPU(1:lmu,0:nmu-1);      INPV_s   => INPV(1:lmv  ,0:nmv-1);
!!$    
!!$    TMPU_u => TMPU(1:lmu+1,:);            TMPV_u => TMPV(1:lmv+1,:);    
!!$    TMPU_v => TMPU(:,1:nmu+1);            TMPV_v => TMPV(:,1:nmv+1);    
!!$    TMPU_n => TMPU(1:lmu,2:nmu+1);        TMPV_n => TMPV(1:lmv,2:nmv+1);
!!$    TMPU_e => TMPU(2:lmu+1,1:nmu);        TMPV_e => TMPV(2:lmv+1,1:nmv);
!!$    TMPU_w => TMPU(1:lmu,1:nmu);          TMPV_w => TMPV(1:lmv,1:nmv);  
!!$    TMPU_s => TMPU(1:lmu,1:nmu);          TMPV_s => TMPV(1:lmv,1:nmv);  
!!$
!!$    
!!$    if (is_u) then    
!!$       if (is_t) then
!!$          OUTU=coeft*INPU
!!$       else 
!!$          OUTU=0_prec
!!$       endif
!!$
!!$       if (is_diff) OUTU_c = OUTU_c - nu*(-coef*INPU_c+invdx2*(INPU_e+INPU_w)+invdy2*(INPU_n+INPU_s))
!!$
!!$       if (is_conv) then
!!$          TMPU_u=U_EN_U(1:lmu+1,:)*(INPU(0:lmu,:)+INPU(1:lmu+1,:)) 
!!$          OUTU_c=OUTU_c+p5invdx*(TMPU_e-TMPU_w)
!!$          TMPU_v=V_EN_U(:,1:nmu+1)*(INPU(:,0:nmu)+INPU(:,1:nmu+1))
!!$          OUTU_c=OUTU_c+p5invdy*(TMPU_n-TMPU_s)
!!$       endif
!!$    endif
!!$
!!$    if (is_v) then    
!!$       if (is_t) then
!!$          OUTV=coeft*INPV
!!$       else
!!$          OUTV=0_prec
!!$       endif
!!$
!!$       if (is_diff) OUTV_c = OUTV_c - nu*(-coef*INPV_c+invdx2*(INPV_e+INPV_w)+invdy2*(INPV_n+INPV_s))
!!$
!!$       if (is_conv) then
!!$          TMPV_u=U_EN_V(1:lmv+1,:)*(INPV(0:lmv,:)+INPV(1:lmv+1,:)) 
!!$          OUTV_c=OUTV_c+p5invdx*(TMPV_e-TMPV_w)
!!$          TMPV_v=V_EN_V(:,1:nmv+1)*(INPV(:,0:nmv)+INPV(:,1:nmv+1))
!!$          OUTV_c=OUTV_c+p5invdy*(TMPV_n-TMPV_s)
!!$       endif
!!$    endif
!!$
!!$    if (is_grdp.or.is_grdv) stop 'NS pas encore implemente en VF2'
!!$
!!$       if (is_grdp.and.is_u)  OUTU_c = OUTU_c + p5invdx*(INPP_ne+INPP_se-INPP_nw-INPP_sw)
!!$
!!$       if (is_grdp.and.is_u)  OUTV_c = OUTV_c + p5invdy*(INPP_ne+INPP_nw-INPP_se-INPP_sw)
!!$
!!$       if (is_div.or.is_grdv) OUTP =            p5invdx*(INPU_ne+INPU_se-INPU_nw-INPU_sw) &
!!$                                              + p5invdy*(INPV_ne+INPV_nw-INPV_se-INPV_sw)
!!$       if (is_grdv) then
!!$          OUTU_c = OUTU_c - rlag*p5invdx*(OUTP_ne+OUTP_se-OUTP_nw-OUTP_sw)
!!$          OUTV_c = OUTV_c - rlag*p5invdy*(OUTP_ne+OUTP_nw-OUTP_se-OUTP_sw)
!!$          OUTP   = 0_prec
!!$       endif
!!$
!!$    call fixe_cl(INPU,INPV,OUTU,OUTV)
!!$
!!$    return
!!$  end subroutine mavecxy_vf2

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine solve_precond(SMU,SMV,SMP,VTU,VTV,PRE,ncgmax)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Appel direct a Precond
    !SK
    !SK================================================================
    !SK
    use uncol, only : nb_tasks
    implicit none

    integer ncgmax
    double precision, dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP
    double precision, dimension(0:size(SMU,1)-1,0:size(SMU,2)-1) :: ru
    double precision, dimension(0:size(SMV,1)-1,0:size(SMV,2)-1) :: rv
    double precision, dimension(0:size(SMP,1)-1,0:size(SMP,2)-1) :: rp

    double precision  :: err0,err
    integer           :: ncg

    call mavecxy(VTU,VTV,PRE,ru,rv,rp)
    if (is_u)  ru=SMU-ru;    if (is_v)  rv=SMV-rv;    if (is_div) rp=SMP-rp

    err0=sqrt(global_ddot(ru,rv,rp,ru,rv,rp))
    err =err0
    ncg =1

    do while (err>epsv*err0.and.err>epsva.and.ncg.lt.abs(ncgmax))

       call precond(SMU,SMV,SMP,VTU,VTV,PRE)
       call mavecxy(VTU,VTV,PRE,ru,rv,rp)
       if (is_u)  ru=SMU-ru;    if (is_v)  rv=SMV-rv;    if (is_div) rp=SMP-rp

       err=sqrt(global_ddot(ru,rv,rp,ru,rv,rp))
       if (ncg==1.and.nb_tasks>1) err0=err
       if (is_print>=3*niv_solve-1) print *,'niveau ',niv_solve,'erreur Precond only : ',ncg,err
       if (is_cv_file.and.niv_solve==1) write (78,*) 'niveau ',niv_solve,'erreur Precond only : ',ncg,err
       if (err<epsv*err0.or.err<epsva) exit

       ncg=ncg+1
    enddo

    if (is_print>=3*niv_solve-2) then
       if (ncg.eq.ncgmax) then
          print *,'Sorry! Precond only n''a pas converge : ncg =',ncg,' err=',err
       else
          print *,'Bingo! Precond only a converge : ncg =',ncg,' err=',err
       endif
    endif
    
    return

  end subroutine solve_precond



  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine inverse(SMU,SMW,SMP,VTU0,VTV0,PRE0,i0,k0)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Resolution du systeme par methode du pivot de Gauss
    !SK
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      implicit none

      real(kind=prec), dimension(0:,0:) :: VTU0,VTV0,PRE0,SMU,SMW,SMP

      real(kind=prec), dimension(0:size(VTU0,1)-1,0:size(VTU0,2)-1) :: VTU,T1X
      real(kind=prec), dimension(0:size(VTV0,1)-1,0:size(VTV0,2)-1) :: VTV,T1Y
      real(kind=prec), dimension(0:size(PRE0,1)-1,0:size(PRE0,2)-1) :: PRE,T1P

      integer :: i0,k0,k,i,i2,iflag=0
      real(kind=prec) :: a, pivot
      real(kind=prec), dimension(9,9) :: mat
      real(kind=prec), dimension(9)   :: b
      

      if (iflag.eq.0) then
         VTU=0._prec
         VTV=0._prec
         PRE=0._prec
!         iflag=1
!         print *,'mise a yero des tableaux'
      endif

      do i=1,9

         select case(i)
         case (1);   VTU(i0,k0)=1._prec
         case (2);   VTU(i0+1,k0)=1._prec
         case (3);   VTU(i0,k0+1)=1._prec
         case (4);   VTU(i0+1,k0+1)=1._prec
         case (5);   VTV(i0,k0)=1._prec
         case (6);   VTV(i0+1,k0)=1._prec
         case (7);   VTV(i0,k0+1)=1._prec
         case (8);   VTV(i0+1,k0+1)=1._prec
         case (9);   PRE(i0,k0)=1._prec
         end select


         call mavecxy_df20(VTU,VTV,PRE,T1X,T1Y,T1P,i0,k0)
         mat(1,i)=T1X(i0,k0)
         mat(2,i)=T1X(i0+1,k0)
         mat(3,i)=T1X(i0,k0+1)
         mat(4,i)=T1X(i0+1,k0+1)
         mat(5,i)=T1Y(i0,k0)
         mat(6,i)=T1Y(i0+1,k0)
         mat(7,i)=T1Y(i0,k0+1)
         mat(8,i)=T1Y(i0+1,k0+1)
         mat(9,i)=T1P(i0,k0)

         select case(i)
         case (1);   VTU(i0,k0)=0._prec
         case (2);   VTU(i0+1,k0)=0._prec
         case (3);   VTU(i0,k0+1)=0._prec
         case (4);   VTU(i0+1,k0+1)=0._prec
         case (5);   VTV(i0,k0)=0._prec
         case (6);   VTV(i0+1,k0)=0._prec
         case (7);   VTV(i0,k0+1)=0._prec
         case (8);   VTV(i0+1,k0+1)=0._prec
         case (9);   PRE(i0,k0)=0._prec
         end select
         
     enddo
         
     b(1)=SMU(i0,k0)
     b(2)=SMU(i0+1,k0)
     b(3)=SMU(i0,k0+1)
     b(4)=SMU(i0+1,k0+1)
     b(5)=SMW(i0,k0)
     b(6)=SMW(i0+1,k0)
     b(7)=SMW(i0,k0+1)
     b(8)=SMW(i0+1,k0+1)
     b(9)=SMP(i0,k0)

     do i=1,8
        if (abs(mat(i,i))<1.e-5) print '(I3,A,E10.3)',i,' pivot --> ',mat(i,i)
        do i2=i+1,9
           pivot=mat(i2,i)/mat(i,i)
           do k=i,9
              mat(i2,k)=mat(i2,k)-pivot*mat(i,k)
           enddo
           b(i2)=b(i2)-pivot*b(i)
        enddo
     enddo

     do i=9,1,-1
        a=0._prec
        do k=i+1,9
           a=a+mat(i,k)*b(k)
        enddo
        b(i)=(b(i)-a)/mat(i,i)
     enddo

     VTU0(i0,k0)     =b(1)
     VTU0(i0+1,k0)   =b(2)
     VTU0(i0,k0+1)   =b(3)
     VTU0(i0+1,k0+1) =b(4)
     VTV0(i0,k0)     =b(5)
     VTV0(i0+1,k0)   =b(6)
     VTV0(i0,k0+1)   =b(7)
     VTV0(i0+1,k0+1) =b(8)
     PRE0(i0,k0)     =b(9)

     return
  end subroutine inverse                                 ! end_out_light

  !************************************************************************

end module ajout

!!$                      e=1._prec/(d*dx2*a*b+c*dx2*a*b+c*dy2*d*b+c*dy2*d*a)
!!$                      VTU(i,k) = -(-c*d*dy2*beta+c*d*dy2*upsilon*b*dx-dy*c*delta*dx*b+mu  &
!!$                           & *dy*d*dx*b-alpha*d*dx2*b-alpha*c*dx2*b-alpha*c*dy2*d)*e
!!$                      VTU(i+1,k) = (alpha*c*dy2*d+c*d*dy2*a*upsilon*dx-dy*c*delta*a*dx+mu*  &
!!$                           & dy*d*dx*a+beta*d*dx2*a+beta*c*dx2*a+c*d*dy2*beta)*e
!!$                      VTV(i,k) = (-dx*d*dy*alpha*b+dx*d*dy*a*beta-d*dy*a*upsilon*b*dx2+dx2  &
!!$                           & *a*b*mu+mu*dy2*d*b+mu*dy2*d*a+delta*a*dx2*b)*e
!!$                      VTV(i,k+1) = (dx*c*dy*alpha*b+delta*c*dy2*a-dx*c*dy*a*beta+delta*c*dy2  &
!!$                           & *b+dx2*a*b*mu+c*dy*a*upsilon*b*dx2+delta*a*dx2*b)*e
!!$                      PRE(i,k) = dy*dx*(c*d*dy*alpha*b-c*d*dy*a*beta+c*d*dy*a*upsilon*b*dx-c  &
!!$                           & *delta*a*dx*b+d*dx*a*b*mu)*e
!!$


!!$                      e=0.5_prec/a**3/(dx2+dy2)
!!$                      VTU(i,k) = -(-c*d*dy**2*beta+c*d*dy**2*upsilon*b*dx-dy*c*delta*dx*b+mu  &
!!$                           & *dy*d*dx*b-alpha*d*dx**2*b-alpha*c*dx**2*b-alpha*c*dy**2*d)*e
!!$                      VTU(i+1,k) = (alpha*c*dy**2*d+c*d*dy**2*a*upsilon*dx-dy*c*delta*a*dx+mu*  &
!!$                           & dy*d*dx*a+beta*d*dx**2*a+beta*c*dx**2*a+c*d*dy**2*beta)*e
!!$                      VTV(i,k) = (-dx*d*dy*alpha*b+dx*d*dy*a*beta-d*dy*a*upsilon*b*dx**2+dx**2  &
!!$                           & *a*b*mu+mu*dy**2*d*b+mu*dy**2*d*a+delta*a*dx**2*b)*e
!!$                      VTV(i,k+1) = (dx*c*dy*alpha*b+delta*c*dy**2*a-dx*c*dy*a*beta+delta*c*dy**2  &
!!$                           & *b+dx**2*a*b*mu+c*dy*a*upsilon*b*dx**2+delta*a*dx**2*b)*e
!!$                      PRE(i,k) = dy*dx*(c*d*dy*alpha*b-c*d*dy*a*beta+c*d*dy*a*upsilon*b*dx-c  &
!!$                           & *delta*a*dx*b+d*dx*a*b*mu)*e

!!$                   alpha  =alpha  -(a*VTU(i,k)  +invdx*PRE(i,k))
!!$                   beta   =beta   -(b*VTU(i+1,k)-invdx*PRE(i,k))
!!$                   mu     =mu     -(c*VTV(i,k)  +invdy*PRE(i,k))
!!$                   delta  =delta  -(d*VTV(i,k+1)-invdy*PRE(i,k))
!!$                   upsilon=upsilon-(-invdx*VTU(i,k)+invdx*VTU(i+1,k) &
!!$                        &              -invdy*VTV(i,k)+invdy*VTV(i,k+1))
!!$                   print '(5e15.3)',alpha,beta,mu,delta,upsilon


form zblas.f90................

    if (1>4.and.is_div) then

       if (ncheck/=df2) stop 'df2 please'

       TMPU=0._prec; OUTU=0._prec
       where (U_en_U<0.) TMPU=1.
       where (V_en_U<0.) OUTU=1.

       TMPV=0._prec; OUTV=0._prec
       where (U_en_V<0.) TMPV=1.
       where (V_en_V<0.) OUTV=1.


       !!  TMPU=0.5_prec; TMPV=0.5_prec; 
       !!  OUTU=0.5_prec; OUTV=0.5_prec; 

       do k=1,nmu
          do i=1,lmu
             OUTU(i,k)=(coef*nu+coeft)*INPU(i,k)&
                  &                  -nu*((INPU(i+1,k)+INPU(i-1,k))*invdx2          & 
                  +(INPU(i,k+1)+INPU(i,k-1))*invdy2)         &
                  + U_EN_U(i,k)*(TMPU(i,k)*(INPU(i+1,k)-INPU(i,k))&
                  &             +(1._prec-TMPU(i,k))*(INPU(i,k)-INPU(i-1,k)))*invdx &
                  + V_EN_U(i,k)*(OUTU(i,k)*(INPU(i,k+1)-INPU(i,k)) &
                  &             +(1._prec-OUTU(i,k))*(INPU(i,k)-INPU(i,k-1)))*invdy &
                  + invdx*(INPP(i,k)-INPP(i-1,k))
          enddo
       enddo


       do k=1,nmv
          do i=1,lmv
             OUTV(i,k)=(coef*nu+coeft)*INPV(i,k)&
                  &                  -nu*((INPV(i+1,k)+INPV(i-1,k))*invdx2          & 
                  +(INPV(i,k+1)+INPV(i,k-1))*invdy2)         &
                  + U_EN_V(i,k)*(TMPV(i,k)*(INPV(i+1,k)-INPV(i,k))&
                  &             +(1._prec-TMPV(i,k))*(INPV(i,k)-INPV(i-1,k)))*invdx &
                  + V_EN_V(i,k)*(OUTV(i,k)*(INPV(i,k+1)-INPV(i,k)) &
                  &             +(1._prec-OUTV(i,k))*(INPV(i,k)-INPV(i,k-1)))*invdy &
                  + invdy*(INPP(i,k)-INPP(i,k-1))
          enddo
       enddo


       do k=0,nmp+1
          do i=0,lmp+1
             OUTP(i,k)=invdx*(INPU(i+1,k)-INPU(i,k))+invdy*(INPV(i,k+1)-INPV(i,k))
          enddo
       enddo

       call fixe_cl(INPU,INPV,INPP,OUTU,OUTV,OUTP)

    else






       if (is_oned) then
          SMU=vit_adv; VTV=0._prec; 
          call prepare_adv_grids(SMU,VTV)
       else if (nexample==17) then
          call exacte
!!$          VTUS=8.*XU/dx+YU/dy
!!$          VTVS=XV/dx+2.*YV/dy
          call prepare_adv_grids(VTUS,VTVS)
       else
          call prepare_adv_grids(VTU,VTV)
       endif

!SK       lmu=size(VTU,1)-2
!SK       nmu=size(VTU,2)-2
!SK       lmv=size(VTV,1)-2
!SK       nmv=size(VTV,2)-2
!SK
!SK       call set_grid(-1)
!SK
!!$       TMPx(1:lmu+1,1:nmu)=8.*(XU(1:lmu+1,1:nmu)-0.5*dx)/dx+YU(1:lmu+1,1:nmu)/dy
!SK       TMPx(1:lmu+1,1:nmu)=10*(t_start-YU(1:lmu+1,1:nmu))&
!SK            & *exp(-1/nu*((t_start-XU(1:lmu+1,1:nmu)+0.5*dx)**2+(t_start-YU(1:lmu+1,1:nmu))**2))
!SK       call outdon(TMPx(1:lmu+1,1:nmu)-U_en_U(1:lmu+1,1:nmu),'diff u_en_u')
!SK
!SK!!$       TMPy(1:lmv,1:nmv+1)=XV(1:lmv,1:nmv+1)/dx+2.*(YV(1:lmv,1:nmv+1)-0.5*dy)/dy
!SK       TMPy(1:lmv,1:nmv+1)=10*(XV(1:lmv,1:nmv+1)-t_start)&
!SK            & *exp(-1/nu*((t_start-XV(1:lmv,1:nmv+1))**2+(t_start-YV(1:lmv,1:nmv+1)+0.5*dy)**2))
!SK       call outdon(TMPy(1:lmv,1:nmv+1)-V_en_V(1:lmv,1:nmv+1),'diff en v_en_v')
!SK
!SK!!$       TMPy(1:lmv,1:nmv)=8.*XU(1:lmv,1:nmv)/dx+(YU(1:lmv,1:nmv)-0.5*dy)/dy
!SK       TMPy(1:lmv,1:nmv)=10*(t_start-YU(1:lmv,1:nmv)+0.5*dy)&
!SK            & *exp(-1/nu*((t_start-XU(1:lmv,1:nmv))**2+(t_start-YU(1:lmv,1:nmv)+0.5*dy)**2))
!SK       call outdon(TMPy(1:lmv,1:nmv)-U_en_V(1:lmv,1:nmv),'diff en u_en_v')
!SK
!SK!!$       TMPx(1:lmu,1:nmu)=(XV(1:lmu,1:nmu)-0.5*dx)/dx+2.*YV(1:lmu,1:nmu)/dy
!SK       TMPx(1:lmu,1:nmu)=10*(XV(1:lmu,1:nmu)-0.5*dx-t_start&
!SK            & )*exp(-1/nu*((t_start-XV(1:lmu,1:nmu)+0.5*dx)**2+(t_start-YV(1:lmu,1:nmu))**2))
!SK       call outdon(TMPx(1:lmu,1:nmu)-v_en_U(1:lmu,1:nmu),'diff en v_en_u')
!SK!!$
!SK!!$       call outdon(VTUS,'VTU')
!SK!!$       call outdon(VTVS,'VTV')
!SK!!$       call outdon(v_en_U(1:lmu,1:nmu),' v_en_u')
!SK!!$       call outdon(TMPx(1:lmu,1:nmu),'ce qu''il faut')
!SK!!$
!SK       call parallel_stop


       call set_grid(-1)
       call set_switch(uv_switch*conv_code*diff_code/is_impl)
       call mavecxy(VTU,VTV,PRE,VTU4,VTV4,PRE)

       call exacte
