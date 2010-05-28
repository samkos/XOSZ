
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  Routines Multigrilles
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mg

  use constante
  use data,     only : is_print,epsv,epsva
  use conv
  use para

contains

  !************************************************************************  ! start_out_light

  !SKffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff       
  subroutine restrict(INPU,INPV,INPP,OUTU,OUTV,OUTP,flag)
    !SKffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    use drapeaux, only : is_u,is_v,is_div,is_decal
    implicit none
    real(kind=prec), dimension(0:,0:) :: INPU,INPV,INPP,OUTU,OUTV,OUTP
    integer :: lmp,nmp,lmpr,nmpr
    logical :: flag

    if (is_decal) then
       lmp=size(INPV,1)-2;  nmp=size(INPU,2)-2
       lmpr=size(OUTV,1)-2; nmpr=size(OUTU,2)-2

!!$       call output('lmp,nmp',lmp,nmp)
!!$       call output('lmpr,nmpr',lmpr,nmpr)
!!$
!!$       call output('size(INPU,1),size(INPU,2)',size(INPU,1),size(INPU,2))
!!$       call output('size(INPV,1),size(INPP,2)',size(INPV,1),size(INPV,2))
!!$       call output('size(INPP,1),size(INPP,2)',size(INPP,1),size(INPP,2))
!!$       call output('size(OUTU,1),size(OUTU,2)',size(OUTU,1),size(OUTU,2))
!!$       call output('size(OUTV,1),size(OUTP,2)',size(OUTV,1),size(OUTV,2))
!!$       call output('size(OUTP,1),size(OUTP,2)',size(OUTP,1),size(OUTP,2))

       if (is_div.and..not.flag) OUTP=INPP(::3,::3)

       if (is_u.or.flag) then
          OUTU(1:lmpr+1,:)=INPU(2:lmp:3,0::3)
          OUTU(0,:)       =INPU(0,0::3)
          OUTU(lmpr+2,:)  =INPU(lmp+2,0::3)
       endif
       if (is_v.or.flag) then
          OUTV(:,1:nmpr+1)=INPV(0::3,2:nmp:3)
          OUTV(:,0)       =INPV(0::3,0)
          OUTV(:,nmpr+2)  =INPV(0::3,nmp+2)
       endif
    else
       if (is_div) stop 'restrict not yet implemented for ns equation'
       if (is_u) OUTU=INPU(::2,::2)
       if (is_v) OUTV=INPV(::2,::2)
    endif

    return
  end subroutine  restrict

  !************************************************************************

  !SKffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
  subroutine prolong(INPU,INPV,INPP,OUTU,OUTV,OUTP)
    !SKffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    use drapeaux, only : is_u,is_v,is_div,is_decal
    implicit none
    real(kind=prec), dimension(0:,0:) :: INPU,INPV,INPP,OUTU,OUTV,OUTP
    integer :: lm,nm,i,k,lmp,nmp,lmpl,nmpl

    if (is_decal) then
       lmp=size(INPV,1)-2; lmpl=size(OUTV,1)-2
       nmp=size(INPU,2)-2; nmpl=size(OUTU,2)-2

       if (is_div) then
          OUTP(0::3,0::3)=INPP;
          OUTP(1::3,0::3)=(2._prec*INPP(0:lmp,:)+INPP(1:lmp+1,:))/3._prec
          OUTP(2::3,0::3)=(INPP(0:lmp,:)+2._prec*INPP(1:lmp+1,:))/3._prec
          OUTP(:,1::3)   =(2._prec*OUTP(:,0:nmpl:3)+OUTP(:,3::3))/3._prec
          OUTP(:,2::3)   =(OUTP(:,0:nmpl:3)+2._prec*OUTP(:,3::3))/3._prec
       endif

       if (is_u) then
          OUTU=0.
          OUTU(2:lmpl:3,0::3)=INPU(1:lmp+1,:)
          OUTU(3:lmpl-2:3,0::3)=(2_prec*INPU(1:lmp,:)+INPU(2:lmp+1,:))/3._prec
          OUTU(4:lmpl-1:3,0::3)=(INPU(1:lmp,:)+2_prec*INPU(2:lmp+1,:))/3._prec
!!$          OUTU(1,0::3)       =2._prec*OUTU(2,0::3)-OUTU(3,0::3)
!!$          OUTU(lmpl+1,0::3)  =2_prec*OUTU(lmpl,0::3)-OUTU(lmpl-1,0::3)
          OUTU(1,0::3)       =0._prec
          OUTU(lmpl+1,0::3)  =0._prec
          OUTU(:,1::3)       =(2._prec*OUTU(:,0:nmpl:3)+OUTU(:,3:nmpl+1:3))/3._prec
          OUTU(:,2::3)       =(OUTU(:,0:nmpl:3)+2._prec*OUTU(:,3:nmpl+1:3))/3._prec
       endif

       if (is_v) then
          OUTV=0.
          OUTV(0::3,2:nmpl:3)=INPV(:,1:nmp+1)
          OUTV(0::3,3:nmpl-1:3)=(2._prec*INPV(:,1:nmp)+INPV(:,2:nmp+1))/3._prec
          OUTV(0::3,4:nmpl-1:3)=(INPV(:,1:nmp)+2._prec*INPV(:,2:nmp+1))/3._prec
!!$          OUTV(0::3,1)       =2._prec*OUTV(0::3,2)-OUTV(0::3,3)       
!!$          OUTV(0::3,nmpl+1)  =2._prec*OUTV(0::3,nmpl)-OUTV(0::3,nmpl-1)       
          OUTV(0::3,1)       = 0._prec
          OUTV(0::3,nmpl+1)  = 0._prec
          OUTV(1::3,:)       =(2._prec*OUTV(0:lmpl:3,:)+OUTV(3:lmpl+1:3,:))/3._prec
          OUTV(2::3,:)       =(OUTV(0:lmpl:3,:)+2._prec*OUTV(3:lmpl+1:3,:))/3._prec
       endif
    else
       lm=size(OUTU,1)-1
       nm=size(OUTU,2)-1

       if (is_div) stop 'restrict not yet implemented for ns equation'

       if (is_u) then
          OUTU(0:lm:2,0:nm:2)=INPU;

          do k=0,nm/2-1
             do i=0,lm/2-1
                OUTU(i*2+1,k*2)=0.5_prec*(INPU(i,k)+INPU(i+1,k))
                OUTU(i*2,k*2+1)=0.5_prec*(INPU(i,k)+INPU(i,k+1))
                OUTU(i*2+1,k*2+1)=0.25_prec*(INPU(i,k)+INPU(i+1,k)+INPU(i,k+1)+INPU(i+1,k+1))
             enddo
          enddo

          do i=0,lm/2-1
             OUTU(i*2+1,nm)=0.5_prec*(INPU(i,nm/2)+INPU(i+1,nm/2))
          enddo
          do k=0,nm/2-1
             OUTU(lm,k*2+1)=0.5_prec*(INPU(lm/2,k)+INPU(lm/2,k+1))
          enddo
       endif

       if (is_v) then
          OUTV(0:lm:2,0:nm:2)=INPV;

          do k=0,nm/2-1
             do i=0,lm/2-1
                OUTV(i*2+1,k*2)=0.5_prec*(INPV(i,k)+INPV(i+1,k))
                OUTV(i*2,k*2+1)=0.5_prec*(INPV(i,k)+INPV(i,k+1))
                OUTV(i*2+1,k*2+1)=0.25_prec*(INPV(i,k)+INPV(i+1,k)+INPV(i,k+1)+INPV(i+1,k+1))
             enddo
          enddo

          do i=0,lm/2-1
             OUTV(i*2+1,nm)=0.5_prec*(INPV(i,nm/2)+INPV(i+1,nm/2))
          enddo
          do k=0,nm/2-1
             OUTV(lm,k*2+1)=0.5_prec*(INPV(lm/2,k)+INPV(lm/2,k+1))
          enddo
       endif

    endif

    return
  end  subroutine prolong



  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine lisse(SMU,SMV,SMP,VTU,VTV,PRE,nb_lis)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   lissage de type Gauss-Seidel
    !SK
    !SK================================================================
    !SK
    use blas
    use data,   only : nu,sor_theta,un_ms_sor_theta, ncheck, ntype_solver,ntype_precond
    use disc,   only : dx,dy,invdx,invdy,invdx2,invdy2,p5invdx,p5invdy &
         ,coeft,U_EN_U, U_EN_V, V_EN_V, V_EN_U
    use drapeaux
    use champs, only : VTUS,VTVS
    implicit none
    real(kind=prec), dimension(0:,0:), target :: VTU,VTV,PRE
    real(kind=prec), dimension(0:,0:)         :: SMU,SMV,SMP
    real(kind=prec) :: a,b,c,d,e,alpha,beta,mu,delta,upsilon,dx2,dy2
    integer :: nb_lis

    real(kind=prec), dimension(0:size(VTU,1)-1,0:size(VTU,2)-1) :: TMPU
    real(kind=prec), dimension(0:size(VTV,1)-1,0:size(VTV,2)-1) :: TMPV
    real(kind=prec), dimension(0:size(PRE,1)-1,0:size(PRE,2)-1) :: TMPP

    real(kind=prec) :: divise, err,err0
    integer :: lmu,nmu,lmv,nmv,lmp,nmp,i,j,k,l,k0,i0,i1,k1


    lmu=size(SMU,1)-2;     nmu=size(SMU,2)-2
    lmv=size(SMV,1)-2;     nmv=size(SMV,2)-2
    lmp=size(SMP,1)-2;     nmp=size(SMP,2)-2

    divise=1._prec/(coeft+nu*(2._prec/dx/dx+2._prec/dy/dy))
    a=1._prec/divise; b=a;  c=a; d=a
    dx2=dx*dx; dy2=dy*dy
    e=0.5_prec/a/(dx2+dy2)

    if (is_decal) then
       VTU(1:lmu,0)=SMU(1:lmu,0); VTU(1:lmu,nmu+1)=SMU(1:lmu,nmu+1); 
       VTV(0,1:nmv)=SMV(0,1:nmv); VTV(lmv+1,1:nmv)=SMV(lmv+1,1:nmv); 
    else
       VTU(0,:)=SMU(0,:); VTU(lmu+1,:)=SMU(lmu+1,:); 
       VTU(:,0)=SMU(:,0); VTU(:,nmu+1)=SMU(:,nmu+1); 
       VTV(0,:)=SMV(0,:); VTV(lmv+1,:)=SMV(lmv+1,:); 
       VTV(:,0)=SMV(:,0); VTV(:,nmv+1)=SMV(:,nmv+1); 
    endif
       
    if (is_ns) then
       if (ncheck/=df4.and.ncheck/=df2) &
            stop 'pour l instant lisse demande DF en NS pour U_en_U'
       
       if (.not.is_decal) print *, 'lisse only implemented for staggered mesh'


       l=0
       err=1.
       err0=1.
    
       do while (err>epsv*err0.and.err>epsva.and.l.lt.abs(nb_lis))
          l=l+1
          if (1>2.and.is_print>=3*niv_solve-1) then
             call mavecxy(VTU,VTV,PRE,TMPU,TMPV,TMPP)
       
             TMPU=SMU-TMPU; TMPV=SMV-TMPV; TMPP=SMP-TMPP;
             err=sqrt(global_ddot(TMPU,TMPV,TMPP,TMPU,TMPV,TMPP))
             if (l==1) err0=err
             if (my_task==0) print *,'niveau ',niv_solve,' erreur lisse : ',l,err
          endif
          do j=-1,1,2
             k0=0; i0=0
             if (j==-1) then; k0=nmp+1; i0=lmp+1; endif
             do k=k0,nmp+1-k0,j
                do i=i0,lmp+1-i0,j
                   if (k+i==0) then
                      VTU(0,0)=2_prec*SMU(0,0)-VTU(1,0)
                      VTV(0,0)=2_prec*SMV(0,0)-VTV(0,1)
                   else if (k==0.and.i==lmp+1) then
                      VTU(lmp+2,0)=2_prec*SMU(lmp+2,0)-VTU(lmp+1,0)
                      VTV(lmp+1,0)=2_prec*SMV(lmp+1,0)-VTV(lmp+1,1)
                   else if (i==0.and.k==nmp+1) then
                      VTU(0,nmp+1)=2_prec*SMU(0,nmp+1)-VTU(1,nmp+1)
                      VTV(0,nmp+2)=2_prec*SMV(0,nmp+2)-VTV(0,nmp+1)
                   else if (i==lmp+1.and.k==nmp+1) then
                      VTU(lmp+2,nmp+1)=2_prec*SMU(lmp+2,nmp+1)-VTU(lmp+1,nmp+1)
                      VTV(lmp+1,nmp+2)=2_prec*SMV(lmp+1,nmp+2)-VTV(lmp+1,nmp+1)
                   else if (k==0) then
                      VTV(i,1)= SMV(i,0)-dy*p5invdx*(VTU(i+1,k)-VTU(i,k)-dx*SMP(i,k))
                      VTV(i,0)= 2._prec*SMV(i,0)-VTV(i,1)
                      k1=1
                      PRE(i,0)=dy*(-SMV(i,k1) +                                      &
                            a*VTV(i,k1)-nu*((VTV(i+1,k1)+VTV(i-1,k1))*invdx2          & 
                           +(VTV(i,k1+1)+VTV(i,k1-1))*invdy2)         &
                           + U_EN_V(i,k1)*(VTV(i+1,k1)-VTV(i-1,k1))*p5invdx &
                           + V_EN_V(i,k1)*(VTV(i,k1+1)-VTV(i,k1-1))*p5invdy &
                           + invdy*PRE(i,k1) )       
                   else if (k==nmp+1) then
                      VTV(i,nmp+1)= SMV(i,nmp+2)+dy*p5invdx*(VTU(i+1,k)-VTU(i,k)-dx*SMP(i,k))
                      VTV(i,nmp+2)= 2._prec*SMV(i,nmp+2)-VTV(i,nmp+1)
                      k1=nmp+1
                      PRE(i,nmp+1)=dy*(SMV(i,k1) -                                      &
                           (a*VTV(i,k1)-nu*((VTV(i+1,k1)+VTV(i-1,k1))*invdx2          & 
                           +(VTV(i,k1+1)+VTV(i,k1-1))*invdy2)         &
                           + U_EN_V(i,k1)*(VTV(i+1,k1)-VTV(i-1,k1))*p5invdx &
                           + V_EN_V(i,k1)*(VTV(i,k1+1)-VTV(i,k1-1))*p5invdy &
                           + invdy*(-PRE(i,nmp) )))        
                   else if (i==0) then
                      VTU(1,k)= SMU(0,k)-dx*p5invdy*(VTV(i,k+1)-VTV(i,k)-dy*SMP(i,k))
                      VTU(0,k)= 2._prec*SMU(0,k)-VTU(1,k)
                      i1=1
                      PRE(0,k)= dx*(-SMU(i1,k) +                                    &
                            a*VTU(i1,k)-nu*((VTU(i1+1,k)+VTU(i1-1,k))*invdx2          & 
                                +(VTU(i1,k+1)+VTU(i1,k-1))*invdy2)         &
                            + U_EN_U(i1,k)*(VTU(i1+1,k)-VTU(i1-1,k))*p5invdx &
                            + V_EN_U(i1,k)*(VTU(i1,k+1)-VTU(i1,k-1))*p5invdy &
                            + invdx*PRE(i1,k))
                   else if (i==lmp+1) then
                      VTU(lmp+1,k)= SMU(lmp+2,k)+dx*p5invdy*(VTV(i,k+1)-VTV(i,k)-dy*SMP(i,k))
                      VTU(lmp+2,k)= 2._prec*SMU(lmp+2,k) - VTU(lmp+1,k)
                      i1=lmp+1
                      PRE(lmp+1,k)=dx*(SMU(i1,k) -                         &
                           (a*VTU(i1,k)-nu*((VTU(i1+1,k)+VTU(i1-1,k))*invdx2          & 
                                +(VTU(i1,k+1)+VTU(i1,k-1))*invdy2)         &
                            + U_EN_U(i1,k)*(VTU(i1+1,k)-VTU(i1-1,k))*p5invdx &
                            + V_EN_U(i1,k)*(VTU(i1,k+1)-VTU(i1,k-1))*p5invdy &
                            + invdx*(-PRE(lmp,k))))
                   else if (1==2.or.i*k==0.or.i==lmp+1.or.k==nmp+1) then
                      VTU(i:i+1,k)=0._prec; VTV(i,k:k+1)=0._prec;  PRE(i,k)=0._prec
                      call mavecxy_df20(VTU,VTV,PRE,TMPU,TMPV,TMPP,i,k)
                      TMPU(i:i+1,k)   = SMU(i:i+1,k) -TMPU(i:i+1,k) 
                      TMPV(i,k:k+1)   = SMV(i,k:k+1) -TMPV(i,k:k+1) 
                      TMPP(i    ,k)   = SMP(i    ,k) -TMPP(i    ,k)    
                      call inverse(TMPU,TMPV,TMPP,VTU,VTV,PRE,i,k)
                   else
                      i1=i
                      alpha  = SMU(i1,k) -                                    &
                           (-nu*((VTU(i1+1,k)+VTU(i1-1,k))*invdx2          & 
                                +(VTU(i1,k+1)+VTU(i1,k-1))*invdy2)         &
                            + U_EN_U(i1,k)*(VTU(i1+1,k)-VTU(i1-1,k))*p5invdx &
                            + V_EN_U(i1,k)*(VTU(i1,k+1)-VTU(i1,k-1))*p5invdy &
                            + invdx*(-PRE(i-1,k)))
                      i1=i+1
                      beta   =SMU(i1,k) -                                      &
                           (-nu*((VTU(i1+1,k)+VTU(i1-1,k))*invdx2          & 
                                +(VTU(i1,k+1)+VTU(i1,k-1))*invdy2)         &
                            + U_EN_U(i1,k)*(VTU(i1+1,k)-VTU(i1-1,k))*p5invdx &
                            + V_EN_U(i1,k)*(VTU(i1,k+1)-VTU(i1,k-1))*p5invdy &
                            + invdx*(PRE(i+1,k) ))
                      k1=k
                      mu     =SMV(i,k1) -                                      &
                           (-nu*((VTV(i+1,k1)+VTV(i-1,k1))*invdx2          & 
                                +(VTV(i,k1+1)+VTV(i,k1-1))*invdy2)         &
                            + U_EN_V(i,k1)*(VTV(i+1,k1)-VTV(i-1,k1))*p5invdx &
                            + V_EN_V(i,k1)*(VTV(i,k1+1)-VTV(i,k1-1))*p5invdy &
                            + invdy*(-PRE(i,k-1) ))
                      k1=k+1
                      delta  =SMV(i,k1) -                                      &
                           (-nu*((VTV(i+1,k1)+VTV(i-1,k1))*invdx2          & 
                                +(VTV(i,k1+1)+VTV(i,k1-1))*invdy2)         &
                            + U_EN_V(i,k1)*(VTV(i+1,k1)-VTV(i-1,k1))*p5invdx &
                            + V_EN_V(i,k1)*(VTV(i,k1+1)-VTV(i,k1-1))*p5invdy &
                            + invdy*(PRE(i,k+1) ))
                      upsilon=SMP(i,k)  

                      b = dy2*upsilon*a*dx-dy*delta*dx+mu*dy*dx
                      VTU(i,k) =    un_ms_sor_theta*VTU(i,k)+sor_theta*&
                           &        (dy2*beta-b+alpha*dx2+alpha*dx2+alpha*dy2)*e
                      VTU(i+1,k) = un_ms_sor_theta*VTU(i+1,k)+sor_theta*&
                           &        (alpha*dy2+b+beta*dx2+beta*dx2+dy2*beta)*e
                      c = dx*dy*alpha+dy*upsilon*a*dx2-dx*dy*beta
                      VTV(i,k) =   un_ms_sor_theta*VTV(i,k)+sor_theta*&
                           &        (-c+dx2*mu+mu*dy2+mu*dy2+delta*dx2)*e
                      VTV(i,k+1) = un_ms_sor_theta*VTV(i,k+1)+sor_theta*&
                           &        (c+delta*dy2+delta*dy2+dx2*mu+delta*dx2)*e
                      PRE(i,k) = un_ms_sor_theta*PRE(i,k)+sor_theta*&
                           &           a*dy*dx*(dy*(alpha-beta+upsilon*a*dx)-delta*dx+dx*mu)*e

                      endif

                enddo
             enddo
          enddo
       enddo

       if (is_decal.and.is_mpp) stop 'no good for mpp and staggered mesh'
       PRE(0,0)=2.*PRE(0,1)-PRE(0,2)
       PRE(0,nmp+1)=2.*PRE(1,nmp+1)-PRE(2,nmp+1)
       PRE(lmp+1,nmp+1)=2.*PRE(lmp,nmp+1)-PRE(lmp-1,nmp+1)
       PRE(lmp+1,0)=2.*PRE(lmp,0)-PRE(lmp-1,0)
       
    else if (is_diff.and.is_conv) then
       do l=1,abs(nb_lis)
          do k=1,nmu
             do i=1,lmu
                VTU(i,k)=un_ms_sor_theta*VTU(i,k)+sor_theta*(      &
                     SMU(i,k)                                      &
                     + nu*((VTU(i+1,k)+VTU(i-1,k))*invdx2          & 
                     +(VTU(i,k+1)+VTU(i,k-1))*invdy2)         &
                     - U_EN_U(i,k)*(VTU(i+1,k)-VTU(i-1,k))*p5invdx &
                     - V_EN_V(i,k)*(VTU(i,k+1)-VTU(i,k-1))*p5invdy &
                     )*divise
             enddo
          enddo
          do k=1,nmv
             do i=1,lmv 
                VTV(i,k)=un_ms_sor_theta*VTV(i,k)+sor_theta*(      &
                     SMV(i,k)                                      &
                     + nu*((VTV(i+1,k)+VTV(i-1,k))*invdx2          &
                     +(VTV(i,k+1)+VTV(i,k-1))*invdy2)         &
                     - U_EN_U(i,k)*(VTV(i+1,k)-VTV(i-1,k))*p5invdx &
                     - V_EN_V(i,k)*(VTV(i,k+1)-VTV(i,k-1))*p5invdy &
                     )*divise
             enddo
          enddo
       enddo
    else if (is_diff) then
       do l=1,abs(nb_lis)
          if (is_decal) then
             VTU(0,:)=2._prec*SMU(0,:)-VTU(1,:); 
             VTV(:,0)=2._prec*SMV(:,0)-VTV(:,1); 
          endif
          do k=1,nmu
             do i=1,lmu
                VTU(i,k)=un_ms_sor_theta*VTU(i,k)+sor_theta*(      &
                     SMU(i,k)                                      &
                     +nu*((VTU(i+1,k)+VTU(i-1,k))*invdx2           &
                     +    (VTU(i,k+1)+VTU(i,k-1))*invdy2)          &
                     )*divise
             enddo
          enddo
          do k=1,nmv
             do i=1,lmv
                VTV(i,k)=un_ms_sor_theta*VTU(i,k)+sor_theta*(      &
                     SMV(i,k)                                      &
                     +nu*((VTV(i+1,k)+VTV(i-1,k))*invdx2           &
                     +(VTV(i,k+1)+VTV(i,k-1))*invdy2)          &
                     )*divise
             enddo
          enddo
          if (is_decal) then
             VTU(lmu+1,:)=2._prec*SMU(lmu+1,:)-VTU(lmu,:); 
             VTV(:,nmv+1)=2._prec*SMV(:,nmv+1)-VTV(:,nmv); 
          endif
       enddo
    else
       stop 'Hum!!! Lisseur pas tres complet!!!'
    endif

    return
  end subroutine lisse

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine mavecxy_df20(INPU,INPV,INPP,OUTU,OUTV,OUTP,i0,k0)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Produit Matrice/Vecteur
    !SK     l'equation de Poisson 2D -/\u=f  est discretisee
    !SK     au 2e ordre en differences finies, soit
    !SK             
    !SK             -in    + 2 in   - in       -in    + 2 in   - in 
    !SK               i,j-1     i,j    i,j+1     i-1,j     i,j    i+1,j
    !SK    out    = ------------------------ + ------------------------
    !SK      i,j                dx2                     dy2
    !SK
    !SK================================================================
    !SK
    use drapeaux
    use para
    use blas,   only : fixe_cl
    use data,   only : nu, rlag
    use disc,   only : dx,dy,invdx,invdy,invdx2,invdy2,p5invdx,p5invdy &
         ,coeft,U_EN_U, U_EN_V, V_EN_V, V_EN_U
    implicit none

    real(kind=prec), dimension(0:,0:), target :: OUTU,OUTV,OUTP,INPU,INPV,INPP

    real(kind=prec) :: coef
    integer lmu,nmu,lmv,nmv,lmp,nmp,i,k,i0,k0

    lmu=size(INPU,1)-2; nmu=size(INPU,2)-2
    lmv=size(INPV,1)-2; nmv=size(INPV,2)-2
    lmp=size(INPP,1)-2; nmp=size(INPP,2)-2
    coef=2._prec/dx/dx+2._prec/dy/dy

    do k=max(k0,1),min(k0+1,nmp)
       do i=max(i0,1),min(i0+1,lmp+1)

          if (is_t) then
             OUTU(i,k)=coeft*INPU(i,k) 
          else
             OUTU(i,k)=0.d0
          endif

          if (is_diff) OUTU(i,k)=OUTU(i,k)&
               - nu*(-coef*INPU(i,k)+invdx2*(INPU(i+1,k)+INPU(i-1,k))+invdy2*(INPU(i,k+1)+INPU(i,k-1)))  
          
          if (is_conv) OUTU(i,k)=OUTU(i,k) &
               + p5invdx*U_EN_U(i,k)*(INPU(i+1,k)-INPU(i-1,k))+p5invdy*V_EN_U(i,k)*(INPU(i,k+1)-INPU(i,k-1))
          
          if (is_div) OUTU(i,k) = OUTU(i,k) + invdx*(INPP(i,k)-INPP(i-1,k))
       enddo
    enddo
             
    do k=max(k0,1),min(k0+1,nmp+1)
       do i=max(i0,1),min(i0+1,lmp)
          
          if (is_t) then
             OUTV(i,k)=coeft*INPV(i,k) 
          else
             OUTV(i,k)=0.d0
          endif
          
          if (is_diff) OUTV(i,k)=OUTV(i,k) &
               - nu*(-coef*INPV(i,k)+invdx2*(INPV(i+1,k)+INPV(i-1,k))+invdy2*(INPV(i,k+1)+INPV(i,k-1))) 

          if (is_conv) OUTV(i,k)=OUTV(i,k) &
               + p5invdx*U_EN_V(i,k)*(INPV(i+1,k)-INPV(i-1,k))+p5invdy*V_EN_V(i,k)*(INPV(i,k+1)-INPV(i,k-1))

          if (is_div)  OUTV(i,k) = OUTV(i,k) + invdy*(INPP(i,k)-INPP(i,k-1))

       enddo
    enddo

    i=i0; k=k0
    OUTP(i,k) = invdx*(INPU(i+1,k)-INPU(i,k)) + invdy*(INPV(i,k+1)-INPV(i,k))

    call fixe_cl(INPU,INPV,INPP,OUTU,OUTV,OUTP)

    return
  end subroutine mavecxy_df20

  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine inverse(SMU,SMW,SMP,VTU0,VTV0,PRE0,i0,k0)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Resolution du systeme par methode du pivot de Gauss
    !SK
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use blas, only : mavecxy
      implicit none

      real(kind=prec), dimension(0:,0:) :: VTU0,VTV0,PRE0,SMU,SMW,SMP

      real(kind=prec), dimension(0:size(VTU0,1)-1,0:size(VTU0,2)-1) :: VTU,T1X,T1Xb
      real(kind=prec), dimension(0:size(VTV0,1)-1,0:size(VTV0,2)-1) :: VTV,T1Y,T1Yb
      real(kind=prec), dimension(0:size(PRE0,1)-1,0:size(PRE0,2)-1) :: PRE,T1P,T1Pb

      integer :: i0,k0,k,i,i2
      real(kind=prec) :: a, pivot
      real(kind=prec), dimension(5 ,5) :: mat,mat2
      real(kind=prec), dimension(9)   :: b,bd

      VTU=0._prec
      VTV=0._prec
      PRE=0._prec


      mat=0._prec
      
      do i=1,5

         select case(i)
         case (1);   VTU(i0,k0)  =1._prec   
         case (2);   VTU(i0+1,k0)=1._prec   
         case (3);   VTV(i0,k0)  =1._prec   
         case (4);   VTV(i0,k0+1)=1._prec   
         case (5);   PRE(i0,k0)  =1._prec   
         end select


         call mavecxy_df20(VTU,VTV,PRE,T1X,T1Y,T1P,i0,k0) 
         call mavecxy(VTU,VTV,PRE,T1Xb,T1Yb,T1Pb)

         mat(1,i)=T1X(i0,k0)    
         mat(2,i)=T1X(i0+1,k0)  
         mat(3,i)=T1Y(i0,k0)    
         mat(4,i)=T1Y(i0,k0+1)  
         mat(5,i)=T1P(i0,k0)    

         mat2(1,i)=T1Xb(i0,k0)   
         mat2(2,i)=T1Xb(i0+1,k0) 
         mat2(3,i)=T1Yb(i0,k0)   
         mat2(4,i)=T1Yb(i0,k0+1) 
         mat2(5,i)=T1Pb(i0,k0)   

         select case(i)
         case (1);   VTU(i0,k0)=0._prec
         case (2);   VTU(i0+1,k0)=0._prec
         case (3);   VTV(i0,k0)=0._prec
         case (4);   VTV(i0,k0+1)=0._prec
         case (5);   PRE(i0,k0)=0._prec
         end select
         
     enddo
         
     b(1)=SMU(i0,k0)
     b(2)=SMU(i0+1,k0)
     b(3)=SMW(i0,k0)
     b(4)=SMW(i0,k0+1)
     b(5)=SMP(i0,k0)
     bd=b
     
     if (sum(abs(mat-mat2))>1.e-10) then
        print *
        print *,i0,k0,sum(abs(mat-mat2))
        call outdon(mat-mat2 ,'mat-mat2')
        call outdon(mat,'mat')
        call outdon(mat,'mat2')
        VTU=0._prec
        VTV=0._prec
        PRE=0._prec
        
        do i=1,5
           
           select case(i)
           case (1);   VTU(i0,k0)  =1._prec   
           case (2);   VTU(i0+1,k0)=1._prec   
           case (3);   VTV(i0,k0)  =1._prec   
           case (4);   VTV(i0,k0+1)=1._prec   
           case (5);   PRE(i0,k0)  =1._prec   
           end select


           call mavecxy(VTU,VTV,PRE,T1Xb,T1Yb,T1Pb)
           call mavecxy_df20(VTU,VTV,PRE,T1X,T1Y,T1P,i0,k0) 
           print '("a ",11X,E11.3,11X,16X,E11.3,11X)',VTV(i0,k0+1),T1Y(i0,k0+1)
           print '("a ",E11.3,E11.3,E11.3,5X,E11.3,E11.3,E11.3)'&
                & ,VTU(i0,k0),PRE(i0,k0),VTU(i0+1,k0),T1X(i0,k0),T1P(i0,k0),T1X(i0+1,k0)
           print '("a ",11X,E11.3,11X,16X,E11.3,11X)',VTV(i0,k0),T1Y(i0,k0)
           print *
           print '("b ",11X,E11.3,11X,16X,E11.3,11X)',VTV(i0,k0+1),T1Yb(i0,k0+1)
           print '("b ",E11.3,E11.3,E11.3,5X,E11.3,E11.3,E11.3)'&
                & ,VTU(i0,k0),PRE(i0,k0),VTU(i0+1,k0),T1Xb(i0,k0),T1Pb(i0,k0),T1Xb(i0+1,k0)
           print '("b ",11X,E11.3,11X,16X,E11.3,11X)',VTV(i0,k0),T1Yb(i0,k0)
           print *
           print *
         
           select case(i)
           case (1);   VTU(i0,k0)=0._prec
           case (2);   VTU(i0+1,k0)=0._prec
           case (3);   VTV(i0,k0)=0._prec
           case (4);   VTV(i0,k0+1)=0._prec
           case (5);   PRE(i0,k0)=0._prec
           end select
     enddo
     endif

!!$     return

!!$     call parallel_stop

     do i=1,4
        if (abs(mat(i,i))<1.e-5) print '(I3,A,E10.3)',i,' pivot --> ',mat(i,i)
!        print '(I3,A,E10.3)',i,' pivot --> ',mat(i,i)
        do i2=i+1,5
           pivot=mat(i2,i)/mat(i,i)
           do k=i,5
              mat(i2,k)=mat(i2,k)-pivot*mat(i,k)
           enddo
           b(i2)=b(i2)-pivot*b(i)
        enddo
     enddo
!!$
!!$     do i=1,5
!!$        print '("mat ",5E11.3)',(mat(i,k),k=1,5)
!!$     enddo
!!$     print *

     do i=5,1,-1
        a=0._prec
        do k=i+1,5
           a=a+mat(i,k)*b(k)
        enddo
        if (abs(mat(i,i))+abs(b(i)-a)>1.e-18)  b(i)=(b(i)-a)/mat(i,i)
     enddo

!!$     print '("solution ",5E11.3)',(b(i)-bd(i),i=1,5)

     VTU0(i0,k0)     =b(1)
     VTU0(i0+1,k0)   =b(2)
     VTV0(i0,k0)     =b(3)
     VTV0(i0,k0+1)   =b(4)
     PRE0(i0,k0)     =b(5)

     return
  end subroutine inverse                                 ! end_out_light

  !************************************************************************

  subroutine test_restrict
   use disc
   use para
   use data, only : ncheck,ntype_solver,ntype_precond
   implicit none

   if (ncheck/=df2.and.(abs(ntype_precond)/=mg_solver.and.abs(ntype_solver)/=mg_solver))&
        & stop 'Verification  restrict necessite DF2 et MG'

   call restrict(YU,YV,YP,grid(ngrid_max)%U2,grid(ngrid_max)%V2,grid(ngrid_max)%P2,.true.)

   call prolong(grid(ngrid_max)%U2,grid(ngrid_max)%V2,grid(ngrid_max)%P2&
        &      ,grid(ngrid_max)%U1,grid(ngrid_max)%V1,grid(ngrid_max)%P1)

   grid(ngrid_max)%U1=grid(ngrid_max)%U1-YU
   grid(ngrid_max)%V1=grid(ngrid_max)%V1-YV
   grid(ngrid_max)%P1=grid(ngrid_max)%P1-YP


   call outdon(grid(ngrid_max)%U1/dx,"difference%U1")
   call outdon(grid(ngrid_max)%V1/dx,"difference%V1")
   call outdon(grid(ngrid_max)%P1/dx,"difference%P1")


   call parallel_stop

    return
  end subroutine test_restrict



  !************************************************************************

end module mg
