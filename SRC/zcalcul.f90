module calcul

  use constante
  use para
  use io,     only : sortie_champ,surface

  use data,     only : is_print,epsv,epsva, ninterne,ndirection
  use drapeaux, only : is_u,is_v,is_div,is_precond,is_decal,is_cv_file,niv_lim
  use blas,     only : mavecxy, global_ddot, fixe_cl
  use disc,     only : it
  use conv
  use champs, only   : VTUS,VTVS,PRES

  integer :: tag_adv_x=543, tag_adv_y=139

contains


  !SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !SK  Routines calculs
  !SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine solve(SMU,SMV,SMP,VTU,VTV,PRE)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Aiguillage solveur
    !SK
    !SK================================================================
    !SK
    use data, only : ntype_solver,npcg,ncheck,nom_file_cv
    use init, only : set_grid,ncheck_first,ncheck_second
    use mg,   only : lisse
    use drapeaux
    implicit none

    real(kind=prec), dimension(:,:) :: VTU,VTV,PRE,SMU,SMV,SMP
    integer :: ok 

    call set_grid(-1)
    call set_switch(global_switch*is_impl*t_code)

    niv_solve=niv_solve+1

    if (is_cv_file.and.niv_solve==1) open(file=nom_file_cv,unit=78,iostat=ok)
    if (is_cv_file.and.ok/=0.and.niv_solve==1) stop 'pb ouverture fichier cv'  

    select case(abs(ntype_solver)) 
    case (1)                           
       call solve_bicgstab  (SMU,SMV,SMP,VTU,VTV,PRE,npcg)         
    case (2:9)                                                   
       call solve_bicgstab_l(SMU,SMV,SMP,VTU,VTV,PRE,npcg,ntype_solver)
    case (10)                                                     ! start_out_light ! out_small out_tiny
      call solve_gmres   (SMU,SMV,SMP,VTU,VTV,PRE,npcg)           ! out_small out_tiny
    case (11)                                                      
       call solve_gmresr    (SMU,SMV,SMP,VTU,VTV,PRE,npcg,.true.)
    case (12) 
       call solve_gmresr    (SMU,SMV,SMP,VTU,VTV,PRE,npcg,.false.) 
    case (13,-13)                                                  ! start_tiny 
       call solve_mg        (SMU,SMV,SMP,VTU,VTV,PRE,npcg,ntype_solver<0)
    case (14) 
       call solve_orthodir  (SMU,SMV,SMP,VTU,VTV,PRE,npcg)        ! start_small   
    case (15) 
       call solve_orthomin  (SMU,SMV,SMP,VTU,VTV,PRE,npcg)         ! end_small 
    case (16) 
       call solve_lisse     (SMU,SMV,SMP,VTU,VTV,PRE,npcg)        ! end_out_light   ! end_tiny
    case (17,18) 
       call solve_shwarz(SMU,SMV,SMP,VTU,VTV,PRE,npcg)   
    case default  
       stop 'solver not available in this configuration'
    end select
 
    if (is_cv_file.and.niv_solve==1)    close(78)
    
    niv_solve=niv_solve-1

    call set_grid(-1)

  end subroutine solve


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine precond(SMU,SMV,SMP,VTU,VTV,PRE)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Aiguillage préconditionneur
    !SK
    !SK================================================================
    !SK
    use champs, only : VTUp,VTVp,PREp,SMUp,SMVp,SMPp
    use data, only   :  ntype_solver,ntype_precond,nprecond,ncheck ,ncheck_precond&
         & ,nrecvddm,nrecvddmx, nrecvddmy
    use disc
    use champs, only : MASK_DIRU,MASK_DIRV
    use drapeaux
    use para, only   : is_mpp,rfr_stn,rfr_ddm_stn,stn_news
    use init,   only : set_grid
    use mg,     only : lisse
    implicit none

    real(kind=prec), dimension(0:,0:)          :: VTU,VTV,PRE,SMU,SMV,SMP
!!$    real(kind=prec), dimension(0:size(SMP,2)-1):: p_raccord_w, p_raccord_e
!!$    real(kind=prec), dimension(0:size(SMP,1)-1):: p_raccord_n, p_raccord_s

    real(kind=prec), dimension(:,:), pointer :: pSMUp,pSMVp,pSMPp&
         & ,pVTUp,pVTVp,pPREp
    integer :: old_ncheck
    logical :: old_is_precond,old_is_mpp
    integer :: if,id,kf,kd,lmp,nmp,nb_t,nb_t_fin

    call set_grid(ngrid_max)

    old_ncheck=ncheck;                ncheck=ncheck_precond
    old_is_precond=is_precond;        is_precond=.false.
    old_is_mpp=is_mpp;                is_mpp=.false.

    niv_solve=niv_solve+1

    nb_t_fin=0
    if (is_schwarz_mult) nb_t_fin=nb_tasks-1

    if (is_u)   SMUp(iudeb:iufin,kudeb:kufin)=SMU; 
    if (is_v)   SMVp(ivdeb:ivfin,kvdeb:kvfin)=SMV; 
    if (is_div) SMPp(ipdeb:ipfin,kudeb:kpfin)=SMP

    if (is_u) call rfr_ddm_stn(SMUp,iudeb,iufin,kudeb,kufin,nrecvddmx,nrecvddmy)
    if (is_v) call rfr_ddm_stn(SMVp,ivdeb,ivfin,kvdeb,kvfin,nrecvddmx,nrecvddmy)
    if (is_div) call rfr_ddm_stn(SMPp,ipdeb,ipfin,kpdeb,kpfin,nrecvddmx,nrecvddmy)


    pVTUp => VTUp(iumgdeb:iumgfin,kumgdeb:kumgfin)
    pVTVp => VTVp(ivmgdeb:ivmgfin,kvmgdeb:kvmgfin)
    pPREp => PREp(ipmgdeb:ipmgfin,kpmgdeb:kpmgfin)
    pSMUp => SMUp(iumgdeb:iumgfin,kumgdeb:kumgfin)
    pSMVp => SMVp(ivmgdeb:ivmgfin,kvmgdeb:kvmgfin)
    pSMPp => SMPp(ipmgdeb:ipmgfin,kpmgdeb:kpmgfin)
 
    if (is_schwarz) then
       if (is_u)   VTUp(iudeb:iufin,kudeb:kufin)=VTU; 
       if (is_v)   VTVp(ivdeb:ivfin,kvdeb:kvfin)=VTV; 
       if (is_div) PREp(ipdeb:ipfin,kudeb:kpfin)=PRE
    else
       if (is_u)   VTUp=0._prec
       if (is_v)   VTVp=0._prec

       if (is_u) where (MASK_DIRU) pSMUp=0._prec
       if (is_v) where (MASK_DIRV) pSMVp=0._prec
    endif

    
    do nb_t=0,0 !nb_t_fin
       
       if (is_schwarz) then
          if (is_u) call rfr_ddm_stn(VTUp,iudeb,iufin,kudeb,kufin,nrecvddmx,nrecvddmy)
          if (is_v) call rfr_ddm_stn(VTVp,ivdeb,ivfin,kvdeb,kvfin,nrecvddmx,nrecvddmy)
          if (is_div) call rfr_ddm_stn(PREp,ipdeb,ipfin,kpdeb,kpfin,nrecvddmx,nrecvddmy)
          call fixe_cl(pVTUp,pVTVp,pPREp,pSMUp,pSMVp,pSMPp)
       endif


       if (.not.is_schwarz_mult.or.my_task==nb_t) then
          select case(abs(ntype_precond))
          case (1) 
             call solve_bicgstab(SMUp,SMVp,SMPp,pVTUp,pVTVp,pPREp,nprecond)
          case (2:9)                                                          
             call solve_bicgstab_l&
                  & (SMUp,SMVp,SMPp,pVTUp,pVTVp,pPREp,nprecond,ntype_precond)
          case (10)                                                           ! start_out_light out_small ! start_tiny
             call solve_gmres(pSMUp,pSMVp,pSMPp,pVTUp,pVTVp,pPREp,nprecond)       ! out_small
          case (13,-13)                                                           
             call solve_mg  (pSMUp,pSMVp,pSMPp,pVTUp,pVTVp,pPREp,nprecond,ntype_precond<0)   
          case (16)                                                           
             call lisse     (pSMUp,pSMVp,pSMPp,pVTUp,pVTVp,pPREp,nprecond)    ! end_out_light  ! end_tiny
          case default  
             stop 'preconditioner not available in this configuration'
          end select
       endif

       if (is_div) then
          if=0; id=0; kf=0; kd=0
          lmp=size(PRE,1)-2;nmp=size(PRE,2)-2;
          if (.not.is_north) kf=nrecvddmy
          if (.not.is_south) kd=nrecvddmy
          if (.not.is_west)  id=nrecvddmx
          if (.not.is_east)  if=nrecvddmx
          PRE(id:lmp+1-if,kd:nmp+1-kf)=PREp(ipdeb+id:ipfin-if,kpdeb+kd:kpfin-kf)
          if (.not.is_schwarz) PRE=PRE-global_minval(PRE)
       endif
!!$    if (is_div) PRE=PREp(ipdeb:ipfin,kpdeb:kpfin)

!!$    p_raccord_e=PREp(ipdeb,kpdeb:kpfin)
!!$    p_raccord_w=PREp(ipfin,kpdeb:kpfin)
!!$    p_raccord_s=PREp(ipdeb:ipfin,kpdeb)
!!$    p_raccord_n=PREp(ipdeb:ipfin,kpfin)
!!$    
!!$    if (is_div) call rfr_ddm_stn(PREp,ipdeb,ipfin,kpdeb,kpfin,nrecvddmx,nrecvddmy)
!!$    
!!$    if (.not.is_east)  PRE=PRE+PREp(ipdeb,kpdeb)-p_raccord_e(0)
!!$    if (.not.is_south) PRE=PRE+PREp(ipdeb,kpdeb)-p_raccord_s(0)

    enddo

    niv_solve=niv_solve-1
    
    if (is_u)   VTU=VTUp(iudeb:iufin,kudeb:kufin); 
    if (is_v)   VTV=VTVp(ivdeb:ivfin,kvdeb:kvfin); 

    is_mpp=old_is_mpp
    ncheck=old_ncheck
    is_precond=old_is_precond
    call set_grid(-1)

  end subroutine precond

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine conv_collect(err0,err,n)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use conv
    implicit none
    integer :: n
    real(kind=prec) :: err0,err,rho

    if (niv_solve<6.and.n>0) then
       flag_niv(niv_solve)=.true.
       rho=(err/err0)**(1._prec/real(n))
       nb_iter(niv_solve)=nb_iter(niv_solve)+1
       total_iter(niv_solve)=total_iter(niv_solve)+n+1
       max_iter(niv_solve)=max(max_iter(niv_solve),n+1)
       min_iter(niv_solve)=min(min_iter(niv_solve),n+1)
       last_iter(niv_solve)=n+1
       total_rho(niv_solve)=total_rho(niv_solve)+rho
       max_rho(niv_solve)=max(max_rho(niv_solve),rho)
       min_rho(niv_solve)=min(min_rho(niv_solve),rho)
       last_rho(niv_solve)=rho
       first_res(niv_solve)=err0
       last_res(niv_solve)=err
    endif

    return
  end subroutine conv_collect


  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive &
       subroutine solve_bicgstab(SMU,SMV,SMP,VTU,VTV,PRE,ncgmax)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Gradient conjugue BiCGSTAB
    !SK
    !SK================================================================
    !SK
    implicit none

    integer ncgmax
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP
    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1) :: &
         vu,ru,pu,su,tu,r0u,yu,zu
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1) :: &
         vv,rv,pv,sv,tv,r0v,yv,zv
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1) :: &
         vp,rp,pp,sp,tp,r0p,yp,zp

    real(kind=prec)  :: rho,rhom1,alpha,beta,omega,err0,err,errdiv=0._prec
    integer           :: ncg

    if (is_u) zu=0._prec;     if (is_v) zv=0._prec;   if (is_div) zp=0._prec;    
    ! car utilise ds precond et non initalise avant!!!!!!!  
    ! yu est initialise a qqch car on appelle mavecxy

    !!call outdon(VTU,'VTu')
    call mavecxy(VTU,VTV,PRE,yu,yv,yp)
    !!call outdon(yu,'yu')
    !!call parallel_stop

    if (is_u)  ru=SMU-yu;    if (is_v)  rv=SMV-yv;    if (is_div) rp=SMP-yp
    if (is_u) r0u=ru;        if (is_v) r0v=rv;        if (is_div) r0p=rp
    if (is_u)  pu=ru;        if (is_v)  pv=rv;        if (is_div)  pp=rp

    if (is_precond) then
       call precond(pu,pv,pp,yu,yv,yp)
       call mavecxy(yu,yv,yp,vu,vv,vp)
    else
       call mavecxy(pu,pv,pp,vu,vv,vp)
    endif

    rho=global_ddot(ru,rv,rp,ru,rv,rp)
    err0=sqrt(rho)
    err=err0
    ncg=0

    do while (err>epsv*err0.and.err>epsva.and.ncg.lt.abs(ncgmax))

!!$       call outdon(vtu,'vtu')
!!$       call outdon(ru,'ru')

       ncg=ncg+1

       alpha=global_ddot(r0u,r0v,r0p,vu,vv,vp)
!!$       call outdon(pu,'pu')
!!$       call outdon(vu,'vu')
!!$       call outdon(r0u,'r0u')
!!$       if (my_task.eq.0) then
!!$          print *,'alpha=vu*r0u ',alpha
!!$       end if
       if (alpha==0.e0)  exit
       alpha=rho/alpha

       if (is_u)   su=ru-alpha*vu;    
       if (is_v)   sv=rv-alpha*vv;  
       if (is_div) sp=rp-alpha*vp

       if (is_precond) then
          call precond(su,sv,sp,zu,zv,zp)
          call mavecxy(zu,zv,zp,tu,tv,tp)
       else
          call mavecxy(su,sv,sp,tu,tv,tp)
       endif
       omega=global_ddot(tu,tv,tp,tu,tv,tp)
       if (omega==0._prec) exit
       omega=global_ddot(su,sv,sp,tu,tv,tp)/omega

       if (is_precond) then
          if (is_u) VTU=VTU+alpha*yu+omega*zu; 
          if (is_v) VTV=VTV+alpha*yv+omega*zv; 
          if (is_div) PRE=PRE+alpha*yp+omega*zp
       else
          if (is_u) VTU=VTU+alpha*pu+omega*su; 
          if (is_v) VTV=VTV+alpha*pv+omega*sv; 
          if (is_div) PRE=PRE+alpha*pp+omega*sp
       endif

       if (is_u) ru=su-omega*tu;            
       if (is_v) rv=sv-omega*tv;  
       if (is_div) rp=sp-omega*tp            

       err=sqrt(global_ddot(ru,rv,rp,ru,rv,rp))
       if (is_print>=3*niv_solve-1.and.is_div) &
            & errdiv=sqrt(global_sum(rp*rp))
       if (ncg==1) err0=err
       if (is_print>=3*niv_solve-1.and.my_task==0)&
            & print '(A,1X,I5,A,1X,I5,2E16.7)','niveau ',niv_solve&
            & ,' erreur BiCGSTAB : ',ncg,err,errdiv
       if (is_cv_file.and.niv_solve==1)&
            & write (78,'(I5,1X,2E16.7)') ncg,err,errdiv
       if (is_div.and.niv_solve==1) call sortie_champ(PRE,'p')
       if (is_div.and.niv_solve==1) call sortie_champ(PRE-PRES,'X',surface)
       if (niv_solve==1) call sortie_champ(VTU-VTUS,'u',surface)
       if (niv_solve==1) call sortie_champ(VTV-VTVS,'v',surface)
       if (niv_solve==1) call sortie_champ(ru,'U',surface)
       if (niv_solve==1) call sortie_champ(rv,'V',surface)
       if (is_div.and.niv_solve==1) call sortie_champ(rp,'P',surface)
       if (err<epsv*err0.or.err<epsva) exit

       rhom1=rho
       rho=global_ddot(ru,rv,rp,r0u,r0v,r0p)
       beta=rho*alpha/rhom1/omega

       if (is_u) pu=ru+beta*(pu-omega*vu);    
       if (is_v) pv=rv+beta*(pv-omega*vv);      
       if (is_div) pp=rp+beta*(pp-omega*vp)

       if (is_precond) then
          call precond(pu,pv,pp,yu,yv,yp)
          call mavecxy(yu,yv,yp,vu,vv,vp)
       else
          call mavecxy(pu,pv,pp,vu,vv,vp)
       endif

    enddo

    call conv_collect(err0,err,ncg-1)

!    if (niv_solve==1) call outdon(SMU,'SMU',1,1)

    if ((niv_solve.eq.1.or.is_print>=3*niv_solve-2).and.my_task==0.and.ncg.eq.ncgmax) then
       print '(A,1X,I5,A,1X,2E16.7)'&
            & ,'Sorry! BiCGSTAB n''a pas converge : ncg =',ncg&
            & ,' err=',err,errdiv
    endif
    if (is_print>=3*niv_solve-2.and.my_task==0.and.ncg.lt.ncgmax) then
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Bingo! BiCGSTAB a converge : ncg =',ncg&
               & ,' err=',err,errdiv
    endif
    
    return

  end subroutine solve_bicgstab


  !************************************************************************    

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss   
  recursive &
       subroutine solve_bicgstab_l(SMU,SMV,SMP,VTU,VTV,PRE,ncgmax,l)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Gradient conjugue BiCGSTAB
    !SK
    !SK================================================================
    !SK
    use data, only : nbmg
    implicit none

    integer ncgmax,l
    real(kind=prec), dimension(:,:) :: VTU,VTV,PRE,SMU,SMV,SMP

    real(kind=prec), dimension(size(VTU,1),size(VTU,2))              :: r0u,vu
    real(kind=prec), dimension(size(VTV,1),size(VTV,2))              :: r0v,vv
    real(kind=prec), dimension(size(PRE,1),size(PRE,2))              :: r0p,vp     
    real(kind=prec), dimension(size(VTU,1),size(VTU,2),0:l) :: ru,uu,xu
    real(kind=prec), dimension(size(VTV,1),size(VTV,2),0:l) :: rv,uv,xv
    real(kind=prec), dimension(size(PRE,1),size(PRE,2),0:l) :: rp,up,xp

    real(kind=prec), dimension(0:l) :: gamma,gamma_1,gamma_2,sigma
    real(kind=prec), dimension(0:l,0:l) :: tau

    real(kind=prec) :: rho_0,rho_1,alpha,beta,omega,err0,err,errdiv=0._prec
    integer :: ncg,i,j

    call mavecxy(VTU,VTV,PRE,r0u,r0v,r0p)

    if (is_u)   r0u=SMU-r0u;
    if (is_v)   r0v=SMV-r0v;
    if (is_div) r0p=SMP-r0p;

    if (is_precond) then
       call precond(r0u,r0v,r0p,vu,vv,vp)
       if (is_u)   r0u=vu
       if (is_v)   r0v=vv
       if (is_div) r0p=vp
    endif

    if (is_u)   then; ru=0._prec; ru(:,:,0)=r0u;  endif
    if (is_v)   then; rv=0._prec; rv(:,:,0)=r0v;  endif
    if (is_div) then; rp=0._prec; rp(:,:,0)=r0p;  endif

    err=sqrt(global_ddot(r0u,r0v,r0p,r0u,r0v,r0p))
    err0=err
    ncg=0
    
    rho_0=1._prec; alpha=0._prec;    omega=1._prec
    if (is_u)   then; uu=0._prec; xu=0._prec; xu(:,:,0)=VTU; endif
    if (is_v)   then; uv=0._prec; xv=0._prec; xv(:,:,0)=VTV; endif
    if (is_div) then; up=0._prec; xp=0._prec; xp(:,:,0)=PRE; endif
    r0u=r0u+1._prec; 

    do while (err>epsv*err0.and.err>epsva.and.ncg.lt.abs(ncgmax))

       ncg=ncg+1

       rho_0=-omega*rho_0

       do j=0,l-1
          rho_1=global_ddot(ru(:,:,j),rv(:,:,j),rp(:,:,j),r0u,r0v,r0p)
          beta=alpha*rho_1/rho_0
          rho_0=rho_1
          if (is_u)   uu(:,:,0:j)=ru(:,:,0:j)-beta*uu(:,:,0:j)
          if (is_v)   uv(:,:,0:j)=rv(:,:,0:j)-beta*uv(:,:,0:j)
          if (is_div) up(:,:,0:j)=rp(:,:,0:j)-beta*up(:,:,0:j)
          if (is_precond) then
             call mavecxy(uu(:,:,j)  ,uv(:,:,j)  ,up(:,:,j),vu,vv,vp)
             call precond(vu,vv,vp,uu(:,:,j+1),uv(:,:,j+1),up(:,:,j+1))
          else
             call mavecxy(uu(:,:,j)  ,uv(:,:,j)  ,up(:,:,j) &
                         ,uu(:,:,j+1),uv(:,:,j+1),up(:,:,j+1))
          endif
          alpha=rho_0/global_ddot(uu(:,:,j+1),uv(:,:,j+1),up(:,:,j+1),r0u,r0v,r0p)
          if (is_u)   ru(:,:,0:j)=ru(:,:,0:j)-alpha*uu(:,:,1:j+1)
          if (is_v)   rv(:,:,0:j)=rv(:,:,0:j)-alpha*uv(:,:,1:j+1)
          if (is_div) rp(:,:,0:j)=rp(:,:,0:j)-alpha*up(:,:,1:j+1)
          if (is_precond) then
             call mavecxy(ru(:,:,j)  ,rv(:,:,j)  ,rp(:,:,j),  vu,vv,vp)
             call precond(vu,vv,vp,ru(:,:,j+1),rv(:,:,j+1),rp(:,:,j+1))
          else
             call mavecxy(ru(:,:,j)  ,rv(:,:,j)  ,rp(:,:,j) &
                         ,ru(:,:,j+1),rv(:,:,j+1),rp(:,:,j+1))
          endif
          if (is_u)   xu(:,:,0)=xu(:,:,0)+alpha*uu(:,:,0)
          if (is_v)   xv(:,:,0)=xv(:,:,0)+alpha*uv(:,:,0)
          if (is_div) xp(:,:,0)=xp(:,:,0)+alpha*up(:,:,0)
       enddo

       do j=1,l
          do i=1,j-1
             tau(i,j)=global_ddot(ru(:,:,j)  ,rv(:,:,j)  ,rp(:,:,j) &
                                 ,ru(:,:,i)  ,rv(:,:,i)  ,rp(:,:,i) )/sigma(i)
             if (is_u)   ru(:,:,j)=ru(:,:,j)-tau(i,j)*ru(:,:,i)
             if (is_v)   rv(:,:,j)=rv(:,:,j)-tau(i,j)*rv(:,:,i)
             if (is_div) rp(:,:,j)=rp(:,:,j)-tau(i,j)*rp(:,:,i)             
          enddo
          sigma(j)=  global_ddot(ru(:,:,j)  ,rv(:,:,j)  ,rp(:,:,j) &
                                ,ru(:,:,j)  ,rv(:,:,j)  ,rp(:,:,j) )
          gamma_1(j)=global_ddot(ru(:,:,0)  ,rv(:,:,0)  ,rp(:,:,0) &
                                ,ru(:,:,j)  ,rv(:,:,j)  ,rp(:,:,j) )/sigma(j)
       enddo
       gamma(l)=gamma_1(l)
       omega=gamma(l)

       do j=l-1,1,-1
          gamma(j)  =gamma_1(j)-dot_product(tau(j,j+1:l),gamma(j+1:l))
       enddo
       do j=1,l-2
          gamma_2(j)=gamma(j+1)+dot_product(tau(j,j+1:l-1),gamma(j+2:l))
       enddo
       gamma_2(l-1)=gamma(l)

       if (is_u)   xu(:,:,0)=xu(:,:,0)+gamma(1)  *ru(:,:,0)
       if (is_v)   xv(:,:,0)=xv(:,:,0)+gamma(1)  *rv(:,:,0)
       if (is_div) xp(:,:,0)=xp(:,:,0)+gamma(1)  *rp(:,:,0)             
                                     
       if (is_u)   ru(:,:,0)=ru(:,:,0)-gamma_1(l)*ru(:,:,l)
       if (is_v)   rv(:,:,0)=rv(:,:,0)-gamma_1(l)*rv(:,:,l)
       if (is_div) rp(:,:,0)=rp(:,:,0)-gamma_1(l)*rp(:,:,l)             

       if (is_u)   uu(:,:,0)=uu(:,:,0)-gamma(l)  *uu(:,:,l)
       if (is_v)   uv(:,:,0)=uv(:,:,0)-gamma(l)  *uv(:,:,l)
       if (is_div) up(:,:,0)=up(:,:,0)-gamma(l)  *up(:,:,l)             

       do j=1,l-1
          if (is_u)   xu(:,:,0)=xu(:,:,0)+gamma_2(j)*ru(:,:,j)
          if (is_v)   xv(:,:,0)=xv(:,:,0)+gamma_2(j)*rv(:,:,j)
          if (is_div) xp(:,:,0)=xp(:,:,0)+gamma_2(j)*rp(:,:,j)             
                                        
          if (is_u)   ru(:,:,0)=ru(:,:,0)-gamma_1(j)*ru(:,:,j)
          if (is_v)   rv(:,:,0)=rv(:,:,0)-gamma_1(j)*rv(:,:,j)
          if (is_div) rp(:,:,0)=rp(:,:,0)-gamma_1(j)*rp(:,:,j)             
   
          if (is_u)   uu(:,:,0)=uu(:,:,0)-gamma(j)  *uu(:,:,j)
          if (is_v)   uv(:,:,0)=uv(:,:,0)-gamma(j)  *uv(:,:,j)
          if (is_div) up(:,:,0)=up(:,:,0)-gamma(j)  *up(:,:,j)             
       enddo

       err=sqrt(global_ddot(ru(:,:,0)  ,rv(:,:,0)  ,rp(:,:,0) &
                           ,ru(:,:,0)  ,rv(:,:,0)  ,rp(:,:,0) ))
       if (is_print>=3*niv_solve-1.and.is_div)&
            & errdiv=sqrt(global_sum(rp(:,:,0)*rp(:,:,0)))

       if (is_print>=3*niv_solve-1.and.my_task==0)&
            & print '(A,1X,I5,A,1X,I5,2E16.7)'&
            & ,'niveau ',niv_solve,' erreur BiCGSTAB(l) : '&
            & ,ncg,err,errdiv
       if (my_task==0.and.is_cv_file.and.niv_solve==1)&
            & write (78,'(I5,1X,2E16.7)')&
            & ncg,err,errdiv
       if (is_div.and.niv_solve==1) call sortie_champ(xp(:,:,0),'p')
       if (is_div.and.niv_solve==1) call sortie_champ(xp(:,:,0)-PRES,'X',surface)
       if (niv_solve==1) call sortie_champ(xu(:,:,0)-VTUS,'u',surface)
       if (niv_solve==1) call sortie_champ(xv(:,:,0)-VTVS,'v',surface)
       if (niv_solve==1) call sortie_champ(ru(:,:,0),'U',surface)
       if (niv_solve==1) call sortie_champ(rv(:,:,0),'V',surface)
       if (is_div.and.niv_solve==1) call sortie_champ(rp(:,:,0),'P',surface)
       if (err<epsv*err0.or.err<epsva) exit

    enddo

    call conv_collect(err0,err,ncg-1)

    if (is_u)   VTU=xu(:,:,0); 
    if (is_v)   VTV=xv(:,:,0); 
    if (is_div) PRE=xp(:,:,0); 
    
    if ((niv_solve.eq.1.or.is_print>=3*niv_solve-2).and.my_task==0.and.ncg.eq.ncgmax) then
       print '(A,1X,I5,A,1X,2E16.7)'&
            & ,'Sorry! BiCGSTAB(l) n''a pas converge : ncg =',ncg&
            & ,' err=',err,errdiv
    endif
    if (is_print>=3*niv_solve-2.and.my_task==0.and.ncg.lt.ncgmax) then
       print '(A,1X,I5,A,1X,2E16.7)'&
            & ,'Bingo! BiCGSTAB(l) a converge : ncg =',ncg&
            & ,' err=',err,errdiv
    endif
    
    return

  end subroutine solve_bicgstab_l

  !************************************************************************    ! start_out_light  start_small ! start_tiny


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive &
       subroutine solve_gmres(SMU,SMV,SMP,VTU,VTV,PRE,ncgmax)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Methode GMRES
    !SK
    !SK================================================================
    !SK
    implicit none

    integer ncgmax
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP

    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1) :: ru,wu,zu
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1) :: rv,wv,zv
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1) :: rp,wp,zp
    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1,1:ninterne+1) :: vu 
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1,1:ninterne+1) :: vv
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1,1:ninterne+1) :: vp

    real(kind=prec), dimension(1:ninterne+1,1:ninterne) :: h
    real(kind=prec), dimension(1:ninterne)                :: cs,sn
    real(kind=prec), dimension(1:ninterne+1)              :: s,y

    real(kind=prec)   :: alpha,err,err0,bnrm2,temp,errdiv=0._prec
    integer, save     :: ncg,i,k

    
    bnrm2 = sqrt(global_ddot(SMU,SMV,SMP,SMU,SMV,SMP))
    if (bnrm2<(0.1_prec)**14) bnrm2=1._prec

    if (is_u) vu=0._prec; if (is_v) vv=0._prec; if (is_div) vp=0._prec; 

    err=1._prec; err0=err
    ncg=0

    do while (err>epsv*err0.and.err>epsva.and.ncg<=ncgmax)

       ncg=ncg+1

       call mavecxy(VTU,VTV,PRE,wu,wv,wp)
       if (is_precond) then
          if (is_u) wu=SMU-wu; if (is_v) wv=SMV-wv; if (is_div) wp=SMP-wp
          call precond(wu,wv,wp,ru,rv,rp)
       else
          if (is_u) ru=SMU-wu; if (is_v) rv=SMV-wv; if (is_div) rp=SMP-wp
       endif
       
       err=sqrt(global_ddot(ru,rv,rp,ru,rv,rp))
       if (is_print>=3*niv_solve-1.and.is_div)&
            & errdiv=sqrt(global_sum(rp*rp))
       if (ncg<=2) err0=err
       if (is_print>=3*niv_solve-1.and.my_task==0)&
            & print '(A,1X,I5,A,1X,I5,2E16.7)'&
            & ,'niveau ',niv_solve,' erreur GMRES : ',ncg,err,errdiv
       if (is_cv_file.and.ncg>=2.and.niv_solve==1)&
            & write (78,'(I5,1X,2E16.7)') err,errdiv
       if (err<epsv*err0.or.err<epsva) exit

       alpha=1._prec/err
       if (is_u) vu(:,:,1)=alpha*ru; if (is_v) vv(:,:,1)=alpha*rv; if (is_div) vp(:,:,1)=alpha*rp
       cs=0._prec; sn=0._prec; h=0._prec

       s=0._prec;  s(1)=err        !SK   s = norm( r )*e1;
       
       !SK   construct orthonormal

       do i=1,ninterne

          if (is_precond) then
             call mavecxy(vu(:,:,i),vv(:,:,i),vp(:,:,i),zu,zv,zp)
             call precond(zu,zv,zp,wu,wv,wp)
          else
             call mavecxy(vu(:,:,i),vv(:,:,i),vp(:,:,i),wu,wv,wp)
          endif

          do k=1,i
             H(k,i)= global_ddot(wu,wv,wp,vu(:,:,k),vv(:,:,k),vp(:,:,k))  
             alpha=H(k,i)
             if (is_u)   wu=wu-alpha*vu(:,:,k); 
             if (is_v)   wv=wv-alpha*vv(:,:,k); 
             if (is_div) wp=wp-alpha*vp(:,:,k);  
          enddo

          H(i+1,i)= sqrt(global_ddot(wu,wv,wp,wu,wv,wp))
          alpha=1._prec/H(i+1,i)

          if (is_u) vu(:,:,i+1)=alpha*wu; if (is_v) vv(:,:,i+1)=alpha*wv;  if (is_div) vp(:,:,i+1)=alpha*wp;  

          do k = 1,i-1        ! apply Glvens rotatlon
             temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i)
             H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i)
             H(k,i)   = temp
          enddo

          call rotmat( H(i,i),H(i+1,i),cs(i),sn(i)) ! form l-th rotatlon matrlx
          temp   = cs(i)*s(i); s(i+1) = -sn(i)*s(i); s(i)   = temp
          H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i)
          H(i+1,i) = 0.0

      enddo

      y(ninterne)=s(ninterne)/H(ninterne,ninterne)  ! call solve_hessenberg(H,s,y,ninterne)
      do i=ninterne-1,1,-1
         temp=0._prec
         do k=i+1,ninterne
            temp=temp+H(i,k)*y(k)
         enddo
         y(i)=(s(i)-temp)/H(i,i)
      enddo

            
      do i=1,ninterne                               !SK   x = x + V(:,1:l)*y
         if (is_u)   VTU=VTU+y(i)*vu(:,:,i);
         if (is_v)   VTV=VTV+y(i)*vv(:,:,i);
         if (is_div) PRE=PRE+y(i)*vp(:,:,i);
      enddo
      
    enddo

    call conv_collect(err0,err,ncg-1)

    if (is_print>=3*niv_solve-2.and.my_task==0) then
       if (ncg.eq.ncgmax) then
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Sorry! GMRES n''a pas converge : ncg =',ncg&
               & ,' err=',err,errdiv
       else
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Bingo! GMRES a converge : ncg =',ncg&
               & ,' err=',err,errdiv
       endif
    endif

    return
  end subroutine solve_gmres

  subroutine rotmat (a,b,c,s)
    !SK
    !SK  Compute the Givens rotation matrix parameters for a and b.
    !SK
    implicit none
    real(kind=prec) :: a,b,c,s,temp
    
    if ( b.le.(0.1_prec)**18 ) then
       c = 1._prec;   s = 0._prec
    else if ( abs(b)> abs(a) ) then 
       temp = a / b
       s = 1._prec / sqrt( 1._prec + temp*temp )
       c = temp * s
    else
       temp = b / a
       c = 1._prec / sqrt( 1._prec + temp*temp  )
       s = temp * c
    endif

    return
  end subroutine rotmat
 
  !************************************************************************          !end_small end_tiny

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive &
       subroutine solve_gmresr(SMU,SMV,SMP,VTU,VTV,PRE,ncgmax,reconj)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Methode GMRESR
    !SK
    !SK================================================================
    !SK
    implicit none

    integer ncgmax
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP
    logical :: reconj
    integer, parameter :: stockmax=5
    real(kind=prec), parameter :: stock_critere=(0.1_prec)**9

    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1) :: vu,ru
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1) :: vv,rv
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1) :: vp,rp
    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1,0:ndirection-1) :: uu,cu
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1,0:ndirection-1) :: uv,cv
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1,0:ndirection-1) :: up,cp

    real(kind=prec), dimension(:,:,:),allocatable, save :: usu,csu
    real(kind=prec), dimension(:,:,:),allocatable, save :: usv,csv
    real(kind=prec), dimension(:,:,:),allocatable, save :: usp,csp
    real(kind=prec), dimension(:)    ,allocatable, save :: coef_stock

    real(kind=prec)  :: min_coef_stock,rho,rho_1,alpha&
         & ,err,err0,errdiv=0._prec
    integer, save     :: old_where,where_min_coef_stock,gmres_deb,ndir_stock 
    integer, save     :: ndir,n,ncg,ok

!SK
!SK=============================================================
!SK   INITIALISATION
!SK=============================================================
!SK


!SK-------------------------------------------------------------
!SK  Mise a yero des directions stockees   
!SK-------------------------------------------------------------

      if (it.eq.1) then
         print *,'it=1'
         if (is_u)   then; uu=0._prec; cu=0._prec; endif
         if (is_v)   then; uv=0._prec; cv=0._prec; endif
         if (is_div) then; up=0._prec; cp=0._prec; endif
         if (reconj) then
            allocate (coef_stock(0:stockmax),stat=ok); 
            if (ok/=0) stop 'coef_stock : error alloc'
            if (is_u) then
               allocate (usu(0:size(SMU,1)-1,0:size(SMU,2)-1,0:stockmax),stat=ok); 
               if (ok/=0) stop 'usu : error alloc'
               allocate (csu(0:size(SMU,1)-1,0:size(SMU,2)-1,0:stockmax),stat=ok); 
               if (ok/=0) stop 'csu : error alloc'
               usu=0._prec; csu=0._prec
            endif
            if (is_v) then
               allocate (usv(0:size(SMV,1)-1,0:size(SMV,2)-1,0:stockmax),stat=ok); 
               if (ok/=0) stop 'usv : error alloc'
               allocate (csv(0:size(SMV,1)-1,0:size(SMV,2)-1,0:stockmax),stat=ok); 
               if (ok/=0) stop 'csv : error alloc'
               usv=0._prec; csv=0._prec
            endif
            if (is_div) then
               allocate (usp(0:size(SMP,1)-1,0:size(SMP,2)-1,0:stockmax),stat=ok); if (ok/=0) stop 'usp : error alloc'
               allocate (csp(0:size(SMP,1)-1,0:size(SMP,2)-1,0:stockmax),stat=ok); if (ok/=0) stop 'csp : error alloc'
               usp=0._prec; csp=0._prec
            endif
            coef_stock=0._prec
            ndir_stock=0
            min_coef_stock=stock_critere
            where_min_coef_stock=0
         endif
      endif



      !SK-------------------------------------------------------------
      !SK   Calcul du reste initial et mise a yero des vecteurs de 
      !SK   travail
      !SK-------------------------------------------------------------

      call mavecxy(VTU,VTV,PRE,vu,vv,vp)

      if (is_u)   ru=SMU-vu
      if (is_v)   rv=SMV-vv
      if (is_div) rp=SMP-vp

      rho=global_ddot(ru,rv,rp,ru,rv,rp)
      err0=sqrt(rho)
      err=err0
      
      !SK-------------------------------------------------------------
      !SK   Calcul du point de depart optimal par reconjugaison avec
      !SK   les directions Stockees
      !SK-------------------------------------------------------------

      !SK      print *,it,ndir_stock

      if (reconj) then
         do n=0,ndir_stock-1
            alpha=global_ddot(csu(:,:,n),csv(:,:,n),csp(:,:,n),ru,rv,rp)
            coef_stock(n)=abs(alpha)
            if (is_u)   then;  VTU=VTU+alpha*usu(:,:,n);    ru=ru-alpha*csu(:,:,n);   endif
            if (is_v)   then;  VTV=VTV+alpha*usv(:,:,n);    rv=rv-alpha*csv(:,:,n);   endif
            if (is_div) then;  PRE=PRE+alpha*usp(:,:,n);    rp=rp-alpha*csp(:,:,n);   endif
            rho_1=global_ddot(ru,rv,rp,ru,rv,rp)
           ! print *,n,alpha,rho_1
         enddo
         if (ndir_stock.gt.0) then
            min_coef_stock=10000._prec
            do n=0,ndir_stock-1
               if (coef_stock(n).lt.min_coef_stock) then
                  min_coef_stock=coef_stock(n)
                  where_min_coef_stock=n
               endif
            enddo
            min_coef_stock=max(min_coef_stock,stock_critere)
         endif
      endif

      ndir=-1
      gmres_deb=0
      alpha=stock_critere/100._prec
      
      !     
      !===============================================================
      !     ITERATIONS GMRESR
      !===============================================================
      !     
      ncg=0
      do while (err>epsv*err0.and.err>epsva.and.ncg.lt.abs(ncgmax))

         ncg=ncg+1
         
       !SK-------------------------------------------------------------
       !SK   On ajoute une nouvelle direction optimale au stock
       !SK-------------------------------------------------------------


         if (reconj.and.abs(alpha)>min_coef_stock) then
            if (is_u)   usu(:,:,where_min_coef_stock)=uu(:,:,ndir)
            if (is_u)   usv(:,:,where_min_coef_stock)=uv(:,:,ndir)
            if (is_div) usp(:,:,where_min_coef_stock)=up(:,:,ndir)
            if (is_u)   csu(:,:,where_min_coef_stock)=cu(:,:,ndir)
            if (is_u)   csv(:,:,where_min_coef_stock)=cv(:,:,ndir)
            if (is_div) csp(:,:,where_min_coef_stock)=cp(:,:,ndir)
            old_where=where_min_coef_stock
            coef_stock(where_min_coef_stock)=abs(alpha)
            ndir_stock=min(ndir_stock+1,stockmax)
            min_coef_stock=10000
            do n=0,ndir_stock
               if (coef_stock(n).lt.min_coef_stock) then
                  min_coef_stock=coef_stock(n)
                  where_min_coef_stock=n
               endif
            enddo
            min_coef_stock=max(min_coef_stock,stock_critere)

!!$          !SK     Orthogonalisation
!!$
!!$          do n=0,ndir_stock-1
!!$             if (n.ne.old_where) then
!!$                alpha=global_ddot(csu(:,:,n),csv(:,:,n),csp(:,:,n) &
!!$                                 ,csu(:,:,old_where),csv(:,:,old_where),csp(:,:,old_where))
!!$                if (is_u)   usu(:,:,old_where)= usu(:,:,old_where)-alpha*usu(:,:,n)
!!$                if (is_v)   usv(:,:,old_where)= usv(:,:,old_where)-alpha*usv(:,:,n)
!!$                if (is_div) usp(:,:,old_where)= usp(:,:,old_where)-alpha*usp(:,:,n)
!!$                if (is_u)   csu(:,:,old_where)= csu(:,:,old_where)-alpha*csu(:,:,n)
!!$                if (is_v)   csv(:,:,old_where)= csv(:,:,old_where)-alpha*csv(:,:,n)
!!$                if (is_div) csp(:,:,old_where)= csp(:,:,old_where)-alpha*csp(:,:,n)
!!$             endif
!!$          enddo
!!$            
!!$          !SK     Normalisation
!!$
!!$          alpha=global_ddot(csu(:,:,old_where),csv(:,:,old_where),csp(:,:,old_where) &
!!$                           ,csu(:,:,old_where),csv(:,:,old_where),csp(:,:,old_where))
!!$          alpha=1_prec/sqrt(alpha)
!!$          if (is_u)   usu(:,:,old_where)= usu(:,:,old_where)*alpha
!!$          if (is_v)   usv(:,:,old_where)= usv(:,:,old_where)*alpha
!!$          if (is_div) usp(:,:,old_where)= usp(:,:,old_where)*alpha
!!$          if (is_u)   csu(:,:,old_where)= csu(:,:,old_where)*alpha
!!$          if (is_v)   csv(:,:,old_where)= csv(:,:,old_where)*alpha
!!$          if (is_div) csp(:,:,old_where)= csp(:,:,old_where)*alpha
!!$            
         endif

         !SK-------------------------------------------------------------           
         !SK     Nouvelle direction pour GMRESR 
         !SK-------------------------------------------------------------
            
         ndir=mod(ndir+1,ndirection-1)
         gmres_deb=max(gmres_deb,ndir)
         
         if (is_precond) then
            call precond(ru,rv,rp,uu(:,:,ndir),uv(:,:,ndir),up(:,:,ndir))
            call mavecxy(uu(:,:,ndir),uv(:,:,ndir),up(:,:,ndir)   &
                        ,cu(:,:,ndir),cv(:,:,ndir),cp(:,:,ndir))
         else
            print *, 'GMRES not implemented  with no precond'
            call parallel_stop
         endif


         !SK-------------------------------------------------------------           
         !SK     Orthogonalisation
         !SK-------------------------------------------------------------

         do n=0,gmres_deb-1
            if (n.ne.ndir) then
               alpha=global_ddot(cu(:,:,n),   cv(:,:,n)   ,cp(:,:,n) &
                                ,cu(:,:,ndir),cv(:,:,ndir),cp(:,:,ndir))
               if (is_u)   uu(:,:,ndir)= uu(:,:,ndir)-alpha*uu(:,:,n)
               if (is_v)   uv(:,:,ndir)= uv(:,:,ndir)-alpha*uv(:,:,n)
               if (is_div) up(:,:,ndir)= up(:,:,ndir)-alpha*up(:,:,n)
               if (is_u)   cu(:,:,ndir)= cu(:,:,ndir)-alpha*cu(:,:,n)
               if (is_v)   cv(:,:,ndir)= cv(:,:,ndir)-alpha*cv(:,:,n)
               if (is_div) cp(:,:,ndir)= cp(:,:,ndir)-alpha*cp(:,:,n)
            endif
         enddo


         !SK-------------------------------------------------------------           
         !SK     Normalisation
         !SK-------------------------------------------------------------

         alpha=1._prec/sqrt(global_ddot(cu(:,:,ndir),cv(:,:,ndir),cp(:,:,ndir) &
                               ,cu(:,:,ndir),cv(:,:,ndir),cp(:,:,ndir)))
         if (is_u)   uu(:,:,ndir)= alpha*uu(:,:,ndir)
         if (is_v)   uv(:,:,ndir)= alpha*uv(:,:,ndir)
         if (is_div) up(:,:,ndir)= alpha*up(:,:,ndir)
         if (is_u)   cu(:,:,ndir)= alpha*cu(:,:,ndir)
         if (is_v)   cv(:,:,ndir)= alpha*cv(:,:,ndir)
         if (is_div) cp(:,:,ndir)= alpha*cp(:,:,ndir)

         !SK-------------------------------------------------------------           
         !SK     Avance d'un pas pour VTU, VTV, R1X, R1Y
         !SK-------------------------------------------------------------

         alpha=global_ddot(cu(:,:,ndir),cv(:,:,ndir),cp(:,:,ndir),ru,rv,rp)
         if (is_u)   then;  VTU=VTU+alpha*uu(:,:,ndir);    ru=ru-alpha*cu(:,:,ndir);   endif
         if (is_v)   then;  VTV=VTV+alpha*uv(:,:,ndir);    rv=rv-alpha*cv(:,:,ndir);   endif
         if (is_div) then;  PRE=PRE+alpha*up(:,:,ndir);    rp=rp-alpha*cp(:,:,ndir);   endif
            
         err=sqrt(global_ddot(ru,rv,rp,ru,rv,rp))
         if (is_print>=3*niv_solve-1.and.is_div)&
            & errdiv=sqrt(global_sum(rp*rp))
         if (is_print>=3*niv_solve-1.and.my_task==0)&
              & print '(A,1X,I5,A,1X,I5,2E16.7)','niveau ',niv_solve&
              & ,' erreur GMRESR : ',ncg,err,errdiv
         if (is_cv_file.and.niv_solve==1)&
              & write (78,'(I5,1X,2E16.7)') ncg,err,errdiv
         if (is_div.and.niv_solve==1) call sortie_champ(PRE,'p')
         if (is_div.and.niv_solve==1) call sortie_champ(PRE-PRES,'X',surface)
         if (niv_solve==1) call sortie_champ(VTU-VTUS,'u',surface)
         if (niv_solve==1) call sortie_champ(VTV-VTVS,'v',surface)
         if (niv_solve==1) call sortie_champ(ru,'U',surface)
         if (niv_solve==1) call sortie_champ(rv,'V',surface)
         if (is_div.and.niv_solve==1) call sortie_champ(rp,'P',surface)
      enddo
     

      call conv_collect(err0,err,ncg-1)

      if (is_print>=3*niv_solve-2.and.my_task==0) then
         if (ncg.eq.ncgmax) then
            print '(A,1X,I5,A,1X,2E16.7)'&
                 & ,'Sorry! GMRESR n''a pas converge : ncg =',ncg&
                 & ,' err=',err,errdiv
         else
            print '(A,1X,I5,A,1X,2E16.7)'&
                 & ,'Bingo! GMRESR a converge : ncg =',ncg&
                 & ,' err=',err,errdiv
         endif
      endif

     return
  end subroutine solve_gmresr

  !************************************************************************ ! start_small  start_tiny

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive &
       subroutine solve_orthodir(SMU,SMV,SMP,VTU,VTV,PRE,ncgmax)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Methode ORTHODIR
    !SK
    !SK================================================================
    !SK
    implicit none

    integer ncgmax
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP

    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1) :: &
         pu,qu,ru,su,tu,uu,vu,wu
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1) :: &
         pv,qv,rv,sv,tv,uv,vv,wv
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1) :: &
         pp,qp,rp,sp,tp,up,vp,wp
    real(kind=prec) ::  alpha,beta,gamma,delta,sigma&
         & ,err0,err,err1,errdiv=0._prec
    integer :: restart,ncg


    call fixe_cl(SMU,SMV,SMP,VTU,VTV,PRE)
       
    if (is_u)   qu=0._prec;
    if (is_v)   qv=0._prec;
    if (is_div) qp=0._prec;

    call mavecxy(VTU,VTV,PRE,vu,vv,vp)

    if (is_u)     ru=SMU-vu
    if (is_v)     rv=SMV-vv
    if (is_div)   rp=SMP-vp

    err=10.
    err1=0.
    restart=0
    do while (err>epsv*err1.and.err>epsva.and.restart<=ncgmax)

       restart=restart+1

       
       if (is_precond) then
          call precond(ru,rv,rp,pu,pv,pp)
       else
          if (is_u)     pu=ru
          if (is_v)     pv=rv
          if (is_div)   pp=rp
       endif

       err=sqrt(global_ddot(ru,rv,rp,ru,rv,rp))
       err0=err
       if (restart==1) err1=err
       !     
       !=======================
       !     ITERATIONS
       !=======================
       !     
       
       sigma=0._prec
       do while (err>epsv*err0.and.err>epsva.and.ncg<=ninterne)

          ncg=ncg+1

          alpha=1._prec/sqrt(global_ddot(pu,pv,pp,pu,pv,pp))
          
          if (is_u)   pu=alpha*pu      ! normalisation de P
          if (is_v)   pv=alpha*pv
          if (is_div) pp=alpha*pp

          if (is_precond) then
             call mavecxy(pu,pv,pp,uu,uv,up)
             call precond(uu,uv,up,vu,vv,vp)
             call mavecxy(vu,vv,vp,wu,wv,wp)

             call precond(ru,rv,rp,su,sv,sp)
             call mavecxy(su,sv,sp,tu,tv,tp)
          else
             call mavecxy(pu,pv,pp,uu,uv,up)
             if (is_u)   vu=uu; 
             if (is_v)   vv=uv; 
             if (is_div) vp=up
             call mavecxy(uu,uv,up,wu,wv,wp)

             call mavecxy(ru,rv,rp,tu,tv,tp)
          endif

          beta = global_ddot(pu,pv,pp,wu,wv,wp)
          alpha= global_ddot(pu,pv,pp,tu,tv,tp)/beta

          if (is_u)    then; VTU=VTU+alpha*pu; ru=ru-alpha*uu; endif
          if (is_v)    then; VTV=VTV+alpha*pv; rv=rv-alpha*uv; endif
          if (is_div)  then; PRE=PRE+alpha*pp; rp=rp-alpha*up; endif

          err=sqrt(global_ddot(ru,rv,rp,ru,rv,rp))
          if (is_print>=3*niv_solve.and.is_div)&
            & errdiv=sqrt(global_sum(rp*rp))
          if (is_print>=3*niv_solve.and.my_task==0)&
               & print '(A,1X,I5,A,1X,I5,2E16.7)','niveau ',niv_solve&
               & ,' erreur ORTHODIR  : ',err,errdiv
          if (err<epsv*err0.or.err<epsva) exit

          if (is_precond) then
             call precond(wu,wv,wp,uu,uv,up)
             call mavecxy(uu,uv,up,wu,wv,wp)
          else
             if (is_u) uu=wu;   if (is_v) uv=wv;  if (is_div) up=wp
             call mavecxy(uu,uv,up,wu,wv,wp)
          endif   

          gamma=global_ddot(pu,pv,pp,wu,wv,wp)/beta
          if (ncg>1) sigma=global_ddot(qu,qv,qp,wu,wv,wp)/delta
          delta=beta
          
          if (is_u)   then; uu=pu; pu=vu-gamma*pu-sigma*qu; qu=uu; endif
          if (is_v)   then; uv=pv; pv=vv-gamma*pv-sigma*qv; qv=uv; endif
          if (is_div) then; up=pp; pp=vp-gamma*pp-sigma*qp; qp=up; endif

       enddo

       if (is_print>=3*niv_solve-1.and.my_task==0)&
            & print '(A,1X,I5,A,1X,I5,2E16.7)'&
            & ,'niveau ',niv_solve,' erreur ORTHODIR : '&
            & ,restart,err,errdiv
       if (is_cv_file.and.niv_solve==1)&
            & write(78,'(I5,1X,2E16.7)') restart,err,errdiv

    enddo


    call conv_collect(err1,err,restart-1)

    if (is_print>=3*niv_solve-2.and.my_task==0)  then
       if (restart.eq.ncgmax+1) then
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Sorry! ORTHODIR n''a pas converge : ncg =',restart&
               & ,' err=',err,errdiv
       else
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Bingo! ORTHODIR a converge : ncg =',ncg&
               & ,' err=',err,errdiv
       endif
    endif

    return
  end subroutine solve_orthodir

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive &
       subroutine solve_orthomin(SMU,SMV,SMP,VTU,VTV,PRE,ncgmax)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Methode ORTHOMIN
    !SK
    !SK================================================================
    !SK
    implicit none

    integer ncgmax
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP
    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1) :: ru,uu,MApu
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1) :: rv,uv,MApv
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1) :: rp,up,MApp
    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1,1:ndirection) :: pu,AMApu,MAMApu
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1,1:ndirection) :: pv,AMApv,MAMApv
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1,1:ndirection) :: pp,AMApp,MAMApp
    real(kind=prec), dimension(1:ndirection) :: delta
    
    real(kind=prec) ::  alpha,beta,err0,err,err1,errdiv=0._prec
    integer :: l,restart,ncg
    
    if (is_u)    then; pu=0._prec; AMApu=0._prec; MAMApu=0._prec; endif
    if (is_v)    then; pv=0._prec; AMApv=0._prec; MAMApv=0._prec; endif
    if (is_div)  then; pp=0._prec; AMApp=0._prec; MAMApp=0._prec; endif
    delta=0.

    call fixe_cl(SMU,SMV,SMP,VTU,VTV,PRE)

    call mavecxy(VTU,VTV,PRE,uu,uv,up)

    if (is_u)     ru=SMU-uu
    if (is_v)     rv=SMV-uv
    if (is_div)   rp=SMP-up

    err=10.
    err1=0.
    restart=0
    do while (err>epsv*err1.and.err>epsva.and.restart<=abs(ncgmax))
       
       restart=restart+1

       if (is_precond) then
          call precond(ru,rv,rp,pu(:,:,ndirection),pv(:,:,ndirection),pp(:,:,ndirection))
       else
          if (is_u)     pu(:,:,ndirection)=ru
          if (is_v)     pv(:,:,ndirection)=rv
          if (is_div)   pp(:,:,ndirection)=rp
       endif

       err=sqrt(global_ddot(ru,rv,rp,ru,rv,rp))
       err0=err
       if (restart==1) err1=err0
      if (is_print>=3*niv_solve-1.and.my_task==0)&
           & print '(A,1X,I5,A,1X,I5,2E16.7)'&
           & ,'niveau ',niv_solve&
           & ,' erreur ORTHOMIN(inside) : ',restart,err,errdiv
      if (is_cv_file.and.niv_solve==1)&
           & write (78,'(I5,1X,2E16.7)') restart,err,errdiv
       if (err<epsv*err1.or.err<epsva) exit
       !     
       !=======================
       !     ITERATIONS
       !=======================
       !     
    
       ncg=1
       do while (err>epsv*err0.and.err>epsva.and.ncg<=ninterne)

          if (is_precond) then
             call mavecxy(pu(:,:,ndirection) ,pv(:,:,ndirection) ,pp(:,:,ndirection), uu,uv,up)
             call precond(uu,uv,up,MApu,MApv,MApp)
             call mavecxy(MApu,MApv,MApp,AMApu(:,:,ndirection),AMApv(:,:,ndirection),AMApp(:,:,ndirection))
             call precond(AMApu (:,:,ndirection),AMApv (:,:,ndirection),AMApp (:,:,ndirection) &
                         ,MAMApu(:,:,ndirection),MAMApv(:,:,ndirection),MAMApp(:,:,ndirection))
             alpha = global_ddot(ru,rv,rp,MApu,MApv,MApp)
          else
             call mavecxy(pu(:,:,ndirection) ,pv(:,:,ndirection) ,pp(:,:,ndirection), uu,uv,up)
             call mavecxy(uu,uv,up,AMApu(:,:,ndirection),AMApv(:,:,ndirection),AMApp(:,:,ndirection))
             alpha = global_ddot(uu,uv,up,ru,rv,rp)
          endif

          delta(ndirection)=&
               global_ddot(pu(:,:,ndirection) ,pv(:,:,ndirection) ,pp(:,:,ndirection) &
                          ,AMApu(:,:,ndirection),AMApv(:,:,ndirection),AMApp(:,:,ndirection))

          alpha=alpha/delta(ndirection)

          if (is_u)    then; VTU=VTU+alpha*pu(:,:,ndirection); ru=ru-alpha*uu; endif
          if (is_v)    then; VTV=VTV+alpha*pv(:,:,ndirection); rv=rv-alpha*uv; endif
          if (is_div)  then; PRE=PRE+alpha*pp(:,:,ndirection); rp=rp-alpha*up; endif

          err=global_ddot(ru,rv,rp,ru,rv,rp)
          if (is_print>=3*niv_solve.and.is_div)&
            & errdiv=sqrt(global_sum(rp*rp))
          if (is_print>=3*niv_solve.and.my_task==0)&
               & print '(A,1X,I5,A,1X,I5,2E16.7)'&
               & ,'niveau ',niv_solve,' erreur ORTHOMIN : ',err,errdiv

          if (is_precond) then
             call precond(ru,rv,rp,uu,uv,rp)
          else
             if (is_u)    uu=ru
             if (is_v)    uv=rv
             if (is_div)  up=rp
          endif

          do l=max(1,ndirection-ncg+1),ndirection
             if (is_precond) then
                beta=-global_ddot(ru,rv,rp,MAMApu(:,:,l),MAMApv(:,:,l),MAMApp(:,:,l))/delta(l)
             else                
                beta=-global_ddot(ru,rv,rp,AMApu(:,:,l),AMApv(:,:,l),AMApp(:,:,l))/delta(l)
             endif
             if (is_u)    uu=uu+beta*pu(:,:,l)
             if (is_v)    uv=uv+beta*pv(:,:,l)
             if (is_div)  up=up+beta*pp(:,:,l)
          enddo

          do l=max(1,ndirection-ncg),ndirection-1
             if (is_u)   pu(:,:,l)   =pu(:,:,l+1)
             if (is_v)   pv(:,:,l)   =pv(:,:,l+1)
             if (is_div) pp(:,:,l)   =pp(:,:,l+1)
             if (is_u)   AMApu(:,:,l)=AMApu(:,:,l+1)
             if (is_v)   AMApv(:,:,l)=AMApv(:,:,l+1)
             if (is_div) AMApp(:,:,l)=AMApp(:,:,l+1)
             if (is_precond) then
                if (is_u)   MAMApu(:,:,l)=MAMApu(:,:,l+1)
                if (is_v)   MAMApv(:,:,l)=MAMApv(:,:,l+1)
                if (is_div) MAMApp(:,:,l)=AMApp(:,:,l+1)
             endif
             delta(l)=delta(l+1)
          enddo

          if (is_u)    pu(:,:,ndirection)= uu
          if (is_v)    pv(:,:,ndirection)= uv
          if (is_div)  pp(:,:,ndirection)= up
          ncg=ncg+1
      enddo

    enddo

    call conv_collect(err0,err,restart-1)


    if (is_print>=3*niv_solve-2.and.my_task==0) then
       if (restart.eq.ncgmax) then
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Sorry! ORTHOMIN n''a pas converge : ncg =',restart&
               & ,' err=',err,errdiv
       else
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Bingo! ORTHOMIN a converge : ncg =',restart&
               & ,' err=',err,errdiv
       endif
    endif

    return
  end subroutine solve_orthomin

  !************************************************************************    ! end_small   

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive &
       subroutine mg_solve(SMU,SMV,SMP,VTU,VTV,PRE,level)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Noyau de resolution multigrille
    !SK
    !SK================================================================
    !SK
    use data,     only : nb_prelis,nb_postlis,nb_cycle,nbmg
    use disc,     only : grid
    use init,     only : set_grid
    use mg
    implicit none
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP
    integer :: level

    real(kind=prec), dimension(:,:), pointer :: U1,U2,U3,V1,V2,V3,P1,P2,P3
    integer :: k

    U1=>grid(level)%U1;     U2=>grid(level)%U2; U3=>grid(level)%U3; 
    V1=>grid(level)%V1;     V2=>grid(level)%V2; V3=>grid(level)%V3; 
    P1=>grid(level)%P1;     P2=>grid(level)%P2; P3=>grid(level)%P3; 

    if (.not.is_decal) call fixe_cl(SMU,SMV,SMP,VTU,VTV,PRE)

    call set_grid(level)

    if (level<=0) then
       niv_solve=niv_solve+1
       if (is_div) then
          call lisse(SMU,SMV,SMP,VTU,VTV,PRE,nbmg);
       else
          call solve_bicgstab(SMU,SMV,SMP,VTU,VTV,PRE,1000)
       endif
!       call solve_bicgstab_l(SMU,SMV,SMP,VTU,VTV,PRE,1000,3)
       niv_solve=niv_solve-1
    else
       call lisse(SMU,SMV,SMP,VTU,VTV,PRE,nb_prelis);

!!$       print *,'avt mavecxy'
!!$       call flush(6)
!!$       print *,"hellllllllllo",size(VTU,1),size(VTV,1),size(PRE,1),size(U1,1),size(V1,1),size(P1,1)
!!$       print *,"hellllllllllo",size(VTU,2),size(VTV,2),size(PRE,2),size(U1,2),size(V1,2),size(P1,2)
!!$       call flush(6)
       call mavecxy(VTU,VTV,PRE,U1,V1,P1);
!!$       print *,'apres mavecxy'
!!$       call flush(6)

       if (is_u) U1=SMU-U1;        if (is_v) V1=SMV-V1;         if (is_div) P1=SMP-P1;
       call restrict(U1,V1,P1,U2,V2,P2,.false.)
       

       do k=1,nb_cycle
          call  mg_solve(U2,V2,P2,U3,V3,P3,level-1);
       enddo

       call set_grid(level)
       call prolong(U3,V3,P3,U1,V1,P1)
       if (is_u) VTU=VTU+U1;       if (is_v) VTV=VTV+V1;       if (is_div) PRE=PRE+P1;
       call lisse(SMU,SMV,SMP,VTU,VTV,PRE,nb_postlis);
    endif

    return
  end subroutine mg_solve

  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive &
       subroutine fmg_start(SMU,SMV,SMP,VTU,VTV,PRE,level)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Noyau de resolution multigrille
    !SK
    !SK================================================================
    !SK
    use data,     only : nb_prelis,nb_postlis,nb_cycle,nbmg
    use disc,     only : grid
    use init,     only : set_grid
    use mg
    implicit none
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP
    integer :: level

    real(kind=prec), dimension(:,:), pointer :: U2,U3,V2,V3,P2,P3

    U2=>grid(level)%U2; U3=>grid(level)%U3; 
    V2=>grid(level)%V2; V3=>grid(level)%V3; 
    P2=>grid(level)%P2; P3=>grid(level)%P3; 

    if (.not.is_decal) call fixe_cl(SMU,SMV,SMP,VTU,VTV,PRE)

    if (level<=0) then
       niv_solve=niv_solve+1
       call set_grid(level)
       if (is_div) then
          call lisse(SMU,SMV,SMP,VTU,VTV,PRE,nbmg);
       else
          call solve_bicgstab(SMU,SMV,SMP,VTU,VTV,PRE,1000)
       endif
       niv_solve=niv_solve-1
    else
       call restrict(SMU,SMV,SMP,U2,V2,P2,.false.)

       call  fmg_start(U2,V2,P2,U3,V3,P3,level-1);
       call prolong(U3,V3,P3,VTU,VTV,PRE)

       call mg_solve(SMU,SMV,SMP,VTU,VTV,PRE,level)

    endif

    return
  end subroutine fmg_start

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive subroutine solve_mg(SMU,SMV,SMP,VTU,VTV,PRE,nmgmax,is_fmg)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Solveur Multigrille
    !SK
    !SK================================================================
    !SK
    use mg
    use data,     only : nb_prelis,nb_postlis,nb_cycle,nbmg
    use data,     only : ncheck
    use disc   
    use init,     only : set_grid
    implicit none

    integer nmgmax
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP
    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1) :: RU
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1) :: RV
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1) :: RP

    real(kind=prec)  :: err0,err,errdiv=0._prec
    integer          :: nmg
    logical          :: is_fmg

    if (ncheck/=df2) stop 'Multigrille uniquement implemente en DF2 '

    if (.not.is_decal) call fixe_cl(SMU,SMV,SMP,VTU,VTV,PRE)

    call set_grid(ngrid_max)
    call mavecxy(VTU,VTV,PRE,RU,RV,RP)

    if (is_u)  RU=SMU-RU;    if (is_v)  RV=SMV-RV;    if (is_div) RP=SMP-RP

    err0=sqrt(global_ddot(RU,RV,RP,RU,RV,RP))
    err=err0
    nmg=0

    if (is_fmg)  call fmg_start(SMU,SMV,SMP,VTU,VTV,PRE,ngrid_max)

    do while (err>epsv*err0.and.err>epsva.and.nmg.lt.abs(nmgmax))

       nmg=nmg+1

       call mg_solve(SMU,SMV,SMP,VTU,VTV,PRE,ngrid_max)

       call mavecxy(VTU,VTV,PRE,RU,RV,RP)
       if (is_u)  RU=SMU-RU;    if (is_v)  RV=SMV-RV;    if (is_div) RP=SMP-RP
       err=sqrt(global_ddot(RU,RV,RP,RU,RV,RP))
       if (is_print>=3*niv_solve-1.and.is_div)&
            & errdiv=sqrt(global_sum(rp*rp))

       if (is_print>=3*niv_solve-1.and.my_task==0)&
            & print '(A,1X,I5,A,1X,I5,2E16.7)','niveau ',niv_solve&
            & ,' erreur MG : ',nmg,err,errdiv
       if (is_cv_file.and.niv_solve==1)&
            & write (78,'(I5,1X,2E16.7)')&
            & nmg,err,errdiv

       if (is_div.and.niv_solve==1) call sortie_champ(PRE,'p')
       if (is_div.and.niv_solve==1) call sortie_champ(PRE-PRES,'X',surface)
       if (niv_solve==1) call sortie_champ(VTU-VTUS,'u',surface)
       if (niv_solve==1) call sortie_champ(VTV-VTVS,'v',surface)
       if (niv_solve==1) call sortie_champ(ru,'U',surface)
       if (niv_solve==1) call sortie_champ(rv,'V',surface)
       if (is_div.and.niv_solve==1) call sortie_champ(rp,'P',surface)

    enddo

    call conv_collect(err0,err,nmg-1)

    if (is_print>=3*niv_solve-2.and.my_task==0) then
       if (nmg.eq.nmgmax) then
          print *,'Sorry! MG n''a pas converge : nmg =',nmg,' err=',err,errdiv
       else
          print *,'Bingo! MG a converge : nmg =',nmg,' err=',err,errdiv
       endif
    endif

    return

  end subroutine solve_mg

  !************************************************************************   

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive subroutine solve_lisse(SMU,SMV,SMP,VTU,VTV,PRE,nmgmax)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Solveur Multigrille
    !SK
    !SK================================================================
    !SK
    use mg
    use data,     only : nb_prelis,nb_postlis,nb_cycle,nbmg
    use data,     only : ncheck
    use disc,     only : ngrid_max
    use init,     only : set_grid
    implicit none

    integer nmgmax
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP
    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1) :: RU
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1) :: RV
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1) :: RP

    real(kind=prec)  :: err0,err,errdiv
    integer           :: nmg

    if (.not.is_decal) call fixe_cl(SMU,SMV,SMP,VTU,VTV,PRE)

    if (ncheck/=df2) stop 'Multigrille uniquement implemente en DF2 '

    call set_grid(ngrid_max)
    call mavecxy(VTU,VTV,PRE,RU,RV,RP)

    if (is_u)  RU=SMU-RU;    if (is_v)  RV=SMV-RV;    if (is_div) RP=SMP-RP

    err0=sqrt(global_ddot(RU,RV,RP,RU,RV,RP))
    err=err0
    nmg=0

    do while (err>epsv*err0.and.err>epsva.and.nmg.lt.abs(nmgmax))

       nmg=nmg+1

       call lisse(SMU,SMV,SMP,VTU,VTV,PRE,-1);

       call mavecxy(VTU,VTV,PRE,RU,RV,RP)
       if (is_u)  RU=SMU-RU;    if (is_v)  RV=SMV-RV;    if (is_div) RP=SMP-RP
       err=sqrt(global_ddot(RU,RV,RP,RU,RV,RP))
       if (is_div) errdiv=sqrt(global_sum(RP*RP))

       if (is_print>=3*niv_solve-1.and.my_task==0)&
            & print '(A,1X,I5,A,1X,I5,2E16.7)'&
            & ,'niveau ',niv_solve,' erreur lisseur : ',nmg&
            & ,err,errdiv
       if (is_cv_file.and.niv_solve==1)&
            & write (78,'(I5,1X,2E16.7)')&
            & nmg,err,errdiv
       if (is_div.and.niv_solve==1) call sortie_champ(PRE,'p')
       if (is_div.and.niv_solve==1) call sortie_champ(PRE-PRES,'X',surface)
       if (niv_solve==1) call sortie_champ(VTU-VTUS,'u',surface)
       if (niv_solve==1) call sortie_champ(VTV-VTVS,'v',surface)
       if (niv_solve==1) call sortie_champ(ru,'U',surface)
       if (niv_solve==1) call sortie_champ(rv,'V',surface)
       if (is_div.and.niv_solve==1) call sortie_champ(rp,'P',surface)

!!$       if (mod(nmg,nbmg)==0) then
!!$          call outdon(RU,'RU')
!!$          call outdon(RV,'RP')
!!$          call outdon(SMP,'SMP')
!!$          call outdon(RP,'RP')
!!$       endif

    enddo

    call conv_collect(err0,err,nmg-1)


    if (is_print>=3*niv_solve-2.and.my_task==0) then
       if (nmg.eq.nmgmax) then
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Sorry! lisseur n''a pas converge : nmg =',nmg&
               & ,' err=',err,errdiv
       else
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Bingo! lisseur a converge : nmg =',nmg&
               & ,' err=',err,errdiv
       endif
    endif

    return

  end subroutine solve_lisse
 

  !************************************************************************              ! end_tiny

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  recursive subroutine solve_shwarz(SMU,SMV,SMP,VTU,VTV,PRE,nmgmax)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Solveur Schwarz
    !SK
    !SK================================================================
    !SK
    use data,     only : ncheck
    implicit none

    integer nmgmax
    real(kind=prec), dimension(0:,0:) :: VTU,VTV,PRE,SMU,SMV,SMP
    real(kind=prec), dimension(0:size(SMU,1)-1,0:size(SMU,2)-1) :: RU
    real(kind=prec), dimension(0:size(SMV,1)-1,0:size(SMV,2)-1) :: RV
    real(kind=prec), dimension(0:size(SMP,1)-1,0:size(SMP,2)-1) :: RP

    real(kind=prec)  :: err0,err,errdiv=0._prec
    integer           :: nmg,lmu,lmv,lmp,nmu,nmv,nmp


    call output('is_mpp av',l1=is_mpp)


    call fixe_cl(SMU,SMV,SMP,VTU,VTV,PRE,.true.);
    if (is_decal) then
       lmu=size(VTU,1)-2; nmu=size(VTU,2)-2; 
       lmv=size(VTV,1)-2; nmv=size(VTV,2)-2; 
       lmp=size(PRE,1)-2; nmp=size(PRE,2)-2; 
       select case (ncheck)
          case (df2,vf2)
             if (is_south) VTV(:,0)=2.*SMV(:,0)-VTV(:,1)
             if (is_north) VTV(:,nmv+1)=2.*SMV(:,nmv+1)-VTV(:,nmv)
             if (is_west)  VTU(0,:)=2.*SMU(0,:)-VTU(1,:)
             if (is_east)  VTU(lmu+1,:)=2.*SMU(lmu+1,:)-VTU(lmu,:)
          case (df4,vf4)
             if (is_south) VTV(:,0)=(SMV(:,0)&
                  & -(0.9375_prec*VTV(:,1)-0.3125_prec*VTV(:,2)+0.0625_prec*VTV(:,3)))/0.3125_prec
             if (is_north) VTV(:,nmv+1)=(SMV(:,nmv+1)&
                  & -(0.9375_prec*VTV(:,nmv)-0.3125_prec*VTV(:,nmv-1)+0.0625_prec*VTV(:,nmv-2)))/0.3125_prec
             if (is_west)  VTU(0,:)=(SMU(0,:)&
                  & -(0.9375_prec*VTU(1,:)-0.3125_prec*VTU(2,:)+0.0625_prec*VTU(3,:)))/0.3125_prec
             if (is_east)  VTU(lmu+1,:)=(SMU(lmu+1,:)&
                  & -(0.9375_prec*VTU(lmu,:)-0.3125_prec*VTU(lmu-1,:)+0.0625_prec*VTU(lmu-2,:)))/0.3125_prec
          end select
    endif

    if (is_u) RU=VTU;  if (is_v) RV=VTV;  if (is_div) RP=PRE
    err0=sqrt(global_ddot(RU,RV,RP,RU,RV,RP))
    err=err0
    nmg=0


    do while (err>epsv*err0.and.err>epsva.and.nmg.lt.abs(nmgmax))

       nmg=nmg+1
       
       if (is_u) RU=VTU;  if (is_v) RV=VTV;  if (is_div) RP=PRE

       call precond(SMU,SMV,SMP,VTU,VTV,PRE);
       
       if (nmg>1) then
          if (is_u) RU=RU-VTU;  if (is_v) RV=RV-VTV;  if (is_div) RP=RP-PRE
          err=sqrt(global_ddot(RU,RV,RP,RU,RV,RP))
          if (is_div) errdiv=sqrt(global_sum(RP*RP))


          if (is_print>=3*niv_solve-1.and.my_task==0)&
               & print '(A,1X,I5,A,1X,I5,2E16.7)'&
               & ,'niveau ',niv_solve,' erreur Schwarz  : ',nmg&
               & ,err,errdiv
          if (is_cv_file.and.niv_solve==1)&
               & write (78,'(I5,1X,2E16.7)')&
               & nmg,err,errdiv
       endif

       if (niv_solve==1) call sortie_champ(VTU-VTUS,'U',surface)
       if (niv_solve==1) call sortie_champ(VTV-VTVS,'V',surface)
       if (is_div) PRE=PRE-PRE(0,0)+PRES(0,0)
!!$       if (is_div) PRE(lmp-2:lmp+1,:)=PRES(lmp-2:lmp+1,:)
       if (is_div.and.niv_solve==1) call sortie_champ(PRE,'P',surface)

    enddo

!    call fixe_cl(VTU,VTV,PRE,VTU,VTV,PRE,.true.)

    call conv_collect(err0,err,nmg-1)

    call output('is_mpp ap',l1=is_mpp)

    if (is_print>=3*niv_solve-2.and.my_task==0) then
       if (nmg.eq.nmgmax) then
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Sorry! Schwarz n''a pas converge : nmg =',nmg&
               & ,' err=',err,errdiv
       else
          print '(A,1X,I5,A,1X,2E16.7)'&
               & ,'Bingo! Schwarz a converge : nmg =',nmg&
               & ,' err=',err,errdiv
       endif
    endif

    return

  end subroutine solve_shwarz
 
  !************************************************************************   ! end_out_light  


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine prepare_adv_grids(VTU,VTV)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc
    use mg
    use init,     only : set_grid
    use para,     only : rfr_stn, stn_news
    use data,     only : ncheck, ncheck_precond, nrecvddmx, nrecvddmy
    use compact,  only : derm0x,derm0y
    use drapeaux, only : is_decal,is_precond,is_conv
    use champs, only : SMUp,SMVp
    implicit none
    real(kind=prec), dimension(0:,0:) :: VTU,VTV
    real(kind=prec), dimension(0:size(VTV,1)+2*nrecvddmx+1,0:size(VTV,2)+2*nrecvddmy+1) :: TMPy   ! out_light
    real(kind=prec), dimension(0:size(VTU,1)+2*nrecvddmx+1,0:size(VTU,2)+2*nrecvddmy+1) :: TMPx   ! out_light
    integer :: lmu,nmu,lmv,nmv,k,i,lmu1,nmu1,lmv1,nmv1


    lmu=size(VTU,1)-2
    nmu=size(VTU,2)-2
    lmv=size(VTV,1)-2
    nmv=size(VTV,2)-2


!!$    VTU=XU/dx+YU/dy
!!$    VTV=XV/dx+YV/dx

!!$    call outdon(VTU,'VTU')
!!$    call outdon(VTV,'VTV')
!!$    

    call set_grid(-1)

    select case (ncheck)
    case(vf4,vf4e)                                                    ! start_out_light
       U_EN_U=derm0x(VTU,lmu_global)
       V_EN_V=derm0y(VTV,nmv_global)
    case(vf2)
       if (is_mpp) call rfr_stn(VTU,stn_news,1,1)
       if (is_mpp) call rfr_stn(VTV,stn_news,1,1)
       if (is_west) U_EN_U(0,:)        = VTU(0,:)
       U_EN_U(1:lmu+1,:)   = 0.5_prec*(VTU(0:lmu,:)+VTU(1:lmu+1,:))
       if (is_east) U_EN_U(lmu+2,:)    = VTU(lmu+1,:)
       if (is_south) V_EN_V(:,0)       = VTV(:,0)
       V_EN_V(:,1:nmv+1)   = 0.5_prec*(VTV(:,0:nmv)+VTV(:,1:nmv+1))
       if (is_north) V_EN_V(:,nmv+2)   = VTV(:,nmv+1)
    case(df4,df4e,df2)
       U_EN_U = VTU
       V_EN_V = VTV
    end select

    if (is_decal) then
       U_en_V=0._prec
       V_en_U=0._prec
       select case(ncheck)
       case(vf4,vf4e)      
          lmu1=lmu+2;   nmu1=nmu+2
          lmv1=lmv+2;   nmv1=nmv+2
          TMPy(0:lmu+1,0:nmu+2)=derm0y(VTU,nmu_global)
          TMPx(0:lmv+2,0:nmv+1)=derm0x(VTV,lmv_global)
          U_EN_V(0:lmv1,0:nmv+1) = TMPy(0:lmv1,0:nmv+1)
          V_EN_U(0:lmu+1,0:nmu1) = TMPx(0:lmu+1,0:nmu1)
          if (.not.is_west) call snd_msg(p_west,tag_adv_x,U_en_V(1,:))
          if (.not.is_east) call rcv_msg(p_east,tag_adv_x,U_en_V(lmv+1,:))
          if (.not.is_south) call snd_msg(p_south,tag_adv_y,V_en_U(:,1))
          if (.not.is_north) call rcv_msg(p_north,tag_adv_y,V_en_U(:,nmu+1))
          tag_adv_x=tag_adv_x+1
          tag_adv_y=tag_adv_y+1
       case(df4,df4e) 
          TMPx(0:lmu+2,0:nmu+1)= derm0x(VTU,lmu_global)
          TMPy(0:lmu+2,0:nmu+2) = derm0y(TMPx(0:lmu+2,0:nmu+1),nmu_global)
          U_EN_V(0:lmv+1,1:nmv) = TMPy(1:lmv+2,1:nmv)
          TMPy(0:lmv+1,0:nmv+2)= derm0y(VTV,nmv_global)
          TMPx(0:lmv+2,0:nmv+2) = derm0x(TMPy(0:lmv+1,0:nmv+2),lmv_global)
          V_EN_U(1:lmu,0:nmu+1) = TMPx(1:lmu,1:nmu+2)
       case(vf2)
          V_EN_U(1:lmu,1:nmu+1)   = 0.5_prec*(VTV(1:lmu,1:nmu+1)+VTV(0:lmu-1,1:nmu+1))
          U_EN_V(1:lmv+1,1:nmv)   = 0.5_prec*(VTU(1:lmv+1,1:nmv)+VTU(1:lmv+1,0:nmv-1))
       case(df2)
          if (is_mpp) call rfr_stn(VTU,stn_news,1,1)
          if (is_mpp) call rfr_stn(VTV,stn_news,1,1)
          V_EN_U(1:lmu,1:nmu) = 0.25_prec*(VTV(1:lmu,  1:nmu)+VTV(1:lmu  ,2:nmu+1)       &
                                          +VTV(0:lmu-1,1:nmu)+VTV(0:lmu-1,2:nmu+1))
          U_EN_V(1:lmv,1:nmv) = 0.25_prec*(VTU(1:lmv  ,1:nmv)+VTU(2:lmv+1,1:nmv)        &
                                          +VTU(1:lmv  ,0:nmv-1)+VTU(2:lmv+1,0:nmv-1))
       end select
    endif

!!$    call outdon(U_en_U,'U_en_U')
!!$    call outdon(V_en_U,'V_en_U')
!!$    call outdon(U_en_V,'U_en_V')
!!$    call outdon(V_en_V,'V_en_V')
!!$    call parallel_stop

    
    if (ngrid_max>=0) then                                         
       SMUp(iudeb:iufin,kudeb:kufin) = VTU
       SMVp(ivdeb:ivfin,kvdeb:kvfin) = VTV

       call rfr_ddm_stn(SMUp,iudeb,iufin,kudeb,kufin,nrecvddmx,nrecvddmy)
       call rfr_ddm_stn(SMVp,ivdeb,ivfin,kvdeb,kvfin,nrecvddmx,nrecvddmy)

       grid(ngrid_max)%U_EN_U=SMUp(iumgdeb:iumgfin,kumgdeb:kumgfin)
       grid(ngrid_max)%V_EN_V=SMVp(ivmgdeb:ivmgfin,kvmgdeb:kvmgfin)

       do k=ngrid_max-1,0,-1                                        ! start_out_light
          call restrict(grid(k+1)%U_EN_U(:,:),grid(k+1)%V_EN_V(:,:),TMPx&
               &       ,grid(k)  %U_EN_U(:,:),grid(k)  %V_EN_V(:,:),TMPx,.true.)
       enddo                                                        ! end_out_light

!!$    do k=ngrid_max,-1,-1
!!$       print *
!!$       print *,'niveau ',k
!!$       call  outdon(grid(k)%U_EN_U,'grid(k)%U_EN_U')
!!$       call  outdon(grid(k)%V_EN_V,'grid(k)%V_EN_V')                                 
!!$    enddo                                                     

       if (ncheck_precond==vf4) stop 'xxxxxxxx vf4 a faire'

       if (ncheck_precond==vf2) then
          do k=ngrid_max,0,-1
             lmu=size(grid(k)%U_EN_U(:,:),1)-3; 
             nmv=size(grid(k)%V_EN_V(:,:),2)-3; 
             grid(k)%U_EN_U(0,:)       = grid(k)%U_EN_U(0,:)
             grid(k)%U_EN_U(lmu+2,:)   = grid(k)%U_EN_U(lmu+1,:)
             grid(k)%U_EN_U(1:lmu+1,:) = 0.5_prec*(grid(k)%U_EN_U(0:lmu,:)+grid(k)%U_EN_U(1:lmu+1,:))
             grid(k)%V_EN_V(:,nmv+2)   = grid(k)%V_EN_V(:,nmv+1)
             grid(k)%V_EN_V(:,0)       = grid(k)%V_EN_V(:,0)
             grid(k)%V_EN_V(:,1:nmv+1) = 0.5_prec*(grid(k)%V_EN_V(:,0:nmv)+grid(k)%V_EN_V(:,1:nmv+1))
          enddo
       endif
          
!!$    call output('iumgdeb,iumgfin,kumgdeb,kumgfin',iumgdeb,iumgfin,kumgdeb,kumgfin)  ! start_out_light
!!$    call output('size(VTU)',size(vtu,1),size(VTU,2))
!!$    call outdon(VTU,'VTU')
!!$    call outdon(grid(ngrid_max)%U_EN_U,'grid(ngrid_max)%U_EN_U')
!!$    call parallel_stop

     endif

    if (is_decal) then                                            ! start_out_light
          do k=ngrid_max,0,-1
             call set_grid(k)
    
             lmu=size(U_EN_U,1)-2
             nmu=size(U_EN_U,2)-2
             lmv=size(V_EN_V,1)-2
             nmv=size(V_EN_V,2)-2

             i=ncheck_precond; is_mpp=.false.
             if (k==-1) then; i=ncheck; is_mpp=nb_tasks>1; endif;
             call output('grid k',k,i)
             call output('size(U_en_U)',size(U_en_U,1),size(U_en_U,2))
             call output('size(U_en_V)',size(U_en_V,1),size(U_en_V,2))            
             call output('size(V_en_U)',size(V_en_U,1),size(V_en_U,2))
             call output('size(V_en_V)',size(V_en_V,1),size(V_en_V,2))            
             call output('size(TMPx)',size(TMPx,1),size(TMPx,2))
             call output('size(TMPy)',size(TMPy,1),size(TMPy,2))
             
             select case(i)
             case(vf4,vf4e)      
                U_en_V=0._prec
                V_en_U=0._prec
                lmu1=lmu+2
                nmu1=nmu+2
                lmv1=lmv+2
                nmv1=nmv+2
                TMPy(0:lmu+1,0:nmu1)=derm0y(SMUp,nmu_global)
                TMPx(0:lmv1,0:nmv+1)=derm0x(SMVp,lmv_global)
                U_EN_V(0:lmv+2,1:nmv) = TMPy(1:lmu+1,1:nmu+1)
                V_EN_U(1:lmu,0:nmu+2) = TMPx(1:lmv+1,1:nmv+1)
             case(df4,df4e) 
                lmu1=lmu+2
                nmu1=nmu+2
                lmv1=lmv+2
                nmv1=nmv+2
                TMPx(0:lmu1,0:nmu+1)=derm0x(U_en_U,lmu_global)
                TMPy(0:lmu1,0:nmu1) = derm0y(TMPx(0:lmu1,0:nmu+1),nmu_global)
                U_EN_V(0:lmv+1,1:nmv) = TMPy(1:lmv+2,1:nmv)
                TMPy(0:lmv+1,0:nmv1)=derm0y(V_en_V,nmv_global)
                TMPx(0:lmv1,0:nmv1) = derm0x(TMPy(0:lmv+1,0:nmv1),lmv_global)
                V_EN_U(1:lmu,0:nmu+1) = TMPx(1:lmu,1:nmu+2)
             case(vf2)
                TMPx(1:lmu,0:nmu+1) = 0.25_prec*(V_en_V(0:lmv,  0:nmv  )+V_en_V(1:lmv+1,0:nmv)          &
                     +V_en_V(0:lmv , 1:nmv+1)+V_en_V(1:lmv+1,1:nmv+1))
                V_EN_U(0,:)             = 0._prec         
                V_EN_U(lmu+1,:)         = 0._prec         
                V_EN_U(1:lmu,0)         = TMPx(1:lmu,0)
                V_EN_U(1:lmu,1:nmu+1)   = 0.5_prec*(TMPx(1:lmu,0:nmu)+TMPx(1:lmu,1:nmu+1))
                V_EN_U(1:lmu,nmu+2)     = TMPx(1:lmu,nmu+1)
                
                TMPy(0:lmv+1,1:nmv) = 0.25_prec*(U_en_U(0:lmu,0:nmu    )+U_en_U(1:lmu+1,0:nmu)          &
                     +U_en_U(0:lmu,1:nmu+1  )+U_en_U(1:lmu+1,1:nmu+1))
                U_EN_V(:,0)             = 0._prec
                U_EN_V(:,nmv+1)         = 0._prec
                U_EN_V(0,1:nmv)         = TMPy(0,1:nmv)
                U_EN_V(1:lmv+1,1:nmv)   = 0.5_prec*(TMPy(0:lmv,1:nmv)+TMPy(1:lmv+1,1:nmv))
                U_EN_V(lmv+2,1:nmv)     = TMPy(lmv+1,1:nmv)
             case(df2)
                V_EN_U(1:lmu,1:nmu) =                    &
                     0.25_prec*                                  &
                     (V_EN_V(1:lmu,  1:nmu)          &
                     +V_EN_V(1:lmu  ,2:nmu+1)       &
                     +V_EN_V(0:lmu-1,1:nmu)          &
                     +V_EN_V(0:lmu-1,2:nmu+1))
                U_EN_V(1:lmv,1:nmv) =                    &
                     0.25_prec*                                  &
                     (U_EN_U(1:lmv  ,1:nmv)            &
                     +U_EN_U(2:lmv+1,1:nmv)        &
                     +U_EN_U(1:lmv  ,0:nmv-1)          &
                     +U_EN_U(2:lmv+1,0:nmv-1))
             end select
          enddo
          is_mpp=nb_tasks>1
    endif                                                           ! end_out_light



    return
  end subroutine prepare_adv_grids

  !************************************************************************   



end module calcul


 
