module blas

  use selectprec
  use constante
  use para

contains

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine mavecxy(INPU,INPV,INPP,OUTU,OUTV,OUTP)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data, only : ncheck
    use drapeaux
    implicit none
    real(kind=prec), dimension(:,:) :: INPU,INPV,INPP,OUTU,OUTV,OUTP

!!$    print *,'dans mavecxy'
!!$    call flush(6)

    select case(ncheck)
    case(vf4,vf4e);    call mavecxy_vf4(INPU,INPV,INPP,OUTU,OUTV,OUTP)    ! out_light odut_tiny
    case(df4,df4e);    call mavecxy_df4(INPU,INPV,INPP,OUTU,OUTV,OUTP)    ! out_light
    case(vf2);         call mavecxy_df2(INPU,INPV,INPP,OUTU,OUTV,OUTP,.false.)
    case(df2);         call mavecxy_df2(INPU,INPV,INPP,OUTU,OUTV,OUTP,.true.)
    case default;      print *,'ncheck=',ncheck; 
                       stop 'shema non encore implemente in mavecxy'
    end select

    return
  end subroutine mavecxy

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine mavecxy_df2(INPU,INPV,INPP,OUTU,OUTV,OUTP,flag)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    use drapeaux
    use para
    use data,   only : nu, rlag, ncheck
    use disc,   only : dx,dy,invdx,invdy,invdx2,invdy2,p5invdx,p5invdy &
         ,coeft,U_EN_U_c, U_EN_V_c, V_EN_V_c, V_EN_U_c &
         ,U_EN_U, U_EN_V, V_EN_V, V_EN_U
    implicit none

    real(kind=prec), dimension(0:,0:), target :: OUTU,OUTV,OUTP,INPU,INPV,INPP
    logical :: flag

    real(kind=prec), dimension(0:size(INPU,1)-1,0:size(INPU,2)-1),target :: TMPU
    real(kind=prec), dimension(0:size(INPV,1)-1,0:size(INPV,2)-1),target :: TMPV

    real(kind=prec), dimension(:,:),pointer ::     &
         INPU_c,  INPV_c,  OUTU_c,  OUTV_c         &
        ,INPU_n,  INPU_e,  INPU_w,  INPU_s         &
        ,INPV_n,  INPV_e,  INPV_w,  INPV_s         &
        ,TMPU_u,TMPU_v,TMPU_n,TMPU_e,TMPU_w,TMPU_s &
        ,TMPV_u,TMPV_v,TMPV_n,TMPV_e,TMPV_w,TMPV_s 

    real(kind=prec), dimension(:,:),pointer ::  &    ! start_out_light
         INPVp_n, INPUp_e, INPUp_w, INPVp_s     &     
        ,INPP_n,  INPP_e,  INPP_w,  INPP_s      &     
        ,OUTP_n,  OUTP_e,  OUTP_w,  OUTP_s      &    
        ,INPU_ne, INPU_se ,INPU_nw, INPU_sw     &  
        ,INPV_ne, INPV_se ,INPV_nw, INPV_sw     &
        ,INPP_ne, INPP_se ,INPP_nw, INPP_sw     &        
        ,OUTP_ne, OUTP_se ,OUTP_nw, OUTP_sw, OUTP_c  ! end_out_light

    real(kind=prec) :: coef
    integer  :: lmu,nmu,lmv,nmv 
    integer :: lmp,nmp,ip0,ip1,kp0,kp1           ! out_light

!!$    print *,"hellllllllllo",size(INPU,1),size(INPV,1),size(INPP,1),size(OUTU,1),size(OUTV,1),size(OUTP,1)
!!$    print *,"hellllllllllo",size(INPU,2),size(INPV,2),size(INPP,2),size(OUTU,2),size(OUTV,2),size(OUTP,2)
!!$    print *,"is_div ",is_div
!!$    call flush(6)

    lmu=size(INPU,1)-2; nmu=size(INPU,2)-2
    lmv=size(INPV,1)-2; nmv=size(INPV,2)-2
    lmp=size(INPP,1)-2; nmp=size(INPP,2)-2
    coef=2._prec*invdx2+2._prec*invdy2

    if (is_u) OUTU_c => OUTU(1:lmu,1:nmu);           if (is_v) OUTV_c => OUTV(1:lmv,1:nmv);    
    if (is_u) INPU_c => INPU(1:lmu,1:nmu);           if (is_v) INPV_c => INPV(1:lmv,1:nmv);    
    if (is_u) INPU_n => INPU(1:lmu,2:nmu+1);         if (is_v) INPV_n => INPV(1:lmv,2:nmv+1);  
    if (is_u) INPU_e => INPU(2:lmu+1,1:nmu);         if (is_v) INPV_e => INPV(2:lmv+1,1:nmv);  
    if (is_u) INPU_w => INPU(0:lmu-1,1:nmu);         if (is_v) INPV_w => INPV(0:lmv-1,1:nmv);  
    if (is_u) INPU_s => INPU(1:lmu,0:nmu-1);         if (is_v) INPV_s => INPV(1:lmv,0:nmv-1);  

    if (.not.flag) then     ! VF2 
       TMPU_u => TMPU(1:lmu+1,:);            TMPV_u => TMPV(1:lmv+1,:);    
       TMPU_v => TMPU(:,1:nmu+1);            TMPV_v => TMPV(:,1:nmv+1);    
       TMPU_n => TMPU(1:lmu,2:nmu+1);        TMPV_n => TMPV(1:lmv,2:nmv+1);
       TMPU_e => TMPU(2:lmu+1,1:nmu);        TMPV_e => TMPV(2:lmv+1,1:nmv);
       TMPU_w => TMPU(1:lmu,1:nmu);          TMPV_w => TMPV(1:lmv,1:nmv);  
       TMPU_s => TMPU(1:lmu,1:nmu);          TMPV_s => TMPV(1:lmv,1:nmv);  
    endif

!!$    print *,"ici"
!!$    call flush(6)

    if (is_grdp.or.is_div.or.is_grdv) then                          ! start_out_light
       lmp=size(INPP,1)-2; nmp=size(INPP,2)-2
       ip0=1; kp0=1; ip1=lmp; kp1=nmp
       if (is_north) kp1=kp1+1
       if (is_south) kp0=kp0-1
       if (is_east ) ip1=ip1+1
       if (is_west ) ip0=ip0-1
       !       print *,'lmp+2,nmp+2',lmp+2,nmp+2

       OUTP_c => OUTP(ip0:ip1,kp0:kp1)

       if (is_decal) then                                           
          INPP_n => INPP(1:lmv,1:nmv)  ;          OUTP_n => OUTP(1:lmv,1:nmv)  
          INPP_e => INPP(1:lmu,1:nmu)  ;          OUTP_e => OUTP(1:lmu,1:nmu)  
          INPP_w => INPP(0:lmu-1,1:nmu);          OUTP_w => OUTP(0:lmu-1,1:nmu)
          INPP_s => INPP(1:lmv,0:nmv-1);          OUTP_s => OUTP(1:lmv,0:nmv-1)
 
          INPUp_e => INPU(ip0+1:ip1+1,kp0:kp1);  INPVp_n => INPV(ip0:ip1,kp0+1:kp1+1)
          INPUp_w => INPU(ip0:ip1    ,kp0:kp1);   INPVp_s => INPV(ip0:ip1,kp0:kp1)
       else                                                         
          INPU_ne=> INPU(ip0+1:ip1+1,kp0+1:kp1+1);INPV_ne=> INPV(ip0+1:ip1+1,kp0+1:kp1+1); 
          INPU_se=> INPU(ip0+1:ip1+1,kp0:kp1);    INPV_se=> INPV(ip0+1:ip1+1,kp0:kp1);     
          INPU_nw=> INPU(ip0:ip1,kp0+1:kp1+1);    INPV_nw=> INPV(ip0:ip1,kp0+1:kp1+1);     
          INPU_sw=> INPU(ip0:ip1,kp0:kp1);        INPV_sw=> INPV(ip0:ip1,kp0:kp1);         

          INPP_ne=> INPP(1:lmu,1:nmu)    ; OUTP_ne=> OUTP(1:lmu,1:nmu)        
          INPP_se=> INPP(1:lmu,0:nmu-1)  ; OUTP_se=> OUTP(1:lmu,0:nmu-1)      
          INPP_nw=> INPP(0:lmu-1,1:nmu)  ; OUTP_nw=> OUTP(0:lmu-1,1:nmu)      
          INPP_sw=> INPP(0:lmu-1,0:nmu-1); OUTP_sw=> OUTP(0:lmu-1,0:nmu-1)    
       endif
    endif                                                                      ! end_out_light

    if (is_mpp) then
       if (is_u.or.(is_div.or.is_grdv)) call rfr_stn(INPU,stn_news,1,1)
       if (is_v.or.(is_div.or.is_grdv)) call rfr_stn(INPV,stn_news,1,1)
       if (is_grdp) call rfr_stn(INPP,stn_news,1,1)                            ! out_light
    endif

    if (is_u) then  
       if (is_t) then
          OUTU=coeft*INPU
       else
          OUTU=0._prec
       endif

       if (is_diff) OUTU_c = OUTU_c - nu*(-coef*INPU_c+invdx2*(INPU_e+INPU_w)&
            &                                         +invdy2*(INPU_n+INPU_s))

       if (is_conv) then
          if (flag) then
             OUTU_c = OUTU_c + p5invdx*U_EN_U_c*(INPU_e-INPU_w)&
                  &           +p5invdy*V_EN_U_c*(INPU_n-INPU_s)
          else
             TMPU_u=U_EN_U(1:lmu+1,:)*(INPU(0:lmu,:)+INPU(1:lmu+1,:)) 
             OUTU_c=OUTU_c+p5invdx*(TMPU_e-TMPU_w)
             TMPU_v=V_EN_U(:,1:nmu+1)*(INPU(:,0:nmu)+INPU(:,1:nmu+1))
             OUTU_c=OUTU_c+p5invdy*(TMPU_n-TMPU_s)
          endif
       endif
    endif


    if (is_v) then
       if (is_t)  then
          OUTV=coeft*INPV
       else
          OUTV=0._prec
       endif

       if (is_diff) OUTV_c = OUTV_c - nu*(-coef*INPV_c+invdx2*(INPV_e+INPV_w)&
            &                                         +invdy2*(INPV_n+INPV_s))

       if (is_conv) then
          if (flag) then
             OUTV_c = OUTV_c + p5invdx*U_EN_V_c*(INPV_e-INPV_w)&
                  &          + p5invdy*V_EN_V_c*(INPV_n-INPV_s)
          else
             TMPV_u=U_EN_V(1:lmv+1,:)*(INPV(0:lmv,:)+INPV(1:lmv+1,:)) 
             OUTV_c=OUTV_c+p5invdx*(TMPV_e-TMPV_w)
             TMPV_v=V_EN_V(:,1:nmv+1)*(INPV(:,0:nmv)+INPV(:,1:nmv+1))
             OUTV_c=OUTV_c+p5invdy*(TMPV_n-TMPV_s)
          endif
       endif
    endif

    if (is_decal) then                                           ! start_out_light
       if (is_grdp.and.is_u)  OUTU_c = OUTU_c + invdx*(INPP_e-INPP_w)

       if (is_grdp.and.is_v)  OUTV_c = OUTV_c + invdy*(INPP_n-INPP_s)

       if (is_div.or.is_grdv) OUTP_c = invdx*(INPUp_e-INPUp_w) + invdy*(INPVp_n-INPVp_s)

       if (is_grdv) then
          OUTU_c = OUTU_c - rlag*invdx*(OUTP_e-OUTP_w)
          OUTV_c = OUTV_c - rlag*invdy*(OUTP_n-OUTP_s)
          OUTP   = 0._prec
       endif
    else          ! maillage non decalle
       if (is_grdp.and.is_u)  OUTU_c = OUTU_c + p5invdx*(INPP_ne+INPP_se-INPP_nw-INPP_sw)

       if (is_grdp.and.is_u)  OUTV_c = OUTV_c + p5invdy*(INPP_ne+INPP_nw-INPP_se-INPP_sw)

       if (is_div.or.is_grdv) OUTP_c =    p5invdx*(INPU_ne+INPU_se-INPU_nw-INPU_sw) &
                                        + p5invdy*(INPV_ne+INPV_nw-INPV_se-INPV_sw)

       if (is_grdv) then
          if (is_grdp) call rfr_stn(OUTP,stn_news,1,1)
          OUTU_c = OUTU_c - rlag*p5invdx*(OUTP_ne+OUTP_se-OUTP_nw-OUTP_sw)
          OUTV_c = OUTV_c - rlag*p5invdy*(OUTP_ne+OUTP_nw-OUTP_se-OUTP_sw)
          OUTP   = 0._prec
       endif
    endif                                              ! end_out_light

!!$    print *,"is_grdp,is_grdv,is_div ",is_grdp,is_grdv,is_div 
!!$    return


    call fixe_cl(INPU,INPV,INPP,OUTU,OUTV,OUTP)

    return
  end subroutine mavecxy_df2

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine mavecxy_df4 (INPU,INPV,INPP,OUTU,OUTV,OUTP)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use compact
    use drapeaux
    use disc
    use data,   only : nu, rlag
    implicit none
    real(kind=prec), dimension(0:,0:), target :: OUTU,OUTV,OUTP,INPU,INPV,INPP
    real(kind=prec), dimension(0:size(INPV,1),0:size(INPV,2)), target   :: TMPdiv
    real(kind=prec), dimension(0:size(INPP,1),0:size(INPP,2)), target   :: TMPP
    real(kind=prec), dimension(0:size(INPU,1),0:size(INPU,2)-1), target   :: TMPu
    real(kind=prec), dimension(0:size(INPV,1)-1,0:size(INPV,2)), target   :: TMPv
    integer :: lmu,nmu,lmv,nmv,lmp,nmp,ip,kp
    real(kind=prec), dimension(:,:),pointer :: &
             OUTU_c,  OUTV_c, TMPP_C, TMPdiv_c, TMPu_c, TMPv_c&
          & ,TMPP_u,TMPP_v

    lmu=size(INPU,1)-2;  nmu=size(INPU,2)-2
    lmv=size(INPV,1)-2;  nmv=size(INPV,2)-2

    if (is_u) OUTU_c => OUTU(1:lmu,1:nmu);     
    if (is_v) OUTV_c => OUTV(1:lmv,1:nmv);   
    if (is_grdp.or.is_div.or.is_grdv) then
       TMPP_c   => TMPP(1:lmu,1:nmv)
       ip=size(INPV,1)-2;  if (.not.is_east) ip=ip+1;  
       kp=size(INPU,2)-2;  if (.not.is_north) kp=kp+1; 
       TMPdiv_c => TMPdiv(1:ip+1,1:kp+1)
       if (is_div) then
          TMPu_c   => TMPu(1:size(OUTP,1),:)
          TMPv_c   => TMPv(:,1:size(OUTP,2))
          TMPP_u   => TMPP(1:lmu,1:nmu)
          TMPP_v   => TMPP(1:lmv,1:nmv)
       endif
    endif

    if (is_u) then    
       if (is_t) then;     OUTU=coeft*INPU; else; OUTU=0._prec; endif;

       if (is_diff) OUTU = OUTU - nu*(der2x(INPU,lmu_global)+der2y(INPU,nmu_global))

       if (is_conv) OUTU = OUTU + U_EN_U*der1x(INPU,lmu_global) + V_EN_U*der1y(INPU,nmu_global)
    endif

    if (is_v) then    
       if (is_t) then;     OUTV=coeft*INPV; else; OUTV=0._prec; endif

       if (is_diff) OUTV = OUTV - nu*(der2x(INPV,lmv_global)+der2y(INPV,nmv_global))

       if (is_conv) OUTV = OUTV + U_EN_V*der1x(INPV,lmv_global) + V_EN_V*der1y(INPV,nmv_global)
    endif

    if (is_decal) then                                                     
       if (is_grdp.and.is_u) then
          lmp=size(INPP,1)-2
          TMPP(0:lmp+2,0:nmu+1)   = derm1x(INPP,lmp_global)
          OUTU_c = OUTU_c+TMPP_u
       endif

       if (is_grdp.and.is_v) then
          nmp=size(INPP,2)-2; 
          TMPP(0:lmv+1,0:nmp+2) = derm1y(INPP,nmp_global)
          OUTV_c = OUTV_c+TMPP_v
       endif
    else
       if (is_grdp.and.is_u) then
          TMPP=derm0y(derm1x(INPP,lmp_global),nmp_global); 
          OUTU_c = OUTU_c + TMPP_c
       endif
       
       if (is_grdp.and.is_v) then
          TMPP=derm0x(derm1y(INPP,nmp_global),lmp_global);  
          OUTV_c = OUTV_c + TMPP_c
       endif
       
    endif
    
    if (is_decal) then
       if (is_div.or.is_grdv) then
          TMPu = derm1x(INPU,lmu_global); TMPv=derm1y(INPV,nmv_global)
          OUTP   = TMPu_c+TMPv_c
       endif

       if (is_grdv) then
          OUTU_c = OUTU_c - rlag*derm1x(OUTP,lmp_global)
          OUTV_c = OUTV_c - rlag*derm1y(OUTP,nmp_global)
          OUTP   = 0._prec
       endif
    else                        
       if (is_div.or.is_grdv) then
          TMPdiv = derm0y(derm1x(INPU,lmu_global),nmu_global) 
          TMPdiv = TMPdiv + derm0x(derm1y(INPV,nmv_global),lmv_global);
          OUTP   = TMPdiv_c
       endif
       
       if (is_grdv) then
          call output('ici,lmp_global,nmp_global',lmp_global,nmp_global)        ! out_small
          TMPP=derm0y(derm1x(OUTP,lmp_global),nmp_global);  OUTU_c = OUTU_c - rlag*TMPP_c
          TMPP=derm0x(derm1y(OUTP,nmp_global),lmp_global);  OUTV_c = OUTV_c - rlag*TMPP_c
          OUTP   = 0._prec
       endif
    endif

    call fixe_cl(INPU,INPV,INPP,OUTU,OUTV,OUTP)


!!$    call outdon(INPU,'INPU',1,1)
!!$    call outdon(INPV,'INPV',1,1)

    return
  end subroutine mavecxy_df4

  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss             ! sdtart_tiny       
  subroutine mavecxy_vf4 (INPU,INPV,INPP,OUTU,OUTV,OUTP)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use compact
    use drapeaux
    use disc
    use data,   only : nu, rlag
    implicit none
    real(kind=prec), dimension(0:,0:),target :: OUTU,OUTV,OUTP
    real(kind=prec), dimension(0:,0:),target ::INPU,INPV,INPP

    real(kind=prec), dimension(0:size(INPU,1),0:size(INPU,2)), target :: TMPU
    real(kind=prec), dimension(0:size(INPV,1),0:size(INPV,2)), target :: TMPV
    real(kind=prec), dimension(:,:), pointer ::                &
          TMPUx,TMPUy, TMPVx,TMPVy, TMPU_c, TMPV_c              &
         ,TMPU_e, TMPU_w, TMPV_n, TMPV_s, OUTU_c,OUTV_c         &
         ,TMPUx_c,TMPUx_w,TMPUx_e,TMPUy_c,TMPUy_s,TMPUy_n       &
         ,TMPVx_c,TMPVx_w,TMPVx_e,TMPVy_c,TMPVy_s,TMPVy_n       &
         ,INPP_n,INPP_e,INPP_w,INPP_s                           &
         ,OUTP_n,OUTP_e,OUTP_w,OUTP_s
    real(kind=prec), dimension(:,:), pointer :: TMPdu_c,TMPdv_c
    real(kind=prec), dimension(0:size(INPU,1),0:size(INPU,2)-1), target   :: TMPdu
    real(kind=prec), dimension(0:size(INPV,1)-1,0:size(INPV,2)), target   :: TMPdv

    integer :: lmu,nmu,lmv,nmv,lmp,nmp

    lmu=size(INPU,1)-2;  nmu=size(INPU,2)-2
    lmv=size(INPV,1)-2;  nmv=size(INPV,2)-2
    

    TMPUx   => TMPU(0:lmu+2,0:nmu+1);      TMPVx   => TMPV(0:lmv+2,0:nmv+1);
    TMPUx_e => TMPU(1:lmu+2,0:nmu+1);      TMPVx_e => TMPV(1:lmv+2,0:nmv+1);
    TMPUx_w => TMPU(0:lmu+1,0:nmu+1);      TMPVx_w => TMPV(0:lmv+1,0:nmv+1);
    TMPUx_c => TMPU(0:lmu+1,0:nmu+1);      TMPVx_c => TMPV(0:lmv+1,0:nmv+1);
    TMPUy   => TMPU(0:lmu+1,0:nmu+2);      TMPVy   => TMPV(0:lmv+1,0:nmv+2);
    TMPUy_n => TMPU(0:lmu+1,1:nmu+2);      TMPVy_n => TMPV(0:lmv+1,1:nmv+2);
    TMPUy_s => TMPU(0:lmu+1,0:nmu+1);      TMPVy_s => TMPV(0:lmv+1,0:nmv+1);
    TMPUy_c => TMPU(0:lmu+1,0:nmu+1);      TMPVy_c => TMPV(0:lmv+1,0:nmv+1);

    if (is_div.or.is_grdp.or.is_grdv) then
       lmp=size(INPP,1)-1; 
       nmp=size(INPP,2)-1; 
       if (is_decal) then
          INPP_w  => INPP(0:lmu-1,0:nmu+1); INPP_e  => INPP(1:lmu,0:nmu+1);
          INPP_s  => INPP(0:lmv+1,0:nmv-1); INPP_n  => INPP(0:lmv+1,1:nmv); 
          OUTP_w  => OUTP(0:lmu-1,0:nmu+1); OUTP_e  => OUTP(1:lmu,0:nmu+1); 
          OUTP_s  => OUTP(0:lmv+1,0:nmv-1); OUTP_n  => OUTP(0:lmv+1,1:nmv); 

          TMPU_e  => TMPU(1:lmp+1,0:nmp  );   TMPV_s  => TMPV(0:lmp  ,0:nmp  )
          TMPU_w  => TMPU(0:lmp  ,0:nmp  );   TMPV_n  => TMPV(0:lmp  ,1:nmp+1);
          OUTU_c  => OUTU(1:lmu,:);           OUTV_c  => OUTV(:  ,1:nmv);
          if (is_div) then
             TMPdu_c   => TMPdu(1:size(OUTP,1),:)
             TMPdv_c   => TMPdv(:,1:size(OUTP,2))
!!$             TMPP_u   => TMPP(1:lmu,1:nmu)
!!$             TMPP_v   => TMPP(1:lmv,1:nmv)
          endif
       else              
          INPP_w  => INPP(0:lmu-1,0:nmu); INPP_e  => INPP(1:lmu,0:nmu);
          INPP_s  => INPP(0:lmv,0:nmv-1); INPP_n  => INPP(0:lmv,1:nmv); 
          OUTP_w  => OUTP(0:lmu-1,0:nmu); OUTP_e  => OUTP(1:lmu,0:nmu);
          OUTP_s  => OUTP(0:lmv,0:nmv-1); OUTP_n  => OUTP(0:lmv,1:nmv); 


          TMPU_c  => TMPU(0:lmu+1,0:nmu);   TMPV_c  => TMPV(0:lmv,0:nmv+1);
          TMPU_e  => TMPU(1:lmp+1,0:nmp);   TMPV_s  => TMPV(0:lmp,0:nmp);
          TMPU_w  => TMPU(0:lmp  ,0:nmp);   TMPV_n  => TMPV(0:lmp,1:nmp+1);
          OUTU_c  => OUTU(1:lmu,1:nmu);     OUTV_c  => OUTV(1:lmv  ,1:nmv);
       endif             
     endif


    if (is_u) then    
       if (is_t) then;  OUTU=coeft*integ_vf(INPU); else; OUTU=0._prec; endif

       if (is_diff) then
          TMPUx=derm1x(INPU,lmu_global);
          TMPUx_c=TMPUx_e-TMPUx_w
          OUTU = OUTU - nu*invdx*integ_y(TMPUx_c)
          TMPUy=derm1y(INPU,nmu_global);
          TMPUy_c=TMPUy_n-TMPUy_s
          OUTU = OUTU - nu*invdy*integ_x(TMPUy_c)
       endif

       if (is_conv) then
          TMPUx=U_EN_U*derm0x(INPU,lmu_global);
          TMPUx_c=TMPUx_e-TMPUx_w
          OUTU = OUTU + invdx*integ_y(TMPUx_c)
          TMPUy=V_EN_U*derm0y(INPU,nmu_global);
          TMPUy_c=TMPUy_n-TMPUy_s
          OUTU = OUTU + invdy*integ_x(TMPUy_c)
       endif
    endif

    if (is_v) then    

       if (is_t) then; OUTV=coeft*integ_vf(INPV); else; OUTV=0._prec; endif

       if (is_diff) then
          TMPVx=derm1x(INPV,lmv_global);
          TMPVx_c=TMPVx_e-TMPVx_w
          OUTV = OUTV -nu*invdx*integ_y(TMPVx_c)
          TMPVy=derm1y(INPV,nmv_global);
          TMPVy_c=TMPVy_n-TMPVy_s
          OUTV = OUTV -nu*invdy*integ_x(TMPVy_c)
       endif

       if (is_conv) then
          TMPVx=U_EN_V*derm0x(INPV,lmv_global)
          TMPVx_c=TMPVx_e-TMPVx_w
          OUTV = OUTV + invdx*integ_y(TMPVx_c)
          TMPVy=V_EN_V*derm0y(INPV,nmv_global)
          TMPVy_c=TMPVy_n-TMPVy_s
          OUTV = OUTV + invdy*integ_x(TMPVy_c)
       endif
    endif

    if (is_grdp)  call rfr_stn(INPP,stn_news,1,1)

    if (is_decal) then
       if (is_grdp.and.is_u)  OUTU_c = OUTU_c + invdx*intega_y(INPP_e-INPP_w)

       if (is_grdp.and.is_v)  OUTV_c = OUTV_c + invdy*intega_x(INPP_n-INPP_s)

!!$       if (is_div.or.is_grdv) then
!!$          TMPdu = derm1x(INPU,lmu_global); TMPdv=derm1y(INPV,nmv_global)
!!$          OUTP   = TMPdu_c+TMPdv_c
!!$       endif
       if (is_div.or.is_grdv) then
             TMPU(0:lmu+1,0:nmu+1) = intega_y(INPU);   
             TMPV(0:lmv+1,0:nmv+1) = intega_x(INPV); 
             call rfr_stn(TMPU(0:lmu+1,0:nmu+1),stn_ew,1,1)
             call rfr_stn(TMPV(0:lmv+1,0:nmv+1),stn_ns,1,1)
             OUTP  = invdx*(TMPU_e-TMPU_w)+invdy*(TMPV_n-TMPV_s)
       endif
       
       if (is_grdv) then
          OUTU_c = OUTU_c - rlag*invdx*integd_y(OUTP_e-OUTP_w)
          OUTV_c = OUTV_c - rlag*invdy*integd_x(OUTP_n-OUTP_s)
          OUTP   = 0._prec
       endif
    else
       if (is_grdp.and.is_u) OUTU_c = OUTU_c + invdx*integd_y(INPP_e-INPP_w)

       if (is_grdp.and.is_v) OUTV_c = OUTV_c + invdy*integd_x(INPP_n-INPP_s)

       if (is_div.or.is_grdv) then
             TMPU_c = integd_y(INPU);   
             TMPV_c = integd_x(INPV); 
             OUTP  = invdx*(TMPU_e-TMPU_w)+invdy*(TMPV_n-TMPV_s)
       endif
       
       if (is_grdv) then
          OUTU_c = OUTU_c - rlag*invdx*integd_y(OUTP_e-OUTP_w)
          OUTV_c = OUTV_c - rlag*invdy*integd_x(OUTP_n-OUTP_s)
          OUTP   = 0._prec
       endif
       
    endif

    call fixe_cl(INPU,INPV,INPP,OUTU,OUTV,OUTP)

    return
  end subroutine mavecxy_vf4
!end_out_light ednd_tiny

  !************************************************************************

  !SKffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
  function  global_ddot(XU,XV,XP,YU,YV,YP) result(r)
    !SKffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    !SK
    !SK   On calcule le produit scalaire des tableaux x,y.
    !SK      d'abord on calcule la contribution a ce produit de chaque proc
    !SK      puis on somme le tout!
    !SK
    !SK================================================================
    !SK
    use para
    use data
    use drapeaux
    use uncol 
    implicit none

    real(kind=prec), dimension(0:,0:) :: XU,XV,XP,YU,YV,YP
    real(kind=prec) :: r
!!$    integer :: lmp,nmp

    r=0._prec

    if (is_u)   r=r+sum(XU*YU)
    if (is_v)   r=r+sum(XV*YV)
    if (is_div) r=r+sum(XP*YP)         ! start_out_light
    if (is_mpp) r=global_add(r)

    return
  end function global_ddot

  !************************************************************************
  
  
  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine fixe_cl(INPU,INPV,INPP,OUTU,OUTV,OUTP,flag)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use para
    use disc,     only : invdx,invdy,p5invdx,p5invdy 
    use drapeaux, only : is_oned,is_u,is_v,is_div,is_decal
    use data,     only : ncheck
    implicit none
    real(kind=prec), dimension(0:,0:) :: INPU,INPV,INPP,OUTU,OUTV,OUTP
    logical, optional :: flag
    integer :: lmu,nmu,lmv,nmv,lmp,nmp

    lmu = size(INPU,1)-2; nmu = size(INPU,2)-2
    lmv = size(INPV,1)-2; nmv = size(INPV,2)-2
    lmp = size(INPP,1)-2; nmp = size(INPP,2)-2
    
    if (is_mpp.and.is_u) then 
       if (.not.is_north) INPU(:,nmu+1)=0._prec;  
       if (.not.is_east ) INPU(lmu+1,:)=0._prec;  
       if (.not.is_west ) INPU(0,:)    =0._prec;  
       if (.not.is_south) INPU(:,0)    =0._prec;  
    endif

    if (is_mpp.and.is_v) then 
       if (.not.is_north) INPV(:,nmv+1)=0._prec;
       if (.not.is_east ) INPV(lmv+1,:)=0._prec;
       if (.not.is_west ) INPV(0,:)    =0._prec;
       if (.not.is_south) INPV(:,0)    =0._prec;
    endif

    if (is_mpp.and.is_div) then                              ! start_out_light
       if (.not.is_north) OUTP(:,nmp+1)=0._prec;
       if (.not.is_east ) OUTP(lmp+1,:)=0._prec;
       if (.not.is_west ) OUTP(0,:)    =0._prec;
       if (.not.is_south) OUTP(:,0)    =0._prec;
    endif                                        ! end_out_light

    if (present(flag).or..not.is_decal) then     ! out_light
       if (is_u) then
          OUTU(0,:)    = INPU(0,:); 
          OUTU(lmu+1,:)= INPU(lmu+1,:)  
          if (is_oned) then
             OUTU(:,nmu+1)= OUTU(:,nmu)
             OUTU(:,0)= OUTU(:,1)
             return
          endif
          OUTU(:,nmu+1)= INPU(:,nmu+1);  
          OUTU(:,0)    = INPU(:,0)    ;  
       endif
       if (is_v) then
          OUTV(:,nmv+1)= INPV(:,nmv+1);  
          OUTV(lmv+1,:)= INPV(lmv+1,:);  
          OUTV(0,:)    = INPV(0,:)    ;  
          OUTV(:,0)    = INPV(:,0)    ;  
       endif
    else                                        ! start_out_light
       select case(ncheck)   
       case(vf2,df2)         
       if (is_u) then
          if (is_north) then; OUTU(1:lmu,nmu+1)=INPU(1:lmu,nmu+1);                    
                        else; OUTU(1:lmu,nmu+1)= 0._prec; endif
          if (is_east ) then; OUTU(lmu+1,:)    =0.5_prec*(INPU(lmu,:)+INPU(lmu+1,:))
                        else; OUTU(lmu+1,:)    = 0._prec; endif
          if (is_west ) then; OUTU(0,:)        =0.5_prec*(INPU(0,:)+INPU(1,:))
                        else; OUTU(0,:)        = 0._prec; endif
          if (is_south) then; OUTU(1:lmu,0)    =INPU(1:lmu,0);   
                        else; OUTU(1:lmu,0)    = 0._prec; endif

       endif            
       if (is_v) then   
          if (is_north) then; OUTV(:,nmv+1)      = 0.5_prec*(INPV(:,nmv)+INPV(:,nmv+1))  
                        else; OUTV(:,nmv+1)      = 0._prec; endif;
          if (is_east ) then; OUTV(lmv+1,1:nmv)  = INPV(lmv+1,1:nmv)                      
                        else; OUTV(lmv+1,1:nmv)  = 0._prec; endif
          if (is_west ) then; OUTV(0,1:nmv)      = INPV(0,1:nmv)                          
                        else; OUTV(0,1:nmv)      = 0._prec; endif;
          if (is_south) then; OUTV(:,0)          = 0.5_prec*(INPV(:,0)+INPV(:,1))        
                        else; OUTV(:,0)          = 0._prec; endif;
       endif            
       case(df4,df4e,vf4) 
          if (is_u) then
             if (.not.is_mpp.or.is_north) then; OUTU(1:lmu,nmu+1) = INPU(1:lmu,nmu+1)
                                          else; OUTU(1:lmu,nmu+1) = 0.; endif
             if (.not.is_mpp.or.is_east ) then
                OUTU(lmu+1,:)     = (0.3125_prec*INPU(lmu+1,:)+0.9375_prec*INPU(lmu,:)      &
                     -0.3125_prec*INPU(lmu-1,:)+0.0625_prec*INPU(lmu-2,:))
                else; OUTU(lmu+1,:) = 0.;  endif
             if (.not.is_mpp.or.is_west ) then
                OUTU(0,:)         = (0.3125_prec*INPU(0,:)    +0.9375_prec*INPU(1,:)        &
                                                 -0.3125_prec*INPU(2,:)    +0.0625_prec*INPU(3,:)    )
                else; OUTU(0,:)=0._prec; endif
             if (.not.is_mpp.or.is_south) then;  OUTU(1:lmu,0)     = INPU(1:lmu,0)    
                                          else;  OUTU(1:lmu,0)     = 0.; endif
          endif
          if (is_v) then
             if (.not.is_mpp.or.is_north) then;
                OUTV(:,nmv+1)     = (0.3125_prec*INPV(:,nmv+1)+0.9375_prec*INPV(:,nmv)      &
                                                 -0.3125_prec*INPV(:,nmv-1)+0.0625_prec*INPV(:,nmv-2))
                else; OUTV(:,nmv+1)=0.; endif
             if (.not.is_mpp.or.is_east )   then; OUTV(lmv+1,1:nmv) = INPV(lmv+1,1:nmv)
                                            else; OUTV(lmv+1,1:nmv) = 0.; endif
             if (.not.is_mpp.or.is_west )   then; OUTV(0,1:nmv)     = INPV(0,1:nmv)    
                                            else; OUTV(0,1:nmv)     = 0.; endif
             if (.not.is_mpp.or.is_south)   then
                OUTV(:,0)         = (0.3125_prec*INPV(:,0)    +0.9375_prec*INPV(:,1)        &
                                                 -0.3125_prec*INPV(:,2)    +0.0625_prec*INPV(:,3)    )
                else; OUTV(:,0)=0.; endif
          endif
       end select                           ! end_out_light

       if (is_div.and.is_decal) then
          select case(ncheck)   
          case (df2,vf2)
             if (is_south.and.is_west) OUTP(0,0)=0.5*(2.*OUTP(0,1)-OUTP(0,2)+2.*OUTP(1,0)-OUTP(2,0))
             if (is_north.and.is_west) OUTP(0,nmp+1)=0.5*(2.*OUTP(1,nmp+1)-OUTP(2,nmp+1)&
                  & +2.*OUTP(0,nmp)-OUTP(0,nmp-1))
             if (is_north.and.is_east) OUTP(lmp+1,nmp+1)=0.5*(2.*OUTP(lmp,nmp+1)-OUTP(lmp-1,nmp+1)&
                  & + 2.*OUTP(lmp+1,nmp)-OUTP(lmp+1,nmp-1))
             if (is_south.and.is_east) OUTP(lmp+1,0)=0.5*(2.*OUTP(lmp,0)-OUTP(lmp-1,0)&
                  & + 2.*OUTP(lmp+1,1)-OUTP(lmp+1,2))
          case(df4,vf4)
             if (is_south.and.is_west) OUTP(0,0)=0.5*(&
                  & 4.*OUTP(0,1)-6.*OUTP(0,2)+4.*OUTP(0,3)-OUTP(0,4)&
                  &+4.*OUTP(1,0)-6.*OUTP(2,0)+4.*OUTP(3,0)-OUTP(4,0))
             if (is_north.and.is_west) OUTP(0,nmp+1)=0.5*(&
                  & +4.*OUTP(1,nmp+1)-6.*OUTP(2,nmp+1)+4.*OUTP(3,nmp+1)-OUTP(4,nmp+1)&
                  & +4.*OUTP(0,nmp)-6.*OUTP(0,nmp-1)+4.*OUTP(0,nmp-2)-OUTP(0,nmp-3))
             if (is_north.and.is_east) OUTP(lmp+1,nmp+1)=0.5*(&
                  & 4.*OUTP(lmp,nmp+1)-6.*OUTP(lmp-1,nmp+1)+4.*OUTP(lmp-2,nmp+1)-OUTP(lmp-3,nmp+1)&
                  &+4.*OUTP(lmp+1,nmp)-6.*OUTP(lmp+1,nmp-1)+4.*OUTP(lmp+1,nmp-2)-OUTP(lmp+1,nmp-3))
             if (is_south.and.is_east) OUTP(lmp+1,0)=0.5*(&
                  & 4.*OUTP(lmp,0)-6.*OUTP(lmp-1,0)+4.*OUTP(lmp-2,0)-OUTP(lmp-3,0)&
                  &+4.*OUTP(lmp+1,1)-6.*OUTP(lmp+1,2)+4.*OUTP(lmp+1,3)-OUTP(lmp+1,4))
!          case (vf4)
          case default
             stop 'fixe_cl a amender'
          end select
!!$          if (is_south.and.is_west) OUTP(0,0)    =p5invdx*(INPU(1,0)+INPU(0,0))+p5invdy*(INPV(0,1)+INPV(0,0))
!!$          if (is_north.and.is_west) OUTP(0,nmp+1)=p5invdx*(INPU(1,nmp+1)+INPU(0,nmp+1))+p5invdy*(INPV(0,nmp+2)+INPV(0,nmp+1))
!!$          if (is_north.and.is_east) OUTP(lmp+1,nmp+1)=p5invdx*(INPU(lmp+2,nmp+1)+INPU(lmp+1,nmp+1))&
!!$               & +p5invdy*(INPV(lmp+1,nmp+2)+INPV(lmp+1,nmp+1))
!!$          if (is_south.and.is_east) OUTP(lmp+1,0)=p5invdx*(INPU(lmp+2,0)+INPU(lmp+1,0))+p5invdy*(INPV(lmp+1,1)+INPV(lmp+1,0))
!!$             if (is_south) OUTP(:,0)=0.
!!$             if (is_north) OUTP(:,nmp+1)=0.
!!$             if (is_east)  OUTP(lmp+1,:)=0.
!!$             if (is_west)  OUTP(0,:)=0.
!!$             if (is_south.and.is_west) OUTP(0,0)=0.
!!$             if (is_north.and.is_west) OUTP(0,nmp+1)=0.
!!$             if (is_north.and.is_east) OUTP(lmp+1,nmp+1)=0.
!!$             if (is_south.and.is_east) OUTP(lmp+1,0)=0.
!!$             if (is_south.and.is_west) INPP(0,0)=0.
!!$             if (is_north.and.is_west) INPP(0,nmp+1)=0.
!!$             if (is_north.and.is_east) INPP(lmp+1,nmp+1)=0.
!!$             if (is_south.and.is_east) INPP(lmp+1,0)=0.
!!$             if (is_south.and.is_west) OUTP(0,0)=INPP(0,0)
       endif


    endif

    return
  end subroutine fixe_cl



  !************************************************************************

end module blas



