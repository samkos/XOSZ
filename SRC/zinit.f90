
module init

  use constante
  use para
  use, intrinsic :: iso_fortran_env , only :  iostat_end
  
  integer, save :: ncheck_first,ncheck_second

  character(len=*), parameter, dimension(1:18) ::       &
       nom_solver_short = (/ 'BiCGSTAB    ' ,'BiCGSTAB_2  ' ,'BiCGSTAB_3  ' &
       ,'BiCGSTAB_4  ' ,'BiCGSTAB_5  ' ,'BiCGSTAB_6  ' ,'BiCGSTAB_7  ' &
       ,'BiCGSTAB_8  ' ,'BiCGSTAB_9  ' ,'GMRES       ' ,'GMRESR_rec  ' &
       ,'GMRESR_norec' ,'Multigrille ' ,'ORTHODIR    ' ,'ORTHOMIN    ' &
       ,'Lisseur     ' ,'Schwarz_add ' ,'Schwarz_mult'/)
contains

  !SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !SK  Routines d'initialisation du calcul
  !SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coinit
  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    use disc
    use data
    use uncol
    use champs
    use drapeaux
    use para
    use compact, only :  initialise_general_coef, initialise_df4_coef,     initialise_vf4_coef 
    use dac
    use plot_flags
    implicit none

    integer :: lmu,lmv,lmp,nmu,nmv,nmp,ok

    it_start=0
    call lecture_input

    ncheck_first=ncheck
    ncheck_second=ncheck_precond

    is_decal         = (is_decale/=0)
    is_output        = (is_print<0)
    is_print         = abs(is_print)
    is_cv_file       = is_cv/=0.and.mod(abs(is_cv),5)==0
    is_cv_global     = is_cv>0
    niv_lim          = abs(is_cv)
    if (mod(niv_lim,5)==0)  niv_lim=niv_lim/5

    nom_file_cv=trim(nom_fic_output)//'cv_'&
         & //trim(nom_solver_short(abs(ntype_solver)))//'.dat'
    call init_para

    if (is_decal) then                                       ! start_out_light
       lmp_global=lm_global+1; nmp_global=nm_global+1; 
       lmu_global=lm_global+2; nmu_global=nm_global+1; 
       lmv_global=lm_global+1; nmv_global=nm_global+2; 

       lmp=int(real(lmp_global-2)/real(nb_i_blocks)+0.99) ! on considere les pts a repartir hors bords  
       nmp=int(real(nmp_global-2)/real(nb_k_blocks)+0.99) ! qui servent soit aux conditions aux limites 
       if (is_east)  lmp=lmp_global-2-lmp*(nb_i_blocks-1) ! soit comme yone de duplication              
       if (is_north) nmp=nmp_global-2-nmp*(nb_k_blocks-1)
       lmv=lmp; nmv=nmp
       lmu=lmp; nmu=nmp
       if (is_east)  lmu=lmu+1
       if (is_north) nmv=nmv+1                         
    else                                                     ! end_out_light
       lmp_global=lm_global+1; nmp_global=nm_global+1; 
       lmu_global=lm_global+2; nmu_global=nm_global+2; 
       lmv_global=lm_global+2; nmv_global=nm_global+2; 

       lmu=int(real(lmu_global-2)/real(nb_i_blocks)+0.99) ! on considere les pts a
       nmu=int(real(nmu_global-2)/real(nb_k_blocks)+0.99) ! repartir hors bords  
       if (is_east)  lmu=lmu_global-2-lmu*(nb_i_blocks-1) ! servant soit aux cond
       if (is_north) nmu=nmu_global-2-nmu*(nb_k_blocks-1) ! aux limites soit comme
       lmv=lmu; nmv=nmu                                   ! zone de duplication              
       
       lmp=lmu; nmp=nmu
       if (is_east)  lmp=lmp-1
       if (is_north) nmp=nmp-1                         
    endif                                                   ! out_light
    
    call output('lm_global ,nm _global        ',lm_global+2 ,nm_global+2)
    call output('lmu_global,nmu_global,lmu,nmu',lmu_global,nmu_global,lmu+2,nmu+2)
    call output('lmv_global,nmv_global,lmv,nmv',lmv_global,nmv_global,lmv+2,nmv+2)
    call output('lmp_global,nmp_global,lmp,nmp',lmp_global,nmp_global,lmp+2,nmp+2)  ! out_light

    ok=min(min(lmu,nmu),min(lmv,nmv))
    ok=min(ok,min(lmp,nmp))
    ok=global_min(ok)
    if (nexample/=1.and.((ncheck==df4.or.ncheck==df4e.or.ncheck==vf4).and.ok<4)) &
         call parallel_stop ('(lm<4 or nm<4) while x-stencil sizes > 4 ')

    allocate(VTU (0:lmu+1,0:nmu+1)); allocate(VTU0(0:lmu+1,0:nmu+1)); 
    allocate(VTU1(0:lmu+1,0:nmu+1)); allocate(VTU2(0:lmu+1,0:nmu+1)); 
    allocate(VTU3(0:lmu+1,0:nmu+1)); allocate(VTU4(0:lmu+1,0:nmu+1)); 
    allocate(VTUS(0:lmu+1,0:nmu+1)); allocate(SMU (0:lmu+1,0:nmu+1),stat=ok); 
    if (ok/=0) stop 'VTU-SMU  : error alloc'
                                                                  
    allocate(VTV (0:lmv+1,0:nmv+1)); allocate(VTV0(0:lmv+1,0:nmv+1)); 
    allocate(VTV1(0:lmv+1,0:nmv+1)); allocate(VTV2(0:lmv+1,0:nmv+1)); 
    allocate(VTV3(0:lmv+1,0:nmv+1)); allocate(VTV4(0:lmv+1,0:nmv+1)); 
    allocate(VTVS(0:lmv+1,0:nmv+1)); allocate(SMV (0:lmv+1,0:nmv+1),stat=ok); 
    if (ok/=0) stop 'VTV-SMV  : error alloc'
                                                                  
    allocate(PRE (0:lmp+1,0:nmp+1)); allocate(PRE0(0:lmp+1,0:nmp+1)); 
    allocate(PRE1(0:lmp+1,0:nmp+1)); allocate(PRE2(0:lmp+1,0:nmp+1)); 
    allocate(PRE3(0:lmp+1,0:nmp+1)); allocate(PRES(0:lmp+1,0:nmp+1)); 
    allocate(SMP (0:lmp+1,0:nmp+1),stat=ok); 
    if (ok/=0) stop 'PRE-SMP  : error alloc'

    allocate(CLDU(0:lmu+1,0:nmu+1),stat=ok); if (ok/=0) stop 'CLDU : error alloc'
    allocate(CLDV(0:lmv+1,0:nmv+1),stat=ok); if (ok/=0) stop 'CLDV : error alloc'


    call output('size(VTU)',size(VTU,1),size(VTU,2))
    call output('size(VTV)',size(VTV,1),size(VTV,2))            
    call output('size(PRE)',size(PRE,1),size(PRE,2))            




    is_impl=diff_code

    select case(nexample)
    case(1)       ! solution test 1D
       global_switch = u_code*oned_code
       uv_switch = u_code
       is_ns=.false.
       is_unsteady=1
       vit_adv=rho_ca
       if (is_kuta>0) then 
          is_impl=is_impl*conv_code
       else
          is_impl=1
       endif
    case(20,22)       ! solution test
       global_switch = u_code
       uv_switch = u_code
       is_ns=.false.
       is_unsteady=0
       tau=1.E+15
       t_start=0._prec; t_all=tau
    case(16)       ! solution test diffusion + temps Bipole
       global_switch = u_code*v_code
       uv_switch = u_code*v_code
       is_ns=.false.
       is_unsteady=1
    case(18)       ! solution test diffusion Bipole
       global_switch = u_code*v_code
       uv_switch = u_code*v_code
       is_ns=.false.
       is_unsteady=0
       tau=1.E+15
       t_all=t_start+0.5*tau
    case(17)       ! solution test bipole advection+diffusion
       global_switch = u_code*v_code
       uv_switch = u_code*v_code
       is_ns=.false.
       is_unsteady=0
       tau=1.E+15
       t_all=t_start+0.5*tau
       if (is_kuta>0) is_impl=is_impl*conv_code
    case(19,21,23,24,25,26,27)  ! Navier-Stokes                           ! start_out_light
       global_switch = u_code*v_code
       uv_switch = u_code*v_code
       if (is_kuta>0) is_impl=is_impl*conv_code
       is_ns=.true.
       select case(nexample)
       case (23,24,27); is_unsteady=0
       case default;    is_unsteady=1
       end select
       select case (ns_solver)
       case(1); global_switch=global_switch                      ! Compressibilite Artificielle
       case(4); global_switch=global_switch*grdv_code            ! Lagrangien Augmente
       case(2); global_switch=global_switch*grdp_code*div_code   ! Resolution Couplee
       case default;  stop 'methode de projection non encore implementee'
       end select                                                ! end_out_light
    case(28,29,30)     ! burgers
       global_switch = u_code*v_code
       uv_switch = u_code*v_code
       if (is_kuta>0) is_impl=is_impl*conv_code
       is_ns = .false. 
       is_unsteady=1
    end select
    is_kuta=abs(is_kuta)

    call set_switch(global_switch*is_impl)
    call output ('global_switch,is_conv,is_diff,is_div'  &
         & ,i1=global_switch,l1=is_conv,l2=is_diff,l3=is_div)
    call output ('is_ns,is_grdp,is_grdv',l1=is_ns,l2=is_grdp,l3=is_grdv)

    un_ms_sor_theta = 1._prec - sor_theta                 ! out_light
    
    is_precond=(ntype_precond/=0)
    if (.not.is_precond) ncheck_precond=0
    is_schwarz_add = ntype_solver==17
    is_schwarz_mult = ntype_solver==18
    is_schwarz = is_schwarz_add.or.is_schwarz_mult

    nt_print=max(int(t_print/tau),1)

    call init_grids


    temps=t_start  
    if (is_unsteady>0.and.it_start==0) then
       temps=temps-3._prec*tau; call exacte; VTU3=VTUS; VTV3=VTVS; PRE3=PRES
       temps=temps+tau;         call exacte; VTU2=VTUS; VTV2=VTVS; PRE2=PRES
       temps=temps+tau;         call exacte; VTU1=VTUS; VTV1=VTVS; PRE1=PRES
       temps=temps+tau;         call exacte; VTU0=VTUS; VTV0=VTVS; PRE0=PRES
                                             VTU =VTUS; VTV =VTVS; PRE =PRES
    else
       VTU3=1._prec;  VTV3=1._prec; PRE3=1._prec
       VTU2=1._prec;  VTV2=1._prec; PRE2=1._prec
       VTU1=1._prec;  VTV1=1._prec; PRE1=1._prec
       VTU0=1._prec;  VTV0=1._prec; PRE0=1._prec
       VTU =1._prec;  VTV =1._prec; PRE =1._prec
    endif
    
    
    call initialise_general_coef                      ! start_out_light
    call initialise_df4_coef(VTU,VTV)
    call initialise_vf4_coef(VTU,VTV)
    if (is_precond.and.(ncheck_precond==df4.or.ncheck_precond==df4e&
         .or.ncheck_precond==vf4)) then
       call initialise_df4_coef(VTUp,VTVp)
       call initialise_vf4_coef(VTUp,VTVp)
    endif
    if (is_mpp) then
       call initialise_df4_part_coefx(lmu_global,lmu)
       call initialise_df4_part_coefy(nmv_global,nmv)
       call initialise_vf4_part_coefx(lmu_global,lmu)
       call initialise_vf4_part_coefy(nmv_global,nmv)
       if (is_ns.or.is_decal) then                       
          call initialise_df4_part_coefx(lmp_global,lmp)
          call initialise_df4_part_coefy(nmp_global,nmp)
          call initialise_vf4_part_coefx(lmp_global,lmp)
          call initialise_vf4_part_coefy(nmp_global,nmp)
       endif
    endif                                               ! end_out_light

!    if (is_restart_save>0)  call restore_champs           ! out_light 

    open(file=trim(nom_fic_output)//'Z.dat',unit=79,iostat=ok)
    write(79,*) controle_string
    write(79,*) 
    close(79)

    return
  end subroutine coinit

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine coend
  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    use drapeaux
    use disc
    use champs
    use data
    use compact, only :  release_general_coef, release_df4_coef, release_vf4_coef 
    use dac

    implicit none
    integer ok

    if (is_mpp) then
       call release_df4_part_coefx(lmu_global)
       call release_df4_part_coefy(nmv_global)
       call release_vf4_part_coefx(lmu_global)
       call release_vf4_part_coefy(nmv_global)
       if (is_ns.or.is_decal) then                       
          call release_df4_part_coefx(lmp_global)
          call release_df4_part_coefy(nmp_global)
          call release_vf4_part_coefx(lmp_global)
          call release_vf4_part_coefy(nmp_global)
       endif
    endif                                               ! end_out_light
 
                                                 ! start_out_light
    if (is_precond.and.(ncheck_precond==df4.or.ncheck_precond==df4e&
         .or.ncheck_precond==vf4)) then
       call release_df4_coef(VTUp,VTVp)
       call release_vf4_coef(VTUp,VTVp)
    endif                                             

    call release_df4_coef(VTU,VTV)
    call release_vf4_coef(VTU,VTV)
    call release_general_coef                      ! end_out_light


    call release_grids

    deallocate(VTU ); deallocate(VTU0); deallocate(VTU1); 
    deallocate(VTU2); deallocate(VTU3); deallocate(VTU4); 
    deallocate(VTUS); deallocate(SMU ,stat=ok); 
    if (ok/=0) stop 'VTU-SMU  : error alloc'
                                                   
    deallocate(VTV ); deallocate(VTV0); deallocate(VTV1); 
    deallocate(VTV2); deallocate(VTV3); deallocate(VTV4); 
    deallocate(VTVS); deallocate(SMV ,stat=ok); 
    if (ok/=0) stop 'VTV-SMV  : error alloc'
                                                   
    deallocate(PRE ); deallocate(PRE0); deallocate(PRE1);
    deallocate(PRE2); deallocate(PRE3); deallocate(PRES);
    deallocate(SMP ,stat=ok); if (ok/=0) stop 'SMP  : error alloc'

    if (ngrid_max.ge.0) then
       deallocate(U_EN_U,stat=ok); if (ok/=0) stop 'U_EN_U : error dealloc'
       deallocate(V_EN_V,stat=ok); if (ok/=0) stop 'V_EN_V : error dealloc'
    endif
    if (is_decal) deallocate(V_EN_U,stat=ok);  ! start_out_light
    if (ok/=0) stop 'V_EN_U : error alloc'  
    if (is_decal) deallocate(U_EN_V,stat=ok); 
    if (ok/=0) stop 'U_EN_V : error alloc'     ! end_out_light


    deallocate(XU); deallocate(YU); if (ok/=0) stop 'YU : error alloc'
    deallocate(XV); deallocate(YV); if (ok/=0) stop 'YV : error alloc'
    deallocate(XP); deallocate(YP,stat=ok); if (ok/=0) stop 'YP : error alloc'

    deallocate(CLDU,stat=ok); if (ok/=0) stop 'CLDU : error alloc'
    deallocate(CLDV,stat=ok); if (ok/=0) stop 'CLDV : error alloc'

    if (is_precond) then
       deallocate(MASK_DIRU,stat=ok); if (ok/=0) stop 'MASK_DIRU : error alloc'
       deallocate(MASK_DIRV,stat=ok); if (ok/=0) stop 'MASK_DIRV : error alloc'
    endif


    return 
  end subroutine coend

  !************************************************************************


  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine init_grids
  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    use data
    use drapeaux
    use champs
    use disc
    use para
    implicit none

    integer :: lmu,lmv,lmp,nmu,nmv,nmp,i,k,ok&
         &    ,lmu1,nmu1,lmv1,nmv1,lmp1,nmp1,lmmg,nmmg&
         &    ,lmu2,nmu2,lmv2,nmv2

    lmu=size(VTUS,1)-2; nmu=size(VTUS,2)-2
    lmv=size(VTVS,1)-2; nmv=size(VTVS,2)-2
    lmp=size(PRES,1)-2; nmp=size(PRES,2)-2

    allocate(XU(0:lmu+1,0:nmu+1)); allocate(YU(0:lmu+1,0:nmu+1)); 
    allocate(XV(0:lmv+1,0:nmv+1)); allocate(YV(0:lmv+1,0:nmv+1)); 
    allocate(XP(0:lmp+1,0:nmp+1)); allocate(YP(0:lmp+1,0:nmp+1),stat=ok); 
    if (ok/=0) stop 'XU-YP : error alloc'

    ! Domaine d'integration
    select case(nexample)             ! start_out_light
    case(1)       ! Solution 1D
       xmin=0.;     ymin=0._prec
       xmax=2.*pi;  ymax=1._prec       
    case(23)       ! Canal Droit            
       xmin=-1._prec; ymin=0._prec
       xmax=1._prec;  ymax=1._prec
    case(24)       ! Kovasznay
       xmin=-0.5_prec; ymin=-0._prec
       xmax= 1._prec;  ymax= 2._prec
    case(19,21,25)       ! Cavite Entrainee    
       xmin=-1._prec; ymin=-1._prec
       xmax= 1._prec; ymax=1._prec        
    case(26)       ! Tourbillon Taylor-Green
       xmin=0.0_prec*pi;  ymin= 0.0_prec*pi
       xmax=1.0_prec*pi;  ymax= 1.0_prec*pi
    case(27)       ! solution analytique NS
       xmin=1.0_prec;  ymin= 1.0_prec*pi
       xmax=2.0_prec;  ymax= 2.0_prec*pi
    case default                      ! end_out_light
       xmin=0._prec; ymin=0._prec
       xmax=1._prec; ymax=1._prec
    end select                        ! out_light

    if (is_decal) then                                              ! start_out_light
       dx=(xmax-xmin)/real(lmp_global-1)
       dy=(ymax-ymin)/real(nmp_global-1)
       xpdl=xmin+icolumn*int(real(lmp_global-2)/real(nb_i_blocks)+0.99)*dx; 
       zpdl=ymin+iline*int(real(nmp_global-2)/real(nb_k_blocks)+0.99)*dy
       xudl=xpdl-0.5_prec*dx;  zudl=zpdl
       xvdl=xpdl;              zvdl=zpdl-0.5_prec*dy
    else                                                            ! end_out_light
       dx=(xmax-xmin)/real(lmu_global-1)
       dy=(ymax-ymin)/real(nmu_global-1)
       xudl=xmin+icolumn*int(real(lm_global)/real(nb_i_blocks)+0.99)*dx; 
       zudl=ymin+iline*int(real(nm_global)/real(nb_k_blocks)+0.99)*dy
       xvdl=xudl;           zvdl=zudl
       xpdl=xudl+0.5_prec*dx;  zpdl=zudl+0.5_prec*dy;
    endif                                                           ! out_light

    invdx   = 1._prec/dx
    invdy   = 1._prec/dy
    invdx2  = 1._prec/dx/dx
    invdy2  = 1._prec/dy/dy
    p5invdx = 0.5_prec/dx
    p5invdy = 0.5_prec/dy
    xwdu=xmin; xedu=xmax
    ysdv=ymin; yndv=ymax
    
    if (is_oned) then
       invdy   = 0._prec
       invdy2  = 0._prec
       p5invdy = 0._prec
    endif




    if (is_decal) then                                         ! start_out_light
       XU = spread( (/ (i, i=0,lmu+1) /)*dx + xudl, 2, nmu+2)
       YU = spread( (/ (k, k=0,nmu+1) /)*dy + zudl, 1, lmu+2)
       XV = spread( (/ (i, i=0,lmv+1) /)*dx + xvdl, 2, nmv+2)
       YV = spread( (/ (k, k=0,nmv+1) /)*dy + zvdl, 1, lmv+2)
       XP = spread( (/ (i, i=0,lmp+1) /)*dx + xpdl, 2, nmp+2)
       YP = spread( (/ (k, k=0,nmp+1) /)*dy + zpdl, 1, lmp+2)
    else                                                       ! end_out_light
       XU = spread( (/ (i, i=0,lmu+1) /)*dx + xudl, 2, nmu+2)
       YU = spread( (/ (k, k=0,nmu+1) /)*dy + zudl, 1, lmu+2)
       XV = spread( (/ (i, i=0,lmv+1) /)*dx + xvdl, 2, nmv+2)
       YV = spread( (/ (k, k=0,nmv+1) /)*dy + zvdl, 1, lmv+2)
       XP = spread( (/ (i, i=0,lmp+1) /)*dx + xpdl, 2, nmp+2)
       YP = spread( (/ (k, k=0,nmp+1) /)*dy + zpdl, 1, lmp+2)
    endif                                                     ! out_light  

    nrecvddmx=abs(nrecvddm)
    nrecvddmy=abs(nrecvddm)
    if (nrecvddm<0.and.nb_k_blocks<nb_i_blocks.and.nb_k_blocks>1)&
         & nrecvddmy=nrecvddmy*2
    nrecvddm=abs(nrecvddm)

    lmu1=lmu
    nmu1=nmu
    if (.not.is_west)  lmu1=lmu1+nrecvddmx-1
    if (.not.is_east)  lmu1=lmu1+nrecvddmx-1
    if (.not.is_north) nmu1=nmu1+nrecvddmy-1
    if (.not.is_south) nmu1=nmu1+nrecvddmy-1
    lmv1=lmv+lmu1-lmu
    nmv1=nmv+nmu1-nmu
    lmp1=lmp+lmu1-lmu
    nmp1=nmp+nmu1-nmu

    iudeb=nrecvddmx-1; kudeb=nrecvddmy-1
    if (is_south) kudeb=0
    if (is_west)  iudeb=0
    ivdeb=iudeb; kvdeb=kudeb
    ipdeb=iudeb; kpdeb=kudeb
    iufin=iudeb+lmu+1; kufin=kudeb+nmu+1
    ivfin=ivdeb+lmv+1; kvfin=kvdeb+nmv+1
    ipfin=ivdeb+lmp+1; kpfin=kvdeb+nmp+1

    ngrid_max=-1
    if (ntype_precond/=0) ngrid_max=0

    if (abs(ntype_solver)==mg_solver.or.abs(ntype_precond)==mg_solver) then                   ! start_out_light
       if ( (is_ns.and..not.is_decal).or. &
            & (nb_tasks>1.and.abs(ntype_solver)==mg_solver) ) &
            & stop 'MG not yet implemented for and not decal meshed' !

       ngrid_max=10

       if (is_decal) then
          lmmg=lmp1
          nmmg=nmp1

          call get_mg_dim(lmp,lmmg,ngrid_max,3)
          call get_mg_dim(nmp,nmmg,ngrid_max,3)

          kpmgdeb=(lmp1-lmmg)/2   
          ipmgdeb=(nmp1-nmmg)/2
          
          if (is_north) kpmgdeb=kpfin-nmmg-1
          if (is_east)  ipmgdeb=ipfin-lmmg-1
          if (is_south) kpmgdeb=0
          if (is_west)  ipmgdeb=0
          ipmgfin=ipmgdeb+lmmg+1; kpmgfin=kpmgdeb+nmmg+1
          iumgdeb=ipmgdeb; iumgfin=ipmgfin+1
          ivmgdeb=ipmgdeb; ivmgfin=ipmgfin
          kumgdeb=kpmgdeb; kumgfin=kpmgfin
          kvmgdeb=kpmgdeb; kvmgfin=kpmgfin+1

          iudeb=ipdeb;   iufin=ipfin+1
          kudeb=ipdeb;   kufin=kpfin
          ivdeb=ipdeb;   ivfin=ipfin
          kvdeb=ipdeb;   kvfin=kpfin+1       
       else
          lmmg=lmu1
          nmmg=nmu1
          
          call get_mg_dim(lmu,lmmg,ngrid_max,2)
          call get_mg_dim(nmu,nmmg,ngrid_max,2)
          
          iumgdeb=(lmu1-lmmg)/2   
          kumgdeb=(nmu1-nmmg)/2
          
          if (is_north) kumgdeb=kufin-nmmg-1
          if (is_east)  iumgdeb=iufin-lmmg-1
          if (is_south) kumgdeb=0
          if (is_west)  iumgdeb=0
          iumgfin=iumgdeb+lmmg+1; kumgfin=kumgdeb+nmmg+1
          
          ivdeb=iudeb;         ivfin=iufin
          kvdeb=kudeb;         kvfin=kufin
          ipdeb=iudeb;         ipfin=iufin-1
          ipdeb=kudeb;         kpfin=kufin-1

          ivmgdeb=iumgdeb;         ivmgfin=iumgfin
          kvmgdeb=kumgdeb;         kvmgfin=kumgfin
          ipmgdeb=iumgdeb;         ipmgfin=iumgfin-1
          ipmgdeb=kumgdeb;         kpmgfin=kufin-1
       endif
    else                                                                              ! end_out_light
       iumgdeb=0; iumgfin=lmu1+1
       kumgdeb=0; kumgfin=nmu1+1
       ivmgdeb=0; ivmgfin=lmv1+1
       kvmgdeb=0; kvmgfin=nmv1+1
       ipmgdeb=0; ipmgfin=lmp1+1
       kpmgdeb=0; kpmgfin=nmp1+1
    endif                                                                             ! out_light

    call output(s='nb grid max / soit  : ',i1=ngrid_max,i2=2**ngrid_max-1)
    call output('lmu,nmu',lmu,nmu)
    call output('lmv,nmv',lmv,nmv)
    call output('lmp,nmp',lmp,nmp)
    call output('lmu1,nmu1',lmu1,nmu1)
    call output('lmv1,nmv1',lmv1,nmv1)
    call output('lmp1,nmp1',lmp1,nmp1)
    call output('lmmg,nmmg',lmmg,nmmg)
    call output('iudeb,iufin,kudeb,kufin',iudeb,iufin,kudeb,kufin)
    call output('ivdeb,ivfin,kvdeb,kvfin',ivdeb,ivfin,kvdeb,kvfin)
    call output('ipdeb,ipfin,kpdeb,kpfin',ipdeb,ipfin,kpdeb,kpfin)
    call output('iumgdeb,iumgfin,kumgdeb,kumgfin',iumgdeb,iumgfin,kumgdeb,kumgfin) 
    call output('ivmgdeb,ivmgfin,kvmgdeb,kvmgfin',ivmgdeb,ivmgfin,kvmgdeb,kvmgfin) 
    call output('ipmgdeb,ipmgfin,kpmgdeb,kpmgfin',ipmgdeb,ipmgfin,kpmgdeb,kpmgfin) 


    allocate(grid(-1:ngrid_max),stat=ok); if (ok/=0) stop 'grid : error alloc'

    grid(-1)%dx      = dx
    grid(-1)%dy      = dy
    grid(-1)%invdx   = invdx  
    grid(-1)%invdy   = invdy  
    grid(-1)%invdx2  = invdx2 
    grid(-1)%invdy2  = invdy2 
    grid(-1)%p5invdx = p5invdx
    grid(-1)%p5invdy = p5invdy
    grid(-1)%epsv    = epsv
    grid(-1)%epsva   = epsva

    if (ncheck==vf4.or.ncheck==vf2) then
       lmu2=lmu+2
       nmv2=nmv+2
!!$       if (ncheck==vf2)  lmu2=lmu+1
!!$       if (ncheck==vf2)  nmv2=nmv+1
       allocate(grid(-1)%U_EN_U(0:lmu2,0:nmu+1),stat=ok); 
       if (ok/=0) stop 'U_EN_U : error alloc'
       allocate(grid(-1)%V_EN_V(0:lmv+1,0:nmv2),stat=ok); 
       if (ok/=0) stop 'V_EN_V : error alloc'
    else
       allocate(grid(-1)%U_EN_U(0:lmu+1,0:nmu+1),stat=ok); 
       if (ok/=0) stop 'U_EN_U : error alloc'
       allocate(grid(-1)%V_EN_V(0:lmv+1,0:nmv+1),stat=ok); 
       if (ok/=0) stop 'V_EN_V : error alloc'
    endif
    
    grid(-1)%U_EN_U_c => grid(-1)%U_EN_U(1:size(grid(-1)%U_en_U,1)-2,1:size(grid(-1)%U_en_U,2)-2)
    grid(-1)%V_EN_V_c => grid(-1)%V_EN_V(1:size(grid(-1)%V_en_V,1)-2,1:size(grid(-1)%V_en_V,2)-2)
    
    if (is_decal) then                    ! start_out_light
       nmu2=nmu+1; if (ncheck==vf4.or.(ncheck==vf2.and.is_north)) nmu2=nmu+2
       lmv2=lmv+1; if (ncheck==vf4.or.(ncheck==vf2.and.is_east))  lmv2=lmv+2
       allocate(grid(-1)%V_EN_U(0:lmu+1,0:nmu2),stat=ok); 
       if (ok/=0) stop 'V_EN_U : error alloc'
       allocate(grid(-1)%U_EN_V(0:lmv2,0:nmv+1),stat=ok); 
       if (ok/=0) stop 'U_EN_V : error alloc'
       grid(-1)%V_EN_U_c => grid(-1)%V_EN_U(1:size(grid(-1)%V_en_U,1)-2,1:size(grid(-1)%V_en_U,2)-2)
       grid(-1)%U_EN_V_c => grid(-1)%U_EN_V(1:size(grid(-1)%U_en_V,1)-2,1:size(grid(-1)%U_en_V,2)-2)
    endif


    lmu=lmu1; nmu=nmu1;
    lmv=lmv1; nmv=nmv1;
    lmp=lmp1; nmp=nmp1

    lmu=iumgfin-iumgdeb-1
    nmu=kumgfin-kumgdeb-1
    lmv=ivmgfin-ivmgdeb-1
    nmv=kvmgfin-kvmgdeb-1
    lmp=ipmgfin-ipmgdeb-1
    nmp=kpmgfin-kpmgdeb-1


    if (ntype_precond/=0.or.(abs(ntype_precond)==mg_solver.or.abs(ntype_solver)==mg_solver)) then
       allocate (VTUp(0:lmu+1,0:nmu+1),stat=ok); 
       allocate (VTVp(0:lmv+1,0:nmv+1),stat=ok); 

       allocate(MASK_DIRU(0:lmu+1,0:nmu+1),stat=ok); 
                 if (ok/=0) stop 'MASK_DIRU : error alloc'
       allocate(MASK_DIRV(0:lmv+1,0:nmv+1),stat=ok); 
                 if (ok/=0) stop 'MASK_DIRV : error alloc'
       

       VTUp=1._prec; VTUp(0,0)=2._prec; VTUp(lmu+1,0)=2._prec;  
       VTUp(0,nmu+1)=2._prec; VTUp(lmu+1,nmu+1)=2._prec; 
       VTVp=1._prec; VTVp(0,0)=2._prec; VTVp(lmv+1,0)=2._prec;  
       VTVp(0,nmv+1)=2._prec; VTVp(lmv+1,nmv+1)=2._prec; 
       if (.not.is_north) VTUp(:,nmu+1) = VTUp(:,nmu+1)-1.;      
       if (.not.is_east)  VTUp(lmu+1,:) = VTUp(lmu+1,:)-1.;      
       if (.not.is_west)  VTUp(0,:)      = VTUp(0,:)-1.;          
       if (.not.is_south) VTUp(:,0)      = VTUp(:,0)-1.;          
       MASK_DIRU = VTUp==0.                       ;         
       if (.not.is_north) VTVp(:,nmv+1) = VTVp(:,nmv+1)-1.  
       if (.not.is_east)  VTVp(lmv+1,:) = VTVp(lmv+1,:)-1.  
       if (.not.is_west)  VTVp(0,:)     = VTVp(0,:)-1.      
       if (.not.is_south) VTVp(:,0)     = VTVp(:,0)-1.      
       MASK_DIRV = VTVp==0.                            

       deallocate (VTUp,stat=ok); 
       deallocate (VTVp,stat=ok); 

       allocate (VTUp(0:lmu1+1,0:nmu1+1),stat=ok); 
       allocate (VTVp(0:lmv1+1,0:nmv1+1),stat=ok); 
       allocate (PREp(0:lmp1+1,0:nmp1+1),stat=ok); 
       allocate (SMUp(0:lmu1+1,0:nmu1+1),stat=ok); 
       allocate (SMVp(0:lmv1+1,0:nmv1+1),stat=ok); 
       allocate (SMPp(0:lmp1+1,0:nmp1+1),stat=ok); 

       VTUp=0._prec; SMUp=0._prec
       VTVp=0._prec; SMVp=0._prec
       PREp=0._prec; SMPp=0._prec
    endif


    do k=ngrid_max,0,-1

       grid(k)%dx      = dx
       grid(k)%dy      = dy
       grid(k)%invdx   = invdx  
       grid(k)%invdy   = invdy  
       grid(k)%invdx2  = invdx2 
       grid(k)%invdy2  = invdy2 
       grid(k)%p5invdx = p5invdx
       grid(k)%p5invdy = p5invdy
       grid(k)%epsv    = epsvc
       grid(k)%epsva   = epsvac

       call output('ici ncheck_precond,vf4,lmu+4,nmu+3',ncheck_precond,vf4,lmu+4,nmu+3)
       if (ncheck_precond==vf4.or.ncheck_precond==vf2) then
          allocate(grid(k)%U_EN_U(0:lmu+2,0:nmu+1),stat=ok); 
          if (ok/=0) stop 'U_EN_U : error alloc'
          allocate(grid(k)%V_EN_V(0:lmv+1,0:nmv+2),stat=ok); 
          if (ok/=0) stop 'V_EN_V : error alloc'
          call output('la size size(grid(k)%U_en_U)',size(grid(k)%U_en_U,1),size(grid(k)%U_en_U,2))
          call output('la size size(grid(k)%V_en_V)',size(grid(k)%V_en_V,1),size(grid(k)%V_en_V,2))
       else
          allocate(grid(k)%U_EN_U(0:lmu+1,0:nmu+1),stat=ok); 
          if (ok/=0) stop 'U_EN_U : error alloc'
          allocate(grid(k)%V_EN_V(0:lmv+1,0:nmv+1),stat=ok); 
          if (ok/=0) stop 'V_EN_V : error alloc'
       endif
       
       grid(k)%U_EN_U_c => grid(k)%U_EN_U(1:size(grid(k)%U_en_U,1)-2,1:size(grid(k)%U_en_U,2)-2)
       grid(k)%V_EN_V_c => grid(k)%V_EN_V(1:size(grid(k)%V_en_V,1)-2,1:size(grid(k)%V_en_V,2)-2)

       if (is_decal) then                    ! start_out_light
          allocate(grid(k)%V_EN_U(0:size(grid(k)%U_en_U,1)-1,0:size(grid(k)%U_en_U,2)-1),stat=ok); 
          if (ok/=0) stop 'V_EN_U : error alloc'
          allocate(grid(k)%U_EN_V(0:size(grid(k)%V_en_V,1)-1,0:size(grid(k)%V_en_V,2)-1),stat=ok); 
          if (ok/=0) stop 'U_EN_V : error alloc'
          grid(k)%V_EN_U_c => grid(k)%V_EN_U(1:size(grid(k)%U_en_U,1)-2,1:size(grid(k)%U_en_U,2)-2)
          grid(k)%U_EN_V_c => grid(k)%U_EN_V(1:size(grid(k)%V_en_V,1)-2,1:size(grid(k)%V_en_V,2)-2)
       endif                                 ! end_out_light
    
       if (is_u.and.(abs(ntype_precond)==mg_solver.or.abs(ntype_solver)==mg_solver)) then   ! start_out_light
          if (is_decal) then
             allocate(grid(k)%U1(0:lmp+2,0:nmp+1),stat=ok);     if (ok/=0) stop 'U1 : error alloc'
             allocate(grid(k)%U2(0:lmp/3+2,0:nmp/3+1),stat=ok); if (ok/=0) stop 'U2 : error alloc'
             allocate(grid(k)%U3(0:lmp/3+2,0:nmp/3+1),stat=ok); if (ok/=0) stop 'U3 : error alloc'
          else
             allocate(grid(k)%U1(0:lmu+1,0:nmu+1),stat=ok);     if (ok/=0) stop 'U1 : error alloc'
             allocate(grid(k)%U2(0:lmu/2+1,0:nmu/2+1),stat=ok); if (ok/=0) stop 'U2 : error alloc'
             allocate(grid(k)%U3(0:lmu/2+1,0:nmu/2+1),stat=ok); if (ok/=0) stop 'U3 : error alloc'
          endif
          grid(k)%U1=0._prec; grid(k)%U2=0._prec; grid(k)%U3=0._prec; 
       endif

       if (is_v.and.(abs(ntype_precond)==mg_solver.or.abs(ntype_solver)==mg_solver)) then
          if (is_decal) then
             allocate(grid(k)%V1(0:lmp+1,0:nmp+2),stat=ok);     if (ok/=0) stop 'V1 : error alloc'
             allocate(grid(k)%V2(0:lmp/3+1,0:nmp/3+2),stat=ok); if (ok/=0) stop 'V2 : error alloc'
             allocate(grid(k)%V3(0:lmp/3+1,0:nmp/3+2),stat=ok); if (ok/=0) stop 'V3 : error alloc'
          else
             allocate(grid(k)%V1(0:lmv+1,0:nmv+1),stat=ok);     if (ok/=0) stop 'V1 : error alloc'
             allocate(grid(k)%V2(0:lmv/2+1,0:nmv/2+1),stat=ok); if (ok/=0) stop 'V2 : error alloc'
             allocate(grid(k)%V3(0:lmv/2+1,0:nmv/2+1),stat=ok); if (ok/=0) stop 'V3 : error alloc'
          endif
          grid(k)%V1=0._prec; grid(k)%V2=0._prec; grid(k)%V3=0._prec; 
       endif

       if ((abs(ntype_precond)==mg_solver.or.abs(ntype_solver)==mg_solver)) then
          if (is_decal) then
             allocate(grid(k)%P1(0:lmp+1,0:nmp+1),stat=ok);     if (ok/=0) stop 'P1 : error alloc'
             allocate(grid(k)%P2(0:lmp/3+1,0:nmp/3+1),stat=ok); if (ok/=0) stop 'P2 : error alloc'
             allocate(grid(k)%P3(0:lmp/3+1,0:nmp/3+1),stat=ok); if (ok/=0) stop 'P3 : error alloc'
          else
             allocate(grid(k)%P1(0:lmp+1,0:nmp+1),stat=ok);     if (ok/=0) stop 'P1 : error alloc'
             allocate(grid(k)%P2(0:lmp/2+1,0:nmp/2+1),stat=ok); if (ok/=0) stop 'P2 : error alloc'
             allocate(grid(k)%P3(0:lmp/2+1,0:nmp/2+1),stat=ok); if (ok/=0) stop 'P3 : error alloc'
          endif
          grid(k)%P1=0._prec; grid(k)%P2=0._prec; grid(k)%P3=0._prec; 
       endif


       call output('grille,lmu+1,nmu+1 ',k,lmu+1,nmu+1)
       call output('grille,lmv+1,nmv+1 ',k,lmv+1,nmv+1)
       call output('grille,lmp+1,nmp+1 ',k,lmp+1,nmp+1)
       call output('dx,invdx,invdx2,p5invdx',r1=dx,r2=invdx,r3=invdx2,r4=p5invdx)
       call output('dy,invdy,invdy2,p5invdy',r1=dy,r2=invdy,r3=invdy2,r4=p5invdy)       ! end_out_light

       if (is_decal) then
          lmp=lmp/3; lmu=lmp+1; lmv=lmp
          nmp=nmp/3; nmu=nmp;   nmv=nmp+1
       else
          lmu=lmu/2;  lmv=lmv/2
          nmu=nmu/2;  nmv=nmv/2
       endif

       if (is_decal) then
          dx=dx*3._prec   ;  dy=dy*3._prec
       else
          dx=dx*2._prec   ;  dy=dy*2._prec
       endif
       invdx   = 1._prec/dx
       invdy   = 1._prec/dy
       invdx2  = 1._prec/dx/dx
       invdy2  = 1._prec/dy/dy
       p5invdx = 0.5_prec/dx
       p5invdy = 0.5_prec/dy

    enddo


    call set_grid(-1)

    return
  end subroutine init_grids

  !************************************************************************  ! start_out_light

  subroutine get_mg_dim(lm,lmmg,ngrid_max,n) 
    implicit none
    integer :: lm,lmmg,ngrid_max,n

    integer :: i,i0,j

    i=0;     i0=0;
    j=lm+1
    do while (j>=lm)
       i=i+1;       
       j=(int((lmmg+1)/n**i))*n**i-1
       if (lmmg-j<=j-lm.and.j/n**i>=2) i0=i
    enddo
    lmmg=(int((lmmg+1)/n**i0))*n**i0-1

    ngrid_max=min(ngrid_max,i0)

    return
  end subroutine get_mg_dim                                               ! end_out_light


  !************************************************************************
    
    subroutine set_grid (level)
      use disc
      use data, only : epsv, epsva
      use drapeaux, only : is_decal
      implicit none
      integer :: level

      dx      =  grid(level)%dx      
      dy      =  grid(level)%dy      
      invdx   =  grid(level)%invdx   
      invdy   =  grid(level)%invdy   
      invdx2  =  grid(level)%invdx2  
      invdy2  =  grid(level)%invdy2  
      p5invdx =  grid(level)%p5invdx 
      p5invdy =  grid(level)%p5invdy 
      epsv    =  grid(level)%epsv
      epsva   =  grid(level)%epsva

      U_EN_U   => grid(level)%U_EN_U
      V_EN_V   => grid(level)%V_EN_V
      U_EN_U_c => grid(level)%U_EN_U_c
      V_EN_V_c => grid(level)%V_EN_V_c

      if (is_decal) then          ! start_out_light
         U_EN_V   => grid(level)%U_EN_V
         V_EN_U   => grid(level)%V_EN_U
         U_EN_V_c => grid(level)%U_EN_V_c
         V_EN_U_c => grid(level)%V_EN_U_c
    else                                      ! end_out_light
         U_EN_V    => grid(level)%U_EN_U
         V_EN_U    => grid(level)%V_EN_V
         U_EN_V_c  => grid(level)%U_EN_U_c
         V_EN_U_c  => grid(level)%V_EN_V_c
      endif                                     ! start_out_light

!      call output('now grid_level ',level)   ! end_out_light

      return
    end subroutine set_grid



  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine release_grids
  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    use data
    use drapeaux
    use champs
    use disc
    use para
    implicit none

    integer level,ok

    do level=0,ngrid_max

       if (is_u.and.(abs(ntype_precond)==mg_solver.or.abs(ntype_solver)==mg_solver)) then   ! start_out_light
          if (is_decal) then
             deallocate(grid(level)%U1,stat=ok);     if (ok/=0) stop 'U1 : error dealloc'
             deallocate(grid(level)%U2,stat=ok); if (ok/=0) stop 'U2 : error dealloc'
             deallocate(grid(level)%U3,stat=ok); if (ok/=0) stop 'U3 : error dealloc'
          else
             deallocate(grid(level)%U1,stat=ok);     if (ok/=0) stop 'U1 : error dealloc'
             deallocate(grid(level)%U2,stat=ok); if (ok/=0) stop 'U2 : error dealloc'
             deallocate(grid(level)%U3,stat=ok); if (ok/=0) stop 'U3 : error dealloc'
          endif
       endif

       if (is_v.and.(abs(ntype_precond)==mg_solver.or.abs(ntype_solver)==mg_solver)) then
          if (is_decal) then
             deallocate(grid(level)%V1,stat=ok); if (ok/=0) stop 'V1 : error dealloc'
             deallocate(grid(level)%V2,stat=ok); if (ok/=0) stop 'V2 : error dealloc'
             deallocate(grid(level)%V3,stat=ok); if (ok/=0) stop 'V3 : error dealloc'
          else
             deallocate(grid(level)%V1,stat=ok); if (ok/=0) stop 'V1 : error dealloc'
             deallocate(grid(level)%V2,stat=ok); if (ok/=0) stop 'V2 : error dealloc'
             deallocate(grid(level)%V3,stat=ok); if (ok/=0) stop 'V3 : error dealloc'
          endif
       endif

       if ((abs(ntype_precond)==mg_solver.or.abs(ntype_solver)==mg_solver)) then
          if (is_decal) then
             deallocate(grid(level)%P1,stat=ok); if (ok/=0) stop 'P1 : error dealloc'
             deallocate(grid(level)%P2,stat=ok); if (ok/=0) stop 'P2 : error dealloc'
             deallocate(grid(level)%P3,stat=ok); if (ok/=0) stop 'P3 : error dealloc'
          else
             deallocate(grid(level)%P1,stat=ok); if (ok/=0) stop 'P1 : error dealloc'
             deallocate(grid(level)%P2,stat=ok); if (ok/=0) stop 'P2 : error dealloc'
             deallocate(grid(level)%P3,stat=ok); if (ok/=0) stop 'P3 : error dealloc'
          endif
       endif
    end do
       
    do level=-1,ngrid_max
       deallocate(grid(level)%U_EN_U,stat=ok); 
       deallocate(grid(level)%V_EN_V,stat=ok); 
       if (is_decal) then                    ! start_out_light
          deallocate(grid(level)%V_EN_U,stat=ok); 
          deallocate(grid(level)%U_EN_V,stat=ok); 
       endif
    end do

    deallocate(grid,stat=ok); if (ok/=0) stop 'grid : error dealloc'


    if (ntype_precond/=0.or.(abs(ntype_precond)==mg_solver.or.abs(ntype_solver)==mg_solver)) then
       deallocate (VTUp,stat=ok); 
       deallocate (VTVp,stat=ok); 
       deallocate (PREp,stat=ok); 
       deallocate (SMUp,stat=ok); 
       deallocate (SMVp,stat=ok); 
       deallocate (SMPp,stat=ok); 
    endif

    return
  end subroutine release_grids

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine lecture_input
  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    use data
    use disc
    use uncol
    implicit none
    integer :: i,i0,ind, iargc, numarg

    integer, parameter :: input_unit = 7, output_unit=6, ok = 0
    integer :: open_status
    character (len=200) :: line, arg, input_file_name
    integer :: dat1, RetCode, nb_input
    real(kind=prec), dimension(96) :: rbuffer=0.
    real(kind=prec), dimension(9) :: nbuffer=0.

    if (my_task==0) then


       !!call timestamp ( )

       numarg = iargc ( )
       
       if (numarg /=1) then
          print *,"usage : zephyr <input file>   by default using ./input"
          input_file_name = "./input"
       else
          call getarg ( 1, input_file_name)
       end if
       


!!       open( unit=output_unit, & 
!!         file='/bgl/FS2/ohess/Codes/zephyr/output', status="new", &
!!!!            iostat=open_status )
       
       ! if( open_status /= ok ) then
       !    print *
       !    print *, "Unable to open file"
       ! end if

       open( unit=input_unit, & 
         file=input_file_name, status="old", &
            iostat=open_status )
       
       if( open_status /= ok ) then
          print *
          print *, "Unable to open input file ",input_file_name
       end if

       read (unit=input_unit, fmt=*) i

       nb_input = 1
       i0=23
       i0=40    ! out_light
       read_loop: do
          read (input_unit, '(A)', iostat=RetCode)  line
          if ( nb_input>i0 +1.or. RetCode == iostat_end)  exit read_loop
           if ( RetCode /= 0 ) then
             print *,"erreur de lecture du fichier input"
             stop
             exit read_loop
           end if
           ind = index (line, "!")
           if (ind  == 1 )  cycle read_loop
           line = line(1:ind-1)
           !!print *,"ligne lue :",line
           if (len_trim(line)==0)  cycle read_loop

           if (nb_input == 1) then
              nom_fic_output = line
              !!print *,"fic_ouput = ",nom_fic_output
           else if (nb_input == i0+1) then
              nom_fic_save = line
              !! print *,"fic_save = ",nom_fic_save
           else
              read (line, fmt=*) rbuffer(nb_input-1)
              !!print *,"data ",nb_input,":",rbuffer(nb_input-1)
           end if
           nb_input = nb_input + 1
       end do read_loop

       ! read (unit=input_unit, fmt=*)   
       ! i0=23
       ! i0=39    ! out_light
       ! do i=1,i0
       !    read (unit=input_unit, fmt=*) rbuffer(i)
       ! enddo
       ! read (unit=input_unit, fmt=*) nom_fic_save
       ! close (unit=input_unit)

!!!       read *,i
!!!       read *,  nom_fic_output
!!!       i0=23
!!!       i0=39    ! out_light
!!!       do i=1,i0
!!!          read *,rbuffer(i)
!!!       enddo
!!!       read *,nom_fic_save

       nom_fic_save = 'SAVE'

       if (rbuffer(39)>0) call restore_param(rbuffer,i0)         ! out_light

       rbuffer(40)    = t_start             ! end_out_light
       rbuffer(41)    = it_start            ! end_out_light
    
       do i=1,nb_tasks-1
          call snd_msg(i,900+i,rbuffer)
       enddo
       
!!$       if (len(nom_fic_output)/=96) stop 'pb de transmission de chaine'
!!$       do i=1,96
!!$          nbuffer(i)=ichar(nom_fic_output(i:i))
!!$       enddo
!!$       print *,nbuffer
!!$       do i=1,nb_tasks-1
!!$          call snd_msg(i,920+i,nbuffer)
!!$       enddo
    else
       call rcv_msg(0,900+my_task,rbuffer)

!!$       call rcv_msg(0,920+my_task,nbuffer)
!!$       write(nom_fic_output,'(96A1)') (char(int(nbuffer(i))),i=1,96)
    endif
    
    nexample        =  rbuffer( 1) 
    t_start         =  rbuffer( 2) 
    t_all           =  rbuffer( 3) 
    tau             =  rbuffer( 4) 
    t_print         =  rbuffer( 5) 
    is_print        =  rbuffer( 6) 
    is_cv           =  rbuffer( 7) 
    is_kuta         =  rbuffer( 8) 
    is_richardson   =  rbuffer( 9) 
    lm_global       =  rbuffer(10) 
    nm_global       =  rbuffer(11) 
    nu              =  rbuffer(12) 
    ncheck          =  rbuffer(13) 
    ncheck_precond  =  rbuffer(14) 
    ntype_solver    =  rbuffer(15) 
    ntype_precond   =  rbuffer(16) 
    npcg            =  rbuffer(17) 
    nprecond        =  rbuffer(18) 
    nrecvddm        =  rbuffer(19) 
    epsv            =  rbuffer(20) 
    epsva           =  rbuffer(21) 
    epsvc           =  rbuffer(22) 
    epsvac          =  rbuffer(23)       ! start_out_light
    ninterne        =  rbuffer(24) 
    ndirection      =  rbuffer(25) 
    ns_solver       =  rbuffer(26) 
    rho_ca          =  rbuffer(27) 
    rlag            =  rbuffer(28) 
    nbmax_ca        =  rbuffer(29) 
    nbp_ca          =  rbuffer(30) 
    nbmg            =  rbuffer(31) 
    epsmg           =  rbuffer(32) 
    epsmga          =  rbuffer(33) 
    nb_prelis       =  rbuffer(34) 
    nb_postlis      =  rbuffer(35) 
    nb_cycle        =  rbuffer(36) 
    sor_theta       =  rbuffer(37)     
    is_decale       =  rbuffer(38) 
    is_restart_save =  rbuffer(39)       ! end_out_light
    !t_start         =  rbuffer(40)       ! end_out_light
    it_start        =  rbuffer(41)       ! end_out_light
    

!!$    print *,my_task,': lm,nm=',lm_global,nm_global

!!$    print *,my_task,':', &                                            ! start_out_light
!!$       is_print,is_kuta,nexample,is_unsteady                &
!!$       ,is_richardson,is_restart_save,is_decale             &
!!$       ,lm_global,nm_global                                 &
!!$       ,ncheck,ncheck_precond,ns_solver,ntype_solver        &
!!$       ,ntype_precond,nrecvddm,nprecond,npcg,ninterne       &
!!$       ,nbmax_ca,nbp_ca                                     &
!!$       ,nbmg,nb_prelis,nb_postlis,nb_cycle
!!$
!!$    print *,my_task,':', &
!!$       t_start,t_all,tau,t_print                        &
!!$       ,nu,epsv,epsvc,epsva,epsvac,epsmg,epsmga         &
!!$       ,rlag,rho_ca,tol_ca,sor_theta                    
    
!!$    call parallel_stop                                                          ! end_out_light

 return
 end subroutine lecture_input

 !************************************************************************          ! start_out_light
 
 !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
 subroutine save_data
   !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
   use data
   use drapeaux
   use disc
   use champs
   use uncol, only : tchrono
   use debug
   implicit none
   integer :: ok
   character (len=200) :: nom_fichier_save

   nom_fichier_save = trim(nom_fic_save)
   if (my_task==0) then
      open(file=nom_fichier_save,unit=69,form='unformatted',iostat=ok)
      if (ok/=0) stop 'pb ouverture fichier sauvegarde'  
      
      is_restart_save=abs(is_restart_save)
      write (69) nexample, t_start, t_all, tau, t_print, is_print, is_cv,&
              & is_kuta, is_richardson, lm_global, nm_global, nu, ncheck,&
              & ncheck_precond, ntype_solver, ntype_precond, npcg, nprecond,&
              & nrecvddm, epsv, epsva, epsvc, epsvac, ninterne, ndirection, ns_solver,&
              & rho_ca, rlag, nbmax_ca, nbp_ca, nbmg, epsmg, epsmga, nb_prelis, nb_postlis,&
              & nb_cycle, sor_theta, is_decale, is_restart_save,temps,it,tchrono
!!$      print *,'write', nexample, t_start, t_all, tau, t_print, is_print, is_cv,&
!!$              & is_kuta, is_richardson, lm_global, nm_global, nu, ncheck,&
!!$              & ncheck_precond, ntype_solver, ntype_precond, npcg, nprecond,&
!!$              & nrecvddm, epsv, epsva, epsvc, epsvac, ninterne, ndirection, ns_solver,&
!!$              & rho_ca, rlag, nbmax_ca, nbp_ca, nbmg, epsmg, epsmga, nb_prelis, nb_postlis,&
!!$              & nb_cycle, sor_theta, is_decale, is_restart_save
      if (debug_save) then
         print *,'parametres globaux sauvés par task 0'
      endif
   endif

   call save_or_retrieve('VTU0',VTU0,1,1,69,p_save); call save_or_retrieve('VTV0',VTV0,1,1,69,p_save); 
   call save_or_retrieve('VTU1',VTU1,1,1,69,p_save); call save_or_retrieve('VTV1',VTV1,1,1,69,p_save); 
   call save_or_retrieve('VTU2',VTU2,1,1,69,p_save); call save_or_retrieve('VTV2',VTV2,1,1,69,p_save); 
   call save_or_retrieve('VTU3',VTU3,1,1,69,p_save); call save_or_retrieve('VTV3',VTV3,1,1,69,p_save); 

   if (is_ns) then
      call save_or_retrieve('PRE0',PRE0,1,1,69,p_save); 
      call save_or_retrieve('PRE1',PRE1,1,1,69,p_save); 
      call save_or_retrieve('PRE2',PRE2,1,1,69,p_save); 
      call save_or_retrieve('PRE3',PRE3,1,1,69,p_save); 
   endif

   if (my_task==0) close(69)
   if (my_task==0) write(*,999) nom_fichier_save,it

999 format ('*** Data Saved in file ',A20,'for it=',I10)

   return
 end subroutine save_data

 !************************************************************************  
 
 !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
 subroutine restore_param(rbuffer,i)
   !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
   use data
   use disc
   use uncol, only : tchrono
   implicit none
   
   integer :: ok
   real(kind=prec) :: t_all0, t_start0
   real(kind=prec), dimension(:) :: rbuffer
   integer :: i
   character (len=200) :: nom_fichier_save

   nom_fichier_save = trim(nom_fic_save)

   if (my_task==0) then   
      open(file=nom_fichier_save,unit=69,status='old',form='unformatted',iostat=ok)
      if (ok/=0) then
         print *,'*** fichier SAVE non encore cree'
         rbuffer(i)=-rbuffer(i)
      else
         read(69) nexample, t_start0, t_all0, tau, t_print, is_print, is_cv,&
              & is_kuta, is_richardson, lm_global, nm_global, nu, ncheck,&
              & ncheck_precond, ntype_solver, ntype_precond, npcg, nprecond,&
              & nrecvddm, epsv, epsva, epsvc, epsvac, ninterne, ndirection, ns_solver,&
              & rho_ca, rlag, nbmax_ca, nbp_ca, nbmg, epsmg, epsmga, nb_prelis, nb_postlis,&
              & nb_cycle, sor_theta, is_decale, is_restart_save,temps,it,tchrono

         print *,'*** fichier SAVE en cours de lecture'
         
!!$      print *,'read', nexample, t_start, t_all, tau, t_print, is_print, is_cv,&
!!$              & is_kuta, is_richardson, lm_global, nm_global, nu, ncheck,&
!!$              & ncheck_precond, ntype_solver, ntype_precond, npcg, nprecond,&
!!$              & nrecvddm, epsv, epsva, epsvc, epsvac, ninterne, ndirection, ns_solver,&
!!$              & rho_ca, rlag, nbmax_ca, nbp_ca, nbmg, epsmg, epsmga, nb_prelis, nb_postlis,&
!!$              & nb_cycle, sor_theta, is_decale, is_restart_save
!!$
         if (temps>=t_all) t_all=temps+max(t_all-t_start,t_all0-t_start0)
         t_start=temps
         it_start=it

         rbuffer( 1) =  nexample        
         rbuffer( 2) =  t_start         
         rbuffer( 3) =  t_all           
         rbuffer( 4) =  tau             
         rbuffer( 5) =  t_print         
         rbuffer( 6) =  is_print        
         rbuffer( 7) =  is_cv           
         rbuffer( 8) =  is_kuta         
         rbuffer( 9) =  is_richardson   
         rbuffer(10) =  lm_global       
         rbuffer(11) =  nm_global       
         rbuffer(12) =  nu              
         rbuffer(13) =  ncheck          
         rbuffer(14) =  ncheck_precond  
         rbuffer(15) =  ntype_solver    
         rbuffer(16) =  ntype_precond   
         rbuffer(17) =  npcg            
         rbuffer(18) =  nprecond        
         rbuffer(19) =  nrecvddm        
         rbuffer(20) =  epsv            
         rbuffer(21) =  epsva           
         rbuffer(22) =  epsvc           
         rbuffer(23) =  epsvac                ! start_out_light
         rbuffer(24) =  ninterne        
         rbuffer(25) =  ndirection      
         rbuffer(26) =  ns_solver       
         rbuffer(27) =  rho_ca          
         rbuffer(28) =  rlag            
         rbuffer(29) =  nbmax_ca        
         rbuffer(30) =  nbp_ca          
         rbuffer(31) =  nbmg            
         rbuffer(32) =  epsmg           
         rbuffer(33) =  epsmga          
         rbuffer(34) =  nb_prelis       
         rbuffer(35) =  nb_postlis      
         rbuffer(36) =  nb_cycle        
         rbuffer(37) =  sor_theta           
         rbuffer(38) =  is_decale       
         rbuffer(39) =  is_restart_save       ! end_out_light
      endif
   endif

   return
 end subroutine restore_param

 !************************************************************************
 
 !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
 subroutine restore_champs
   !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
   use data
   use drapeaux
   use champs
   implicit none
   
!!$   read (69) VTU0,VTU1,VTU2,VTU3,VTUS,VTV0,VTV1,VTV2,VTV3,VTVS
!!$   if (is_ns) read (69) PRE0,PRE1,PRE2,PRE3,PRES
!!$
   call save_or_retrieve('VTU0',VTU0,1,1,69,p_retrieve); call save_or_retrieve('VTV0',VTV0,1,1,69,p_retrieve); 
   call save_or_retrieve('VTU1',VTU1,1,1,69,p_retrieve); call save_or_retrieve('VTV1',VTV1,1,1,69,p_retrieve); 
   call save_or_retrieve('VTU2',VTU2,1,1,69,p_retrieve); call save_or_retrieve('VTV2',VTV2,1,1,69,p_retrieve); 
   call save_or_retrieve('VTU3',VTU3,1,1,69,p_retrieve); call save_or_retrieve('VTV3',VTV3,1,1,69,p_retrieve); 
   VTU=VTU0; VTV=VTV0

   if (is_ns) then
      call save_or_retrieve('PRE0',PRE0,1,1,69,p_retrieve); 
      call save_or_retrieve('PRE1',PRE1,1,1,69,p_retrieve); 
      call save_or_retrieve('PRE2',PRE2,1,1,69,p_retrieve); 
      call save_or_retrieve('PRE3',PRE3,1,1,69,p_retrieve); 
      PRE=PRE0
   endif

   if (my_task==0) close(69)

   print *,'*** data Restored... at t_start=',t_start,' calculation to ',t_all

   return
 end subroutine restore_champs                                     ! end_out_light

 !************************************************************************
 
 !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
 subroutine exacte
   !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    !SK
    !SK   Calcule la solution pour chaque point dans 'tab'
    !SK================================================================
    !SK      
    use data
    use disc
    use garbage
    use champs
    use para
    implicit none

    integer   nmu,lmu,nmv,lmv
    real(kind=prec) :: s1

    lmu=size(VTUS,1)-2; nmu=size(VTUS,2)-2;
    lmv=size(VTVS,1)-2; nmv=size(VTVS,2)-2;

    select case(nexample)
    case(1)       ! solution test 1D
       VTUS = sin(XU-vit_adv*temps)*exp(-nu*temps)
       VTVS = 0._prec
       PRES = 0._prec
       SMU  = 0._prec
       SMV  = 0._prec
       SMP  = 0._prec

       CLDU=0._prec
       CLDU(0,0:nmu+1)     = sin(xwdu-vit_adv*temps)*exp(-nu*temps) 
       CLDU(lmu+1,0:nmu+1) = sin(xedu-vit_adv*temps)*exp(-nu*temps) 

       CLDV=0._prec

    case(16)      ! solution test laplacien+temps Bipole 
       VTUS = 10*(temps-YU)*exp(-1/nu*((temps-XU)**2+(temps-YU)**2))
       VTVS = 10*(XV-temps)*exp(-1/nu*((temps-XV)**2+(temps-YV)**2))
       SMU = -10._prec*exp(-(2._prec*temps**2-2._prec*temps*XU+XU**2-2._prec*temps*YU+YU**2)/nu)/nu&
            & *(-nu+4._prec*temps**2-2._prec*temps*XU-6._prec*temps*YU+2._prec*YU*XU+2._prec*YU**2-8._prec&
            & *temps*nu+8._prec*temps**3-8._prec*temps**2*XU+4._prec*temps*XU**2-16._prec*temps**2*YU+12._prec&
            & *temps*YU**2+8._prec*YU*nu+8._prec*YU*temps*XU-4._prec*YU*XU**2-4._prec*YU**3)
      SMV = 10._prec*exp(-(2._prec*temps**2-2._prec*temps*XV+XV**2-2._prec*temps*YV+YV**2)/nu)/nu&
           & *(-nu+4._prec*temps**2-6._prec*temps*XV-2._prec*temps*YV+2._prec*XV**2+2._prec*YV*XV-8._prec&
           & *temps*nu+8._prec*temps**3-16._prec*temps**2*XV+12._prec*temps*XV**2-8._prec*temps**2*YV&
           & +4._prec*temps*YV**2+8._prec*XV*nu-4._prec*XV**3+8._prec*YV*temps*XV-4._prec*XV*YV**2)


       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = 10*(temps-YU(0,0:nmu+1))&
            & *exp(-1/nu*((temps-xwdu)**2+(temps-YU(0,0:nmu+1))**2))
       CLDU(lmu+1,0:nmu+1) = 10*(temps-YU(0,0:nmu+1))&
            & *exp(-1/nu*((temps-xedu)**2+(temps-YU(0,0:nmu+1))**2))

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     = 10*(XV(0:lmv+1,0)-temps)&
            & *exp(-1/nu*((temps-XV(0:lmv+1,0))**2+(temps-ysdv)**2))
       CLDV(0:lmv+1,nmv+1) = 10*(XV(0:lmv+1,0)-temps)&
            & *exp(-1/nu*((temps-XV(0:lmv+1,0))**2+(temps-yndv)**2))

    case(17)      ! solution test advection diffusion coef variable ->  Bipole 
       VTUS = 10*(t_start-YU)*exp(-1/nu*((t_start-XU)**2+(t_start-YU)**2))
       VTVS = 10*(XV-t_start)*exp(-1/nu*((t_start-XV)**2+(t_start-YV)**2))
       SMU = second_membre_u(XU,YU)
       SMV = second_membre_v(XV,YV)


       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = 10*(t_start-YU(0,0:nmu+1))&
            & *exp(-1/nu*((t_start-xwdu)**2+(t_start-YU(0,0:nmu+1))**2))
       CLDU(lmu+1,0:nmu+1) = 10*(t_start-YU(0,0:nmu+1))&
            & *exp(-1/nu*((t_start-xedu)**2+(t_start-YU(0,0:nmu+1))**2))

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     = 10*(XV(0:lmv+1,0)-t_start)&
            & *exp(-1/nu*((t_start-XV(0:lmv+1,0))**2+(t_start-ysdv)**2))
       CLDV(0:lmv+1,nmv+1) = 10*(XV(0:lmv+1,0)-t_start)&
            & *exp(-1/nu*((t_start-XV(0:lmv+1,0))**2+(t_start-yndv)**2))

    case(18)      ! solution test laplacien Bipole 
       VTUS = 10*(t_start-YU)*exp(-1/nu*((t_start-XU)**2+(t_start-YU)**2))
       VTVS = 10*(XV-t_start)*exp(-1/nu*((t_start-XV)**2+(t_start-YV)**2))
       SMU  = 40._prec/nu*(-t_start+YU)*exp(-(2._prec*t_start**2-2._prec*t_start*XU&
            & +XU**2-2._prec*t_start*YU+YU**2)/nu)*(-2._prec*nu+2._prec*t_start**2&
            & -2._prec*t_start*XU+XU**2-2._prec*t_start*YU+YU**2)
       SMV  =  -40._prec*nu*(XV-t_start)*exp(-(2._prec*t_start**2-2._prec*t_start*XV+XV**2&
            & -2._prec*t_start*YV+YV**2)/nu)*(-2._prec*nu+2._prec*t_start**2-2._prec&
            & *t_start*XV+XV**2-2._prec*t_start*YV+YV**2)/nu**2

       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = 10*(t_start-YU(0,0:nmu+1))&
            & *exp(-1/nu*((t_start-xwdu)**2+(t_start-YU(0,0:nmu+1))**2))
       CLDU(lmu+1,0:nmu+1) = 10*(t_start-YU(0,0:nmu+1))&
            & *exp(-1/nu*((t_start-xedu)**2+(t_start-YU(0,0:nmu+1))**2))

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     = 10*(XV(0:lmv+1,0)-t_start)&
            & *exp(-1/nu*((t_start-XV(0:lmv+1,0))**2+(t_start-ysdv)**2))
       CLDV(0:lmv+1,nmv+1) = 10*(XV(0:lmv+1,0)-t_start)&
            & *exp(-1/nu*((t_start-XV(0:lmv+1,0))**2+(t_start-yndv)**2))

    case(20)       ! solution test
       VTUS = XU**2
       VTVS = YV**2
       PRES = 0._prec
       SMU  = -nu*2._prec
       SMV  = -nu*2._prec
       SMP  = 0._prec

       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = xwdu**2
       CLDU(lmu+1,0:nmu+1) = xedu**2

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     = ysdv**2
       CLDV(0:lmv+1,nmv+1) = yndv**2


    case(22)       ! solution test
       VTUS = exp(2._prec*XU)*exp(2._prec*YU)
       VTVS = exp(2._prec*XV)*exp(2._prec*YV)
       PRES = 0._prec
       SMU  = -8._prec*nu*VTUS
       SMV  = -8._prec*nu*VTVS
       SMP  = 0._prec

       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = exp(2._prec*xwdu)*exp(2._prec*YU(0 ,0:nmu+1))
       CLDU(lmu+1,0:nmu+1) = exp(2._prec*xedu)*exp(2._prec*YU(lmu+1,0:nmu+1))

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)   = exp(2._prec*XV(0:lmv+1,0  ))*exp(2._prec*ysdv)
       CLDV(0:lmv+1,nmv+1) = exp(2._prec*XV(0:lmv+1,nmv+1))*exp(2._prec*yndv)


    case(23)        ! solution canal droit                  ! start_out_light
       VTUS = 4._prec*YU*(1._prec-YU)
       VTVS = 0._prec
       PRES = -8._prec*nu*XP
       SMU  = 0._prec
       SMV  = 0._prec
       SMP  = 0._prec

       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)    = 4._prec*YU(0  ,0:nmu+1)*(1._prec-YU(0 ,0:nmu+1))
       CLDU(lmu+1,0:nmu+1)  = 4._prec*YU(lmu+1,0:nmu+1)*(1._prec-YU(lmu+1,0:nmu+1))

       CLDV=0._prec

    case(19,21,25)     ! solution cavite entrainee
       VTUS = 0._prec; 
       VTVS = 0._prec
       PRES = 0._prec
       SMU  = 0._prec
       SMV  = 0._prec
       SMP  = 0._prec

       CLDU=0._prec
       CLDU(:,nmu+1) = ((1._prec-XU(:,nmu+1)**2)**2)
       if (nexample==21) CLDU(:,nmu+1)=1._prec
       if (nexample==19) CLDU(:,nmu+1)=-CLDU(:,nmu+1)
       
       CLDV=0._prec

    case(24)        ! solution Kovasznay
       s1   = 0.5_prec/nu-sqrt(0.25_prec/nu/nu+4._prec*pi*pi)
       VTUS = 1._prec-exp(s1*XU)*cos(2._prec*pi*YU)
       VTVS = s1/2._prec/pi*exp(s1*XV)*sin(2._prec*pi*YV)
       PRES = 0.5_prec*(1._prec-exp(2._prec*s1*XP))
       SMU  = 0._prec
       SMV  = 0._prec
       SMP  = 0._prec

       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = 1._prec-exp(s1*xwdu)*cos(2._prec*pi*YU(0  ,0:nmu+1))
       CLDU(lmu+1,0:nmu+1) = 1._prec-exp(s1*xedu)*cos(2._prec*pi*YU(lmu+1,0:nmu+1))


       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     = s1/2._prec/pi*exp(s1*XV(0:lmv+1,0))*sin(2._prec*pi*ysdv)
       CLDV(0:lmv+1,nmv+1) = s1/2._prec/pi*exp(s1*XV(0:lmv+1,nmv+1))*sin(2._prec*pi*yndv)

    case(26)       ! toubillon 2D de Green-Taylor

       VTUS = -cos(XU)*sin(YU)*exp(-2._prec*temps*nu)
       VTVS =  sin(XV)*cos(YV)*exp(-2._prec*temps*nu)
       PRES =  -.25_prec*(cos(2._prec*XP)+cos(2._prec*YP))*exp(-4._prec*temps*nu)
       SMU  = 0._prec
       SMV  = 0._prec
       SMP  = 0._prec   

       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = -cos(xwdu)*sin(YU(0 ,0:nmu+1))*exp(-2._prec*temps*nu)
       CLDU(lmu+1,0:nmu+1) = -cos(xedu)*sin(YU(lmu+1,0:nmu+1))*exp(-2._prec*temps*nu)

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     =  sin(XV(0:lmv+1,0  ))*cos(ysdv)*exp(-2._prec*temps*nu)
       CLDV(0:lmv+1,nmv+1) =  sin(XV(0:lmv+1,nmv+1))*cos(yndv)*exp(-2._prec*temps*nu)
                                               
    case(27)        ! solution Analytique polynomiale NS avec div/=0
       VTUS = XU**6-YU**6
       VTVS = XV**3*YV**3
       PRES = XP**4*YP**2
       SMU  = -30*nu*XU**4+30*nu*YU**4+6*XU**11-6*XU**5*YU**6 &
            & -6*XU**3*YU**8+4*XU**3*YU**2
       SMV  = -6*nu*XV*YV**3-6*nu*XV**3*YV+3*XV**8*YV**3 &
            & -3*XV**2*YV**9+3*XV**6*YV**5+2*XV**4*YV
!!$  en d(uu)/dx + d(vu)/dy
!!$       SMU  = -30*nu*XU**4+30*nu*YU**4+6*XU**11&
!!$            & -6*XU**5*YU**6-6*XU**3*YU**8+4*XU**3*YU**2
!!$       SMV  = -6*nu*XV*YV**3-6*nu*XV**3*YV+3*XV**8*YV**3&
!!$            & -3*XV**2*YV**9+3*XV**6*YV**5+2*XV**4*YV

       SMP  = 6._prec*XP**5+3._prec*XP**3*YP**2

       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = xwdu**6-YU(0,0:nmu+1)**6
       CLDU(lmu+1,0:nmu+1) = xedu**6-YU(lmu+1,0:nmu+1)**6

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     = XV(0:lmv+1,0)**3*ysdv**3
       CLDV(0:lmv+1,nmv+1) = XV(0:lmv+1,nmv+1)**3*yndv**3


    case(28)      ! Bipole se deplacant sur diagonale
       VTUS = 10*(temps-YU)*exp(-1/nu*((temps-XU)**2+(temps-YU)**2))
       VTVS = 10*(XV-temps)*exp(-1/nu*((temps-XV)**2+(temps-YV)**2))
       SMU  = second_membre_u(XU,YU)
       SMV  = second_membre_v(XV,YV)

       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = 10*(temps-YU(0,0:nmu+1))&
            & *exp(-1/nu*((temps-xwdu)**2+(temps-YU(0,0:nmu+1))**2))
       CLDU(lmu+1,0:nmu+1) = 10*(temps-YU(0,0:nmu+1))&
            & *exp(-1/nu*((temps-xedu)**2+(temps-YU(0,0:nmu+1))**2))

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     = 10*(XV(0:lmv+1,0)-temps)&
            & *exp(-1/nu*((temps-XV(0:lmv+1,0))**2+(temps-ysdv)**2))
       CLDV(0:lmv+1,nmv+1) = 10*(XV(0:lmv+1,0)-temps)&
            & *exp(-1/nu*((temps-XV(0:lmv+1,0))**2+(temps-yndv)**2))

    case(29)      ! Gaussienne se deplacant sur diagonale
       VTUS = exp(-1/nu*((temps-XU)**2+(temps-YU)**2))
       VTVS = sqrt(pi)*exp(-(temps-XV)**2/nu)*derf0((temps-YV)/sqrt(nu))&
            & *(temps-XV)/sqrt(nu)
       SMU  = second_membre_u(XU,YU)
       SMV  = second_membre_v(XV,YV)

       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = exp(-1/nu*((temps-xwdu)**2+(temps-YU(0,0:nmu+1))**2))
       CLDU(lmu+1,0:nmu+1) = exp(-1/nu*((temps-xedu)**2+(temps-YU(0,0:nmu+1))**2))

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     = sqrt(pi)*exp(-(temps-XV(0:lmv+1,0))**2/nu)&
            & *derf0((temps-ysdv)/sqrt(nu)) *(temps-XV(0:lmv+1,0))/sqrt(nu)
       CLDV(0:lmv+1,nmv+1) = sqrt(pi)*exp(-(temps-XV(0:lmv+1,0))**2/nu)&
            & *derf0((temps-yndv)/sqrt(nu)) *(temps-XV(0:lmv+1,0))/sqrt(nu)
    case(30)      ! Solution de Burgers stationaire
       VTUS = 2*nu*(2+25*cos(25*YU)*exp(25*XU-25)-25*cos(25*YU)*exp(-25*XU+25))&
            & /(1+2*XU+cos(25*YU)*exp(25*XU-25)+cos(25*YU)*exp(-25*XU+25))
       VTVS = -50*nu*(exp(25*XV-25)+exp(-25*XV+25))*sin(25*YV)/(1+2*XV+cos(25*YV)&
            & *exp(25*XV-25)+cos(25*YV)*exp(-25*XV+25))

       SMU  = 0._prec
       SMV  = 0._prec

       CLDU=0._prec
       CLDU(:,0)           = VTUS(:,0)
       CLDU(:,nmu+1)       = VTUS(:,nmu+1)
       CLDU(0,0:nmu+1)     = 2*nu*(2+25*cos(25*YU(0,0:nmu+1))*exp(25*xwdu-25)&
            & -25*cos(25*YU(0,0:nmu+1))*exp(-25*xwdu+25))&
            & /(1+2*xwdu+cos(25*YU(0,0:nmu+1))*exp(25*xwdu-25)&
            & +cos(25*YU(0,0:nmu+1))*exp(-25*xwdu+25))
       CLDU(lmu+1,0:nmu+1) = 2*nu*(2+25*cos(25*YU(0,0:nmu+1))*exp(25*xedu-25)&
            & -25*cos(25*YU(0,0:nmu+1))*exp(-25*xedu+25))&
            & /(1+2*xedu+cos(25*YU(0,0:nmu+1))*exp(25*xedu-25)&
            & +cos(25*YU(0,0:nmu+1))*exp(-25*xedu+25))

       CLDV=0._prec
       CLDV(0,:)           = VTVS(0,:)
       CLDV(lmv+1,:)       = VTVS(lmv+1,:)
       CLDV(0:lmv+1,0)     = -50*nu*(exp(25*XV(0:lmv+1,0)-25)+exp(-25*XV(0:lmv+1,0)+25))&
            & *sin(25*ysdv)/(1+2*XV(0:lmv+1,0)+cos(25*ysdv)*exp(25*XV(0:lmv+1,0)-25)&
            & +cos(25*ysdv)*exp(-25*XV(0:lmv+1,0)+25))
       CLDV(0:lmv+1,nmv+1) = -50*nu*(exp(25*XV(0:lmv+1,0)-25)+exp(-25*XV(0:lmv+1,0)+25))&
            & *sin(25*yndv)/(1+2*XV(0:lmv+1,0)+cos(25*yndv)*exp(25*XV(0:lmv+1,0)-25)&
            & +cos(25*yndv)*exp(-25*XV(0:lmv+1,0)+25))

    end select

    if (.not.is_north) then; VTUS(:,nmu+1)=0._prec; VTVS(:,nmv+1)=0._prec; 
                             SMU(:,nmu+1) =0._prec; SMV(:,nmv+1) =0._prec; endif
    if (.not.is_east ) then; VTUS(lmu+1,:)=0._prec; VTVS(lmv+1,:)=0._prec; 
                             SMU(lmu+1,:) =0._prec; SMV(lmv+1,:) =0._prec; endif
    if (.not.is_west ) then; VTUS(0,:)    =0._prec; VTVS(0,:)    =0._prec;
                              SMU(0,:)    =0._prec; SMV(0,:)     =0._prec; endif
    if (.not.is_south) then; VTUS(:,0)    =0._prec; VTVS(:,0)    =0._prec; 
                              SMU(:,0)    =0._prec; SMV(:,0)     =0._prec; endif

    return 
  end subroutine exacte

  !************************************************************************

end module init
