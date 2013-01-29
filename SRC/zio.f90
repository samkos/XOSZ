module io

  use constante
  use para
  use garbage, only : type_machine
  use plot_flags
 
  character(len=100) :: nom_solution, nom_solver,nom_precond

  character(len=*), parameter, dimension(0:2) ::       &
       nom_guess    = (/ '1         '    ,'Current   '    ,'Richardson' /) &
      ,nom_impl     = (/ 'implicitement' ,'explicitement' ,'explicitement' /)

  character(len=*), parameter, dimension(0:1) ::       &           ! start_out_light
       nom_maillage = (/ 'grille B  ', 'grille MAC' /)
  character(len=*), parameter, dimension(1:4) ::       &
       nom_ns_solver = (/ 'Compressibilite Artificielle ' &
       ,'Resolution Couplee           ' &
       ,'Methode de Projection        ' &
       ,'Lagrangien Augmente          '/)      
  character(len=1), parameter, dimension(4) ::&
       & nom_mg_cycle=(/'V','W','X','Y'/)                                  ! end_out_light

  integer, parameter :: ncs2fw=88
  
  
contains

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine sortie_entete
    use data
    use drapeaux
    use disc
    use uncol
    use champs, only : VTUS,VTVS,PRES, vit_adv
    implicit none
!!$    integer, parameter : task_code=0,lmu_code=1,nmu_code=2&
!!$         & ,lmv_code=1,nmv_code=2, lm_code=1,nmu_code=2 &
    integer :: i,k,lm,nm,lmu,nmu,lmv,nmv,lmp,nmp
    integer, dimension(0:2*nb_tasks-1) ::&
         & all_lmu,all_nmu,all_lmv,all_nmv&
         & ,all_lmp,all_nmp,all_nb_grid &
         & ,all_ddm_lmu,all_ddm_nmu,all_ddm_lmv,all_ddm_nmv&
         & ,all_ddm_lmp,all_ddm_nmp


    select case(ntype_solver)
       case (1)      ; write (nom_solver,'(A,I1)')        'BiCGSTAB_',ntype_solver
       case (2:9)    ; write (nom_solver,'(A,I1)')        'BiCGSTAB_',ntype_solver        
       case (10)     ; write (nom_solver,'(A,I1)')        'GMRES_',ninterne               ! start_out_light 
       case (11)     ; write (nom_solver,'(A,I1)')        'GMRESR_rec_',ninterne
       case (12)     ; write (nom_solver,'(A,I1)')        'GMRESR_norec_',ninterne
       case (13)     ; write (nom_solver,'(A,A1,"_",I1,"/",I1)')&                                   ! out_light
            &                            'MG_',nom_mg_cycle(nb_cycle),nb_prelis,nb_postlis          ! out_light 
       case (-13)    ; write (nom_solver,'(A,A1,"_",I1,"/",I1)')&                                   ! out_light
            &                            'FMG_',nom_mg_cycle(nb_cycle),nb_prelis,nb_postlis          ! out_light 
       case (14)     ; write (nom_solver,'(A,I1)')        'ORHTHODIR_',ninterne
       case (15)     ; write (nom_solver,'(A,I1)')        'ORTHOMIN_',ninterne
       case (16)     ; write (nom_solver,'(A,I1)')        'lisse'                         ! end_out_light   
       case (17)     ; write (nom_solver,'(A,I1)')        'Schwarz_add' 
       case (18)     ; write (nom_solver,'(A,I1)')        'Schwarz_mult' 
    end select

    select case(ntype_precond)
       case (1)      ; write (nom_precond,'(A,I1)')        'BiCGSTAB_',ntype_precond
       case (2:9)    ; write (nom_precond,'(A,I1)')        'BiCGSTAB_',ntype_precond                 
       case (10)     ; write (nom_precond,'(A,I1)')        'GMRES_',ninterne                         ! out_light
       case (13)     ; write (nom_precond,'(A,A1,"_",I1,"/",I1)')&                                   ! out_light
            &                             'MG_',nom_mg_cycle(nb_cycle),nb_prelis,nb_postlis          ! out_light 
       case (-13)    ; write (nom_precond,'(A,A1,"_",I1,"/",I1)')&                                   ! out_light
            &                             'FMG_',nom_mg_cycle(nb_cycle),nb_prelis,nb_postlis          ! out_light 
       case (16)     ; write (nom_precond,'(A,I1)')        'lisse'                                   ! out_light 
    end select

    select case(nexample)
       case(1)      ; nom_solution='Solution Test 1D       '
       case(19)     ; nom_solution='Cavite entrainee inv   '
       case(16)     ; nom_solution='Solution Test temps+diffusion Bipole'
       case(17)     ; nom_solution='Solution Test advection diffusion Bipole'
       case(18)     ; nom_solution='Solution Test Poisson Bipole'
       case(20)     ; nom_solution='Solution Test Poisson 1'
       case(21)     ; nom_solution='Cavite entrainee       '  ! out_light 
       case(22)     ; nom_solution='Solution Test Poisson 2'
       case(23)     ; nom_solution='Ecoulement de Poiseuil '  ! start_out_light 
       case(24)     ; nom_solution='Fluide de Kovasynay    '
       case(25)     ; nom_solution='Cavite entrainee reg   '
       case(26)     ; nom_solution='Taylor-Green           '
       case(27)     ; nom_solution='Solution Polynomiale NS'  ! end_out_light 
       case(28)     ; nom_solution='Bipole                 '
       case(29)     ; nom_solution='Gaussienne             '
       case(30)     ; nom_solution='Burgers 2D stationaire '
       case default ; stop 'no example defined'
    end select
    
    i=0
    if (is_precond) then
       lm=size(VTUS,1)-2
       nm=size(VTVS,1)-2
       if (((is_east.or.is_west)  .and.nrecvddmx>lm)&
            &.or.((is_north.or.is_south).and.nrecvddmy>nm)) i=my_task
       i=global_max(i)
       if (i>0) then 
          print *, 'Recouvrement Superieur a lm/nm du sous-domaine ',i
          call parallel_stop
       endif
    endif

    if (my_task.eq.master_task) then
       print *,'***********************************'
       print *,'*                                 *'
       print *,'*  Zephyr v 0.6 (Janvier 2013)    *'
       print *,'*                                 *'
       print *,'***********************************'
       print *,'>> '
       print *,'>> Execution sur : ', type_machine 
       print *,'>> '
       print *,'>> Solution test :   ',trim(nom_solution)
       print *,'>>'
       print *,'>> Viscosite     :   ',nu
       print *,'>>'
       print *,'>> tau         = ',tau
       print *,'>> dx          = ',dx,'       ou lm = ',lm_global
       print *,'>> dy          = ',dy,'       ou nm = ',nm_global
       if (is_oned) then
          print *,'>> CFL     nb  = ',tau/dx/dx*nu
          print *,'>> Courant nb  = ',tau/dx*vit_adv
       endif
       print *,'>>'
       print *,'>> Maillage     --> ',nom_maillage(is_decale)          ! out_light
       print *,'>>'                                                    ! out_light 
       print *,'>> Termes diffusifs  traites ',nom_impl(mod(is_impl,2))
       print *,'>> Termes convectifs traites ',nom_impl(mod(is_impl,3))
       print *,'>> Guess Initial --> ',nom_guess(is_richardson)
       print *,'>>'
       print *,'>> integration temporelle d''ordre  : ' ,is_kuta
       print *,'>>'
       print '(1X,6A)','>> Solveur                         : '&
            & ,trim(nom_solver),'(',trim(nom_disc(ncheck)),')'
       if (is_precond) &
            print '(1X,6A)','>> Preconditione par               : '&
            & ,trim(nom_precond),'(',trim(nom_disc(ncheck_precond)),')'
       if (is_precond) then
          print '(1X,A,I2)','>> Demi-recouvrement delta x       : '&
               & ,nrecvddmx
          print '(1X,A,I2)','>>                   delta y       : '&
               & ,nrecvddmy
       endif
       print *,'>>'
       if (is_ns)  &                                      ! start_out_light      
            print *,'>> Methode Navier-Stokes           : ',nom_ns_solver(ns_solver)
       print *,'>> '                                      ! end_out_light
       print *,'>> Convergence externe             : ',grid(-1)%epsv
       print *,'>> Convergence externe abs         : ',grid(-1)%epsva
       if (is_precond) then
          print *,'>> Convergence interne             : ',grid(ngrid_max)%epsv
          print *,'>> Convergence interne abs         : ',grid(ngrid_max)%epsva
       endif

       if (abs(is_restart_save)>0) then                    ! start_out_light
          print *,'>> '                                      
          print '(1X,A,A)'&
               & ,'>> Sauvegarde globale dans         : ',trim(nom_fic_save)
          print *,'>> Tous les <n> pas de temps       : ',is_restart_save
       endif                                               ! end_out_light

    endif

    all_lmu=0; all_lmv=0; all_lmp=0; all_nb_grid=0
    all_nmu=0; all_nmv=0; all_nmp=0; 

    all_ddm_lmu=0; all_ddm_lmv=0; all_ddm_lmp=0; 
    all_ddm_nmu=0; all_ddm_nmv=0; all_ddm_nmp=0; 

    lmu=size(VTUS,1)-2; nmu=size(VTUS,2)-2
    lmv=size(VTVS,1)-2; nmv=size(VTVS,2)-2
    lmp=size(PRES,1)-2; nmp=size(PRES,2)-2


    all_lmu(my_task)=lmu;     all_nmu(my_task)=nmu
    all_lmv(my_task)=lmv;     all_nmv(my_task)=nmv
    all_lmp(my_task)=lmp;     all_nmp(my_task)=nmp 
    all_nb_grid(my_task)  =ngrid_max;


    all_nb_grid=global_add(all_nb_grid);
    all_lmu=global_add(all_lmu);    all_nmu=global_add(all_nmu);
    all_lmv=global_add(all_lmv);    all_nmv=global_add(all_nmv);
    all_lmp=global_add(all_lmp);    all_nmp=global_add(all_nmp);

    if (ngrid_max.ge.0) then
       lmu=size(grid(ngrid_max)%U1,1)-2; nmu=size(grid(ngrid_max)%U1,2)-2
       lmv=size(grid(ngrid_max)%V1,1)-2; nmv=size(grid(ngrid_max)%V1,2)-2
       lmp=size(grid(ngrid_max)%P1,1)-2; nmp=size(grid(ngrid_max)%P1,2)-2
       

       all_ddm_lmu(my_task)=iumgfin-iumgdeb;     all_ddm_nmu(my_task)=kumgfin-kumgdeb
       all_ddm_lmv(my_task)=ivmgfin-ivmgdeb;     all_ddm_nmv(my_task)=kvmgfin-kvmgdeb
       all_ddm_lmp(my_task)=ipmgfin-ipmgdeb;     all_ddm_nmp(my_task)=kpmgfin-kpmgdeb

       all_ddm_lmu=global_add(all_ddm_lmu);    all_ddm_nmu=global_add(all_ddm_nmu);
       all_ddm_lmv=global_add(all_ddm_lmv);    all_ddm_nmv=global_add(all_ddm_nmv);
       all_ddm_lmp=global_add(all_ddm_lmp);    all_ddm_nmp=global_add(all_ddm_nmp);
    endif


    if (my_task.eq.master_task) then
       print *,'>>'
       write(*,'(" >> ",100A1)') ('=',i=0,78)
       print *,'>>'
       write(*,'(A)',advance="no") " >> task/nbgrid : "
       do k=nb_k_blocks-1,0,-1
          write (*,'(16I6)',advance="no") (i,i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          write (*,'(A)',advance="no")    "           "
          write (*,'(16I6)',advance="no") (all_nb_grid(i)+1,i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          if (k>0) then; print *; write(*,'(A)',advance="no") " >>               "; endif
       enddo
       print *; print *,'>>'
       write(*,'(A)',advance="no") " >> lmu/lmu_ddm : "
       do k=nb_k_blocks-1,0,-1
          write (*,'(16I6)',advance="no") (all_lmu(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          write (*,'(A)',advance="no")    "           "
          write (*,'(16I6)',advance="no") (all_ddm_lmu(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          if (k>0) then; print *; write(*,'(A)',advance="no") " >>               "; endif
       enddo
       print *; print *,'>>'
       write(*,'(A)',advance="no") " >> nmu/nmu_ddm : "
       do k=nb_k_blocks-1,0,-1
          write (*,'(16I6)',advance="no") (all_nmu(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          write (*,'(A)',advance="no")    "           "
          write (*,'(16I6)',advance="no") (all_ddm_nmu(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          if (k>0) then; print *; write(*,'(A)',advance="no") " >>               "; endif
       enddo
       print *; print *,'>>'
       print *,'>>'
       write(*,'(A)',advance="no") " >> lmv/lmv_ddm : "
       do k=nb_k_blocks-1,0,-1
          write (*,'(16I6)',advance="no") (all_lmv(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          write (*,'(A)',advance="no")    "           "
          write (*,'(16I6)',advance="no") (all_ddm_lmv(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          if (k>0) then; print *; write(*,'(A)',advance="no") " >>               "; endif
       enddo
       print *; print *,'>>'
       write(*,'(A)',advance="no") " >> nmv/nmv_ddm : "
       do k=nb_k_blocks-1,0,-1
          write (*,'(16I6)',advance="no") (all_nmv(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          write (*,'(A)',advance="no")    "           "
          write (*,'(16I6)',advance="no") (all_ddm_nmv(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
          if (k>0) then; print *; write(*,'(A)',advance="no") " >>               "; endif
       enddo
       if (is_ns) then
          print *; print *,'>>'
          print *,'>>'
          write(*,'(A)',advance="no") " >> lmp/lmp_ddm : "
          do k=nb_k_blocks-1,0,-1
             write (*,'(16I6)',advance="no") (all_lmp(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
             write (*,'(A)',advance="no")    "           "
             write (*,'(16I6)',advance="no") (all_ddm_lmp(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
             if (k>0) then; print *; write(*,'(A)',advance="no") " >>               "; endif
          enddo
          print *; print *,'>>'
          write(*,'(A)',advance="no") " >> nmp/nmp_ddm : "
          do k=nb_k_blocks-1,0,-1
             write (*,'(16I6)',advance="no") (all_nmp(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
             write (*,'(A)',advance="no")    "           "
             write (*,'(16I6)',advance="no") (all_ddm_nmp(i),i=k*nb_i_blocks,(k+1)*nb_i_blocks-1)
             if (k>0) then; print *; write(*,'(A)',advance="no") " >>               "; endif
          enddo
       endif
       print *; print *,'>>'
       write(*,'(" >> ",100A1)') ('=',i=0,78)
       print *,'>>'
    endif

!!    call parallel_stop

    return
  end subroutine sortie_entete


  !************************************************************************

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine erreur (VTU,VTV,PRE,VTUS,VTVS,PRES)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use disc
    use data
    use drapeaux
    use para
    use blas
    use init
    use conv
    use uncol
    implicit none

    real(kind=prec), dimension(:,:) :: VTU,VTV,PRE,VTUS,VTVS,PRES
    real(kind=prec), dimension(size(VTU,1),size(VTU,2)) :: ERRSX
    real(kind=prec), dimension(size(VTV,1),size(VTV,2)) :: ERRSY
    real(kind=prec), dimension(size(PRE,1),size(PRE,2)) :: ERRSP

    real(kind=prec) ::  err,errx,erry, errp, solmax,solumax,solwmax, solpmax &
         , solmax2,solumax2,solwmax2,solpmax2,err2,errx2,erry2,errp2, div
    integer :: n


    call timer_stop(0)

    if (it-it_start<=2.or.mod(it,nt_print)==0) then

!       call fixe_cl(VTU,VTV,PRE,VTU,VTV,PRE,.true.)  ! met les ghost-cells  a 0
       
       solumax=global_maxval(abs(VTU))
       solwmax=global_maxval(abs(VTV))
       if (is_ns) solpmax=global_maxval(abs(PRE))
       solmax =max(solumax,solwmax)

       solumax2 = global_sum(VTU*VTU)
       solwmax2 = global_sum(VTV*VTV)
       if (is_ns) solpmax2=global_maxval(PRE*PRE)
       solmax2  = solumax2+solwmax2
       
       ERRSX=abs(VTUS-VTU)
       ERRSY=abs(VTVS-VTV)
       if (is_ns) errp=global_minval(PRES)-global_minval(PRE)
       if (is_ns) ERRSP=abs(PRES-PRE-errp)
       call fixe_cl(ERRSX,ERRSY,ERRSP,ERRSX,ERRSY,ERRSP,flag=.true.)

       if (is_decal) then                                   ! start_out_light
          if (is_west)  ERRSX(1,:)=0._prec
          if (is_east)  ERRSX(size(ERRSX,1),:)=0._prec
          if (is_south) ERRSY(:,1)=0._prec
          if (is_north) ERRSY(:,size(ERRSY,2))=0._prec
!!$          ERRSP(1,1)=0.
!!$          ERRSP(size(ERRSP,1),1)=0.
!!$          ERRSP(1,size(ERRSP,2))=0.
!!$          ERRSP(size(ERRSP,1),size(ERRSP,2))=0.
       endif                                                ! end_out_light

       errx2=global_sum(ERRSX*ERRSX)
       erry2=global_sum(ERRSY*ERRSY)
       if (is_ns) errp2=global_sum(ERRSP*ERRSP)
       err2 =errx2+erry2

       errx=global_maxval(ERRSX)
       erry=global_maxval(ERRSY)
       if (is_ns) errp=global_maxval(ERRSP)                    ! out_light
       err=max(errx,erry)

       solmax2=sqrt(solmax2*dx*dy)
       solumax2=sqrt(solumax2*dx*dy)
       solwmax2=sqrt(solwmax2*dx*dy)
       if (is_ns) solpmax2=sqrt(solpmax2*dx*dy)                ! out_light
       err2=sqrt(err2*dx*dy)
       errx2=sqrt(errx2*dx*dy)
       erry2=sqrt(erry2*dx*dy)
       if (is_ns) errp2=sqrt(errp2*dx*dy)


       if (is_ns) then                                        ! start_out_light     
          call set_grid(-1)
          call set_switch(div_code*is_impl*t_code)
          call mavecxy(VTU,VTV,PRE,ERRSX,ERRSY,ERRSP)
          div=global_maxval(abs(ERRSP))
       endif

       !SK      call calcul_div_df(VTU,VTV,maxdivdf,maxdivdf2,lm,nm)
       !SK      call calcul_div_vf(VTU,VTV,maxdivvf,maxdivvf2,lm,nm)   ! end_out_light

       if (my_task==master_task) then
          print *
          print '(A,A,A)','============================='       &
                      ,'========================='              &
                      ,'============================='
          print *,'                t=',temps,' (it=',it,') '   
          print '(A,A,A)','=============================='      &
                      ,'========================'               &
                      ,'============================='

          print 123,'Erreur Linfini      -->  err    = '        &
                      ,err,err/max(epsilon,solmax)              &
                      ,err2,err2/max(epsilon,solmax2)
          print 123,'                    -->  errvtu = '        &
                      ,errx,errx/max(epsilon,solumax)           &
                      ,errx2,errx2/max(epsilon,solumax2)
          print 123,'                    -->  errvtv = '        &
                      ,erry,erry/max(epsilon,solwmax)           &
                      ,erry2,erry2/max(epsilon,solwmax2)
          if (is_ns) &                                              ! start_out_light
          print 123,'                    -->  errpre = '        &
                      ,errp,errp/max(epsilon,solpmax)           &
                      ,errp2,errp2/max(epsilon,solpmax2)
          if (is_ns) &
          print 123,'                    -->  Diverg = '        &
                      ,div                                          ! end_out_light
          print 123,'Reynolds de Maille  -->  Rh     = '        &    
                      ,max(solumax*dx/nu,solwmax*dy/nu)
          print 123,'                    -->  Rh(x)  = '        &
                      ,solumax*dx/nu
          print 123,'                    -->  Rh(y)  = '        &
                      ,solwmax*dy/nu

          print '(A,A,A)','============================='       &
                      ,'========================='              &
                      ,'============================='

          if (is_cv_global.and.any(flag_niv)) then
             print *, '>>                                  CONVERGENCE'
             print '(A,A,A)','=============================='      &
                  ,'========================'               &
                  ,'============================='
             do n=0,niv_lim
                if (flag_niv(n)) then
                   print '(A3,A16,i2,a10,a3,a5,a3,4a10)',' >>'&
                        & ,'niveau ',n,' : total',' / ','nb_ap',' = ' &
                        ,' moy','max','min','last'
                   print '(A3,A21,i7,a3,i5,a3,4i10)',' >>'&
                        & ,'iteration : ',total_iter(n),' / ',nb_iter(n),' = ' &
                        ,total_iter(n)/max(nb_iter(n),1)&
                        & ,max_iter(n),min_iter(n),last_iter(n)
                   
                   print '(A,A,A)','------------------------------'      &
                        ,'------------------------'               &
                        ,'-----------------------------'
                   
                   print '(A3,A20,6A9)',' >>','         :'&
                        & ,'moy','max','min','last','resf','resl'
                   
                   print '(A3,A21,1X,6E9.2)',' >>','rayon cv  : '&
                        & ,total_rho(n)/max(nb_iter(n),1),max_rho(n) &
                        ,min_rho(n),last_rho(n),first_res(n),last_res(n)
                   
                   print '(A,A,A)','=============================='      &
                        ,'========================'               &
                        ,'============================='
                endif
             enddo
           endif
           
           print *


       endif


       call timer_print("reduction",2)
       call timer_print('messages',1)
       call timer_print('total',0)

    endif
    
    

    call timer_start(0)



123 format (A35,10E12.4)

    return
  end subroutine erreur

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine infothemis
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

    use uncol 
    use disc
    use data
    implicit none

    if (my_task==master_task) then

       open(file="CHECKPOINT_NOW", unit=92, status='OLD', action='READ', err=900)
       close(92)
       is_checkpoint_forced = 10

900 continue

       !write(*,3012) it,t_start,t_all,tau,1+int((t_all-t_start)/tau),timer_get(3)   ! THEMIS
       write(ncs2fw,3011) it,1+int((t_all-t_start)/tau),timer_get(3)   ! THEMIS


       call flush(ncs2fw)

    endif
!THEMIS
 3011 format(' CPU TIME FOR THE TIME STEP  ',I7,' FROM ',I10,' : ',E14.5)
 3012 format(' CPU TIME FOR THE TIME STEP  ',I7,' t_start,t_all,tau ',3E14.5,' FROM ',I10,' : ',E14.5)
    return
 3023 format(' On the order of THEMIS Framework...')

 3020 format(/,/,                                                 &
 ' Sortie intermediaire de fichiers suite',/,                     &
 '   Sauvegarde a l''iteration ', I10, ', Temps physique ',E14.5,/,/)

  end subroutine infothemis

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine ioinit
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

    open (file="THEMIS_Code2Themis", unit=ncs2fw, form='formatted', status='unknown', &
            & err=900)
900 continue
  
    is_checkpoint_forced = 0

    return
  end subroutine ioinit

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine ioend
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

    use data
    implicit none

    if (is_checkpoint_forced/=0)   write(ncs2fw,3024)

3024 format(                                                     &
          ' Code stops now after having succesfully checkpointed ',/,       &
          'on the order of THEMIS Framework')

    close(ncs2fw)



  end subroutine ioend
  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine read_plot_flags
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

    use data,  only : nom_fic_output
    integer :: ok

    return

    do
       open(file=trim(nom_fic_output)//'Z.dat',unit=79,iostat=ok)
       if (ok==0) exit
    enddo
    read(79,*) controle_string
    close(79)
    is_plot_psi=index(controle_string,'h',.false.)>0
    is_plot_w  =index(controle_string,'w',.false.)>0

    return
  end subroutine read_plot_flags

  !************************************************************************

  !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine sortie_champ(VTU,string,flag)
    !SKssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use drapeaux
    use para,  only : is_mpp
    use data,  only : nom_fic_output
    use uncol, only : timer_stop,timer_start
    implicit none

    real(kind=prec), dimension(:,:)              :: VTU
    character(len=*) :: string
    logical, optional :: flag
    real(kind=prec) :: valmin,valmax

    integer, save :: nb_plot=0
    integer :: splot




    nb_plot=1+nb_plot

    call read_plot_flags

     if (index(controle_string,string,.false.)<=0) then
       call plot( string ,1,100,valmin,valmax)       
       return
    endif

    valmax=global_maxval(VTU)
    valmin=global_minval(VTU)
    if (valmax-valmin<1.e-8) then
       print *,'champs ',string,' petit everywhere'
       return
    endif

    splot=2
    if (present(flag).and..not.flag) splot=1
    if (present(flag).and.flag) splot=10

    call timer_stop(0)

    if (my_task==0) open(file=trim(nom_fic_output)//trim(string)//'.dat',unit=79)
    call save_or_retrieve('VTU',VTU,1,1,79,p_save_txt)
    if (my_task==0) close(79)
    if (my_task==0) call plot( &
         & trim(nom_fic_output)//trim(string)//'.dat' &
         &,len_trim(trim(nom_fic_output)//trim(string)//'.dat'),splot&
         &,valmin,valmax)

    call timer_start(0)

    return
  end subroutine sortie_champ


  !************************************************************************



end module io
