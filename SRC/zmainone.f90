!SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SK!!                                                             !! 
!SK!!                      ZEPHYR 0.6                             !! 
!SK!!                                                             !!
!SK!!                           Samuel KORTAS,   le 16/03/97      !!  
!SK!!                      (kortas@marius.univ-mrs.fr)            !!
!SK!!                                                             !! 
!SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SK
!SK   Credits :
!SK
!SK     Robert Cimrman 
!SK       Stencils a recouvrement differents en x et y  : Avril 97
!SK       Criteres de convergence differents par grille : Avril 97
!SK
!SK
!SK
!SK
!SK
!SK
!SK
!SK
!SK
!SK
!SK
!SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!SP2 module test
!SP2 end module test



 !SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SK  Programme Principal
!SK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program principal
  use constante
  use blas 
  use calcul 
  use champs 
  use compact
  use data   
  use disc   
  use drapeaux
  use garbage
  use init   
  use io     
  use mg     
  use navier 
  use para   
  use plot_flags
  use test
  use uncol
  implicit none

  call initialize_processors

  !*********************************
  !     call your_program here     *
  !*********************************
  call mainone
  !*********************************

  call release_processors

contains

  subroutine mainone
    use para
    use data
    use disc
    use champs
    implicit none

    integer          :: nitca                        ! out_light
    real(kind=prec) :: maxdiv=1.,diff                ! out_light

    call timer_clear(0); call timer_start(0); call timer_stop(0)

    !! timer communications send receive
    call timer_clear(1); call timer_start(1); call timer_stop(1)
    !! timer communications globales
    call timer_clear(2); call timer_start(2); call timer_stop(2)

    call ioinit
    call coinit

    call sortie_entete

    
    it=it_start
    temps=t_start

    
!!$    call test_df4

    call exacte
    call erreur(VTU0,VTV0,PRE0,VTUS,VTVS,PRES)
    call sortie_fichier(VTU0,VTV0,PRE0)                 ! out_light

    call timer_start(0)

    do while (temps<=t_all .and. (is_checkpoint_forced==0))

       ! timing one iteration
       call timer_clear(3); call timer_start(3);

       temps=temps+tau
       it=it+1

       select case(is_kuta)
       case (1)
       case (2)
          VTU  = 2._prec*VTU0-VTU1
          VTV  = 2._prec*VTV0-VTV1
       case (3)
          VTU  = 3._prec*VTU0-3._prec*VTU1+VTU2
          VTV  = 3._prec*VTV0-3._prec*VTV1+VTV2
       case (4)
          VTU  = 4._prec*VTU0-6._prec*VTU1+4._prec*VTU2-VTU3
          VTV  = 4._prec*VTV0-6._prec*VTV1+4._prec*VTV2-VTV3
       case default
          stop 'schema temporel non valide'
       end select

       if (is_oned) then
          SMU=vit_adv; VTV=0._prec; 
          call prepare_adv_grids(SMU,VTV)
       else if (nexample==17) then
          call exacte
          call prepare_adv_grids(VTUS,VTVS)
       else
          call prepare_adv_grids(VTU,VTV)
       endif


       call set_grid(-1)
       call set_switch(uv_switch*conv_code*diff_code/is_impl)
       call mavecxy(VTU,VTV,PRE,VTU4,VTV4,PRE)

       call exacte

       select case(is_kuta)
       case (1)
          VTU   =SMU+VTU0/tau
          VTV   =SMV+VTV0/tau
       case (2)
          VTU   =SMU+(4._prec*VTU0-VTU1)/2._prec/tau
          VTV   =SMV+(4._prec*VTV0-VTV1)/2._prec/tau
       case (3)
          VTU   =SMU+(18._prec*VTU0-9._prec*VTU1+2._prec*VTU2)/6._prec/tau
          VTV   =SMV+(18._prec*VTV0-9._prec*VTV1+2._prec*VTV2)/6._prec/tau
       case (4)
          VTU   =SMU+(48._prec*VTU0-36._prec*VTU1+16._prec*VTU2&
               &                    -3._prec*VTU3)/12._prec/tau 
          VTV   =SMV+(48._prec*VTV0-36._prec*VTV1+16._prec*VTV2&
               &                    -3._prec*VTV3)/12._prec/tau 
       case default
          stop 'schema temporel non valide'
       end select


       coeft =1._prec; call set_switch(uv_switch*t_code)
       call mavecxy(VTU,VTV,PRE,SMU,SMV,PRE)

       select case(nexample)
          case (16,18,20,22)   ! on ne rajoute pas de terme convectif!
          case default
          SMU=SMU-VTU4; SMV=SMV-VTV4
       end select

       call set_switch(global_switch*is_impl*t_code)

      select case(is_kuta)
       case (1)
          coeft = 1._prec/tau
          VTU   = VTU0
          VTV   = VTV0
          if (is_div) PRE   = PRE0
       case (2)
          coeft = 1.5_prec/tau
          VTU   = 2._prec*VTU0-VTU1
          VTV   = 2._prec*VTV0-VTV1
          if (is_div) PRE   = 2._prec*PRE0-PRE1
       case (3)
          coeft =11._prec/6._prec/tau
          VTU   = 3._prec*VTU0-3._prec*VTU1+VTU2
          VTV   = 3._prec*VTV0-3._prec*VTV1+VTV2
          if (is_div) PRE   = 3._prec*PRE0-3._prec*PRE1+PRE2
       case (4)
          coeft =25._prec/12._prec/tau
          VTU   = 4._prec*VTU0-6._prec*VTU1+4._prec*VTU2-VTU3
          VTV   = 4._prec*VTV0-6._prec*VTV1+4._prec*VTV2-VTV3
          if (is_div) PRE   = 4._prec*PRE0-6._prec*PRE1+4._prec*PRE2-PRE3
       case default
          stop 'schema temporel non valide'
       end select

       select case(is_richardson)
       case(0); VTU=1._prec;   VTV=1._prec
       case(1); VTU=VTU0;   VTV=VTV0; if (is_div) PRE=PRE0
       end select
       PRE=PRE0

       if (is_ns.and..not.is_div) then                              ! start_out_light
          nitca=1
          do while (nitca<=nbmax_ca)
             call ajout_grad    (PRE,SMU,SMV,VTU3,VTV3)

             call fixe_cl(CLDU,CLDV,PRE,VTU3,VTV3,SMP,flag=.true.)   
             call set_switch(global_switch*is_impl*t_code)
             call solve(VTU3,VTV3,SMP,VTU,VTV,PRE)

             call ajout_compress(PRE,VTU,VTV,rho_ca,maxdiv)


             if (mod(nitca,nbp_ca).eq.0) &
                  print '(A,I3,A,A,A,5E15.7)','++++++ ',nitca   &
                  ,' iterations de ',nom_ns_solver(ns_solver),'  div = ',maxdiv

             nitca=nitca+1
          enddo
       else                                                       ! end_out_light
          call fixe_cl(CLDU,CLDV,PRE,SMU,SMV,SMP,flag=.true.)     ! bords int = 0
          
!!$          call set_switch(global_switch*div_code)
!!$          call mavecxy(VTUS,VTVS,PRES,VTU3,VTV3,PRE3)
!!$          
!!$          call outdon(VTUS,'VTUS',1,1)
!!$          call outdon(VTVS,'VTVS',1,1)
!!$          call outdon(PRES,'PRES',1,1)
!!$          call outdon(VTU3,'VTU3')
!!$
!!$          call outdon(VTU3,'VTU3',1,1)
!!$          call outdon(VTV3,'VTV3',1,1)
!!$          call outdon(PRE3,'PRE3',1,1)
!!$
!!$          call outdon(PRE3,'PRE3',1,1)
!!$
!!$
!!$          diff=max(global_maxval(abs(VTU2)),diff)
!!$          print *,'erreur ',diff
!!$
!!$          if (my_task==0) then
!!$             print *,'*** fin ?  **** '
!!$             coeft=0
!!$             do while (coeft==0)
!!$             enddo
!!$          endif
!!$
!!$          call  parallel_stop

!!$          VTU=VTUS; VTV=VTVS; PRE=PRES
!!$          call outdon(vtus,'vtus')
!!$          call outdon(smu,'smu')

          call mavecxy(VTU,VTV,PRE,cldu,cldv,pre)
!!$          call outdon(cldu,'cldu')          
!!!          stop

          call solve(SMU,SMV,SMP,VTU,VTV,PRE)
       endif                                                      ! out_light

       if (is_ns) PRE=PRE+global_minval(PRES)-global_minval(PRE)  ! out_light  
       
       VTU3=VTU2;   VTV3=VTV2;  if (is_div) PRE3=PRE2;
       VTU2=VTU1;   VTV2=VTV1;  if (is_div) PRE2=PRE1;
       VTU1=VTU0;   VTV1=VTV0;  if (is_div) PRE1=PRE0;
       VTU0=VTU;    VTV0=VTV;   if (is_div) PRE0=PRE; 

       ! timing one iteration
       call timer_stop(3)       


       call erreur(VTU,VTV,PRE,VTUS,VTVS,PRES)
       call sortie_fichier(VTU,VTV,PRE)                           ! start_out_light

       call infothemis

       !print *,my_task,'restart_save,it,temps',is_restart_save,it,temps
       if (is_restart_save/=0.or.is_checkpoint_forced/=0) then
          if (is_checkpoint_forced/=0&
               & .or.mod(it,abs(is_restart_save))==0&
               & .or.temps+0.9*tau>t_all)   then
             print *, my_task,'ici avant save'
             call flush(6)
             call save_data
             endif
       endif                                                      ! end_out_light

       if (is_ns) then                                           ! start_out_light
          SMU=VTU-VTU1; SMV=VTV-VTV1; SMP=PRE-PRE1
          call fixe_cl(SMU,SMV,SMP,SMU,SMV,SMP,flag=.true.)   ! met les bords internes a 0
!          diff=global_ddot(SMU,SMV,SMP,SMU,SMV,SMP)
          diff=global_maxval(abs(SMU))
          diff=max(global_maxval(abs(SMV)),diff)
          diff=max(global_maxval(abs(SMP)),diff)
          print *,'diff',diff
          call sortie_champ(SMU,'a',surface)
          call sortie_champ(SMV,'b',surface)
          call sortie_champ(SMP,'c',surface)
          if (is_div.and.diff<=1.E-6) then
             print *,"(is_div.and.diff<=1.E-6)"
             exit
          endif
          if (.not.is_div.and.maxdiv<=5.E-4) then
             print *,"(.not.is_div.and.maxdiv<=5.E-4)"
             exit
          endif
       endif                                                     ! end_out_light
       
    enddo


1000 continue

    PRE=0._prec


!!$    call fixe_cl(VTUS,VTVS,PRES,VTU,VTV,PRE,.true.)
!!$    print *,'max erreur   u  ',global_maxval(VTU-VTUS)
!!$    print *,'erreur eucl u  ',sqrt(global_sum((VTU-VTUS)*(VTU-VTUS)))
!!$    print *,'max errevr   v  ',global_maxval(VTV-VTVS)
!!$    print *,'errevr evcl v  ',sqrt(global_sum((VTV-VTVS)*(VTV-VTVS)))
!!$
!!$    call sortie_champ(VTU-VTUS,'y',surface)
!!$    call sortie_champ(VTU,'Y',surface)
!!$    call sortie_champ(VTV,'Z',surface)
!!$    call sortie_champ(VTV-VTVS,'z',surface)
!!$    call outdon(VTUS-VTU,'diff_u')
!!$    call outdon(VTVS-VTV,'diff_v')



    call coend

    call ioend


!SP2  return

!!$    if (my_task==0) then
!!$       print *,'*** fin ? **** '
!!$       coeft=0
!!$       do while (coeft==0)
!!$       enddo
!!$    endif

  end subroutine mainone

end program principal


