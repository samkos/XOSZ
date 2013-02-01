module selectprec
  integer, parameter :: prec=kind(1.d0)
end module selectprec


module constante
  use selectprec

  integer, parameter :: df2=2,vf2=1,df4=5,vf4=3,vf4e=4,df4e=6,mg_solver=13

  character(len=*), parameter, dimension(1:7) ::       &
       nom_disc     = (/ 'VF2','DF2','VF4','V4e','DF4','D4e','CHB'/)

  real(kind=prec), parameter :: pi=3.141592653589793_prec&
       &                       ,epsilon=(0.1_prec)**16

end module constante

module data
  
  use selectprec

  integer, save :: is_print,is_kuta,nexample,is_unsteady                   &
       ,is_richardson,is_restart_save,is_decale                            &
       ,lm_global,nm_global,ncheck,ncheck_precond,ns_solver,ntype_solver   &
       ,ntype_precond,nrecvddm,nrecvddmx,nrecvddmy,nprecond,npcg,ninterne                      &
       ,ndirection,nbmax_ca,nbp_ca,nbmg,nb_prelis,nb_postlis,nb_cycle,is_cv&
       &, is_checkpoint_forced

  real(kind=prec)  :: t_start,t_all,tau,t_print                       &
       ,nu,epsv,epsvc,epsva,epsvac,epsmg,epsmga                            &
       ,rlag,rho_ca,tol_ca,sor_theta ,xmin,ymin,xmax,ymax,un_ms_sor_theta  

  character(len=96),save :: nom_fic_output,nom_file_cv,nom_fic_save

end module data

module disc
  use selectprec
  real(kind=prec),  save      :: xdl,ydl                &
       ,dx,dy,invdx,invdy,invdx2,invdy2,p5invdx,p5invdy  &
       ,temps,coeft,xudl,xvdl,xpdl,zudl,zvdl,zpdl

  integer, save ::  it,it_start,nt_print &
       ,iudeb, iufin, kudeb, kufin &
       ,iumgdeb, iumgfin, kumgdeb, kumgfin &
       ,ivdeb, ivfin, kvdeb, kvfin &
       ,ivmgdeb, ivmgfin, kvmgdeb, kvmgfin &
       ,ipdeb, ipfin, kpdeb, kpfin &
       ,ipmgdeb, ipmgfin, kpmgdeb, kpmgfin 

  integer, save :: lmu_global,nmu_global, lmv_global,nmv_global &
       ,lmp_global,nmp_global

  real(kind=prec), save       :: xwdu,xedu,ysdv,yndv

  real(kind=prec), dimension (:,:), allocatable, save :: &
       XU,YU,XV,YV,XP,YP

  integer, save :: ngrid_max

  type grid_const
     real(kind=prec) :: &
          dx,dy,invdx,invdy,invdx2,invdy2,p5invdx,p5invdy&
          & ,epsv, epsva
     real(kind=prec), dimension (:,:), pointer :: &  
          U_EN_U,   V_EN_U,   U_EN_V,   V_EN_V    &
         ,fois_U,   plus_U,   fois_V,   plus_V    &
         ,U_EN_U_c, V_EN_U_c, U_EN_V_c, V_EN_V_c  &
         ,fois_U_c, plus_U_c, fois_V_c, plus_V_c  &
         ,U1,U2,U3, V1,V2,V3, P1,P2,P3
  end type grid_const

  type (grid_const), dimension(:), save, allocatable :: grid
  
  real(kind=prec), dimension (:,:), pointer :: &
       U_EN_U,   V_EN_U,   U_EN_V,   V_EN_V   &
      ,U_EN_U_c, V_EN_U_c, U_EN_V_c, V_EN_V_c &
      ,fois_U, plus_U, fois_V, plus_V &
      ,fois_U_c, plus_U_c, fois_V_c, plus_V_c
  
end module disc

module champs

  use selectprec

  real(kind=prec), dimension (:,:), allocatable, save :: &
       VTUS,VTVS,PRES , VTU ,VTV ,PRE  &
       , VTU0,VTV0,PRE0 , VTU1,VTV1,PRE1 , VTU2,VTV2,PRE2 , VTU3,VTV3,PRE3 &
       , SMU ,SMV, SMP ,  CLDU,CLDV ,      VTU4,VTV4 

  real(kind=prec), dimension (:,:), allocatable, save, target ::&
       & VTUp,VTVp,PREp,  SMUp,SMVp,SMPp  

  real(kind=prec), save :: vit_adv=1._prec

  logical, dimension (:,:), allocatable, save :: MASK_DIRU,MASK_DIRV

end module champs


module drapeaux

  integer, parameter ::&
       & t_code    = 7 , u_code    = 11 ,v_code    = 13   &
        ,grdv_code = 17 ,grdp_code = 19 ,conv_code = 3    &
        ,diff_code = 2 , div_code  = 23 ,int_code  = 5&
       &,oned_code = 29 

  integer, save      :: global_switch, is_impl,niv_lim,current_switch&
       & , uv_switch

  logical, save      :: is_t,    is_u,    is_v,    is_grdp &
                       ,is_grdv, is_conv, is_diff, is_div, is_int &
                       ,is_decal, is_ns,is_precond,is_cv_global,is_cv_file&
                      &,is_oned , is_schwarz, is_schwarz_add, is_schwarz_mult
                     
  contains
    
    subroutine set_switch (switch)
      implicit none
      integer :: switch

      is_t    =(mod(switch,t_code   )==0)
      is_u    =(mod(switch,u_code   )==0)
      is_v    =(mod(switch,v_code   )==0)
      is_grdv =(mod(switch,grdv_code)==0)
      is_grdp =(mod(switch,grdp_code)==0)
      is_conv =(mod(switch,conv_code)==0)
      is_diff =(mod(switch,diff_code)==0)
      is_div  =(mod(switch,div_code )==0)
      is_int  =(mod(switch,int_code )==0)
      is_oned =(mod(switch,oned_code )==0)

      current_switch=switch

      return
    end subroutine set_switch

    subroutine print_switch                        ! start_out_light
      implicit none

      print *
      print *,'is_t       :',is_t    
      print *,'is_u       :',is_u    
      print *,'is_v       :',is_v    
      print *,'is_grdv    :',is_grdv 
      print *,'is_grdp    :',is_grdp 
      print *,'is_conv    :',is_conv 
      print *,'is_diff    :',is_diff 
      print *,'is_div     :',is_div  
      print *,'is_int     :',is_int  
      print *,'is_ns      :',is_ns
      print *,'is_decal   :',is_decal
      print *,'is_precond :',is_precond

      return
    end subroutine print_switch                   ! end_out_light

  end module drapeaux

module conv
  
  use selectprec

  integer, save :: niv_solve=0
  
  logical, dimension(0:6), save         :: flag_niv=.false.
  integer, dimension(0:6), save         :: nb_iter=0 ,max_iter=-1e8&
       & ,min_iter=1e8 ,last_iter ,total_iter=0
  real(kind=prec), dimension(0:6), save :: max_rho=-1.e8,min_rho=1.e8&
       & ,last_rho,total_rho=0.,first_res,last_res

end module conv


module plot_flags
  logical, parameter :: surface=.true.,one_d=.false.

  character(len=255),save :: controle_string='1'  ! on met 1 pour que le fichier soit tjrs plein
  logical, save :: is_plot_psi,is_plot_w
end module plot_flags


module debug

  logical, save      :: debug_snd=.false.,  debug_rcv=.false., debug_save=.false.
                     
  contains
    

end module debug


