module garbage

  use selectprec
  use constante

  interface derf0
     module procedure derf_tab,derf_real
  end interface

 character(len=*), parameter :: type_machine='Blue GENE/L'
!LIN character(len=*), parameter :: type_machine='PC Linux'
!SGI character(len=*), parameter :: type_machine='Silicon Graphics'
!SP2 character(len=*), parameter :: type_machine='IBM SP2'
!T3D character(len=*), parameter :: type_machine='Cray T3E'

contains

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function  derf_real(x) result(y)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
 
    real(kind=prec) :: x,y
    
    y=0.

    return
  end function derf_real


  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function derf_tab(x) result (smu)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none
    real(kind=prec), dimension(:,:) :: x
    real(kind=prec), dimension(size(x,1),size(x,2)) :: smu
    integer :: i,k

    do k=1,size(x,2)
       do i=1,size(x,1)
          smu(i,k)=derf_real(x(i,k))
       enddo
    enddo

    return
  end function derf_tab

  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function second_membre_u(x,y) result (smu)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    implicit none

    real(kind=prec), dimension(:,:) :: x,y
    real(kind=prec), dimension(size(x,1),size(x,2)) :: smu
    integer :: i,k

    do k=1,size(x,2)
       do i=1,size(x,1)
          call sm_u(smu(i,k),x(i,k),y(i,k))
       end do
    end do

    return
  end function second_membre_u


  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  function second_membre_v(x,y) result (smv)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

    implicit none

    real(kind=prec), dimension(:,:) :: x,y
    real(kind=prec), dimension(size(x,1),size(x,2)) :: smv
    integer :: i,k

    do k=1,size(x,2)
       do i=1,size(x,1)
          call sm_v(smv(i,k),x(i,k),y(i,k))
       end do
    end do

    return
  end function second_membre_v




  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine sm_u(u,x,y)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data
    use disc
    implicit none
    real(kind=prec) x,y,u,t,s1,s2,s3,s4,s5

    t=temps

    select case(nexample)
    case(29)
       s1 = -0.2_prec*(0.1_prec)**8
       s4=2000000000._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*& 
           &t*y+y**2)/nu)*sqrt(nu)*t-1000000000._prec*exp(-(2._prec*t**2-2._p& 
           &rec*t*x+x**2-2._prec*t*y+y**2)/nu)*sqrt(nu)*x-1000000000._prec*ex& 
           &p(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*sqrt(nu)*& 
           &y-2000000000._prec*sqrt(nu**3)*exp(-(2._prec*t**2-2._prec*t*x+x**& 
           &2-2._prec*t*y+y**2)/nu)+2000000000._prec*sqrt(nu)*exp(-(2._prec*t& 
           &**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x**2-4000000000._prec*& 
           &sqrt(nu)*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu& 
           &)*t*x+4000000000._prec*sqrt(nu)*exp(-(2._prec*t**2-2._prec*t*x+x*& 
           &*2-2._prec*t*y+y**2)/nu)*t**2
       s5=s4+2000000000._prec*sqrt(nu)*exp(-(2._prec*t**2-2._prec*t*x+x& 
           &**2-2._PREC*t*y+y**2)/nu)*y**2-4000000000._prec*sqrt(nu)*exp(-(2.& 
           &_prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*y+1000000000.& 
           &_prec*exp(-2._prec*(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**& 
           &2)/nu)*sqrt(nu)*x
       s3=s5-1000000000._prec*exp(-2._prec*(2._prec*t**2-2._prec*t*x+x*& 
           &*2-2._prec*t*y+y**2)/nu)*sqrt(nu)*t+1772453851._prec*exp(-(2._pre& 
           &c*x**2-4._prec*t*x+3._PREC*t**2-2._prec*t*y+y**2)/nu)*derf_real((& 
           &y-t)/sqrt(nu))*x*y-1772453851._prec*exp(-(2._prec*x**2-4._prec*t*& 
           &x+3._prec*t**2-2._prec*t*y+y**2)/nu)*derf_real((y-t)/sqrt(nu))*t*& 
           &x-1772453851._prec*exp(-(2._prec*x**2-4._prec*t*x+3._prec*t**2-2.& 
           &_PREC*t*y+y**2)/nu)*derf_real((y-t)/sqrt(nu))*t*y+1772453851._pre& 
           &c*exp(-(2._prec*x**2-4._prec*t*x+3._prec*t**2-2._prec*t*y+y**2)/n& 
           &u)*derf_real((y-t)/sqrt(nu))*t**2
       s4 = 1/sqrt(nu**3)
       s2 = s3*s4
       u= s1*s2
    case(28)
       s1 = -10._prec
       s3 = 1/nu
       s5=-6._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2& 
            &)/nu)*t*y-10._prec*exp(-2._prec*(2._prec*t**2-2._prec*t*x+x**2-2.& 
            &_prec*t*y+y**2)/nu)*t*nu+10._prec*exp(-2._prec*(2._prec*t**2-2._p& 
            &rec*t*x+x**2-2._prec*t*y+y**2)/nu)*x*nu-2._prec*exp(-(2._prec*t**& 
            &2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*x+2._prec*exp(-(2._pre& 
            &c*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y*x-8._prec*exp(-(2& 
            &._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*nu-8._prec*e& 
            &xp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t**2*x+4& 
            &._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*& 
            &t*x**2-16._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y& 
            &**2)/nu)*t**2*y
       s4=s5+12._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y& 
            &**2)/nu)*t*y**2+8._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._pr& 
            &ec*t*y+y**2)/nu)*y*nu-4._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2& 
            &-2._prec*t*y+y**2)/nu)*y*x**2-exp(-(2._prec*t**2-2._prec*t*x+x**2& 
            &-2._prec*t*y+y**2)/nu)*nu+4._prec*exp(-(2._prec*t**2-2._prec*t*x+& 
            &x**2-2._prec*t*y+y**2)/nu)*t**2+2._prec*exp(-(2._prec*t**2-2._pre& 
            &c*t*x+x**2-2._prec*t*y+y**2)/nu)*y**2+8._prec*exp(-(2._prec*t**2-& 
            &2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t**3-4._prec*exp(-(2._prec& 
            &*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y**3+8._prec*exp(-(2& 
            &._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y*t*x
       s2 = s3*s4
       u = s1*s2
    case(17)
       t=t_start
       s1 = -20._prec
       s3 = 1/nu
       s5 = -4._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*nu+4&
            &._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t**3-4._prec*exp&
            &(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t**2*x+2._prec*exp(-(2.&
            &E0*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*x**2-8._prec*exp(-(2._prec*t*&
            &*2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t**2*y+6._prec*exp(-(2._prec*t**2-2.&
            &E0*t*x+x**2-2._prec*t*y+y**2)/nu)*t*y**2
       s4 = s5+4._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y*nu&
            &+4._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y*t*x-2._prec*&
            &exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y*x**2-2._prec*exp(-&
            &(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y**3-5._prec*exp(-2._prec*(2&
            &._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*nu+5._prec*exp(-2._prec*(2.E&
            &0*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x*nu
       s2 = s3*s4
       u = s1*s2
    case default 
       stop 'danger : second membre non defini'
    end select


    return
  end subroutine sm_u




  !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine sm_v(v,x,y)
    !SKsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    use data
    use disc
    implicit none

    real(kind=prec) x,y,v,t,s1,s2,s3,s4,s5

    t=temps

    select case(nexample)
    case (29)
       s1 = (0.1_prec)**9
       s4=3544907702._prec*exp(-(x-t)**2/nu)*derf_real((y-t)/sqrt(nu))*& 
           &x**2-7089815404._prec*exp(-(x-t)**2/nu)*derf_real((y-t)/sqrt(nu))& 
           &*t*x+3544907702._prec*exp(-(x-t)**2/nu)*derf_real((y-t)/sqrt(nu))& 
           &*t**2-2000000000._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._pre& 
           &c*t*y+y**2)/nu)*sqrt(nu)*x+2000000000._prec*exp(-(2._prec*t**2-2.& 
           &_prec*t*x+x**2-2._prec*t*y+y**2)/nu)*sqrt(nu)*t-1772453851._prec*& 
           &exp(-(x-t)**2/nu)*derf_real((y-t)/sqrt(nu))*nu+10634723110._prec*& 
           &x*exp(-(x-t)**2/nu)*derf_real((y-t)/sqrt(nu))*nu-7089815404._prec& 
           &*exp(-(x-t)**2/nu)*derf_real((y-t)/sqrt(nu))*x**3
       s5=s4+21269446210._prec*exp(-(x-t)**2/nu)*derf_real((y-t)/sqrt(n& 
           &u))*t*x**2-21269446210._prec*x*exp(-(x-t)**2/nu)*derf_real((y-t)/& 
           &sqrt(nu))*t**2+4000000000._prec*x*exp(-(2._prec*t**2-2._prec*t*x+& 
           &x**2-2._prec*t*y+y**2)/nu)*sqrt(nu)*y-4000000000._prec*sqrt(nu)*e& 
           &xp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*x
       s3=s5-10634723110._prec*t*exp(-(x-t)**2/nu)*derf_real((y-t)/sqrt& 
           &(nu))*nu+7089815404._prec*exp(-(x-t)**2/nu)*derf_real((y-t)/sqrt(& 
           &nu))*t**3-4000000000._prec*sqrt(nu)*exp(-(2._prec*t**2-2._prec*t*& 
           &x+x**2-2._prec*t*y+y**2)/nu)*t*y+4000000000._prec*sqrt(nu)*exp(-(& 
           &2._prec*t**2-2._prec*t*x+x**2-2._PREC*t*y+y**2)/nu)*t**2+17724538& 
           &51._prec*derf_real((y-t)/sqrt(nu))*nu*exp(-(2._prec*x**2-4._prec*& 
           &t*x+3._prec*t**2-2._prec*t*y+y**2)/nu)
       s4 = 1/sqrt(nu**3)
       s2 = s3*s4
       v = s1*s2
    case (28)
       s1 = 10._prec
       s4=-10._prec*exp(-2._prec*(2._prec*t**2-2._prec*t*x+x**2-2._prec& 
            &*t*y+y**2)/nu)*y*nu-2._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2& 
            &._prec*t*y+y**2)/nu)*t*y+10._prec*exp(-2._prec*(2._prec*t**2-2._p& 
            &rec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*nu+8._prec*exp(-(2._prec*t**& 
            &2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x*nu-4._prec*exp(-(2._pr& 
            &ec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x*y**2-6._prec*exp& 
            &(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*x+2._pre& 
            &c*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y*x-8& 
            &._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*& 
            &t*nu-16._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**& 
            &2)/nu)*t**2*x
       s3=s4+12._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y& 
            &**2)/nu)*t*x**2-8._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._pr& 
            &ec*t*y+y**2)/nu)*t**2*y+4._prec*exp(-(2._prec*t**2-2._prec*t*x+x*& 
            &*2-2._prec*t*y+y**2)/nu)*t*y**2-exp(-(2._prec*t**2-2._prec*t*x+x*& 
            &*2-2._prec*t*y+y**2)/nu)*nu+4._prec*exp(-(2._prec*t**2-2._prec*t*& 
            &x+x**2-2._prec*t*y+y**2)/nu)*t**2+8._prec*exp(-(2._prec*t**2-2._p& 
            &rec*t*x+x**2-2._prec*t*y+y**2)/nu)*t**3+8._prec*exp(-(2._prec*t**& 
            &2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y*t*x+2._prec*exp(-(2._p& 
            &rec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x**2-4._prec*exp(& 
            &-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x**3
       s4 = 1/nu
       s2 = s3*s4
       v = s1*s2
    case (17)
       t=t_start
      s1 = 20._prec
      s4 = -4._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*nu+4&
     &._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t**3-8._prec*exp&
     &(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t**2*x+6._prec*exp(-(2.&
     &E0*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*x**2-4._prec*exp(-(2._prec*t*&
     &*2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t**2*y+2._prec*exp(-(2._prec*t**2-2.&
     &E0*t*x+x**2-2._prec*t*y+y**2)/nu)*t*y**2
      s3 = s4+4._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x*nu&
     &-2._prec*exp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x**3+4._prec*e&
     &xp(-(2._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x*t*y-2._prec*exp(-(2&
     &._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*x*y**2+5._prec*exp(-2._prec*(2&
     &._prec*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*t*nu-5._prec*exp(-2._prec*(2.E&
     &0*t**2-2._prec*t*x+x**2-2._prec*t*y+y**2)/nu)*y*nu
      s4 = 1/nu
      s2 = s3*s4
      v = s1*s2

    case default
       stop 'danger : second membre non defini'
    end select


    return
  end subroutine sm_v





end module garbage



!***********************************************************************

!SP2      subroutine flush(n)
!SP2      return
!SP2      end subroutine flush

!SGI      subroutine get_time(s)
!SGI      use selectprec
!SGI      implicit none
!SGI      real(kind=prec), intent(out) :: s
!SGI      integer ::  mclock
!SGI      external mclock 
!SGI      
!SGI      s=0.01*mclock()
!SGI      
!SGI      end subroutine get_time


