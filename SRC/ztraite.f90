
program traite

    implicit none

    double precision, dimension(:,:), allocatable :: VTX,VTZ,PRE,PSI,XU,ZU,XV,ZV,XP,ZP
    integer :: i,k,i0,k0,ok,lmu,nmu,lmv,nmv,lmp,nmp
    integer, dimension(2) :: pos
    integer, dimension(1) :: pos1d
    double precision :: dx,dz
    character(len=100) :: nom_fic_out='/home/panisse/kortas/COMPACT/ALIZEB/SRC/CONVERT/PCG/RES/'

    open(file=trim(nom_fic_output)//'disc.dat' &
         ,unit=79,form='unformatted',iostat=ok)
    if (ok/=0) stop 'pb ouverture fichier XU, XV, XP'  
    read (79) lmu,nmu,lmv,nmv,lmp,nmp

    allocate(XU(0:lmu+1,0:nmu+1),stat=ok); if (ok/=0) stop 'XU : error alloc'
    allocate(ZU(0:lmu+1,0:nmu+1),stat=ok); if (ok/=0) stop 'ZU : error alloc'
    allocate(XV(0:lmv+1,0:nmv+1),stat=ok); if (ok/=0) stop 'XV : error alloc'
    allocate(ZV(0:lmv+1,0:nmv+1),stat=ok); if (ok/=0) stop 'ZV : error alloc'
    allocate(XP(0:lmp  ,0:nmp)  ,stat=ok); if (ok/=0) stop 'XP : error alloc'
    allocate(ZP(0:lmp  ,0:nmp)  ,stat=ok); if (ok/=0) stop 'ZP : error alloc'

    allocate (VTX(0:lmu+1,0:nmu+1),stat=ok);   if (ok/=0) stop 'VTX : error alloc'
    allocate (VTZ(0:lmv+1,0:nmv+1),stat=ok);   if (ok/=0) stop 'VTZ : error alloc'
    allocate (PRE(0:nmp  ,0:nmp),stat=ok);     if (ok/=0) stop 'PRE : error alloc'
    allocate (PSI(0:lmu  ,0:nmv),stat=ok);     if (ok/=0) stop 'PSI : error alloc'
  
    read (79) XU,ZU,XV,ZV,XP,ZP
    close(79)



    open(file=trim(nom_fic_output)//'all.dat' &
         ,unit=80,form='unformatted',iostat=ok)
    if (ok/=0) stop 'pb ouverture fichier XU, XV, XP'  
    read (80) VTX,VTZ,PRE,PSI
    close(80)

    dx=XU(1,1)-XU(0,1)
    dz=ZU(1,1)-ZU(1,0)



    pos=minloc(abs(XU))
    i0=pos(1)

    pos=minloc(abs(ZU))
    k0=pos(2)

    write(*,127) VTX(i0,:)

    write(*,127) VTZ(:,k0)

    print *,'extrema'

    pos1d=minloc(VTX(i0,:))
    k=pos1d(1)
    print '(A,e12.5,A,e12.5)','umin ',VTX(i0,k), ' en y=',(ZU(i0,k)+1.d0)*0.5d0

    pos1d=maxloc(VTZ(:,k0))
    i=pos1d(1)
    print '(A,e12.5,A,e12.5)','vmax ',VTZ(i,k0), ' en x=',(XU(i,k0)+1.d0)*0.5d0

    pos1d=minloc(VTZ(:,k0))
    i=pos1d(1)
    print '(A,e12.5,A,e12.5)','vmin ',VTZ(i,k0), ' en x=',(XU(i,k0)+1.d0)*0.5d0

127 format(255e16.7)
    


end program traite

