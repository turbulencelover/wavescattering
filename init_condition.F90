!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module generates Initial Conditions (IC) for vorticity and temperature (buoyancy) fields  !!
!! It can read a NetCDF as IC                                                                     !!
!! or generated IC based on a description in physical or Fourier space                            !!
!! or IC can be a superposition of the above                                                      !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module init_condition
  use param
  use param_fftw
  use IO_netcdf
  use velvorproj
  use nm_decomp
  implicit none
  ! -----------------------------------------------------------------  
  ! PARAMETERS (for initial conditions)
  ! set the mode for generating/reading IC
  ! icmode = 0 ---> generate IC only based on descriptions in Fourier and Real Spaces
  ! icmode = 1 ---> read IC from a netCDF file
  ! icmode = 2 ---> a superpostion of the netCDF file and describing functions
  integer, parameter :: icmode =   2
  ! ki is the index of wavenumber at which IC peaks (in energy spectrum)
  real,    parameter :: ki     =   24.0    
  
  ! Make internal variables and functions private
  PRIVATE :: ranno, ran1
  !PRIVATE :: iuSPCH,iuSPCZ,iuSPC,iuTRNH,iuTRNZ,iuTRN

CONTAINS
  
subroutine init_cond(zx,zy,zz,tt,uu,vv,ww,ge,g1,g2,zxr,zyr,zzr,ur,vr,wr,tr,ts)
! Define the initial condition
  implicit none
  include 'mpif.h'
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,uu,vv,ww,ge,g1,g2
  real,    intent(out), dimension(n1d,n3d,n2dp)    :: zxr,zyr,zzr,tr,ur,vr,wr
  real,    intent(out) :: ts
  complex, dimension(:,:,:), allocatable :: zx0, zy0, zz0, tt0
  real     :: wvfr,kxwv,kywv,khwv,kzwv,wven

  
  if (icmode==2) then
     ! superpostion of input NetCDF file and user defined fields
     ! for user defined fields you can use the functions in this module like single_planewave ...
     allocate(zx0(iktx,ikty,iktzp))
     allocate(zy0(iktx,ikty,iktzp))
     allocate(zz0(iktx,ikty,iktzp))
     allocate(tt0(iktx,ikty,iktzp))
     call ncreadin(zx,zy,zz,tt,ts)

     ! setting the initial wave field (that is going to be added to background QG
     wvfr  = sqrt(2.0)*cor                    ! frequency of the wave initially shooted
     ! kxwv  = 12.0                           ! kx of the initial wave  
     ! kywv  = 0.0                            ! ky of the initial wave
     khwv  = 12.0   !sqrt(kxwv**2+kywv**2)    ! kh of the initial wave
     kzwv  = khwv * sqrt((bf2 - wvfr**2)/(wvfr**2 - cor**2)) ! derive kz based on dispersion rel.
     wven  = 0.1 
     ! call single_planewave(zx0,zy0,zz0,tt0,uu,vv,ww,ge,g1,g2,wven,kxwv,kywv,kzwv)
     call waves_ring(zx0,zy0,zz0,tt0,uu,vv,ww,ge,g1,g2,wven,khwv,kzwv)
     
     ! superposing the wave field and background QG
     zx = zx + zx0
     zy = zy + zy0
     zz = zz + zz0
     tt = tt + tt0

     ! free the memory for auxilliary variables 
     deallocate(zx0,zy0,zz0,tt0)

  elseif (icmode==1) then
     call ncreadin(zx,zy,zz,tt,ts)
  else ! (icmode==0)
     call exp_peak(zx,zy,zz,tt)
  endif

  call wtoab(zx,zy,zz,tt,ge,g1,g2,uu,vv,ww)
  call proj(zx,zy,zz)

! -----------------------------------------------------------------
! Print out Initial Conditions
! -----------------------------------------------------------------
  if (mype.eq.0) then
     print*,'                '
     print*,'Initial Condition: ---------------------------------------'
     print*,'Starting from Zk.in.ncf + wave field'
     !print*,'Single plane wave at :' 
     !print*,'         kxwv,kywv,kzwv = ', kxwv, kywv, kzwv
     print*,'    Ring of waves  at :' 
     print*,'             khwv,kzwv = ', khwv, kzwv
     print*,'    Initial wave energy = ',wven
     print*,'                '
  endif
  
  return
end subroutine init_cond

subroutine single_planewave(zx,zy,zz,tt,uu,vv,ww,ge,g1,g2,wven,kxwv,kywv,kzwv)
!! Generate a single wave at (kxwv,kywv,kzwv) with the energy = wven 
  implicit none
  include 'mpif.h'
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,uu,vv,ww,ge,g1,g2
  integer :: ikx,iky,ikz,ikza
  real    :: kx,ky,kz,kh
  real    :: kxwv,kywv,khwv,kzwv,wven

  ge = cmplx(0.,0.)
  g1 = cmplx(0.,0.)
  g2 = cmplx(0.,0.)
 
  do ikx = 1,iktx
     kx = kxa(ikx)
     if (abs(kx - kxwv).lt.(twopi/2/L1)) then
        do iky = 1,ikty
           ky = kya(iky)
           kh = sqrt(kx*kx + ky*ky)
           if (abs(ky - kywv).lt.(twopi/2/L2)) then
              do ikz = 1,iktzp
                 ikza = mype*iktzp+ikz
                 kz = kza(ikza)
                 if (abs(kz - kzwv).lt.(twopi/2/L3)) then
                    g1(ikx,iky,ikz) = kh*sqrt(wven)
                   !g2(ikx,iky,ikz) = kh*sqrt(wven)
                 endif
              enddo
           endif
        enddo
     endif
  enddo
   
  call atowb(ge,g1,g2,zx,zy,zz,tt,uu,vv,ww)
  call proj(zx,zy,zz)
  
end subroutine single_planewave


subroutine waves_ring(zx,zy,zz,tt,uu,vv,ww,ge,g1,g2,wven,khwv,kzwv)
!! Generate a ring of waves at kh = constant and kz = another constant with the energy = wven 
  implicit none
  include 'mpif.h'
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,uu,vv,ww,ge,g1,g2
  integer :: ikx,iky,ikz,ikza
  real    :: kx,ky,kz,kh,nmd,nmdtot
  real    :: khwv,kzwv,wven

  ge = cmplx(0.,0.)
  g1 = cmplx(0.,0.)
  g2 = cmplx(0.,0.)

  nmd = 0.0

 
  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     if (abs(kz - kzwv).le.(twopi/L3/2.0)) then
        do iky = 1,ikty
           ky = kya(iky)
           do ikx = 1,iktx
              kx = kxa(ikx)
              kh = sqrt(kx*kx + ky*ky)
              if (abs(kh - khwv).le.(twopi/L1/2.0)) then
                 g1(ikx,iky,ikz) = kh*sqrt(wven/2.0)
                 g2(ikx,iky,ikz) = kh*sqrt(wven/2.0)
                 nmd = nmd + 1
              endif
           enddo
        enddo
     endif
  enddo

  call mpi_allreduce(nmd,nmdtot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,istatus);
  if (istatus.ne.0) print*,'mpi_allreduce for nmdtot is screwed up '
  
  if (nmdtot.eq.0) then
     print*, ' ------> something is screwed up in I.C. (ring of waves)'
     call MPI_Abort(MPI_COMM_WORLD)
  else
     ! print*,'mype,nmdtot =', mype,nmdtot
     g1 = g1/sqrt(nmdtot)
     g2 = g2/sqrt(nmdtot)
  endif
  
  call atowb(ge,g1,g2,zx,zy,zz,tt,uu,vv,ww)
  call proj(zx,zy,zz)
  
end subroutine waves_ring


subroutine exp_peak(zx,zy,zz,tt)
!!! Initialize in Fourier space with Gausian peak at ki (defined above)  
  implicit none
  include 'mpif.h'
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  integer :: ikx,iky,ikz,ikza
  real    :: kx,ky,kz,wk,kh,khn,wkn,kzn
  real    ::  phase,EK
  
  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx = 1,iktx
           kx = kxa(ikx)
           wk  = sqrt(kx*kx + ky*ky + kz*kz)
           khn = sqrt((kx*L1/twopi)**2 + (ky*L2/twopi)**2)
           wkn = sqrt((kx*L1/twopi)**2 + (ky*L2/twopi)**2 + (kz*L3/twopi)**2)
           kzn = kz * L3/twopi
           if (L(ikx,iky,ikz).eq.1) then
              EK = exp(-((wkn-ki)/1.)**2)
              phase = ranno(0)
              zx(ikx,iky,ikz) = wk*sqrt(EK)*cexp(ZI*phase)
              phase = ranno(0)
              zy(ikx,iky,ikz) = wk*sqrt(EK)*cexp(ZI*phase)
              phase = ranno(0)
              zz(ikx,iky,ikz) = wk*sqrt(EK)*cexp(ZI*phase)
              phase = ranno(0)
              tt(ikx,iky,ikz) =    sqrt(EK)*cexp(ZI*phase)
              phase = ranno(0)
              if (kx.eq.0.) then
                 zx(ikx,iky,ikz) = cmplx(real(zx(ikx,iky,ikz)),0.)
                 zy(ikx,iky,ikz) = cmplx(real(zy(ikx,iky,ikz)),0.)
                 zz(ikx,iky,ikz) = cmplx(real(zz(ikx,iky,ikz)),0.)
              endif
           endif
        enddo
     enddo
  enddo
     
end subroutine exp_peak

function ranno (i)

! Controls random number generator.
!-----------------------------------------
! - If argument i.ne.0 it performs initialization with i=seed no.
! - If argument i.eq.0 it draws a random no.
!-----------------------------------------
  implicit none
  real :: ranno
  integer :: i
  integer :: junk,ihold
  real :: twopi,ran1
  save junk

  twopi = 4.*asin(1.)
  if (i.ne.0) then
    if (i.gt.0) i = - i
    junk  = i
    ranno = (ran1(i)-0.5)*twopi
  else
    junk  = junk - 1
    ihold = junk
    ranno = (ran1(ihold)-0.5)*twopi
  endif
  return
end function ranno

function ran1(idum)
  implicit none
  real :: ran1
  integer :: idum

  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ntab=32,ndiv=1+(im-1)/ntab
  real, parameter :: am=1./im,eps=1.2e-7,rnmx=1.-eps
  integer :: j,k,iv(32),iy
  save iv,iy
  data iv /ntab*0/, iy /0/
  if (idum.le.0.or.iy.eq.0) then
    idum=max(-idum,1)
    do j=ntab+8,1,-1
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      if (j.le.ntab) iv(j)=idum
    enddo
    iy=iv(1)
  endif
  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  if (idum.lt.0) idum=idum+im
  j=1+iy/ndiv
  iy=iv(j)
  iv(j)=idum
  ran1=min(am*iy,rnmx)
  return
end function ran1
   
end module init_condition

   
