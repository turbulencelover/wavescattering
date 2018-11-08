!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module calculates 2D energy spectra as a function of (kh,kz) & stores them in a NetCDF    !!
!! Different types of energy spectrum can be selected to be stored:                               !!
!! Potential, Kinetic, Geostrophic, Ageostrophic and Buoyancy Flux                                !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module ncf2Dspc
  use param
  use netcdf
  use velvorproj
  use nm_decomp
  implicit none
  
  !! %% PARAMETERs %% 
  ! Here we determine which energy spectra should be stored:
  integer, parameter :: putE = 0, putKE = 0,putPE = 0, putGE = 1, putAE = 0, putBF = 0
  !! ---> putE: total energy, putKE:kinetic energy ,putPE: poterntial energy
  !! ---> putGE: geostrophic energy, putAE: ageostrophic energy, putBF: buoyancy flux 
  !! >>NOTE<< set putE = 0 if either putKE=putPE=1 or putGE=putAE=1 to save storage space
  integer, parameter :: putKEWV = 1, putPEWV = 1
  !! putKEWV, putPEWV potential and kinetic energy of wave (ageostrophic) modes ONLY
  !! >>NOTE<< calulating putKEWV, putPEWV are costly! set them zero if not necessary in analysis
  !! >>NOTE<< set putAE = 0 if either putKEWV=putPEWV=1 as putAE = putKEWV + putPEWV
  integer, parameter :: putKEMF = 1, putPEMF = 1
  !! putKEMF, putPEMF potential and kinetic energy of mean flow aka (geostrophic) modes ONLY
  !! >>NOTE<< calulating putKEMF, putPEMF are costly! set them zero if not necessary in analysis
  !! >>NOTE<< set putGE = 0 if either putKEMF=putPEMF=1 as putGE = putKEMF + putPEMF
  !number of modes kept (starts from 0 to numkh-1)
  integer, parameter :: numkh =floor(2*ktx/3.0), numkz =floor(2*ktz/3.0)
  
  !! the frequency of outputting
  integer, parameter :: nsp2dout = 20                    ! number of time 2D spectrum is dumped in output file
  integer, parameter :: nsp2devery=floor(nstop/float(nsp2dout))! dump 2D spectra every nsp2devery timesteps
  
  
  !! %% VARIABLES %%
  integer, save  :: idspc2d,idkh2,idkz2,idtm2,idkhkh2,idkzkz2,idnmd2 ! NetCDF IDs
  integer, save  :: iden2,idke2,idpe2,idge2,idae2,idbf2,idkewv2,idpewv2,idpemf2,idkemf2 
  integer, save  :: ispc2d = 0

  !! Make internal variables and functions private
  PRIVATE        :: putE, putKE, putPE, putGE, putAE, putBF, numkh, numkz
  PRIVATE        :: check_ncf2
  
CONTAINS

subroutine check_ncf2(status)
  implicit none
  include 'mpif.h'  
  integer, intent (in) :: status
  if(status /= nf90_noerr) then
     print *, 'ERROR in ncf2Dspc.F90'
     print *, trim(nf90_strerror(status))
     call MPI_Abort(MPI_COMM_WORLD)
  end if
end subroutine check_ncf2
  
subroutine prep_ncf2Dspc()
  ! creates a netcdf file for 2D energy spectrum.
  implicit none
  include 'mpif.h'
    
  !! >>NOTE<< only mype=0 sees the NetCDF file. It's unnecessary for other processes
  !! unlike IO NetCDF files, 2Dspc data is smaller and can be dumped by just one process

  istatus =  nf90_create("spc2D.ncf",ior(NF90_NETCDF4,NF90_MPIIO),idspc2d,comm= MPI_COMM_WORLD,info = MPI_INFO_NULL)
  if (istatus.ne.0) print*,'Yo! error in creating spc2De.ncf! Darn!'
  ! Define the dimensions of 2D spectrum: kh, kz, and time
  call check_ncf2( nf90_def_dim(idspc2d,"kh",numkh,idkh2) )
  call check_ncf2( nf90_def_dim(idspc2d,"kz",numkz,idkz2) )
  call check_ncf2( nf90_def_dim(idspc2d, "ttt" ,nsp2dout+1,idtm2) ) ! nsp2dout+ >>1<< for dumping spc of IC
  ! Define the time variable
  !call check_ncf2( nf90_def_var(idspc2d,"times",NF90_FLOAT,(/ idtm2 /),idtsp2))
  ! Define number of Modes variable for each (kh,kz)
  call check_ncf2( nf90_def_var(idspc2d,"nModes",NF90_INT,(/idkh2,idkz2/),idnmd2))
  ! Define wavenumber arrays (dimensional)
  call check_ncf2( nf90_def_var(idspc2d,"kkh",NF90_FLOAT,(/idkh2,idkz2/),idkhkh2))
  call check_ncf2( nf90_def_var(idspc2d,"kkz",NF90_FLOAT,(/idkh2,idkz2/),idkzkz2))
  ! Define the variables which are the fields that going to be stored
  ! i.e. zx, zy, zz and tt
  if (putE .eq.1) call check_ncf2( nf90_def_var(idspc2d,"E"  ,NF90_FLOAT,(/idkh2,idkz2,idtm2/),iden2))
  if (putKE.eq.1) call check_ncf2( nf90_def_var(idspc2d,"KE" ,NF90_FLOAT,(/idkh2,idkz2,idtm2/),idke2))
  if (putPE.eq.1) call check_ncf2( nf90_def_var(idspc2d,"PE" ,NF90_FLOAT,(/idkh2,idkz2,idtm2/),idpe2))
  if (putGE.eq.1) call check_ncf2( nf90_def_var(idspc2d,"GE" ,NF90_FLOAT,(/idkh2,idkz2,idtm2/),idge2))
  if (putAE.eq.1) call check_ncf2( nf90_def_var(idspc2d,"AE" ,NF90_FLOAT,(/idkh2,idkz2,idtm2/),idae2))
  if (putBF.eq.1) call check_ncf2( nf90_def_var(idspc2d,"BFlux",NF90_FLOAT,(/idkh2,idkz2,idtm2/),idbf2))
  if (putKEWV.eq.1) call check_ncf2( nf90_def_var(idspc2d,"KEWV" ,NF90_FLOAT,(/idkh2,idkz2,idtm2/),idkewv2))
  if (putPEWV.eq.1) call check_ncf2( nf90_def_var(idspc2d,"PEWV" ,NF90_FLOAT,(/idkh2,idkz2,idtm2/),idpewv2))
  if (putKEMF.eq.1) call check_ncf2( nf90_def_var(idspc2d,"KEMF" ,NF90_FLOAT,(/idkh2,idkz2,idtm2/),idkemf2))
  if (putPEMF.eq.1) call check_ncf2( nf90_def_var(idspc2d,"PEMF" ,NF90_FLOAT,(/idkh2,idkz2,idtm2/),idpemf2))
  ! End define mode. This tells netCDF we are done defining metadata.
  call check_ncf2( nf90_enddef(idspc2d) )

  return
end subroutine prep_ncf2Dspc

subroutine dump_ncf2Dspc(zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,zxwv,zywv,zzwv,ttwv)
  
  implicit none 
  include 'mpif.h'
  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,ge,g1,g2
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ux,uy,uz,zxwv,zywv,zzwv,ttwv
  integer :: ikx,iky,ikz,ikza
  integer :: jh,jv,ierr
  integer :: ncoll,ncoll6,ncoll10    
  real    :: kx,ky,kz,wk,wkh2,wkh,wkh2n,kzn,kz2
  real    :: vz,vt,vzx,vzy,vzz
  complex :: vcx,vcy,vcz,vct
  integer, dimension(3) :: nccount2,ncstart2
  real, dimension(1:numkh,1:numkz) :: Ejunk
  complex, dimension(:,:,:), allocatable :: gewv
  real, dimension(:,:,:), allocatable    :: spc2d,spctot2d
  integer, dimension(0:numkh,0:numkz)    :: nm2d,nmtot2d
  real, dimension(1:numkh,1:numkz)       :: khkh,kzkz 

  ispc2d = ispc2d + 1
  if (ispc2d==1) then
     call prep_ncf2Dspc()
  endif
  
  if (putKEWV+putPEWV+putKEMF+putPEMF.ne.0) then
     allocate(gewv(iktx,ikty,iktzp))
     allocate(spc2d(0:numkh,0:numkz,10),spctot2d(0:numkh,0:numkz,10))  
  else
     allocate(spc2d(0:numkh,0:numkz,6),spctot2d(0:numkh,0:numkz,6))
  endif
  
  ncoll  = (numkh+1)*(numkz+1)
  ncoll6 = (numkh+1)*(numkz+1)*6
  ncoll10 = (numkh+1)*(numkz+1)*10
  nm2d    = 0
  spc2d   = 0.0

  call velo (zx,zy,zz,ux,uy,uz)

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kz2 = kz*kz
     kzn = kz*L3/twopi
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           wkh2  = kx*kx+ky*ky
           wkh2n = wkh2 * (L1/twopi)**2
           wkh   = sqrt(wkh2)
           wk    = sqrt(kx*kx+ky*ky+kz*kz)	              
           jh = ifix(wkh*L1/twopi+0.5)
           jv = ifix(abs(kz)*L3/twopi+0.5)          

           if (((jv.le.numkz).and.(jh.le.numkh)).and.(L(ikx,iky,ikz).eq.1)) then           
              wk   = max(wk,  1.e-15)
              wkh2 = max(wkh2,1.e-15)

! Kinetic and potential energy.
              vzx      = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
              spc2d(jh,jv,1) = spc2d(jh,jv,1) + vzx/wk**2
              vzy      = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
              spc2d(jh,jv,1) = spc2d(jh,jv,1) + vzy/wk**2
              vzz      = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
              spc2d(jh,jv,1) = spc2d(jh,jv,1) + vzz/wk**2
              vt       = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )
              spc2d(jh,jv,2) = spc2d(jh,jv,2) + vt*aj/bj
! Geo, ageo decompostition.
! k \in R_k
              if (wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 vz       = real( ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)) )
                 spc2d(jh,jv,4) = spc2d(jh,jv,4) + vz/wkh2
                 vz       = real( g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)) )
                 spc2d(jh,jv,5) = spc2d(jh,jv,5) + vz/wkh2
                 vz       = real( g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)) )
                 spc2d(jh,jv,5) = spc2d(jh,jv,5) + vz/wkh2
              endif
! Special cases: i) k_h=0, ii) k_z=0.
              vzx=real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              vzy=real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              vzz=real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )
! k \in V_k
              if(wkh2n.lt.1.e-10.and.abs(kzn).gt.1.e-10) then
                 spc2d(jh,jv,4) = spc2d(jh,jv,4) + vzz + vt*aj/bj
                 spc2d(jh,jv,5) = spc2d(jh,jv,5) + vzx + vzy
              endif
! k \in B_k
              if(abs(kzn).lt.1.e-10.and.wkh2n.gt.1.e-10) then
                 spc2d(jh,jv,4) = spc2d(jh,jv,4) + vzx + vzy
                 spc2d(jh,jv,5) = spc2d(jh,jv,5) + vzz + vt*aj/bj
              endif
! Buoyancy Flux
              if (aj.gt.0.) then 
                 spc2d(jh,jv,6) = spc2d(jh,jv,6) + aj*real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
              else
                 spc2d(jh,jv,6) = spc2d(jh,jv,6) +    real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
              endif
              nm2d(jh,jv) = nm2d(jh,jv) + 2 
           endif
        enddo
     enddo
  enddo

  spc2d(:,:,3) = spc2d(:,:,1) + spc2d(:,:,2)

  if (putKEWV+putPEWV+putKEMF+putPEMF.ne.0) then 
     gewv=cmplx(0.,0.)
     call atowb(gewv,g1,g2,zxwv,zywv,zzwv,ttwv,ux,uy,uz)
     do ikz=1,iktzp
        ikza = mype*iktzp+ikz
        kz = kza(ikza)
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              wkh   = sqrt(kx*kx+ky*ky)
              wk    = sqrt(kx*kx+ky*ky+kz*kz)	              
              jh = ifix(wkh*L1/twopi+0.5)
              jv = ifix(abs(kz)*L3/twopi+0.5)          
              if (((jv.le.numkz).and.(jh.le.numkh)).and.(L(ikx,iky,ikz).eq.1)) then           
                 wk   = max(wk,  1.e-15)
                 ! Kinetic and potential energy of the waves
                 vzx      = real( zxwv(ikx,iky,ikz)*conjg(zxwv(ikx,iky,ikz)) )
                 vzy      = real( zywv(ikx,iky,ikz)*conjg(zywv(ikx,iky,ikz)) )
                 vzz      = real( zzwv(ikx,iky,ikz)*conjg(zzwv(ikx,iky,ikz)) )
                 vt       = real( ttwv(ikx,iky,ikz)*conjg(ttwv(ikx,iky,ikz)) )
                 spc2d(jh,jv,7) = spc2d(jh,jv,7) + vzx/wk**2 + vzy/wk**2 + vzz/wk**2
                 spc2d(jh,jv,8) = spc2d(jh,jv,8) + vt*aj/bj
                 vcx      = zx(ikx,iky,ikz)-zxwv(ikx,iky,ikz)    
                 vcy      = zy(ikx,iky,ikz)-zywv(ikx,iky,ikz)    
                 vcz      = zz(ikx,iky,ikz)-zzwv(ikx,iky,ikz)    
                 vct       = tt(ikx,iky,ikz)-ttwv(ikx,iky,ikz)    
                 spc2d(jh,jv,9)  = spc2d(jh,jv,9)+real(vcx*conjg(vcx))/wk**2 &
                       + real(vcy*conjg(vcy))/wk**2 + real(vcz*conjg(vcz))/wk**2
                 spc2d(jh,jv,10) = spc2d(jh,jv,10) + real(vct*conjg(vct))*aj/bj
              endif
           enddo
        enddo
     enddo
  endif

  if (putKEWV+putPEWV+putKEMF+putPEMF.ne.0) then
     deallocate(gewv)
     call mpi_reduce(spc2d, spctot2d, ncoll10,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierr)
  else
     call mpi_reduce(spc2d, spctot2d, ncoll6,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierr)
  endif
  call mpi_reduce( nm2d,  nmtot2d, ncoll ,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  if (mype.eq.0) then
  ! At the first dumping of spectra, store wavenumber field and number of modes summed over
     if (ispc2d == 1) then
        call check_ncf2( nf90_put_var(idspc2d,idnmd2,nmtot2d(0:numkh-1,0:numkz-1)) )
        do jv=0,numkz-1
           do jh=0,numkh-1
               khkh(jh+1,jv+1) = jh * twopi/L1
               kzkz(jh+1,jv+1) = jv * twopi/L3
            enddo
        enddo
        call check_ncf2( nf90_put_var(idspc2d,idkhkh2,khkh) )
        call check_ncf2( nf90_put_var(idspc2d,idkzkz2,kzkz) )
     endif
     ncstart2(1:2)=1
     ncstart2(3)=ispc2d
     nccount2(1)=numkh
     nccount2(2)=numkz
     nccount2(3)=1
     if (putE.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,3) 
        call check_ncf2( nf90_put_var(idspc2d,iden2,Ejunk,ncstart2,nccount2) )
     endif
     if (putKE.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,1) 
        call check_ncf2( nf90_put_var(idspc2d,idke2,Ejunk,ncstart2,nccount2) )
     endif
     if (putPE.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,2)
        call check_ncf2( nf90_put_var(idspc2d,idpe2,Ejunk,ncstart2,nccount2) )
     endif
     if (putGE.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,4)
        call check_ncf2( nf90_put_var(idspc2d,idge2,Ejunk,ncstart2,nccount2) )
     endif
     if (putAE.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,5)
        call check_ncf2( nf90_put_var(idspc2d,idae2,Ejunk,ncstart2,nccount2) )
     endif
     if (putBF.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,6)
        call check_ncf2( nf90_put_var(idspc2d,idbf2,Ejunk,ncstart2,nccount2) )
     endif
     if (putKEWV.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,7) 
        call check_ncf2( nf90_put_var(idspc2d,idkewv2,Ejunk,ncstart2,nccount2) )
     endif
     if (putPEWV.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,8)
        call check_ncf2( nf90_put_var(idspc2d,idpewv2,Ejunk,ncstart2,nccount2) )
     endif
     if (putKEMF.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,9) 
        call check_ncf2( nf90_put_var(idspc2d,idkemf2,Ejunk,ncstart2,nccount2) )
     endif
     if (putPEMF.eq.1)  then
        Ejunk(1:numkh,1:numkz) = spctot2d(0:numkh-1,0:numkz-1,10)
        call check_ncf2( nf90_put_var(idspc2d,idpemf2,Ejunk,ncstart2,nccount2) )
     endif
  endif
  deallocate(spc2d,spctot2d)

   if (ispc2d == nsp2dout+1) then
     call close_ncf2Dspc()
  endif
  
  return
end subroutine dump_ncf2Dspc

subroutine close_ncf2Dspc()
  ! close creates a netcdf file for 2D energy spectrum.
  implicit none
  include 'mpif.h'
  call check_ncf2( nf90_close(idspc2d) )
  return
end subroutine close_ncf2Dspc

end module ncf2Dspc
 
