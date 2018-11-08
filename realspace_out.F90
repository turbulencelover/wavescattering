!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module  stores them in a NetCDF    !!
!! Different types of energy spectrum can be selected to be stored:                               !!
!! Potential, Kinetic, Geostrophic, Ageostrophic and Buoyancy Flux                                !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module realspace_out
  use param
  use param_fftw
  use netcdf
  use velvorproj
  use nm_decomp
  implicit none
  
  !! %% PARAMETERs %%
  !! ============== 2D CUTS (HORIZONTAL AND VERTICAL SLICES) ==============
  ! Here we determine which fields should be stored:
  ! 1) total fields 
  integer, parameter :: rsd_zx = 0,rsd_zy = 0,rsd_zz = 1,rsd_tt = 0,rsd_u = 0,rsd_v = 0,rsd_w = 0
  ! 2) just wave component
  integer, parameter :: rsd_zxwv=0,rsd_zywv=0,rsd_zzwv=0,rsd_ttwv=0,rsd_uwv=1,rsd_vvw=1,rsd_wvw=1
  ! @@@@@@@ Horizontal Cut @@@@@@@ 
  integer, parameter :: ihor   = 1    ! index of horizontal slice. ihor=0 no hor. slices stored
  integer, parameter :: nskipH = 1    ! skip (horizontal) physical points every nskipH (nskip=1 keep all) 
  integer, parameter :: nHmax  = n1   ! starting from index=1 keep up to nHmax in the horizontal dir
  ! @@@@@@@  Vertical Cut  @@@@@@@ 
  integer, parameter :: iver   = 1    ! index of vertical slice. ihor=0 no ver. slices stored
  integer, parameter :: nskipV = 1    ! skip (vertical) physical points every nskipV 
  integer, parameter :: nVmax  = n3

  ! the frequency of outputting
  integer, parameter :: nrs2out  = 20             ! number of time slices are dumped in the output file
  ! dump the real space slices every nrs2evry timesteps
  integer, parameter :: nrs2evry = floor(nstop/float(nrs2out))
  
  !! %% VARIABLES %%
  integer, save  :: idrs2,idxx,idyy,idzz,idti2 ! NetCDF IDs
  
  integer, save  :: irs2 = 0   ! counter: count when fields are dumped. increase each time you dump them

  !! Make internal variables and functions private
  PRIVATE        :: putE, putKE, putPE, putGE, putAE, putBF, numkh, numkz
  PRIVATE        :: check_ncf2
  
CONTAINS 

subroutine check_rs(status)
  implicit none
  include 'mpif.h'  
  integer, intent (in) :: status
  if(status /= nf90_noerr) then
     print *, 'ERROR in realspace_out.F90'
     print *, trim(nf90_strerror(status))
     call MPI_Abort(MPI_COMM_WORLD)
  end if
end subroutine check_rs
  
subroutine prep_realspace()
  ! creates a netcdf file for 2D energy spectrum.
  implicit none
  include 'mpif.h'
  integer :: numxx,numyy,numzz
  integer, dimension(3) :: ncdimh,ncdimv

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF,
  ! just to be safe "nf90_64bit_offset" is used for large data, but it seems unnecessary for netCDF4 version
  ! similar IO_NetCDF files, all the processor have access to the output files (unlike ncf2Dspc)
  
  istatus =  nf90_create("realspace_slices.ncf",ior(NF90_NETCDF4,NF90_MPIIO),idrs2,comm= MPI_COMM_WORLD,info = MPI_INFO_NULL)
  if (istatus.ne.0) print*,'Yo! error in creating rs2De.ncf! Darn!'
  ! Define the dimensions of 2D spectrum: kh, kz, and time
  numxx = floor((nHmax+1.0)/real(nskipH))
  numyy = numxx
  numzz = floor((nVmax+1.0)/real(nskipV))
  
  call check_rs( nf90_def_dim(idrs2,"xx",numxx,idxx) )
  call check_rs( nf90_def_dim(idrs2,"yy",numyy,idyy) )
  call check_rs( nf90_def_dim(idrs2,"zz",numzz,idzz) )
  call check_rs( nf90_def_dim(idrs2,"time", nrs2out+1,idti2) ) ! nrs2out + >>1<< for dumping IC
  
  ! Define the time variable
  call check_rs( nf90_def_var(idrs2,"times",NF90_FLOAT,(/idti2/),idtime2))

  ncdimh = (/idxx,idyy,idti2/)
  ncdimv = (/idxx,idzz,idti2/)
  
  ! Define the variables which are the fields that going to be stored
  ! i.e. zx, zy, zz and tt
  if (ihor.ne.0) then
     if (rsd_zx.eq.1) call check_rs( nf90_def_var(idrs2,"ZXh",NF90_FLOAT,ncdimh,idzxh))
     if (rsd_zy.eq.1) call check_rs( nf90_def_var(idrs2,"ZYh",NF90_FLOAT,ncdimh,idzyh))
     if (rsd_zz.eq.1) call check_rs( nf90_def_var(idrs2,"ZZh",NF90_FLOAT,ncdimh,idzzh))
     if (rsd_tt.eq.1) call check_rs( nf90_def_var(idrs2,"TTh",NF90_FLOAT,ncdimh,idtth))
     if (rsd_u .eq.1) call check_rs( nf90_def_var(idrs2,"Uh" ,NF90_FLOAT,ncdimh,iduh ))
     if (rsd_v .eq.1) call check_rs( nf90_def_var(idrs2,"Vh" ,NF90_FLOAT,ncdimh,idvh ))
     if (rsd_w .eq.1) call check_rs( nf90_def_var(idrs2,"Wh" ,NF90_FLOAT,ncdimh,idwh ))
     ! wave component of variables
     if (rsd_zxvw.eq.1) call check_rs( nf90_def_var(idrs2,"ZXh_wv",NF90_FLOAT,ncdimh,idzxhwv))
     if (rsd_zyvw.eq.1) call check_rs( nf90_def_var(idrs2,"ZYh_wv",NF90_FLOAT,ncdimh,idzyhwv))
     if (rsd_zzvw.eq.1) call check_rs( nf90_def_var(idrs2,"ZZh_wv",NF90_FLOAT,ncdimh,idzzhwv))
     if (rsd_ttvw.eq.1) call check_rs( nf90_def_var(idrs2,"TTh_wv",NF90_FLOAT,ncdimh,idtthwv))
     if (rsd_uvw .eq.1) call check_rs( nf90_def_var(idrs2,"Uh_wv" ,NF90_FLOAT,ncdimh,iduhwv ))
     if (rsd_vvw .eq.1) call check_rs( nf90_def_var(idrs2,"Vh_wv" ,NF90_FLOAT,ncdimh,idvhwv ))
     if (rsd_wvw .eq.1) call check_rs( nf90_def_var(idrs2,"Wh_wv" ,NF90_FLOAT,ncdimh,idwhwv ))
  endif

  if (iver.ne.0) then
     if (rsd_zx.eq.1) call check_rs( nf90_def_var(idrs2,"ZXv",NF90_FLOAT,ncdimv,idzxv))
     if (rsd_zy.eq.1) call check_rs( nf90_def_var(idrs2,"ZYv",NF90_FLOAT,ncdimv,idzyv))
     if (rsd_zz.eq.1) call check_rs( nf90_def_var(idrs2,"ZZv",NF90_FLOAT,ncdimv,idzzv))
     if (rsd_tt.eq.1) call check_rs( nf90_def_var(idrs2,"TTv",NF90_FLOAT,ncdimv,idttv))
     if (rsd_u .eq.1) call check_rs( nf90_def_var(idrs2,"Uv" ,NF90_FLOAT,ncdimv,iduv ))
     if (rsd_v .eq.1) call check_rs( nf90_def_var(idrs2,"Vv" ,NF90_FLOAT,ncdimv,idvv ))
     if (rsd_w .eq.1) call check_rs( nf90_def_var(idrs2,"Wv" ,NF90_FLOAT,ncdimv,idwv ))
     ! wave component of variables
     if (rsd_zxvw.eq.1) call check_rs( nf90_def_var(idrs2,"ZXv_wv",NF90_FLOAT,ncdimv,idzxvwv))
     if (rsd_zyvw.eq.1) call check_rs( nf90_def_var(idrs2,"ZYv_wv",NF90_FLOAT,ncdimv,idzyvwv))
     if (rsd_zzvw.eq.1) call check_rs( nf90_def_var(idrs2,"ZZv_wv",NF90_FLOAT,ncdimv,idzzvwv))
     if (rsd_ttvw.eq.1) call check_rs( nf90_def_var(idrs2,"TTv_wv",NF90_FLOAT,ncdimv,idttvwv))
     if (rsd_uvw .eq.1) call check_rs( nf90_def_var(idrs2,"Uv_wv" ,NF90_FLOAT,ncdimv,iduvwv ))
     if (rsd_vvw .eq.1) call check_rs( nf90_def_var(idrs2,"Vv_wv" ,NF90_FLOAT,ncdimv,idvvwv ))
     if (rsd_wvw .eq.1) call check_rs( nf90_def_var(idrs2,"Wv_wv" ,NF90_FLOAT,ncdimv,idwvwv ))
  endif
  
  ! End define mode. This tells netCDF we are done defining metadata.
  call check_rs( nf90_enddef(idrs2) )

  return
end subroutine prep_realspace

subroutine dump_realspace2(zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,zxwv,zywv,zzwv,ttwv)
  
  implicit none 
  include 'mpif.h'
  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,ge,g1,g2
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ux,uy,uz,zxwv,zywv,zzwv,ttwv
 

  ispc2d = ispc2d + 1
  if (ispc2d==1) then
     call prep_ncf2Dspc()
  endif
  
  if ((putKEWV+putPEWV+putKEMF+putPEMF).ne.0) then
     allocate(gewv(iktx,ikty,iktzp))
  endif

  if ((rsd_u+rsd_v+rsd_w).ne.0)  call velo (zx,zy,zz,ux,uy,uz)
  
  ncoll  = (numkh+1)*(numkz+1)
  ncoll6 = (numkh+1)*(numkz+1)*6
  ncoll10 = (numkh+1)*(numkz+1)*10
  nm2d    = 0
  spc2d   = 0.0

 

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
           jh = floor(wkh*L1/twopi+0.5)
           jv = floor(kzn+0.5)+numkz/2.0          

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
        call check_rs( nf90_put_var(idspc2d,idnmd2,nmtot2d(0:numkh-1,0:numkz-1)) )
        do jv=0,numkz-1
           do jh=0,numkh
               khkh(jh+1,jv+1) = jh * twopi/L1
               kzkz(jh+1,jv+1) = (jv - numkz/2.0) * twopi/L3
            enddo
        enddo
        call check_rs( nf90_put_var(idspc2d,idkhkh2,khkh) )
        call check_rs( nf90_put_var(idspc2d,idkzkz2,kzkz) )
     endif
     ncstart2(1:2)=1
     ncstart2(3)=ispc2d
     nccount2(1)=numkh
     nccount2(2)=numkz
     nccount2(3)=1
     if (putE.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,3) 
        call check_rs( nf90_put_var(idspc2d,iden2,Ejunk,ncstart2,nccount2) )
     endif
     if (putKE.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,1) 
        call check_rs( nf90_put_var(idspc2d,idke2,Ejunk,ncstart2,nccount2) )
     endif
     if (putPE.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,2)
        call check_rs( nf90_put_var(idspc2d,idpe2,Ejunk,ncstart2,nccount2) )
     endif
     if (putGE.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,4)
        call check_rs( nf90_put_var(idspc2d,idge2,Ejunk,ncstart2,nccount2) )
     endif
     if (putAE.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,5)
        call check_rs( nf90_put_var(idspc2d,idae2,Ejunk,ncstart2,nccount2) )
     endif
     if (putBF.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,6)
        call check_rs( nf90_put_var(idspc2d,idbf2,Ejunk,ncstart2,nccount2) )
     endif
     if (putKEWV.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,7) 
        call check_rs( nf90_put_var(idspc2d,idkewv2,Ejunk,ncstart2,nccount2) )
     endif
     if (putPEWV.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,8)
        call check_rs( nf90_put_var(idspc2d,idpewv2,Ejunk,ncstart2,nccount2) )
     endif
     if (putKEMF.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,9) 
        call check_rs( nf90_put_var(idspc2d,idkemf2,Ejunk,ncstart2,nccount2) )
     endif
     if (putPEMF.eq.1)  then
        Ejunk(1:numkh,1:numkz+1) = spctot2d(0:numkh-1,0:numkz,10)
        call check_rs( nf90_put_var(idspc2d,idpemf2,Ejunk,ncstart2,nccount2) )
     endif
  endif
  deallocate(spc2d,spctot2d)

   if (ispc2d == nsp2dout+1) then
     call close_ncf2Dspc()
  endif
  
  return
end subroutine dump_realspace2

subroutine close_ncf2Dspc()
  ! close creates a netcdf file for 2D energy spectrum.
  implicit none
  include 'mpif.h'
  call check_rs( nf90_close(idspc2d) )
  return
end subroutine close_ncf2Dspc

end module ncf2Dspc2
 
