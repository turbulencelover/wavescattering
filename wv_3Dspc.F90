!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module calculates and stores 3D spectrum of wave energy in a netcdf file                  !!
!! allows to truncate the spectrum to make the stored file lighter in high-res simulations        !!
!!                                                                                                !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module wv_3Dspc
  use param
  use netcdf
  implicit none
  
  !! %% PARAMETERs %% 
  ! truncation wavenumber (non-dimensional i.e. just the index)
  ! The difference between iktxsp & ktxsp:
  ! iktxws shows the number of wavenumbers up to truncation, both positive and negative
  ! ktxws shows the number of all the same signed wavenumbers kept up to truncation
  ! k t w s --> k:(wavenumber) t:truncation w:wave s:spectrum
  integer, parameter :: ktxws = n1/8,     ktyws =n2/8,       ktzws = n3/8 
  integer, parameter :: iktxws= ktxws+1 , iktyws=2*ktyws ,   iktzws=2*ktzws
  ! numner of times the wave spc is dumped:
  integer, parameter :: nwsout = 30
  ! output the wave spectrum every nwsout
  integer, parameter :: nwsevery= floor(nstop/float(nwsout))
  
  !! %% VARIABLES %%
  integer, save  :: iwavespcdump = 0  ! count how many 3D wave spectrum are stored
  ! NetCDF IDs
  integer, save  :: idwv3dspc,idkwx,idkwy,idkwz,idwtm,idkxws,idkyws,idkzws,idew1,idew2,idtws
  ! Start and end indices (and the number of wavenumbers from start to end)
  integer, save  :: ikxstr,iikxstr,ikystrp,iikystrp,ikystrm,iikystrm
  ! PAY ATTENTION to differences:
  integer, save  :: ikzstr !(starting) index for variables in each process (between 1:iktzp)
  integer, save  :: ikzastr !(starting) index for kza (between 1:iktz)
  integer, save  :: iikzstr !(starting) index for kzws (between 1:iktzsw-1) 
  ! (N)umber of (x-, y-, z-) wavenumbers (IN) (W)avce (S)pectrum
  integer, save  :: nxinws,nyinwsp,nyinwsm,nzinws
  
  !! Make internal variables and functions private
 
  
CONTAINS
  
subroutine check_ncfwv(status)
  implicit none
  include 'mpif.h'  
  integer, intent (in) :: status
  if(status /= nf90_noerr) then
     print *, 'ERROR in wv_3Dspc.F90'
     print *, trim(nf90_strerror(status))
     call MPI_Abort(MPI_COMM_WORLD)
  end if
end subroutine check_ncfwv

subroutine prep_ncf_wv3Dspc()
  ! creates and defines a netcdf file for 3D wave energy spectrum
  implicit none
  include 'mpif.h'
  integer, dimension(4) :: ncdimwv
  integer, dimension(4) :: ncstart,nccount

  ! if (mype.eq.0)  print*,'Yo! creating netCDF file for 3D wave spectrum! Yayyyy!'
  istatus =  nf90_create("WVspc3D.ncf",ior(NF90_NETCDF4,NF90_MPIIO),idwv3dspc,comm= MPI_COMM_WORLD,info = MPI_INFO_NULL)
  if (istatus.ne.0) print*,'Yo! error in creating WVspc3D.ncf! Darn!'
  ! Define the dimensions of 3D spectrum: kwx, kwy, kwz, and time
  call check_ncfwv( nf90_def_dim(idwv3dspc,"ikx",int(iktxws-1,4),idkwx))
  call check_ncfwv( nf90_def_dim(idwv3dspc,"iky",int(iktyws-1,4),idkwy))
  call check_ncfwv( nf90_def_dim(idwv3dspc,"ikz",int(iktzws-1,4),idkwz))
  call check_ncfwv( nf90_def_dim(idwv3dspc,"itime" ,  nwsout+1 ,idwtm)) ! nwsout + >>1<< for dumping spc of IC

  ncdimwv = (/idkwx,idkwy,idkwz,idwtm/)
  ! Define the time variable
  call check_ncfwv( nf90_def_var(idwv3dspc,"time",NF90_FLOAT,(/idwtm/),idtws))
  ! Define the wavenumber variables (dimensional)
  call check_ncfwv( nf90_def_var(idwv3dspc,"kxw",NF90_FLOAT,idkwx,idkxws))
  call check_ncfwv( nf90_def_var(idwv3dspc,"kyw",NF90_FLOAT,idkwy,idkyws))
  call check_ncfwv( nf90_def_var(idwv3dspc,"kzw",NF90_FLOAT,idkwz,idkzws))
  ! Define the wave energy variables 
  call check_ncfwv( nf90_def_var(idwv3dspc,"EWV1", NF90_FLOAT,ncdimwv,idew1))
  call check_ncfwv( nf90_def_var(idwv3dspc,"EWV2", NF90_FLOAT,ncdimwv,idew2))
  ! End define mode. This tells netCDF we are done defining metadata.
  call check_ncfwv( nf90_enddef(idwv3dspc) )


  return
end subroutine prep_ncf_wv3Dspc

subroutine initkws()
  
  implicit none
  include 'mpif.h'
  
  integer :: ikz,ikza,iikz,stflag
  integer, dimension(iktzws-1) :: kzws

  
  ! Checks for the number of modes 
  if ( mod(npe,2).ne.0 ) then
     print*, 'Dude! the wave spectrum is screwed!!! npe gotta be even!'
     call MPI_Abort(MPI_COMM_WORLD)
  endif

  ! Find the wavenumber indices
  ikxstr   = 1                ! index from which we start storing FOR THE FULL VARIABLES 
  iikxstr  = ikxstr           ! index in the output file where we store
  ikystrp  = 1                ! index from which we start storing in NetCDF file for ky>=0
  ikystrm  = ikty-(ktyws-1)+1 ! index from which we start storing in NetCDF file for ky<0
  iikystrm = 1 
  iikystrp = (iktyws-1)-ktyws+1        
  
  ! count number of wavenumbers in each direction that are stored:
  nxinws  = iktxws-1 
  nyinwsp = ktyws
  nyinwsm = ktyws-1

  ! Derive (dimensional) vertical wavenumbers used for 3D spectra
  do ikz = 1,iktzws - 1
     kzws(ikz) = float(ikz - ktzws) * twopi/L3
  enddo

  iikzstr = -1 ! starting index of kzws 
  ikzstr  = -1 ! starting index between 1:iktzp that communicates with kzws
  ikzastr = -1 ! starting index between 1:iktz
  ! to better understand ---> ikzastr = mype*iktzp + ikzstr
  nzinws  = 0  ! number of vert. modes in this thread/proc that are used for 3D spec
  stflag = 0  ! a flag to show if we started taking from this proc

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     if (ikza.le.iktz/2) then
        do iikz = iktzws/2,iktzws-1
           if (kzws(iikz).eq.kza(ikza)) then
              if (stflag.eq.0) then       
                 nzinws = nzinws + 1
                 ikzstr   = ikz
                 ikzastr  = ikza
                 iikzstr  = iikz
                 stflag = 1
                 exit
              else
                 nzinws = nzinws + 1
                 exit
              endif
           endif
        enddo
     else
        do iikz = 1,iktzws/2-1
           if (kzws(iikz).eq.kza(ikza)) then
              if (stflag.eq.0) then       
                 nzinws = nzinws + 1
                 ikzstr   = ikz
                 ikzastr  = ikza
                 iikzstr  = iikz
                 stflag = 1
                 exit
              else
                 nzinws = nzinws + 1
                 exit
              endif
           endif
        enddo
     endif
  enddo
  

  if (mype.eq.0) then
     call check_ncfwv( nf90_put_var(idwv3dspc,idkxws,kxa(ikxstr:ikxstr+nxinws-1) &
          ,start=(/ iikxstr /),count=(/ nxinws /)) )
      call check_ncfwv( nf90_put_var(idwv3dspc,idkyws,kya(ikystrm:ikystrm+nyinwsm-1) &
          ,start=(/ iikystrm /),count=(/ nyinwsm /)) )
     call check_ncfwv( nf90_put_var(idwv3dspc,idkyws,kya(ikystrp:ikystrp+nyinwsp-1) &
          ,start=(/ iikystrp /),count=(/ nyinwsp /)) )
  endif
  if (nzinws.ne.0) then
!!$     print*, 'mype, kza(start), kza(end),ikzstr, iikzstr, nzinws = ', &
!!$          mype, kza(ikzastr)/32.0, kza(ikzastr+nzinws-1)/32.0, ikzastr, iikzstr, nzinws
     call check_ncfwv( nf90_put_var(idwv3dspc,idkzws,kza(ikzastr:ikzastr+nzinws-1) &
          ,start=(/ iikzstr /),count=(/ nzinws /)) )
  endif
  
end subroutine initkws

subroutine dump_wv3Dspc(g1,g2)
  implicit none
  include 'mpif.h'
  complex, intent(in), dimension(iktx,ikty,iktzp) :: g1,g2
  real, dimension(:,:,:), allocatable :: ew1,ew2
  integer  :: ix, iy, iz
  real     :: kx, ky, wkh2
  complex  :: gg1, gg2
  
  iwavespcdump = iwavespcdump + 1
  if (iwavespcdump==1) then
     call prep_ncf_wv3Dspc()
     call initkws()
  endif

  if (nzinws.ne.0) then      
     allocate(ew1(nxinws,nyinwsm+nyinwsp,nzinws))
     allocate(ew2(nxinws,nyinwsm+nyinwsp,nzinws))
     ew1 = 0.0
     ew2 = 0.0 
     do iz = 0, nzinws - 1    
        do iy = 0, nyinwsm - 1
           ky = kya(ikystrm+iy)
           do ix = 0, nxinws - 1
              kx = kxa(ikxstr+ix)
              wkh2  = kx*kx+ky*ky
              if (wkh2.gt.1.e-10) then
                 gg1 = g1(ikxstr+ix,ikystrm+iy,ikzstr+iz)
                 ew1(iikxstr+ix,iikystrm+iy,1+iz)=real(gg1*conjg(gg1))/wkh2
                 gg2 = g2(ikxstr+ix,ikystrm+iy,ikzstr+iz)
                 ew2(iikxstr+ix,iikystrm+iy,1+iz)=real(gg2*conjg(gg2))/wkh2
              endif
           enddo
        enddo
     enddo
     
     do iz = 0, nzinws - 1    
        do iy = 0, nyinwsp - 1
           ky = kya(ikystrp+iy)
           do ix = 0, nxinws - 1
              kx = kxa(ikxstr+ix)
              wkh2  = kx*kx+ky*ky
              if (wkh2.gt.1.e-10) then
                 gg1 = g1(ikxstr+ix,ikystrp+iy,ikzstr+iz)
                 ew1(iikxstr+ix,iikystrp+iy,1+iz)=real(gg1*conjg(gg1))/wkh2
                 gg2 = g2(ikxstr+ix,ikystrp+iy,ikzstr+iz)
                 ew2(iikxstr+ix,iikystrp+iy,1+iz)=real(gg2*conjg(gg2))/wkh2
              endif
           enddo
        enddo
     enddo


     call check_ncfwv( nf90_put_var(idwv3dspc,idtws,time,start=(/ iwavespcdump /)) )
         
     call check_ncfwv( nf90_put_var(idwv3dspc,idew1,ew1,start=(/ 1, 1, iikzstr, iwavespcdump /) &
          ,count=(/ nxinws, nyinwsm+nyinwsp, nzinws, 1 /)) )
     call check_ncfwv( nf90_put_var(idwv3dspc,idew2,ew2,start=(/ 1, 1, iikzstr, iwavespcdump /) &
          ,count=(/ nxinws, nyinwsm+nyinwsp, nzinws, 1 /)) )

     deallocate(ew1)
     deallocate(ew2)
  endif
   
  if (iwavespcdump==nwsout+1) then
     if (mype.eq.0) print*, 'Dude! closing WVspc3D.ncf'
     call close_ncf_wv3Dspc()
  endif
     
  
end subroutine dump_wv3Dspc

      
subroutine close_ncf_wv3Dspc()
  ! close the netcdf file for 3D wave energy spectrum.
  implicit none
  include 'mpif.h'
  call check_ncfwv( nf90_close(idwv3dspc) )
  return
end subroutine close_ncf_wv3Dspc




end module wv_3Dspc
 
