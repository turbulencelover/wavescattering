module IO_netcdf
  use param
  use netcdf
  implicit none
  
  integer :: idink,idoutk,idkx,idky,idkz,idkri,idktm,istat
  integer :: idzxk,idzyk,idzzk,idttk,idtimek 

CONTAINS
  
subroutine ncdumpout(zx,zy,zz,tt,ts)
  !! Create netcdf files and write the fields for restart (output)
  !! Only zx,zy,zz and tt are dumped in a netcdf file
  implicit none
  include 'mpif.h'
  complex, intent(in),    dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  real,    intent(in)   :: ts
  integer, dimension(4) :: ncdimsk
  integer, dimension(4) :: nccountk,ncstartk

!---------------------------------------------------------------------
!     DEFINING THE OUTPUT NETCDF FILE 

!!! prep restart file (output) (define the variables and dims)
  if (mype.eq.0)  print*,'Yo! dumping in netcdf restart file'
  
  ! Create the netCDF file. The nf90_clobber parameter tells netCDF,
  ! just to be safe "nf90_64bit_offset" is used for large data, but it seems unnecessary for netCDF4 version

!!$  print*, 'my rank, zx, zy ',mype, zx, zy
  ! print*, 'my rank, tt(1,1,1) ',mype, tt(1,1,1)
  call check( nf90_create("Zk.out.ncf",ior(NF90_NETCDF4,NF90_MPIIO),idoutk, &
       comm = MPI_COMM_WORLD,info = MPI_INFO_NULL) )
  
  ! Define the dimension: kx, ky, kz. time. RI used as another dim to distinguish between real and imag parts
  call check( nf90_def_dim(idoutk, "kx",int(iktx,4),idkx) )
  call check( nf90_def_dim(idoutk, "ky",int(ikty,4),idky) )
  call check( nf90_def_dim(idoutk, "kz",int(iktz,4),idkz) )
  call check( nf90_def_dim(idoutk, "RI",2 , idkri) )
  call check( nf90_def_dim(idoutk, "t" ,1 , idktm) )
     
  ! ncdimsk array used to pass the dimension IDs to the variables
  ncdimsk = (/idkx,idky,idkz,idkri/)
  
  ! Define the variables which are the fields that going to be stored
  ! i.e. zx, zy, zz and tt
  call check( nf90_def_var(idoutk,"zx" ,NF90_FLOAT,ncdimsk,idzxk))
  call check( nf90_def_var(idoutk,"zy" ,NF90_FLOAT,ncdimsk,idzyk))
  call check( nf90_def_var(idoutk,"zz" ,NF90_FLOAT,ncdimsk,idzzk))
  call check( nf90_def_var(idoutk,"tt" ,NF90_FLOAT,ncdimsk,idttk))
  call check( nf90_def_var(idoutk,"time",NF90_FLOAT,idktm,idtimek))
     
  
  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(idoutk) )

!---------------------------------------------------------------------
  !     WRITING VARIABLES IN NETCDF FILE
  
  ! how many to count in each dimension when writing files
  nccountk(1) = iktx
  nccountk(2) = ikty
  nccountk(3) = iktzp
  nccountk(4) = 1
  ! where to start on the output file
  ncstartk(1) = 1
  ncstartk(2) = 1
  
  call check( nf90_put_var(idoutk, idzxk, real(zx)) )

  ! wrtie time (time is written only one time so just root process is used
  if (mype.eq.0) then 
     call check( nf90_put_var(idoutk, idtimek, ts) )
  endif

  ! in the z-direction mype=0 is in the first place
  ncstartk(3) = mype*iktzp+1
  ncstartk(4) = 1 ! (real part)
  call check( nf90_put_var(idoutk, idzxk, real(zx),ncstartk,nccountk))    
  call check( nf90_put_var(idoutk, idzyk, real(zy),ncstartk,nccountk))
  call check( nf90_put_var(idoutk, idzzk, real(zz),ncstartk,nccountk))
  call check( nf90_put_var(idoutk, idttk, real(tt),ncstartk,nccountk))
  ncstartk(4) = 2 ! (imag part)
  call check( nf90_put_var(idoutk, idzxk, aimag(zx),ncstartk,nccountk))    
  call check( nf90_put_var(idoutk, idzyk, aimag(zy),ncstartk,nccountk))
  call check( nf90_put_var(idoutk, idzzk, aimag(zz),ncstartk,nccountk))
  call check( nf90_put_var(idoutk, idttk, aimag(tt),ncstartk,nccountk))

  
!---------------------------------------------------------------------
!     CLOSING the  NETCDF FILE 
  call check (nf90_close(idoutk))
  
contains
  subroutine check(status)
    integer, intent (in) :: status
    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped!!!"
    end if
  end subroutine check
  
end subroutine ncdumpout

subroutine ncreadin(zx,zy,zz,tt,ts)
! Read the netcdf data from an input file (Zk.in.ncf)
  implicit none
  include 'mpif.h'
  
  complex, intent(out),    dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  real, dimension(iktx,ikty,iktzp) :: wr,wi
  real,    intent(out) :: ts
  integer :: iktx1,ikty1,iktz1
  integer, dimension(5) :: ncstartr,ncstarti,nccount
  
  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access
  if (mype == 0) print*,'Yo! reading from netcdf restart file'
  call check( nf90_open('Zk.in.ncf', NF90_NOWRITE,idink) )

  ! Get the dimensions IDs based on their name
  call check( nf90_inq_dimid(idink, "kx", idkx) )
  call check( nf90_inq_dimid(idink, "ky", idky) )
  call check( nf90_inq_dimid(idink, "kz", idkz) )
  call check( nf90_inq_dimid(idink, "RI",idkri) )
  call check( nf90_inq_dimid(idink, "t", idktm) )

  ! Get the dimension length and check if the grid resolution matches
  call check( nf90_inquire_dimension(idink,idkx,len=iktx1))  
  call check( nf90_inquire_dimension(idink,idky,len=ikty1))
  call check( nf90_inquire_dimension(idink,idkz,len=iktz1))  

  if (iktx1.ne.iktx .or. ikty1.ne.ikty .or. iktz1.ne.iktz) then
    print*,'Sorry, do not know how to change resolution.'
    stop
  endif

  ! Get the variables IDs
  call check( nf90_inq_varid(idink,"time",idtimek))
  call check( nf90_inq_varid(idink, "zx" ,idzxk) )
  call check( nf90_inq_varid(idink, "zy" ,idzyk) )
  call check( nf90_inq_varid(idink, "zz" ,idzzk) )
  call check( nf90_inq_varid(idink, "tt" ,idttk) )

  ! prepare for reading variables
  ncstartr    = 1
  ncstarti    = 1
  ncstarti(4) = 2 ! for imag part start from 2 for dim = 4
  nccount(1) = iktx
  nccount(2) = ikty
  nccount(3) = iktzp
  nccount(4) = 1

  ncstartr(3)= int(mype*iktzp+1)
  ncstarti(3)= int(mype*iktzp+1)

  call check(  nf90_get_var(idink,idtimek,ts))

  call check( nf90_get_var(idink,idzxk,wr,ncstartr,nccount))
  call check( nf90_get_var(idink,idzxk,wi,ncstarti,nccount))
  zx = wr + zi * wi
  call check( nf90_get_var(idink,idzyk,wr,ncstartr,nccount))
  call check( nf90_get_var(idink,idzyk,wi,ncstarti,nccount))
  zy = wr + zi * wi
  call check( nf90_get_var(idink,idzzk,wr,ncstartr,nccount))
  call check( nf90_get_var(idink,idzzk,wi,ncstarti,nccount))
  zz = wr + zi * wi
  call check( nf90_get_var(idink,idttk,wr,ncstartr,nccount))
  call check( nf90_get_var(idink,idttk,wi,ncstarti,nccount))
  tt = wr + zi * wi

  call check (nf90_close(idink))
  
contains
  subroutine check(status)
    integer, intent (in) :: status
    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped!!!"
    end if
  end subroutine check
  
end subroutine ncreadin

end module IO_netcdf
 
