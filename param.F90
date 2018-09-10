module param
  use, intrinsic :: iso_c_binding
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters and global varables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! -----------------------------------------------------------------
! Set Model Resolution
  integer(C_INTPTR_T), parameter :: n1=384, n2=384, n3=384
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Set Number of Processors (must divide n2 and iktz)
  integer, save         :: npe=192
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Math constants
  real, parameter            :: twopi = 4.*asin(1.)
  real, parameter            :: sqrt2 = sqrt(2.)
  complex, parameter :: zi    = cmplx(0.,1.)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Time
  real, save         :: time                             ! global variable showing time 
  real,    parameter :: delt  = 0.2/float(max(n1,n2,n3)) ! timestep
  real,    parameter :: tstop = 8.0                      ! length of integration
  integer, parameter :: nstop  = int(tstop/delt)         ! number of timesteps
  integer, parameter :: noutput = 50                    ! number of times that data dumped in spc/trn/eng/run.list
  integer, parameter :: nout=floor(nstop/float(noutput)) ! do diagnostics every nout timesteps
! ----------------------------------------------------------------- 

! -----------------------------------------------------------------
! Stratification and Rotation
  real, parameter :: PrRa =   32                          ! Prandtle Ratio is set = N/f
  real, parameter :: aj   =   640                        ! thermal expansivity. NOTE: N^2 = aj*bj
  real, parameter :: bj   =   aj                          ! background theta gradient
  real, parameter :: bf2  =   aj*bj                       ! Brunt-Vaisalla frequency squared
  real, parameter :: bf   =   sqrt(bf2)                   ! Brunt-Vaisalla frequency
  real, parameter :: cor  =   bf/PrRa                     ! Coriolis parameter
  real, parameter :: cor2 =   cor*cor                     ! Coriolis parameter squared
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Set Domain Size
  real, parameter :: L1=twopi, L2=twopi
! real, parameter :: L3=twopi*float(n3)/float(n1)        ! This scaling of L3 make the grid uniform, dx=dz
  real, parameter :: L3=twopi/PrRa                       ! makes the grid uniform in charney-scaled coord., dx=N/fdz
! -----------------------------------------------------------------  

! -----------------------------------------------------------------
! Physical and Fourier Array Dimensions
  integer(C_INTPTR_T), parameter :: n1d = n1+2, n2d=n2, n3d=n3         ! size of physical arrays (padded for in-place transforms)
  integer(C_INTPTR_T), save      :: n2dp, n2p                          ! size of local arrays with mpi, n2dp=n2d/npe, n2p=n2/npe
  integer(C_INTPTR_T), parameter :: ktx = n1/2,  kty =n2/2, ktz =n3/2  ! max integer wavenumber
  integer(C_INTPTR_T), parameter :: iktx= ktx+1, ikty=n2,   iktz=n3    ! number of wavenumbers
  integer(C_INTPTR_T), save      :: iktzp                              ! number of local wavenumbers with mpi, iktzp=iktz/npe
  integer(C_INTPTR_T), parameter :: kts = n1                           ! for spectra; should be max(ktx,ktz)
  real, parameter                :: ktrunc_x = twopi/L1*float(n1)/3.   ! dimensional truncation wavenumber (x) --> not index
  real, parameter                :: ktrunc_y = twopi/L2*float(n2)/3.   ! dimensional truncation wavenumber (y)
  real, parameter                :: ktrunc_z = twopi/L3*float(n3)/3.   ! dimensional truncation wavenumber (z)
! -----------------------------------------------------------------  
  
! -----------------------------------------------------------------
! Dissipation
  integer, parameter :: ilap = 4, ilap2 = 2*ilap                       ! hyperviscosity order (ilap=1 is regular viscosity)
  real, parameter    :: visch   =   5e-15                              ! viscosity coeff horizontal
  real, parameter    :: viscz   =   visch*(1/PrRa)**(ilap2)            ! viscosity coeff vertical
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Wavenumbers 
  integer, dimension(:,:,:), allocatable,save :: L               ! Mask for wavenumber truncation
  real, dimension(:), allocatable, save       :: kxa, kya, kza   ! Wavenumbers
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Input/Output related
  integer, parameter :: nread = 0                        ! flag for reading I.C. from a NetCDF input file
  integer, parameter :: ndump    = nstop                 ! time step at which fields are dumped in a restart netcdf
                                                         ! ncwrite = 0 nothing will be dumped in any output files
! -----------------------------------------------------------------

! -----------------------------------------------------------------  
! MPI Stuff
  integer, save :: mype
  integer(C_INTPTR_T), save :: locz, loczstart             ! size and position of block in z   
  integer(C_INTPTR_T), save :: alloc_local                 ! local data size for malloc
  integer, save :: istatus                                 ! used for mpi initiation 
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Misc Stuff
  integer, parameter :: truncmode = 0                      ! 1 for spherical truncation; 0 for cylindrical
  real,    parameter :: fftnorm = float(n1*n2*n3)          ! used to normal fft/fftinverse
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Horizontal and vertical cuts parameters       (push them to a module later)
  integer, parameter :: dumpslices = 1          ! flag to indicate whether or not the real space fields will be dumped
  integer, parameter :: nrsoutput= 400          ! number of times that real space fields are dumped
  integer, parameter :: nrsout=floor(nstop/float(nrsoutput)) ! dump every nrsout timesteps
  integer, save      :: iRSslice = 0            ! count when fields are dumped. increase each time you dump them
  integer, parameter :: nHmax  =  n1            ! starting from index =1 keep up to nHmax in the horizontal dir
                                                ! nHmax  =  n1 keep all the point along x
  integer, parameter :: nVmax  =  n3   
  integer, parameter :: nskipH =  1             ! skip (horizontal) physical points every nskipH (nskip=1 keep all) 
  integer, parameter :: nskipV =  1             ! vertical
! -----------------------------------------------------------------

end module param
