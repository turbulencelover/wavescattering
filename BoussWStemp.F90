! Triply Periodic Boussinesq
! Original Code by Peter Bartello
! Modified by Michael Waite: mwaite@uwaterloo.ca
! f90, mpi: Winter 2008-2009
! Further modifications and personalization by Hossein Kafiabad
! Winter 2016
  
! Set preprocessor flags:
! SPHTRUNC: 1 for spherical truncation1; 0 for cylindrical
! SCHEME: 1 for AB3; 0 for leapfrog
! NETCDF: 1 to use netcdf i/o.  NOTE: at the moment, this is required.
! MPI: 1 to compile with MPI; 0 for serial
! BTbalance: This turns-on the baer-Tribbia variables 
! UnbalEvol : Turns on the calculation of dAb/dt and interaction between balance and imbalance
! KeepLM : keeping timeseries of limited modes for frequency spectra  
! ----> NOTE to save memory BTinit should be turned off if Baer-Tribbia is not used
#define SPHTRUNC    0  
#define SCHEME      1 
#define NETCDF      1    
#define MPI         1   !(without MPI is disabled for NOW)
#define BTbalance   0   ! Make sure this is on when kpGAb2G = 1
#define UnbalEvol   0  
#define KeepLM      0  
#define HCVC4movie  0


!     NOTES: 
!     kx,ky,kz are wavenumbers. 
!     ikx,y,z  are corresponding array indices.
!     KTX,Y,Z  are truncation wavenumbers.  
!     IKTX,Y,Z are corresponding indices.
!     zx,y,z is vorticity. u,v,w is velocity. 
!     nzx,y,z is nonlinear term in vorticity equation.
  

PROGRAM MAIN
  implicit none 
  include 'param.inc'

! ======== Flow Field Variables ==================
  complex, dimension(iktx,ikty,iktzp) :: zxok,zyok,zzok,ttok,zxnk,zynk,zznk,ttnk
  real,    dimension(n1d,n3d,n2dp)    :: zxor,zyor,zzor,ttor,zxnr,zynr,zznr,ttnr
  complex, dimension(iktx,ikty,iktzp) :: nzxk,nzyk,nzzk,nttk,uk,  vk,  wk
  real,    dimension(n1d,n3d,n2dp)    :: nzxr,nzyr,nzzr,nttr,ur,  vr,  wr
  complex, dimension(iktx,ikty,iktzp) :: rhzx,rhzy,rhzz,rhtt,geok,gw1k,gw2k

#if (SCHEME == 1) /* adding an extra field variable for Adam-Bashforth time-stepping */
  complex, dimension(iktx,ikty,iktzp) :: zx1k,zy1k,zz1k,tt1k
#endif


  integer, parameter :: nfmax = 1200  ! max number of forced modes
  complex :: avzx,avzy,avzz,avtt
  complex :: termzx,termzy,termzz,termtt,tzx,tzy,tzz,ttt
  complex :: u,v,w,c1,c2,c3,zi
  complex :: gtau(nfmax,4)
  real  :: aj,bj,bf2,cor,PrRa
  real  :: visch,viscz,v2h,v2z,dmz,dpz,r1,r2
  real  :: ki,Etot,ke,pe
  real  :: iwfr,kxwv,kywv,khwv,kzwv,iweng
  real  :: kf,alpha,tau,ampv,ampw,eps
  real  :: tstop,time,delt,d2,d12,robert,ek,kh_g
  real  :: kx,ky,kz,kh,khn,k2h,k2z,wk2,k,ktrunc_x,ktrunc_y,ktrunc_z
  real  :: kxa(iktx),kya(ikty),kza(iktz)
  real  :: L1,L2,L3
  real  :: ranno,fseed
  real  :: time1(2),time2(2),time3
  real  :: ts=0.,twopi
  real  :: etime
      
  integer :: ikx,iky,ikz,ikza
  integer :: L(iktx,ikty,iktzp)
  integer :: linear,rsflag,noutput
  integer :: irest,myrest,nstop,nout,ndump,nbig,nrsp,nmod,nrst,ntdump,nbigmod
  integer :: ilap,ilap2,iseed,kpTimeScaleSPC
  integer :: nt=0,io=0,nt0

  ! clean this up later please!
  integer :: ikxsp(iktxsp),ikysp(iktysp-1),ikzsp(iktzsp-1)
  real    :: kxsp(iktxsp),kysp(iktysp-1),kzsp(iktzsp-1)
  integer :: ikz1,ikzsp1,nkzsp
  integer :: dumpsp(4),nspdump,totnzsp,nspwv,nspwvevry

#if (BTbalance == 1)
  ! Baer-Tribbia balance model variables  
  complex :: Ord1gw1k,Ord1gw2k,Ord2gw1k,Ord2gw2k
  complex :: dApdt1k,dApdt2k,NbGAp1k,NbGAp2k
  real    :: wkh,wkh2,wkhn,freq
  ! flags for optional BaerTribbia-related diagnosis features
  integer :: BalUnbalTrnsFlag, BalUnbalSpcFlag, BalUnbalEngFlag, nBTout
#endif

#if (UnbalEvol == 1)
  complex, dimension(iktx,ikty,iktzp) :: Mtotk,Ntot1k,Ntot2k,NaMtMt1k,NaMtMt2k,NaGpMt1k,NaGpMt2k
  complex, dimension(iktx,ikty,iktzp) :: dalfapdt1,dalfapdt2,NMtAp1k,NMtAp2k,NGdAp1k,NGdAp2k
  complex, dimension(iktx,ikty,iktzp) :: NcdApdAp1k,NcdApdAp2k,dMdtk,Mtotoldk,NadMdM1k,NadMdM2k
  complex, dimension(iktx,ikty,iktzp) :: NaGpdM1k,NaGpdM2k,dAbdt1,dAbdt2,dAudt1,dAudt2,Ab1,Ab2
  complex :: NbMtAp1k,NbMtAp2k,FirstTrm1,FirstTrm2,NbGdAp1k,NbGdAp2k,SecndTrm1,SecndTrm2
  complex :: ThirdTrm1,ThirdTrm2,NaGdM1k,NaGdM2k,FourthTrm1,FourthTrm2
  integer :: nMtotOld,nMtotNew
  equivalence(dalfapdt1,dAudt1)
  equivalence(dalfapdt2,dAudt2)
#endif

#if (HCVC4movie == 1)
  complex, dimension(iktx,ikty,iktzp) :: wntdFldk
  real,    dimension(n1d,n3d,n2dp)    :: wntdFldr
  integer :: nskipH,nskipV,nHmax,nVmax,nRSout,nRSevery,justHor,dumprmsNGs,tran2gOn
  integer :: nRSGAb,nRSGAu,nRSGG,nRSAuAu,nRSAbAb,nRSAbAu,nRSdiv,nRSzz,nRSzzg,nRSzzab,nRSdau
  integer :: kpGAb,kpGAu,kpGG,kpAuAu,kpAbAb,kpAbAu,kpdiv,kpzz,kpzzg,kpzzab,kpdAu
  equivalence(wntdFldk,wntdFldr)
#endif

#if (UnbalEvol == 1) || (HCVC4movie == 1)
  complex, dimension(iktx,ikty,iktzp) :: zxhelpk,zyhelpk,zzhelpk,tthelpk
  real,    dimension(n1d,n3d,n2dp)    :: zxhelpr,zyhelpr,zzhelpr,tthelpr
  equivalence(zxhelpk,zxhelpr)
  equivalence(zyhelpk,zyhelpr)
  equivalence(zzhelpk,zzhelpr)
  equivalence(tthelpk,tthelpr)
#endif

#if (BTbalance == 1) || (HCVC4movie == 1)
  complex, dimension(iktx,ikty,iktzp) :: zxbalk,zybalk,zzbalk,ttbalk
  real,    dimension(n1d,n3d,n2dp)    :: zxbalr,zybalr,zzbalr,ttbalr
  complex, dimension(iktx,ikty,iktzp) :: zxunbalk,zyunbalk,zzunbalk,ttunbalk
  real,    dimension(n1d,n3d,n2dp)    :: zxunbalr,zyunbalr,zzunbalr,ttunbalr
  complex, dimension(iktx,ikty,iktzp) :: MaGGk,NaGG1k,NaGG2k,NcApAp1k,NcApAp2k,NGAp1k,NGAp2k
  complex, dimension(iktx,ikty,iktzp) :: NaMM1k,NaMM2k,NaGpM1k,NaGpM2k
  complex, dimension(iktx,ikty,iktzp) :: NtotGk,NGGk,NGAk,NAAk,NtotBGk,NGAbk,NAbAbk,NGAuk,NAuAuk,NAbAuk
  equivalence(zxbalk,zxbalr)
  equivalence(zybalk,zybalr)
  equivalence(zzbalk,zzbalr)
  equivalence(ttbalk,ttbalr)
  equivalence(zxunbalk,zxunbalr)
  equivalence(zyunbalk,zyunbalr)
  equivalence(zzunbalk,zzunbalr)
  equivalence(ttunbalk,ttunbalr)
  equivalence(NtotGk,MaGGk)
  equivalence(NGGk,NaGG1k)
  equivalence(NGAk,NaGG2k)
  equivalence(NAAk,NcApAp1k)
  equivalence(NtotBGk,NcApAp2k)
  equivalence(NGAbk,NGAp1k)
  equivalence(NAbAbk,NGAp2k)
  equivalence(NGAuk,NaMM1k)
  equivalence(NAuAuk,NaMM2k)
  equivalence(NAbAuk,NaGpM1k)
#endif


  
  external :: initialize,out,convol,spec,constr,velo,transf,ranno,forset,force,killwaves,killgeo,nltrmsgrms
 
! FFTW stuff
  integer, parameter :: fftw_forward         = -1, fftw_backward                  = 1
  integer, parameter :: fftw_real_to_complex = -1, fftw_complex_to_real           = 1
  integer, parameter :: fftw_estimate        =  0, fftw_measure=1,fftw_use_wisdom = 16
  integer, parameter :: fftw_in_place        =  8
  integer*8 :: plan2_rk,plan2_kr,plan1_rk,plan1_kr
  common/FTW2D/plan2_rk,plan2_kr
  common/FTW1D/plan1_rk,plan1_kr

#if (NETCDF == 1)
  include 'netcdf.inc'
! netcdf stuff for real space output:
  integer :: idnc,idtimes,idu,idv,idw,idth,idzx,idzy,idzz
  integer :: iputu,iputv,iputw,iputth,iputzx,iputzy,iputzz
  integer, dimension(4) :: ncstart,nccount
  common/netcdfrsp/idnc,ncstart,nccount,&
        iputu,iputv,iputw,iputth,iputzx,iputzy,iputzz, &
        idtimes,idu,idv,idw,idth,idzx,idzy,idzz

! netcdf stuff for restart file
  integer :: idnck,idtimesk,idzxk,idzyk,idzzk,idttk
  integer, dimension(5) :: ncstartk,nccountk
  common/netcdfrst/idnck,ncstartk,nccountk,idtimesk,idzxk,idzyk,idzzk,idttk

!! 3D spectrum of waves (and potentially other variables)
  integer :: idnsp,idkxsp,idkysp,idkzsp,idtsp
  integer, dimension(4) :: ncdimsp,ncstartsp,nccountsp
  integer, dimension(3) :: ncdimkk,ncstartkk,nccountkk
  integer :: idtimesp,idkkxx,idkkyy,idkkzz,idwvp,idwvm
  
  common/netcdfwvsp/idnsp,idkxsp,idkysp,idkzsp,idtsp,ncstartsp,nccountsp, &
       ncstartkk,nccountkk,idtimesp,idkkxx,idkkyy,idkkzz,idwvp,idwvm

!! limited number of modes vars
  integer :: idncm
  integer :: idtimesm,idzxm,idzym,idzzm,idttm,idwnxm,idwnym,idwnzm
  integer :: idmx,idmy,idmz,idmri,idmt
  integer, dimension(5) :: ncstartm,nccountm,ncdimsm
  integer, dimension(3) :: ncdim4wn
  common/netcdfmod/idncm,ncstartm,nccountm,idtimesm,idzxm,idzym,idzzm,idttm

#endif

! mpi stuff
#if (MPI == 1)
  include 'mpif.h'
#endif
  integer :: mype,size,ierror
  common/mpi/mype

  integer, parameter :: iscr=iktx*ikty*iktzp
  complex, dimension(iscr) :: zscr1,zscr2,zscr3  
  common/scratch/zscr1,zscr2,zscr3

! real/complex equivalences (in-place transforms)
  equivalence(zxok,zxor)
  equivalence(zyok,zyor)
  equivalence(zzok,zzor)
  equivalence(ttok,ttor)
  equivalence(zxnk,zxnr)
  equivalence(zynk,zynr)
  equivalence(zznk,zznr)
  equivalence(ttnk,ttnr)
  equivalence(nzxk,nzxr)
  equivalence(nzyk,nzyr)
  equivalence(nzzk,nzzr)
  equivalence(nttk,nttr)
  equivalence(uk,ur)
  equivalence(vk,vr)
  equivalence(wk,wr)

#if (KeepLM  == 1 )  
  ! storing limited modes for time frequency
  integer :: nFreqout
  real    :: khmode(nkhm)         ! the horizontal radious of each ring 
  real    :: kxmode(nmodeh,nkhm), kymode(nmodeh,nkhm), kzmode(nkhm,nkzm) !picked wavenumbers
  integer :: ikxmo(nmodeh,nkhm), ikymo(nmodeh,nkhm), ikzmo(nkhm,nkzm) !index of picked wavenumbers
  integer :: imyp(nkhm,nkzm)           !the number of processor for each vertical level 
#endif


#if (MPI == 1 )
! initialize mpi
  call mpi_init(ierror)
  call mpi_comm_size(mpi_comm_world,size,ierror)
  call mpi_comm_rank(mpi_comm_world,mype,ierror)
  if (mype.eq.0) print*,'size = ',size
  if (size.ne.npe) then
     print*,size,npe,'Wrong number of processors!'
     stop
  endif
#else
  mype = 0
  if (npe.ne.1) then
     print*,'Must have npe=1 for serial code'
     stop
  endif
#endif

! -----------------------------------------------------------------
! universal constants
  twopi    =  4.*asin(1.)
  zi       =  cmplx(0.,1.)
  
! -----------------------------------------------------------------
! Parameters
! -----------------------------------------------------------------


  ! ============= Aspect/Prandtl Ratio =============================================================
  PrRa     =  32.0
  
  ! ============= DOMAIN SIZE ======================================================================
  L1       = twopi                         ! domain size in (x) direction
  L2       = twopi                         !                (y)
  L3       = twopi/PrRa                    !                (z)
  
  ! ============ TRUNCATION WAVENUMBERS (automatically adjust based on 1/3 rule) ===================
  ktrunc_x = twopi/L1 * float(n1)/3.       ! dimensional truncation wavenumber (x) --> not index
  ktrunc_y = twopi/L2 * float(n2)/3.       !                                   (y)
  ktrunc_z = twopi/L3 * float(n3)/3.       !                                   (z)

  ! ============ Time Varaibles ====================================================================
  delt     = .1 / float(max(n1,n3))         ! timestep
  tstop    = 4.                             ! length of integration
  nstop    = ifix(tstop/delt)               ! number of timesteps

  ! ============ Outputing =========================================================================
  ! These variables set the frequency of outputting diagnostics in spc/trns/eng files AND in out routine (aka run.list file)
  noutput  = 100                           ! the total number of times that the results are ouputted in spc/trns/eng and run.list
  nout     = nstop/noutput                 ! every "nout" timesteps the results are ouputted in spc/trns/eng and run.list
  ! The frequency of outputting in the restart file
  nout = max(nout,1)
  ndump    = nstop/1                       ! dump restart file every ndump timesteps 
                                           !---> if ndump=nstop then only the last time step is stored in Zk?.out
  nrst     = nstop/ndump                   ! number of restart fields to dump
  
  ! ============ DISSIPATION AND DAMPING  ==========================================================
  ilap     = 4                             ! hyperviscosity order (1=regular viscosity)
  visch    = 2e-14                         ! horizontal viscosity coefficient
  viscz    = visch*(1/PrRa)**(2*ilap)      ! vertical viscosity coefficient
  robert   = 5.e-3                         ! Robert filter parameter
  ek       = 1./(3600.*24.)                ! vshf damping time scale

  ! ============ CORIOLIS AND BUOYANCY  ============================================================
  !aj       = 339.53466/3                   ! thermal expansivity
  aj       = 1440 
  bj       = aj                            ! background theta gradient
  bf2      = aj*bj                         ! Brunt-Vaisalla frequency squared
  cor      = sqrt(aj*bj)/PrRa              ! Coriolis parameter

  ! ============ INITIALIZATION FEATURES  ==========================================================
  irest    = 1                             ! flag for starting from restart file
  iwfr     = sqrt(2.0)*cor                   ! frequency of the wave initially shooted
  kxwv     = 3.0                           ! kx of the initial wave  
  kywv     = 0.0                           ! ky of the initial wave
  khwv     = sqrt(kxwv**2+kywv**2)         ! kh of the initial wave
  kzwv     = khwv * sqrt((bf2 - iwfr**2)/(iwfr**2 - cor**2)) ! derive kz of the initial wave based on dispersion rel.
  iweng    = 0.1                           ! energy of the initial wave
  nspwv    = 200                           ! number of time the 3D wave spectrum is dumped
  nspwvevry= nstop/nspwv
  ki       = (L1/twopi) * (20.)            ! wavenumber for the peak of  initial conditions
  etot     = 1.e-0                         ! initial energy
  pe       = 0.                            ! initial potential energy 
  ke       = etot - pe                     ! initial kinetic energy
  myrest   = irest                         ! to use different random seed (not sure what it is)
  
  ! ============ FORCING PARAMETERS  ===============================================================
  kf       = ki                            ! wavenumber for forcing
  tau      = 25.*delt                      ! time scale for forcing
  eps      = 6.e-5                         ! dissipation rate for forcing amplitude
  ! -----> if ampw and ampv = 0 forget about the rest (no forcing case)
  ampw     = 0.                            ! random wave forcing amplitude
  ampv     = 0.0*sqrt(eps/tau)             ! random vort forcing amplitude

  ! ============ This part is messy for now ---> will be cleaned-up later 
  nbig     = 0                             ! dump real space fields every nbig timesteps
  linear   = 0                             ! flag for omitting nonlinear terms
  rsflag   = 0                             ! flag for outputting real space fields
  nrsp     = rsflag*(nstop/nbig)           ! number of real space fields to dump
  nbigmod  = 0                             ! keep limited no. of modes every nbigmod timesteps
  nmod     = 0 !nstop/nbigmod              ! number of times that limited number of modes are kept


  ! ============ Storing Limited Number of Modes For Frequency Spectra =============================
#if (KeepLM == 1)
  nFreqout = 2                    ! Limited modes are stored every "nFreqout" time steps
  khmode = (/ 8, 20, 40, 60, 100, 150/)  !make sure the number of elements = nkhm (in param.inc)
  kzmode(1,:) = (/ 0, 1, 8, 30 /)          !make sure the number of elements = nkzm 
  kzmode(2,:) = (/ 0, 1, 20, 30 /)
  kzmode(3,:) = (/ 0, 1, 20, 30 /)
  kzmode(4,:) = (/ 0, 1, 20, 30 /)
  kzmode(5,:) = (/ 0, 1, 20, 30 /)
  kzmode(6,:) = (/ 0, 1, 20, 30 /)
#endif
  
  ! =========== OPTIONAL DIAGNOSIS FEATURES THAT CAN BE OUTPUT =====================================
  ! ------ NOTE THESE FEATURE NEED EXTRA OF MEMORY AND INCREASE THE COMPUATION TIME  
  ! ------ SHOULD BE TURN ON IF IT IS NECESSARY


#if (BTbalance == 1)  
  ! ======= Spectrum of/Transfer to balanced ageostrophic and unbalanced agoestrophic separately 
 
  BalUnbalTrnsFlag  =  1          ! Flag for turning on the calculation of the transfer to 
                                   ! balanced and unlanced ageostrophic modes separately
  BalUnbalSpcFlag   =  1          ! Flag for turning on the calculation of  
                                   ! balanced and unlanced ageostrophic spectra separately
  BalUnbalEngFlag   =  1          ! Flag for turning on the calculation of (un)balanced energies

  nBTout = nout                    ! Baer-Tribbia balance model executed every "nout" time steps
#endif
  
#if (HCVC4movie == 1)
  ! ============ Horizontal and Vertical Cuts in Real Space ========================================
  ! This is for making movies. Just cuts not the total field 
  nRSout    = 100                          ! the total number of times that the verical and horizontal cuts are dumped
  nRSevery  = nstop/nRSout                 ! every "nRSout" dump vertical and horizontal real space cuts   
  nHmax     = n1                           ! maximum of horizontal physical points that are kept 
  nVmax     = n3                           !            vertical
  justHor   = 0                            ! Just keep the horizontal slices if  = 1, 
  nskipH    = 2                            ! every nskipH is kept in horizontal direction
  nskipV    = 1                            ! every nskipH is kept in vertical
  kpGAb     = 0                            ! A flag to keep GxA_b's in the RHS of dG/dt 
  kpGAu     = 0                            ! A flag to keep GxA_u's in the RHS of dG/dt 
  kpGG      = 0                            ! A flag to keep GxG's in the RHS of dG/dt
  kpAuAu    = 0                            ! A flag to keep A_uxA_u's in the RHS of dG/dt
  kpAbAb    = 0                            ! A flag to keep A_bxA_b's in the RHS of dG/dt
  kpAbAu    = 0                            ! A flag to keep AbxAu's in the RHS of dG/dt 
  kpdiv     = 0                            ! A flag to keep Horizontal Divergencet
  kpzz      = 1                            ! A flag to keep the total vertical vorticity 
  kpzzg     = 0                            ! A flag to keep vertical vorticity just due to G
  kpzzab    = 0                            ! A flag to keep vertical vorticity just due to A_b
  kpdau     = 1
#if (UnbalEvol == 0)
  kpdau     = 0
#endif  
  Dumprmsngs= 1                            ! A flag to dump the rms NGG, NGAu ...
  tran2gOn  = 1                            ! A flag to dump transfer stuff to different nonlinear terms in dG/dt
#endif
  
  kpTimeScaleSPC = 1 
#if (UnbalEvol == 0)
  kpTimeScaleSPC = 0
#endif

! -----------------------------------------------------------------  
! Open files for diagnostics
! -----------------------------------------------------------------
  if (mype.eq.0) then
!     open (57,file='spc2d.dat',form='formatted')  
     open (51,file='spcz.dat', form='formatted')
     open (52,file='spch.dat', form='formatted')
     open (46,file='eng.dat',  form='formatted')
     open (48,file='trnh.dat', form='formatted')
     open (49,file='trnz.dat', form='formatted')
     open (79,file='eps.dat',  form='formatted')
     open (393,file='kk3D.dat',  form='formatted')
     if (kpTimeScaleSPC .eq. 1) open (89,file='TSspch.dat',  form='formatted')
     if (kpTimeScaleSPC .eq. 1) open (99,file='TSspcz.dat',  form='formatted')

#if (HCVC4movie == 1)
     if (kpzz  .eq.1)  open (111,file='HCzz.out',     form='formatted')
     if (kpzzg .eq.1)  open (112,file='HCzzg.out',    form='formatted')
     if (kpzzab.eq.1)  open (113,file='HCzzab.out',   form='formatted')
     if (kpGAb .eq.1)  open (114,file='HCgab.out',    form='formatted')
     if (kpGAu .eq.1)  open (115,file='HCgau.out',    form='formatted')
     if (kpGG  .eq.1)  open (116,file='HCgg.out',     form='formatted')
     if (kpAbAb.eq.1)  open (117,file='HCabab.out',   form='formatted')
     if (kpAbAu.eq.1)  open (118,file='HCabau.out',   form='formatted')
     if (kpAuAu.eq.1)  open (119,file='HCauau.out',   form='formatted')
     if (kpdiv .eq.1)  open (110,file='HCdiv.out',    form='formatted')
     if (kpdau .eq.1) then 
        open (151,file='HCdau1.out',    form='formatted')
        open (152,file='HCdau2.out',    form='formatted')
     endif
     if (justHor.ne.1) then
        if (kpzz  .eq.1)  open (121,file='VCzz.out',     form='formatted')
        if (kpzzg .eq.1)  open (122,file='VCzzg.out',    form='formatted')
        if (kpzzab.eq.1)  open (123,file='VCzzab.out',   form='formatted')
        if (kpGAb .eq.1)  open (124,file='VCgab.out',    form='formatted')
        if (kpGAu .eq.1)  open (125,file='VCgau.out',    form='formatted')
        if (kpGG  .eq.1)  open (126,file='VCgg.out',     form='formatted')
        if (kpAbAb.eq.1)  open (127,file='VCabab.out',   form='formatted')
        if (kpAbAu.eq.1)  open (128,file='VCabau.out',   form='formatted')
        if (kpAuAu.eq.1)  open (129,file='VCauau.out',   form='formatted')
        if (kpdiv .eq.1)  open (120,file='VCdiv.out',    form='formatted')
        if (kpdau .eq.1)  then
           open (161,file='VCdau1.out',    form='formatted')
           open (162,file='VCdau2.out',    form='formatted')
        endif
     endif
     if (dumprmsNGs.eq.1) open (131,file='rmsNG.dat',    form='formatted')
     if (tran2gOn.eq.1)   open (141,file='trnh2g.dat',   form='formatted')
#endif


#if (BTbalance == 1)  
     if (BalUnbalTrnsFlag.eq.1) then
        open (61,file='trnhBal.dat', form='formatted')
        open (62,file='trnzBal.dat', form='formatted')
        open (63,file='trnhUnb.dat', form='formatted')
        open (64,file='trnzUnb.dat', form='formatted')
     endif
     if (BalUnbalSpcFlag.eq.1) then
        open (71,file='spchBal.dat', form='formatted')
        open (72,file='spczBal.dat', form='formatted')
        open (73,file='spchUnb.dat', form='formatted')
        open (74,file='spczUnb.dat', form='formatted')
     endif
     if (BalUnbalEngFlag.eq.1) then
        open (101,file='engBal.dat', form='formatted')
        open (102,file='engUnb.dat', form='formatted')
        open (91,file='epsBal.dat',  form='formatted')
        open (92,file='epsUnb.dat',  form='formatted')
     endif
#endif
#if (UnbalEvol == 1)
     open (211,file='trnh2unb.dat', form='formatted')
     open (212,file='trnz2unb.dat', form='formatted')
#endif
#if (KeepLM == 1)
     open (41,file='timeseries.dat', form='formatted')
#endif
  endif
    
! -----------------------------------------------------------------
! Set up a few last things
! -----------------------------------------------------------------

! Initialize FFTW
  if (mype.eq.0) print*,'Initializing FFTW'
  call rfftw2d_f77_create_plan(plan2_rk,n1,n3,fftw_real_to_complex,fftw_estimate+fftw_in_place)
  call rfftw2d_f77_create_plan(plan2_kr,n1,n3,fftw_complex_to_real,fftw_estimate+fftw_in_place)
  call    fftw_f77_create_plan(plan1_rk,n2,   fftw_forward,        fftw_estimate+fftw_in_place)
  call    fftw_f77_create_plan(plan1_kr,n2,   fftw_backward,       fftw_estimate+fftw_in_place)

! Initialize NETCDF
#if (NETCDF ==1)
  iputu  = 0
  iputv  = 0
  iputw  = 0
  iputth = 0
  iputzx = 0
  iputzy = 0
  iputzz = 0
  call ncprep(nrsp,nrst,nspwv)

#endif

! Set random seed
  fseed = ranno(1946 + myrest + mype)
  fseed = ranno(0)/6.29 + 0.5
  iseed = ifix(10000.*fseed) + 1
  fseed = ranno(iseed)

! Misc defs
  
  ilap2 = 2*ilap
  d2    = 2.*delt
  d12   = delt/12.
  alpha = exp(-delt/tau)
  k2h   = visch*delt
  k2z   = viscz*delt
  v2h   = visch*delt
  v2z   = viscz*delt
  ek    = ek*delt

! -----------------------------------------------------------------
! Initialize fields, L, etc.
! -----------------------------------------------------------------

  call initialize(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k,         &
            zxor,zyor,zzor,ur,vr,wr,ttor,delt,ke,pe,L,irest,     &
            kxa,kya,kza,ktrunc_x,ktrunc_y,ktrunc_z,aj,bj,cor,ts,L1,L2,L3,iwfr,kxwv,kywv,khwv,kzwv,iweng)
#if (KeepLM == 1)
  call initkmodes(kxa,kya,kza,kxmode,kymode,kzmode,khmode,imyp,ikxmo,ikymo,ikzmo)
#endif
  
#if (HCVC4movie == 1)
  ! counter of # of Real Space slices that is dumped. gets increased each time you dump!
  if (kpzz  .eq.1)  nRSzz   = 0
  if (kpzzg .eq.1)  nRSzzg  = 0
  if (kpzzab.eq.1)  nRSzzab = 0
  if (kpGAb .eq.1)  nRSGAb  = 0
  if (kpGAu .eq.1)  nRSGAu  = 0
  if (kpGG  .eq.1)  nRSGG   = 0
  if (kpAbAb.eq.1)  nRSAbAb = 0
  if (kpAbAu.eq.1)  nRSAbAu = 0
  if (kpAuAu.eq.1)  nRSAuAu = 0
  if (kpdiv .eq.1)  nRSdiv  = 0
  if (kpdAu .eq.1)  nRSdau  = 0
#endif


  ! initialize a few last things
  ! ts = 5.0
  time = ts + nt*delt
  !print*, 'time after init', time
  gtau = cmplx(0.,0.)

  call initksp(kxa,kya,kza,kxsp,kysp,kzsp,ikxsp,ikysp,ikzsp,L1,L2,L3,dumpsp,totnzsp)
  nspdump = 1
  call wvspc3D(gw1k,gw2k,time,nspdump,kxa,kya,kza,ikxsp,ikysp,ikzsp,dumpsp,totnzsp)


#if (SCHEME == 1)
  zx1k=cmplx(0.,0.)
  zy1k=cmplx(0.,0.)
  zz1k=cmplx(0.,0.)
  tt1k=cmplx(0.,0.)
#endif

! -----------------------------------------------------------------
! Print outall Parameters
! -----------------------------------------------------------------
  if (mype.eq.0) then
     print*,'                '
     print*,'                '
     print*,'Parameters: --------------------------------'
     print*,'              N1,N2,N3 =  ', n1,n2,n3
     print*,'              L1,L2,L3 =  ', L1,L2,L3
     print*,'              dx,dy,dz =  ', L1/n1,L2/n2,L3/n3
     print*,'              IKTX,Y,Z =  ', iktx,ikty,iktz
     print*,'                     KT = ', ktrunc_x,ktrunc_y,ktrunc_z
     print*,'Order of Laplacian diss.=  ', ilap 
     print*,'                  VISCH = ', visch
     print*,'                  VISCZ = ', viscz
     print*,'              tau_VISCH = ', (visch*ktrunc_x**ilap2)**(-1)
     print*,'              tau_VISCZ = ', (viscz*ktrunc_z**ilap2)**(-1)
     if (linear.eq.1) print*,'Nonlinear terms switched off.' 
     print*,'               Timestep = ', delt
     print*,'     Integration length = ', nstop*delt,' = ', nstop,' DT.'
     print*,'       Output frequency = ', nout*delt, ' = ', nout,' dt.'
     print*,'  Rspace dump frequency = ', nbig*delt, ' = ', nbig,' dt.'
     print*,'  '
     if (irest.eq.0) then
        print*,' Starting from ICs '
        print*,' Initial kinetic energy = ', ke
        print*,'   "    potential  "    = ', pe
        print*,'   "      total    "    = ', ke + pe
     else
        print*,'  Starting from restart file'
        print*,'  Restart from rec.     = ', irest
        print*,'  Restart from time     = ', ts
     endif
     print*,'    '
     print*,'    Thermal expansivity = ', aj
     print*,'    Vertical T gradient = ', bj
     print*,'    Brunt-Vaisala freq. = ', sqrt(bf2)
     print*,'     Coriolis parameter = ', cor
     print*,'Robert filter parameter = ', robert
     print*,'   '
     if (max(ampv,ampw).ne.0.) then
        print*,' Random forcing: ' 
        print*,'      Amplitude (vortical) = ', ampv
        print*,'      Amplitude (wave)     = ', ampw
        print*,'                Wavenumber = ', kf
        print*,'                    Memory = ', tau,' = ', ifix(tau/delt),' DT.'
     else
        print*,'      No forcing.'
     endif
     print*,'                '
  endif
 
! -----------------------------------------------------------------
! Diagnostics on ICs
! -----------------------------------------------------------------
  if (irest.ge.0) then
     call wtoab(zxok,zyok,zzok,ttok,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
     call spec(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,2,52,kh_g)
     call spec(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,3,51,kh_g)
 
!    call spec2d(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,57)
     call out(zxok,zyok,zzok,ttok,uk,vk,wk,geok,gw1k,gw2k,   &
          nt,delt,ts,nstop,io,L,kxa,kya,kza,aj,bj,cor,visch,viscz,ilap,L1,L2,L3,kh_g) 
!!$#if (HCVC4movie == 1)
!!$     call dumprealHCVCzz(zznk,zznr,nRSslice,noutput,time,58,59,nskipH,nskipV,nHmax,nVmax)
!!$     call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L,  &
!!$       uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
!!$     call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NNN1,NNN2,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
!!$#endif       

  endif ! irest

! -----------------------------------------------------------------
! Beginning of the first timestep.
! -----------------------------------------------------------------
  nt   = 1
  time = ts + nt*delt

  if (linear.ne.1) then
     call constr (zxok,zyok,zzok,ttok,nzxk,nzyk,nzzk,nttk,L,  &
          uk,vk,wk,ur,vr,wr,zxor,zyor,zzor,ttor,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
  else
     nzxk = cmplx(0.,0.)
     nzyk = cmplx(0.,0.)
     nzzk = cmplx(0.,0.)
     nttk = cmplx(0.,0.)
  endif

  if (max(ampv,ampw).ne.0.) then
     call force(nzxk,nzyk,nzzk,nttk,ampv,ampw,gtau,L,nfmax,alpha,kf,kxa,kya,kza,aj,bj,cor,time,L1,L2,L3)
  endif

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx=1,iktx

           if ( L(ikx,iky,ikz).ne.1 ) then
              zxnk(ikx,iky,ikz) = cmplx(0.,0.)
              zynk(ikx,iky,ikz) = cmplx(0.,0.)
              zznk(ikx,iky,ikz) = cmplx(0.,0.)
              ttnk(ikx,iky,ikz) = cmplx(0.,0.)
           else
              kx  = kxa(ikx)
              kh  = sqrt( kx*kx + ky*ky )
              khn = kh * L1/twopi
              wk2 = kx*kx + ky*ky + kz*kz
              k  = sqrt(wk2)

              c1 = +  ky * zzok(ikx,iky,ikz) -  kz * zyok(ikx,iky,ikz)
              c2 = +  kz * zxok(ikx,iky,ikz) -  kx * zzok(ikx,iky,ikz)
              c3 = +  kx * zyok(ikx,iky,ikz) -  ky * zxok(ikx,iky,ikz)
              u = zi * c1 / wk2
              v = zi * c2 / wk2
              w = zi * c3 / wk2

#if (SPHTRUNC == 1)
              r1 = v2h/delt*wk2**ilap
              r2 = k2h/delt*wk2**ilap
#else
              r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
              r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif

              if (khn.lt.0.1) then
                r1 = r1 + ek/delt
                r2 = r2 + ek/delt
              endif

              termzx = nzxk(ikx,iky,ikz) + aj*zi*ky*ttok(ikx,iky,ikz) + cor*zi*kz*u - zxok(ikx,iky,ikz)*r1
              termzy = nzyk(ikx,iky,ikz) - aj*zi*kx*ttok(ikx,iky,ikz) + cor*zi*kz*v - zyok(ikx,iky,ikz)*r1
              termzz = nzzk(ikx,iky,ikz)                              + cor*zi*kz*w - zzok(ikx,iky,ikz)*r1
              termtt = nttk(ikx,iky,ikz)                              - bj*w        - ttok(ikx,iky,ikz)*r2

              rhzx(ikx,iky,ikz) = termzx
              zxnk(ikx,iky,ikz) = zxok(ikx,iky,ikz) + delt*termzx
              rhzy(ikx,iky,ikz) = termzy
              zynk(ikx,iky,ikz) = zyok(ikx,iky,ikz) + delt*termzy
              rhzz(ikx,iky,ikz) = termzz
              zznk(ikx,iky,ikz) = zzok(ikx,iky,ikz) + delt*termzz
              rhtt(ikx,iky,ikz) = termtt
              ttnk(ikx,iky,ikz) = ttok(ikx,iky,ikz) + delt*termtt

           endif ! L
        enddo
     enddo
  enddo

  if (linear.ne.1) then
     call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L, &
          uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
  else
     nzxk = cmplx(0.,0.)
     nzyk = cmplx(0.,0.)
     nzzk = cmplx(0.,0.)
     nttk = cmplx(0.,0.)
  endif

  if (max(ampv,ampw).ne.0.) then
     call force(nzxk,nzyk,nzzk,nttk,ampv,ampw,gtau,L,nfmax,alpha,kf,kxa,kya,kza,aj,bj,cor,time,L1,L2,L3)
  endif

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,ikty
        ky = kya(iky)             
        do ikx = 1,iktx 
           if ( L(ikx,iky,ikz).ne.1 ) then
              zxnk(ikx,iky,ikz) = cmplx(0.,0.)
              zynk(ikx,iky,ikz) = cmplx(0.,0.)
              zznk(ikx,iky,ikz) = cmplx(0.,0.)
              ttnk(ikx,iky,ikz) = cmplx(0.,0.)
           else
              kx  = kxa(ikx)             
              kh  = sqrt( kx*kx + ky*ky )
              khn = kh * L1/twopi
              wk2 = kx*kx + ky*ky + kz*kz
              k   = sqrt( wk2 )

              c1 = +  ky*zznk(ikx,iky,ikz) -  kz*zynk(ikx,iky,ikz)
              c2 = +  kz*zxnk(ikx,iky,ikz) -  kx*zznk(ikx,iky,ikz)
              c3 = +  kx*zynk(ikx,iky,ikz) -  ky*zxnk(ikx,iky,ikz)
              u  = zi * c1 / wk2
              v  = zi * c2 / wk2
              w  = zi * c3 / wk2

#if (SPHTRUNC == 1)
              r1 = v2h/delt*wk2**ilap
              r2 = k2h/delt*wk2**ilap
#else 
              r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
              r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif 

              if (khn.lt.0.1) then
                r1 = r1 + ek/delt
                r2 = r2 + ek/delt
              endif

              termzx = nzxk(ikx,iky,ikz) + aj*zi*ky*ttnk(ikx,iky,ikz) + cor*zi*kz*u - zxok(ikx,iky,ikz)*r1
              termzy = nzyk(ikx,iky,ikz) - aj*zi*kx*ttnk(ikx,iky,ikz) + cor*zi*kz*v - zyok(ikx,iky,ikz)*r1
              termzz = nzzk(ikx,iky,ikz)                              + cor*zi*kz*w - zzok(ikx,iky,ikz)*r1
              termtt = nttk(ikx,iky,ikz)                              - bj*w        - ttok(ikx,iky,ikz)*r2

#if (SCHEME == 0)
              zxnk(ikx,iky,ikz) = zxok(ikx,iky,ikz)  +   delt*(termzx + rhzx(ikx,iky,ikz))/2.
              zynk(ikx,iky,ikz) = zyok(ikx,iky,ikz)  +   delt*(termzy + rhzy(ikx,iky,ikz))/2.
              zznk(ikx,iky,ikz) = zzok(ikx,iky,ikz)  +   delt*(termzz + rhzz(ikx,iky,ikz))/2.
              ttnk(ikx,iky,ikz) = ttok(ikx,iky,ikz)  +   delt*(termtt + rhtt(ikx,iky,ikz))/2.
#endif
#if (SCHEME == 1) 
              zxok(ikx,iky,ikz) = zxok(ikx,iky,ikz)  +   delt*(termzx + rhzx(ikx,iky,ikz))/2.
              zyok(ikx,iky,ikz) = zyok(ikx,iky,ikz)  +   delt*(termzy + rhzy(ikx,iky,ikz))/2.
              zzok(ikx,iky,ikz) = zzok(ikx,iky,ikz)  +   delt*(termzz + rhzz(ikx,iky,ikz))/2.
              ttok(ikx,iky,ikz) = ttok(ikx,iky,ikz)  +   delt*(termtt + rhtt(ikx,iky,ikz))/2.

              zx1k(ikx,iky,ikz) = rhzx(ikx,iky,ikz) + zxok(IKX,IKY,IKZ)*r1
              zy1k(ikx,iky,ikz) = rhzy(ikx,iky,ikz) + zyok(IKX,IKY,IKZ)*r1
              zz1k(ikx,iky,ikz) = rhzz(ikx,iky,ikz) + zzok(IKX,IKY,IKZ)*r1
              tt1k(ikx,iky,ikz) = rhtt(ikx,iky,ikz) + ttok(IKX,IKY,IKZ)*r2

#endif
           endif ! L
        enddo
     enddo
  enddo



#if (SCHEME == 1) 
! -----------------------------------------------------------------
! Beginning of the second timestep (only necessary for AB3).
! -----------------------------------------------------------------
  nt   = 2
  time = ts + nt*delt

  if (LINEAR.ne.1) then
     call constr (zxok,zyok,zzok,ttok,nzxk,nzyk,nzzk,nttk,L, &
          uk,vk,wk,ur,vr,wr,zxor,zyor,zzor,ttor,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
  else
     nzxk = cmplx(0.,0.)
     nzyk = cmplx(0.,0.)
     nzzk = cmplx(0.,0.)
     nttk = cmplx(0.,0.)
  endif

  if (max(ampv,ampw).ne.0.) then
     call force(nzxk,nzyk,nzzk,nttk,AMPV,AMPW,gtau,L,nfmax,alpha,kf,kxa,kya,kza,aj,bj,cor,time,L1,L2,L3)
  endif

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx=1,iktx

           if ( L(ikx,iky,ikz).ne.1 ) then
              zxnk(ikx,iky,ikz) = cmplx(0.,0.)
              zynk(ikx,iky,ikz) = cmplx(0.,0.)
              zznk(ikx,iky,ikz) = cmplx(0.,0.)
              ttnk(ikx,iky,ikz) = cmplx(0.,0.)
           else
              kx  = kxa(ikx)
              kh  = sqrt( kx*kx + ky*ky )
              khn = kh * L1/twopi
              wk2 = kx*kx + ky*ky + kz*kz
              k  = sqrt(wk2)

              c1 = +  ky * zzok(ikx,iky,ikz) -  kz * zyok(ikx,iky,ikz)
              c2 = +  kz * zxok(ikx,iky,ikz) -  kx * zzok(ikx,iky,ikz)
              c3 = +  kx * zyok(ikx,iky,ikz) -  ky * zxok(ikx,iky,ikz)
              u = ZI * c1 / wk2
              v = ZI * c2 / wk2
              w = ZI * c3 / wk2

#if (SPHTRUNC == 1)
              r1 = V2H/delt*wk2**ILAP
              r2 = K2H/delt*wk2**ILAP
#else
              r1 = V2H/delt*kh**ILAP2 + V2Z/delt*kz**ILAP2
              r2 = K2H/delt*kh**ILAP2 + K2Z/delt*kz**ILAP2
#endif

              if (khn.lt.0.1) then
                r1 = r1 + ek/delt
                r2 = r2 + ek/delt
              endif

              termzx = nzxk(ikx,iky,ikz) + aj*ZI*ky*ttok(ikx,iky,ikz) + cor*ZI*kz*u - zxok(IKX,IKY,IKZ)*r1
              termzy = nzyk(ikx,iky,ikz) - aj*ZI*kx*ttok(ikx,iky,ikz) + cor*ZI*kz*v - zyok(IKX,IKY,IKZ)*r1
              termzz = nzzk(ikx,iky,ikz)                              + cor*ZI*kz*w - zzok(ikx,iky,ikz)*r1
              termtt = nttk(ikx,iky,ikz)                              - bj*w        - ttok(ikx,iky,ikz)*r2

              rhzx(ikx,iky,ikz) = termzx
              zxnk(ikx,iky,ikz) = zxok(ikx,iky,ikz) + delt*termzx
              rhzy(ikx,iky,ikz) = termzy
              zynk(ikx,iky,ikz) = zyok(ikx,iky,ikz) + delt*termzy
              rhzz(ikx,iky,ikz) = termzz
              zznk(ikx,iky,ikz) = zzok(ikx,iky,ikz) + delt*termzz
              rhtt(ikx,iky,ikz) = termtt
              ttnk(ikx,iky,ikz) = ttok(ikx,iky,ikz) + delt*termtt
           endif ! L
        enddo
     enddo
  enddo

  if (linear.ne.1) then
     call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L,  &
          uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
  else
     nzxk = cmplx(0.,0.)
     nzyk = cmplx(0.,0.)
     nzzk = cmplx(0.,0.)
     nttk = cmplx(0.,0.)
  endif

  if (max(ampv,ampw).ne.0.) then
     call force(nzxk,nzyk,nzzk,nttk,AMPV,AMPW,gtau,L,nfmax,alpha,kf,kxa,kya,kza,aj,bj,cor,time,L1,L2,L3)
  endif

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,ikty
        ky = kya(iky)             
        do ikx = 1,iktx 
           if ( L(ikx,iky,ikz).ne.1 ) then
              zxnk(ikx,iky,ikz) = cmplx(0.,0.)
              zynk(ikx,iky,ikz) = cmplx(0.,0.)
              zznk(ikx,iky,ikz) = cmplx(0.,0.)
              ttnk(ikx,iky,ikz) = cmplx(0.,0.)
           else
              kx  = kxa(ikx)             
              kh  = sqrt( kx*kx + ky*ky )
              khn = kh * L1/twopi
              wk2 = kx*kx + ky*ky + kz*kz
              k   = sqrt( wk2 )

              c1 = +  ky*zznk(ikx,iky,ikz) -  kz*zynk(ikx,iky,ikz)
              c2 = +  kz*zxnk(ikx,iky,ikz) -  kx*zznk(ikx,iky,ikz)
              c3 = +  kx*zynk(ikx,iky,ikz) -  ky*zxnk(ikx,iky,ikz)
              u  = zi * c1 / wk2
              v  = zi * c2 / wk2
              w  = zi * c3 / wk2

#if (SPHTRUNC == 1)
              r1 = v2h/delt*wk2**ilap
              r2 = k2h/delt*wk2**ilap
#else 
              r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
              r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif 

              if (khn.lt.0.1) then
                r1 = r1 + ek/delt
                r2 = r2 + ek/delt
              endif
              
              termzx = nzxk(ikx,iky,ikz) + aj*zi*ky*ttnk(ikx,iky,ikz) + cor*zi*kz*u - zxok(ikx,iky,ikz)*r1
              termzy = nzyk(ikx,iky,ikz) - aj*zi*kx*ttnk(ikx,iky,ikz) + cor*zi*kz*v - zyok(ikx,iky,ikz)*r1
              termzz = nzzk(ikx,iky,ikz)                              + cor*zi*kz*w - zzok(ikx,iky,ikz)*r1
              termtt = nttk(ikx,iky,ikz)                              - bj*w        - ttok(ikx,iky,ikz)*r2

              zxnk(ikx,iky,ikz) = zxok(ikx,iky,ikz)  +   delt*(termzx + rhzx(ikx,iky,ikz))/2.
              zynk(ikx,iky,ikz) = zyok(ikx,iky,ikz)  +   delt*(termzy + rhzy(ikx,iky,ikz))/2.
              zznk(ikx,iky,ikz) = zzok(ikx,iky,ikz)  +   delt*(termzz + rhzz(ikx,iky,ikz))/2.
              ttnk(ikx,iky,ikz) = ttok(ikx,iky,ikz)  +   delt*(termtt + rhtt(ikx,iky,ikz))/2.

              zxok(ikx,iky,ikz) = rhzx(ikx,iky,ikz)  + zxok(IKX,IKY,IKZ)*r1
              zyok(ikx,iky,ikz) = rhzy(ikx,iky,ikz)  + zyok(IKX,IKY,IKZ)*r1
              zzok(ikx,iky,ikz) = rhzz(ikx,iky,ikz)  + zzok(IKX,IKY,IKZ)*r1
              ttok(ikx,iky,ikz) = rhtt(ikx,iky,ikz)  + ttok(IKX,IKY,IKZ)*r2
                  
           endif ! L
        enddo
     enddo
  enddo
#endif


  time3 = etime(time1)

#if (UnbalEvol == 1)
!------------ creating Mtotoldk ----------------------------------
  call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L, &
          uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)     
  call wtoab(nzxk,nzyk,nzzk,nttk,Mtotoldk,Ntot1k,Ntot2k,rhzx,rhzy,rhzz,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

#if (SCHEME == 0) 
  nMtotOld=1   ! leapfrog: start at nt=2
#endif
#if (SCHEME == 1) 
  nMtotOld=2   ! AB3: start at nt=3
#endif

#endif



! -----------------------------------------------------------------
!                        Subsequent Timesteps
! -----------------------------------------------------------------
#if (SCHEME == 0) 
  nt0=2   ! leapfrog: start at nt=2
#endif
#if (SCHEME == 1) 
  nt0=3   ! AB3: start at nt=3
#endif

  do nt = nt0,NSTOP
     time = ts + nt*delt

     if (LINEAR.ne.1) then
        call convol (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L,uk,vk,wk,ur,vr,wr, &
             rhzx,rhzy,rhzz,rhtt,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
     else
        nzxk = cmplx(0.,0.)
        nzyk = cmplx(0.,0.)
        nzzk = cmplx(0.,0.)
        nttk = cmplx(0.,0.)
     endif

     if (max(ampv,ampw).ne.0.) then
        call force(nzxk,nzyk,nzzk,nttk,ampv,ampw,gtau,L,nfmax,alpha,kf,kxa,kya,kza,aj,bj,cor,time,L1,L2,L3)
     endif

     call velo (zxnk,zynk,zznk,rhzx,rhzy,rhzz,L,kxa,kya,kza)

     do ikz = 1,iktzp
        ikza = mype*iktzp+ikz
        kz = kza(ikza)
        do iky = 1,ikty
           ky = kya(iky)
           do ikx = 1,iktx
              kx  = kxa(ikx)
              kh  = sqrt( kx*kx + ky*ky )
              khn = kh * L1/twopi
              wk2 = kx*kx + ky*ky + kz*kz
              k   = sqrt( wk2 )
! Equations:
              termzx = nzxk(ikx,iky,ikz) +  aj*ZI*ky*ttnk(ikx,iky,ikz) + cor*ZI*kz*rhzx(ikx,iky,ikz)
              termzy = nzyk(ikx,iky,ikz) -  aj*ZI*kx*ttnk(ikx,iky,ikz) + cor*ZI*kz*rhzy(ikx,iky,ikz)
              termzz = nzzk(ikx,iky,ikz)                               + cor*ZI*kz*rhzz(ikx,iky,ikz)
              termtt = nttk(ikx,iky,ikz)                               -        bj*rhzz(ikx,iky,ikz)


#if (SCHEME == 0) 
! leapfrog timestepping
#if (SPHTRUNC == 1)
              r1 = v2h*wk2**ilap
              r2 = k2h*wk2**ilap
#else
              r1 = v2h*kh**ilap2 + v2z*kz**ilap2 
              r2 = k2h*kh**ilap2 + k2z*kz**ilap2
#endif

              if (khn.lt.0.1) then
                r1 = r1 + ek
                r2 = r2 + ek
              endif

              dpz = 1. + r1
              dmz = 1. - r1
              tzx = ( zxok(ikx,iky,ikz)*dmz + d2*termzx )/dpz
              tzy = ( zyok(ikx,iky,ikz)*dmz + d2*termzy )/dpz
              tzz = ( zzok(ikx,iky,ikz)*dmz + d2*termzz )/dpz

              dpz = 1. + r2
              dmz = 1. - r2
              ttt = ( ttok(ikx,iky,ikz)*dmz + D2*termtt )/dpz

! The Robert filter (see Asselin,R.A.,1972,Mon.Wea.Rev.,100,p.487-490)
              avzx = zxnk(ikx,iky,ikz) + robert * (tzx - 2.*zxnk(ikx,iky,ikz) + zxok(ikx,iky,ikz))
              avzy = zynk(ikx,iky,ikz) + robert * (tzy - 2.*zynk(ikx,iky,ikz) + zyok(ikx,iky,ikz))
              avzz = zznk(ikx,iky,ikz) + robert * (tzz - 2.*zznk(ikx,iky,ikz) + zzok(ikx,iky,ikz))
              avtt = ttnk(ikx,iky,ikz) + robert * (ttt - 2.*ttnk(ikx,iky,ikz) + ttok(ikx,iky,ikz))
              zxok(ikx,iky,ikz) = avzx*L(ikx,iky,ikz)
              zxnk(ikx,iky,ikz) =  tzx*L(ikx,iky,ikz)
              zyok(ikx,iky,ikz) = avzy*L(ikx,iky,ikz)
              zynk(ikx,iky,ikz) =  tzy*L(ikx,iky,ikz)
              zzok(ikx,iky,ikz) = avzz*L(ikx,iky,ikz)
              zznk(ikx,iky,ikz) =  tzz*L(ikx,iky,ikz)
              ttok(ikx,iky,ikz) = avtt*L(ikx,iky,ikz)
              ttnk(ikx,iky,ikz) =  ttt*L(ikx,iky,ikz)

#endif /* leapfrog */ 

#if (SCHEME == 1)  
! Adams-Bashforth timestepping
#if (SPHTRUNC == 1)
              r1 = v2h*wk2**ilap
              r2 = k2h*wk2**ilap
#else
              r1 = v2h*kh**ilap2 + v2z*kz**ilap2 
              r2 = k2h*kh**ilap2 + k2z*kz**ilap2
#endif

              if (khn.lt.0.1) then
                 r1 = r1 + ek
                 r2 = r2 + ek
              endif

              dpz = 1. + r1/2.
              dmz = 1. - r1/2.
              tzx = ( zxnk(ikx,iky,ikz)*dmz + D12*( 23.*termzx - 16.*zxok(ikx,iky,ikz) & 
                  + 5.*zx1k(ikx,iky,ikz) ) )/dpz
              tzy = ( zynk(ikx,iky,ikz)*dmz + D12*( 23.*termzy - 16.*zyok(ikx,iky,ikz) & 
                  + 5.*zy1k(ikx,iky,ikz) ) )/dpz
              tzz = ( zznk(ikx,iky,ikz)*dmz + D12*( 23.*termzz - 16.*zzok(ikx,iky,ikz) & 
                  + 5.*zz1k(ikx,iky,ikz) ) )/dpz

              dpz = 1. + r2/2.
              dmz = 1. - r2/2.
              ttt = ( ttnk(ikx,iky,ikz)*dmz + D12*( 23.*termtt - 16.*ttok(ikx,iky,ikz) & 
                  + 5.*tt1k(ikx,iky,ikz) ) )/dpz

              zx1k(ikx,iky,ikz) = zxok(ikx,iky,ikz)*L(ikx,iky,ikz)
              zy1k(ikx,iky,ikz) = zyok(ikx,iky,ikz)*L(ikx,iky,ikz)
              zz1k(ikx,iky,ikz) = zzok(ikx,iky,ikz)*L(ikx,iky,ikz)
              tt1k(ikx,iky,ikz) = ttok(ikx,iky,ikz)*L(ikx,iky,ikz)

              zxok(ikx,iky,ikz) = termzx*L(ikx,iky,ikz)
              zyok(ikx,iky,ikz) = termzy*L(ikx,iky,ikz)
              zzok(ikx,iky,ikz) = termzz*L(ikx,iky,ikz)
              ttok(ikx,iky,ikz) = termtt*L(ikx,iky,ikz)

              zxnk(ikx,iky,ikz) =    tzx*L(ikx,iky,ikz)
              zynk(ikx,iky,ikz) =    tzy*L(ikx,iky,ikz)
              zznk(ikx,iky,ikz) =    tzz*L(ikx,iky,ikz)
              ttnk(ikx,iky,ikz) =    ttt*L(ikx,iky,ikz)

#endif /* adams bashforth */

           enddo
        enddo
     enddo
     
! Write to restart file
     if ( mod(nt,NDUMP).eq. 0 ) then
        ntdump=nt/ndump
#if (NETCDF == 1)
        call ncdumprst(zxnk,zynk,zznk,ttnk,uk,ntdump,time)
#else
      ! call bindumprst(zxnk,zynk,zznk,ttnk,uk,ntdump,time)  
      
#endif
     endif ! ndump

#if (KeepLM == 1)    
     if ( mod(nt,nFreqout).eq. 0 ) then
!     if ( nt .eq. 3 ) then
        call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call vartimeseries(zznk,ttnk,geok,gw1k,L,time,kxmode,kymode,kzmode,khmode,imyp,ikxmo,ikymo,ikzmo,41)
     endif
#endif


     ! Compute diagnostics
     if ( mod(nt,nspwvevry).eq. 0) then
        nspdump = nspdump + 1
        call wvspc3D(gw1k,gw2k,time,nspdump,kxa,kya,kza,ikxsp,ikysp,ikzsp,dumpsp,totnzsp)
     endif
     if ( mod(nt,nout).eq. 0 ) then
        call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call spec(zxnk,zynk,zznk,ttnk,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,2,52,kh_g)
        call spec(zxnk,zynk,zznk,ttnk,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,3,51,kh_g)
!        call spec2d(zxnk,zynk,zznk,ttnk,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,57)
        call out(zxnk,zynk,zznk,ttnk,uk,vk,wk,geok,gw1k,gw2k, &
             nt,delt,ts,nstop,io,L,kxa,kya,kza,aj,bj,cor,visch,viscz,ilap,L1,L2,L3,kh_g) 
        !print*, 'check dAudt1(200,200,20)',dAudt1(200,200,20)
        
        if (linear.ne.1) then
           ! Compute basic transfer   
           call convol(zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L,uk,vk,wk,ur,vr,wr, &
                rhzx,rhzy,rhzz,rhtt,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call transf(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,      & 
                uk,vk,wk,L,time,kxa,kya,kza,aj,bj,cor,L1,L2,L3,2,48)
           call transf(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,      & 
                uk,vk,wk,L,time,kxa,kya,kza,aj,bj,cor,L1,L2,L3,3,49)
        endif ! linear



     endif ! nout 

#if (BTbalance == 1)
     if ( (mod(nt,nBTout).eq.0).or.(nt.eq.nt0) ) then

        ! Baer-Tribbia
        zxbalk = zxnk
        zybalk = zynk
        zzbalk = zznk  
        ttbalk = ttnk

        call killwaves(zxbalk,zybalk,zzbAlk,ttbalk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        ! ---------------- finding NGG
        call constr (zxbalk ,zybalk,zzbalk,ttbalk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxbalr,zybalr,zzbalr,ttbalr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(zxbalk ,zybalk,zzbalk,ttbalk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call wtoab(nzxk,nzyk,nzzk,nttk,MaGGk,NaGG1k,NaGG2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

        ! ---------------- finding NcApAp
        do ikz = 1,iktzp
           ikza = mype*iktzp+ikz
           kz = kza(ikza)
           do iky = 1,IKTY
              ky = kya(iky)
              do ikx=1,IKTX
                 if ( L(ikx,iky,ikz).eq.1 ) then
                    kx  = kxa(ikx)
                    kh  = sqrt( kx*kx + ky*ky )
                    khn = kh * L1/twopi
                    wk2 = kx*kx + ky*ky + kz*kz
                    k  = sqrt(wk2)
                    freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k          
#if (SPHTRUNC == 1)
                    r1 = v2h/delt*wk2**ilap
                    r2 = k2h/delt*wk2**ilap
#else
                    r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
                    r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif
                    
                    gw1k(ikx,iky,ikz)=NaGG1k(ikx,iky,ikz)/(-zi*freq+r1)  != A' (first order balance == Machenhauer
                    gw2k(ikx,iky,ikz)=NaGG2k(ikx,iky,ikz)/(+zi*freq+r1)
                 else
                    gw1k(ikx,iky,ikz)=cmplx(0.,0.)
                    gw2k(ikx,iky,ikz)=cmplx(0.,0.)
                 endif ! Ly
              enddo
           enddo
        enddo
        
        geok=cmplx(0.,0.)
        
        call atowb(geok,gw1k,gw2k,zxbalk ,zybalk,zzbalk,ttbalk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxbalk ,zybalk,zzbalk,ttbalk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxbalr,zybalr,zzbalr,ttbalr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NcApAp1k,NcApAp2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
 
        
        ! ------------ finding NGAp=NaGG+NbGAp+NcApAp
        call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        
        do ikz = 1,iktzp
           ikza = mype*iktzp+ikz
           kz = kza(ikza)
           do iky = 1,IKTY
              ky = kya(iky)
              do ikx=1,IKTX
                 if ( L(ikx,iky,ikz).eq.1 ) then
                    kx  = kxa(ikx)
                    kh  = sqrt( kx*kx + ky*ky )
                    khn = kh * L1/twopi
                    wk2 = kx*kx + ky*ky + kz*kz
                    k  = sqrt(wk2)
                    freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k
#if (SPHTRUNC == 1)
                    r1 = v2h/delt*wk2**ilap
                    r2 = k2h/delt*wk2**ilap
#else
                    r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
                    r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif
                    
                    gw1k(ikx,iky,ikz)=NaGG1k(ikx,iky,ikz)/(-zi*freq+r1)
                    gw2k(ikx,iky,ikz)=NaGG2k(ikx,iky,ikz)/(+zi*freq+r1)
                 else
                    gw1k(ikx,iky,ikz)=cmplx(0.,0.)
                    gw2k(ikx,iky,ikz)=cmplx(0.,0.)
                 endif ! Ly
              enddo
           enddo
        enddo
        
        call atowb(geok,gw1k,gw2k,zxbalk ,zybalk,zzbalk,ttbalk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxbalk ,zybalk,zzbalk,ttbalk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxbalr,zybalr,zzbalr,ttbalr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NGAp1k,NGAp2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
               
        
        
        ! ---------------- Calculating dApdT = 2 GdG/dt = Na(dG/dt)(G+dG/dt) -NaGG -NdG/dtdG/t
        
        ! ---------------- STEP#1 calculating NaNdG/dtdG/t =NaMM (M=dG/dt)
        
        geok=MaGGk
        gw1k=cmplx(0.,0.)
        gw2k=cmplx(0.,0.)
        
        call atowb(geok,gw1k,gw2k,zxbalk ,zybalk,zzbalk,ttbalk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxbalk ,zybalk,zzbalk,ttbalk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxbalr,zybalr,zzbalr,ttbalr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NaMM1k,NaMM2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        
        ! ---------------- STEP#2 calculating NaGpM=N(G+dG)(G+dG)
        
        call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        geok=MaGGk+geok
        gw1k=cmplx(0.,0.)
        gw2k=cmplx(0.,0.) 
        call atowb(geok,gw1k,gw2k,zxbalk ,zybalk,zzbalk,ttbalk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxbalk ,zybalk,zzbalk,ttbalk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxbalr,zybalr,zzbalr,ttbalr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NaGpM1k,NaGpM2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        
       
        ! ---------------- FINAL STEP calculating A_b
        call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        do ikz = 1,iktzp
           ikza = mype*iktzp+ikz
           kz = kza(ikza)
           do iky = 1,IKTY
              ky = kya(iky)
              do ikx=1,IKTX
                 if ( L(ikx,iky,ikz).eq.1 ) then
                    kx  = kxa(ikx)
                    kh  = sqrt( kx*kx + ky*ky )
                    khn = kh * L1/twopi
                    wk2 = kx*kx + ky*ky + kz*kz
                    k  = sqrt(wk2)
                    freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k
                    
#if (SPHTRUNC == 1)
                    r1 = v2h/delt*wk2**ilap
                    r2 = k2h/delt*wk2**ilap
#else
                    r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
                    r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif
                    
                    ! ----------------- This is just Machenhauer ----------------
                    Ord1gw1k=NaGG1k(ikx,iky,ikz)/(-zi*freq+r1)
                    Ord1gw2k=NaGG2k(ikx,iky,ikz)/(+zi*freq+r1)
                    
                    ! ----------------- deriving dA'/dt and NbGA' and adding them together
                    NbGAp1k=NGAp1k(ikx,iky,ikz)-NaGG1k(ikx,iky,ikz)-NcApAp1k(ikx,iky,ikz)
                    NbGAp2k=NGAp2k(ikx,iky,ikz)-NaGG2k(ikx,iky,ikz)-NcApAp2k(ikx,iky,ikz)
                    
                    dApdT1k=NaGpM1k(ikx,iky,ikz)-NaMM1k(ikx,iky,ikz)-NaGG1k(ikx,iky,ikz)
                    dApdT2k=NaGpM2k(ikx,iky,ikz)-NaMM2k(ikx,iky,ikz)-NaGG2k(ikx,iky,ikz)
                    
                    Ord2gw1k=NbGAp1k/(-zi*freq+r1)-(dApdT1k)/(-zi*freq+r1)**2
                    Ord2gw2k=NbGAp2k/(+zi*freq+r1)-(dApdT2k)/(+zi*freq+r1)**2
                    
                    gw1k(ikx,iky,ikz)=Ord1gw1k+Ord2gw1k
                    gw2k(ikx,iky,ikz)=Ord1gw2k+Ord2gw2k
                    Ab1(ikx,iky,ikz)=gw1k(ikx,iky,ikz)
                    Ab2(ikx,iky,ikz)=gw2k(ikx,iky,ikz)
                 else
                    gw1k(ikx,iky,ikz)=cmplx(0.,0.)
                    gw2k(ikx,iky,ikz)=cmplx(0.,0.)
                    Ab1(ikx,iky,ikz)=gw1k(ikx,iky,ikz)
                    Ab2(ikx,iky,ikz)=gw2k(ikx,iky,ikz)
                 endif ! Ly
              enddo
           enddo
        enddo    
        call atowb(geok,gw1k,gw2k,zxbalk,zybalk,zzbalk,ttbalk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

        if (BalUnbalTrnsFlag.eq.1) then
           call convol(zxbalk,zybalk,zzbalk,ttbalk,nzxk,nzyk,nzzk,nttk,L,uk,vk,wk,ur,vr,wr, &
                rhzx,rhzy,rhzz,rhtt,zxbalr,zybalr,zzbalr,ttbalr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call transf(zxbalk,zybalk,zzbalk,ttbalk,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,      & 
                uk,vk,wk,L,time,kxa,kya,kza,aj,bj,cor,L1,L2,L3,2,61)
           call transf(zxbalk,zybalk,zzbalk,ttbalk,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,      & 
                uk,vk,wk,L,time,kxa,kya,kza,aj,bj,cor,L1,L2,L3,3,62)
        endif
        if (BalUnbalSpcFlag.eq.1) then
           call spec(zxbalk,zybalk,zzbalk,ttbalk,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,2,71,kh_g)
           call spec(zxbalk,zybalk,zzbalk,ttbalk,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,3,72,kh_g)
        endif
        if (BalUnbalEngFlag.eq.1) then
           call engeps(zxbalk,zybalk,zzbalk,ttbalk,uk,vk,wk,geok,gw1k,gw2k,nt,delt,ts,  &
                     L,kxa,kya,kza,aj,bj,cor,visch,viscz,ilap,L1,L2,L3,101,91,kh_g)
        endif

        zxunbalk = zxnk-zxbalk
        zyunbalk = zynk-zybalk
        zzunbalk = zznk-zzbalk
        ttunbalk = ttnk-ttbalk
        call wtoab(zxunbalk,zyunbalk,zzunbalk,ttunbalk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        
        if (BalUnbalTrnsFlag.eq.1) then
           call convol(zxunbalk,zyunbalk,zzunbalk,ttunbalk,nzxk,nzyk,nzzk,nttk,L,uk,vk,wk,ur,vr,wr, &
                rhzx,rhzy,rhzz,rhtt,zxunbalr,zyunbalr,zzunbalr,ttunbalr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call transf(zxunbalk,zyunbalk,zzunbalk,ttunbalk,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,      & 
                uk,vk,wk,L,time,kxa,kya,kza,aj,bj,cor,L1,L2,L3,2,63)
           call transf(zxunbalk,zyunbalk,zzunbalk,ttunbalk,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,      & 
                uk,vk,wk,L,time,kxa,kya,kza,aj,bj,cor,L1,L2,L3,3,64)
        endif
        if (BalUnbalSpcFlag.eq.1) then
           call spec(zxunbalk,zyunbalk,zzunbalk,ttunbalk,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,2,73,kh_g)
           call spec(zxunbalk,zyunbalk,zzunbalk,ttunbalk,uk,vk,wk,geok,gw1k,gw2k,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,3,74,kh_g)
        endif
        if (BalUnbalEngFlag.eq.1) then
           call engeps(zxunbalk,zyunbalk,zzunbalk,ttunbalk,uk,vk,wk,geok,gw1k,gw2k,nt,delt,ts,  &
                     L,kxa,kya,kza,aj,bj,cor,visch,viscz,ilap,L1,L2,L3,102,92,kh_g)
        endif

#if (UnbalEvol == 1) /* start of executing dA_b/dt */
        ! ------------------ Finding the evolution of A_b = dA_b/dt (not dA_b/dT) ---------------------

        ! ------------------ Derive first approximation dalf'/dt --------------------------------------
        ! ---------------- finding Mtot = MaGG + MbGA + McAA ------------------------------------------
        zxhelpk = zxnk
        zyhelpk = zynk
        zzhelpk = zznk  
        tthelpk = ttnk 
        call constr (zxhelpk ,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,Mtotk,Ntot1k,Ntot2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        ! ----------------- Calculate NaMtMt ---------------------------------------------------------  
        gw1k = cmplx(0.,0.)
        gw2k = cmplx(0.,0.)
        geok = Mtotk 
        call atowb(geok,gw1k,gw2k,zxhelpk ,zyhelpk,zzhelpk,tthelpk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxhelpk ,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NaMtMt1k,NaMtMt2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        
        !------------------ NaGpMt ------------------------------------------------------------------
        call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        geok=Mtotk+geok
        gw1k=cmplx(0.,0.)
        gw2k=cmplx(0.,0.)     
        call atowb(geok,gw1k,gw2k,zxhelpk ,zyhelpk,zzhelpk,tthelpk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxhelpk ,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NaGpMt1k,NaGpMt2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

         do ikz = 1,iktzp
           ikza = mype*iktzp+ikz
           kz = kza(ikza)
           do iky = 1,IKTY
              ky = kya(iky)
              do ikx=1,IKTX
                 if ( L(ikx,iky,ikz).eq.1 ) then
                    kx  = kxa(ikx)
                    kh  = sqrt( kx*kx + ky*ky )
                    khn = kh * L1/twopi
                    wk2 = kx*kx + ky*ky + kz*kz
                    k  = sqrt(wk2)
                    freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k   
#if (SPHTRUNC == 1)
                    r1 = v2h/delt*wk2**ilap
                    r2 = k2h/delt*wk2**ilap
#else
                    r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
                    r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif         
                    dalfapdt1(ikx,iky,ikz)=(NaGpMt1k(ikx,iky,ikz)-NaMtMt1k(ikx,iky,ikz)-NaGG1k(ikx,iky,ikz))/(-zi*freq+r1)**2
                    dalfapdt2(ikx,iky,ikz)=(NaGpMt2k(ikx,iky,ikz)-NaMtMt2k(ikx,iky,ikz)-NaGG2k(ikx,iky,ikz))/(+zi*freq+r1)**2 
                 else
                    gw1k(ikx,iky,ikz)=cmplx(0.,0.)
                    gw2k(ikx,iky,ikz)=cmplx(0.,0.)
                 endif ! Ly
              enddo
           enddo
        enddo    

        ! ---------- FINDING dalf"/dt = first term +second term + third term +fourth term-------------
        ! ---------------- STEP#1a finding first term = NbMtAp ---------------------------------------
        ! We already have NaMtMt and NcApAp, we just need to derive NMtAp = NaMtMt + NbMtAp + NcApAp 
        do ikz = 1,iktzp
           ikza = mype*iktzp+ikz
           kz = kza(ikza)
           do iky = 1,IKTY
              ky = kya(iky)
              do ikx=1,IKTX
                 if ( L(ikx,iky,ikz).eq.1 ) then
                    kx  = kxa(ikx)
                    kh  = sqrt( kx*kx + ky*ky )
                    khn = kh * L1/twopi
                    wk2 = kx*kx + ky*ky + kz*kz
                    k  = sqrt(wk2)
                    freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k          
#if (SPHTRUNC == 1)
                    r1 = v2h/delt*wk2**ilap
                    r2 = k2h/delt*wk2**ilap
#else
                    r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
                    r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif                    
                    gw1k(ikx,iky,ikz)=NaGG1k(ikx,iky,ikz)/(-zi*freq+r1)
                    gw2k(ikx,iky,ikz)=NaGG2k(ikx,iky,ikz)/(+zi*freq+r1)
                 else
                    gw1k(ikx,iky,ikz)=cmplx(0.,0.)
                    gw2k(ikx,iky,ikz)=cmplx(0.,0.)
                 endif ! Ly
              enddo
           enddo
        enddo
        geok= Mtotk 
        call atowb(geok,gw1k,gw2k,zxhelpk ,zyhelpk,zzhelpk,tthelpk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxhelpk ,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NMtAp1k,NMtAp2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

        ! ---------------- STEP#2 finding second term = NbGdAp ---------------------------------------
        ! ---------------- step#2a finding second term = NGdAp = NaGG + NbGdAp + NcdApdAp ------------
        call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        gw1k=dalfapdt1
        gw2k=dalfapdt2        
        call atowb(geok,gw1k,gw2k,zxhelpk ,zyhelpk,zzhelpk,tthelpk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxhelpk ,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NGdAp1k,NGdAp2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        ! ---------------- step#2b finding second term = NcdApdAp ------------------------------------
        geok=cmplx(0.,0.)
        gw1k=dalfapdt1
        gw2k=dalfapdt2        
        call atowb(geok,gw1k,gw2k,zxhelpk ,zyhelpk,zzhelpk,tthelpk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxhelpk ,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NcdApdAp1k,NcdApdAp2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        ! ---------------- STEP#4 finding fourth term = NaGdM ---------------------------------------
        ! ---------------- step#4a derive dM/dt -----------------------------------------------------
        nMtotNew=nt
        dMdtk = (Mtotk- Mtotoldk)/ delt 
        !if (mype.eq.0) print*, 'nMtotNew & nMtotOld',nMtotNew, nMtotOld
        ! ---------------- step#4b calculating NadMdM
        geok=dMdtk
        gw1k=cmplx(0.,0.)
        gw2k=cmplx(0.,0.)
        call atowb(geok,gw1k,gw2k,zxhelpk ,zyhelpk,zzhelpk,tthelpk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxhelpk ,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NadMdM1k,NadMdM2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        ! ---------------- step#4c calculating NaGpdM 
        call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        geok=dMdtk+geok
        gw1k=cmplx(0.,0.)
        gw2k=cmplx(0.,0.)
        call atowb(geok,gw1k,gw2k,zxhelpk ,zyhelpk,zzhelpk,tthelpk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        call constr (zxhelpk ,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
             uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
        call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NaGpdM1k,NaGpdM2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        

        ! ---------------- FINAL STEP of calculating dalf"/dt ------------------
        do ikz = 1,iktzp
           ikza = mype*iktzp+ikz
           kz = kza(ikza)
           do iky = 1,IKTY
              ky = kya(iky)
              do ikx=1,IKTX
                 if ( L(ikx,iky,ikz).eq.1 ) then
                    kx  = kxa(ikx)
                    kh  = sqrt( kx*kx + ky*ky )
                    khn = kh * L1/twopi
                    wk2 = kx*kx + ky*ky + kz*kz
                    k  = sqrt(wk2)
                    freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k
                    
#if (SPHTRUNC == 1)
                    r1 = v2h/delt*wk2**ilap
                    r2 = k2h/delt*wk2**ilap
#else
                    r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
                    r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif                   
                    NbMtAp1k=NMtAp1k(ikx,iky,ikz)-NaMtMt1k(ikx,iky,ikz)-NcApAp1k(ikx,iky,ikz)
                    NbMtAp2k=NMtAp2k(ikx,iky,ikz)-NaMtMt2k(ikx,iky,ikz)-NcApAp2k(ikx,iky,ikz)
                    FirstTrm1=NbMtAp1k/(-zi*freq+r1)
                    FirstTrm2=NbMtAp2k/(+zi*freq+r1)
                    NbGdAp1k=NGdAp1k(ikx,iky,ikz)-NcdApdAp1k(ikx,iky,ikz)-NaGG1k(ikx,iky,ikz)
                    NbGdAp2k=NGdAp2k(ikx,iky,ikz)-NcdApdAp2k(ikx,iky,ikz)-NaGG2k(ikx,iky,ikz)
                    SecndTrm1=NbGdAp1k/(-zi*freq+r1)
                    SecndTrm2=NbGdAp2k/(+zi*freq+r1)
                    ThirdTrm1=-2*NaMtMt1k(ikx,iky,ikz)/(-zi*freq+r1)**2
                    ThirdTrm2=-2*NaMtMt2k(ikx,iky,ikz)/(+zi*freq+r1)**2
                    NaGdM1k=NaGpdM1k(ikx,iky,ikz)-NadMdM1k(ikx,iky,ikz)-NaGG1k(ikx,iky,ikz)
                    NaGdM2k=NaGpdM2k(ikx,iky,ikz)-NadMdM2k(ikx,iky,ikz)-NaGG2k(ikx,iky,ikz)
                    FourthTrm1=-NaGdM1k/(-zi*freq+r1)**2
                    FourthTrm2=-NaGdM2k/(+zi*freq+r1)**2
                    dAbdt1(ikx,iky,ikz)=dalfapdt1(ikx,iky,ikz)+FirstTrm1+SecndTrm1+ThirdTrm1+FourthTrm1
                    dAbdt2(ikx,iky,ikz)=dalfapdt2(ikx,iky,ikz)+FirstTrm2+SecndTrm2+ThirdTrm2+FourthTrm2
                 else
                    dAbdt1(ikx,iky,ikz)=cmplx(0.,0.)
                    dAbdt2(ikx,iky,ikz)=cmplx(0.,0.)
                 endif ! Ly
              enddo
           enddo
        enddo    
        
        call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
        
        dAudt1 = cmplx(0.,0.)
        dAudt2 = cmplx(0.,0.)

        do ikz = 1,iktzp
           ikza = mype*iktzp+ikz
           kz = kza(ikza)
           do iky = 1,IKTY
              ky = kya(iky)
              do ikx=1,IKTX
                 if ( L(ikx,iky,ikz).eq.1 ) then
                    kx  = kxa(ikx)
                    kh  = sqrt( kx*kx + ky*ky )
                    khn = kh * L1/twopi
                    wk2 = kx*kx + ky*ky + kz*kz
                    k  = sqrt(wk2)
                    freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k
                    
#if (SPHTRUNC == 1)
                    r1 = v2h/delt*wk2**ilap
                    r2 = k2h/delt*wk2**ilap
#else
                    r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
                    r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif                   
                   dAudt1(ikx,iky,ikz)=Ntot1k(ikx,iky,ikz)+(+zi*freq-r1)*gw1k(ikx,iky,ikz)-dAbdt1(ikx,iky,ikz)
                   dAudt2(ikx,iky,ikz)=Ntot2k(ikx,iky,ikz)+(-zi*freq-r1)*gw2k(ikx,iky,ikz)-dAbdt2(ikx,iky,ikz)
                else
                   dAudt1(ikx,iky,ikz)=cmplx(0.,0.)
                   dAudt2(ikx,iky,ikz)=cmplx(0.,0.)
                 endif ! Ly
              enddo
           enddo
        enddo   
        
        call transf2unbal(zxunbalk,zyunbalk,zzunbalk,ttunbalk,geok,gw1k,gw2k,dAudt1,dAudt2, &
                  uk,vk,wk,L,time,kxa,kya,kza,aj,bj,cor,delt,v2h,k2h,v2z,k2z,ilap,ilap2,L1,L2,L3,2,211)
        call transf2unbal(zxunbalk,zyunbalk,zzunbalk,ttunbalk,geok,gw1k,gw2k,dAudt1,dAudt2, &
                  uk,vk,wk,L,time,kxa,kya,kza,aj,bj,cor,delt,v2h,k2h,v2z,k2z,ilap,ilap2,L1,L2,L3,3,212)

#endif  /* end of executing dA_b/dt */

 endif ! end of Baer-Tribbia 

 if (kpTimeScaleSPC .eq. 1) then
    call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
    call specTimeScales(Mtotk,geok,gw1k,gw2k,dAbdt1,Ab1,dAbdt2,Ab2,dAudt1,dAudt2, &
         L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,2,89)
    call specTimeScales(Mtotk,geok,gw1k,gw2k,dAbdt1,Ab1,dAbdt2,Ab2,dAudt1,dAudt2,  &
         L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,3,99)
    
 endif

#endif /* end of executing Baer-Tribbia model */


#if (UnbalEvol == 1) /* start of executing dA_b/dt */
    if ((mod(nt+1,nBTout).eq.0)) then
       !-------------creating Mtotoldk
       nMtotOld=nt
       call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L, &
          uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)  
       call wtoab(nzxk,nzyk,nzzk,nttk,Mtotoldk,Ntot1k,Ntot2k,rhzx,rhzy,rhzz,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
    endif
#endif  /* end of executing dA_b/dt */


#if (HCVC4movie == 1)
     if ( mod(nt,nRSevery).eq. 0 ) then

        if (kpzz  .eq.1)  call dumprealHCVC(zznk,zznr,nRSzz,nRSout,time,111,121,nskipH,nskipV,nHmax,nVmax,justHor)
        if (kpdiv .eq.1)  then
           call velo (zxnk,zynk,zznk,uk,vk,wk,L,kxa,kya,kza)
           do iky=1,ikty
              ky = kya(iky)
              do ikx=1,iktx
                 kx = kxa(ikx)
                 do ikz=1,iktzp
                    wntdFldk(ikx,iky,ikz) = zi*(kx*uk(ikx,iky,ikz)+ky*vk(ikx,iky,ikz))
                 enddo
              enddo
           enddo
           call dumprealHCVC(wntdFldk,wntdFldr,nRSdiv,nRSout,time,110,120,nskipH,nskipV,nHmax,nVmax,justHor)
        endif
        
        
           
#if (UnbalEvol == 1)
        if (kpdAu .eq.1) then
           wntdFldk = dAudt1
           call dumprealHCVC(wntdFldk,wntdFldr,nRSdau,nRSout,time,151,161,nskipH,nskipV,nHmax,nVmax,justHor)
           wntdFldk = dAudt2
           call dumprealHCVC(wntdFldk,wntdFldr,nRSdau,nRSout,time,152,162,nskipH,nskipV,nHmax,nVmax,justHor)
        endif
#endif
           
        if (kpGAb+kpGAu+kpGG+kpAuAu+kpAbAb+kpAbAu+kpzzg+kpzzab.gt.0) then
           
           zxhelpk = zxnk
           zyhelpk = zynk
           zzhelpk = zznk  
           tthelpk = ttnk
           call constr (zxhelpk,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
                uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call wtoab(nzxk,nzyk,nzzk,nttk,NtotGk,rhzx,rhzy,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           zxhelpk = zxnk
           zyhelpk = zynk
           zzhelpk = zznk  
           tthelpk = ttnk
           call killwaves(zxhelpk,zyhelpk,zzhelpk,tthelpk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           if (kpzzg .eq.1)  then
              wntdFldk=zzhelpk
              call dumprealHCVC(wntdFldk,wntdFldr,nRSzzg,nRSout,time,112,122,nskipH,nskipV,nHmax,nVmax,justHor)
           endif
          
           call constr (zxhelpk,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
                uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call wtoab(nzxk,nzyk,nzzk,nttk,NGGk,rhzx,rhzy,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           if (kpgg .eq.1)  then
              wntdFldk=NGGk
              call dumprealHCVC(wntdFldk,wntdFldr,nRSgg,nRSout,time,116,126,nskipH,nskipV,nHmax,nVmax,justHor)
           endif
           zxhelpk = zxnk
           zyhelpk = zynk
           zzhelpk = zznk
           tthelpk = ttnk
           call killgeo(zxhelpk,zyhelpk,zzhelpk,tthelpk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           call constr (zxhelpk,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
                uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call wtoab(nzxk,nzyk,nzzk,nttk,NAAk,rhzx,rhzy,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           NGAk = NtotGk - NGGk - NAAk        
           call constr (zxbalk,zybalk,zzbalk,ttbalk,nzxk,nzyk,nzzk,nttk,L,  &
                uk,vk,wk,ur,vr,wr,zxbalr,zybalr,zzbalr,ttbalr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call wtoab(nzxk,nzyk,nzzk,nttk,NtotBGk,rhzx,rhzy,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           zxhelpk = zxbalk 
           zyhelpk = zybalk
           zzhelpk = zzbalk  
           tthelpk = ttbalk 
           call killgeo(zxhelpk,zyhelpk,zzhelpk,tthelpk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           if (kpzzab .eq.1)  then
              wntdFldk=zzhelpk
              call dumprealHCVC(wntdFldk,wntdFldr,nRSzzab,nRSout,time,113,123,nskipH,nskipV,nHmax,nVmax,justHor)
           endif
           call constr (zxhelpk,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
                uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call wtoab(nzxk,nzyk,nzzk,nttk,NAbAbk,rhzx,rhzy,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           if (kpAbAb .eq.1)  then
              wntdFldk = NAbAbk
              call dumprealHCVC(wntdFldk,wntdFldr,nRSabab,nRSout,time,117,127,nskipH,nskipV,nHmax,nVmax,justHor)
           endif
           NGAbk = NtotBGk - NGGk - NAbAbk
           if (kpGAb .eq.1)  then
              wntdFldk = NGAbk
              call dumprealHCVC(wntdFldk,wntdFldr,nRSGAb,nRSout,time,114,124,nskipH,nskipV,nHmax,nVmax,justHor)
           endif       
           call constr (zxunbalk,zyunbalk,zzunbalk,ttunbalk,nzxk,nzyk,nzzk,nttk,L,  &
                uk,vk,wk,ur,vr,wr,zxunbalr,zyunbalr,zzunbalr,ttunbalr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call wtoab(nzxk,nzyk,nzzk,nttk,NtotGk,rhzx,rhzy,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           zxhelpk = zxunbalk 
           zyhelpk = zyunbalk
           zzhelpk = zzunbalk   
           tthelpk = ttunbalk 
           call killgeo(zxhelpk,zyhelpk,zzhelpk,tthelpk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           call constr (zxhelpk,zyhelpk,zzhelpk,tthelpk,nzxk,nzyk,nzzk,nttk,L,  &
                uk,vk,wk,ur,vr,wr,zxhelpr,zyhelpr,zzhelpr,tthelpr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
           call wtoab(nzxk,nzyk,nzzk,nttk,NAuAuk,rhzx,rhzy,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
           if (kpAuAu .eq.1)  then
              wntdFldk = NAuAuk
              call dumprealHCVC(wntdFldk,wntdFldr,nRSAuAu,nRSout,time,119,129,nskipH,nskipV,nHmax,nVmax,justHor)
           endif
           NGAuk = NtotGk - NAuAuk
           if (kpGAu .eq.1)  then
              wntdFldk = NGAuk
              call dumprealHCVC(wntdFldk,wntdFldr,nRSGAu,nRSout,time,115,125,nskipH,nskipV,nHmax,nVmax,justHor)
           endif
           NAbAuk = NAAk - NAbAbk-NAuAuk
           if (kpAbAu .eq.1)  then
              wntdFldk = NAbAuk
              call dumprealHCVC(wntdFldk,wntdFldr,nRSAbAu,nRSout,time,118,128,nskipH,nskipV,nHmax,nVmax,justHor)
           endif
           if (dumprmsNGs.eq.1)  call nltrmsgrms(NGGk,NGAuk,NGAbk,NAbAbk,NAbAuk,NAuAuk,nt,delt,ts,L,131)
           if (tran2gOn.eq.1) then 
              call transtrms2g(NGGk,NGAbk,NGAuk,NAbAbk,NAbAuk,NAuAuk,zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k, &
                   uk,vk,wk,L,time,kxa,kya,kza,aj,bj,cor,L1,L2,L3,2,141)
           endif

        endif

     endif
#endif    

     

  enddo ! nt
  
  time3 = etime(time2)
  time3 = time2(1) - time1(1)
  
#if (NETCDF == 1)
  if (mype.eq.0) then
     ierror = nf_close(idnc)
!     if (ierror.ne.0) print*,'Error closing out.ncf'
     ierror = nf_close(idnck)
!     if (ierror.ne.0) print*,'Error closing Zk.out.ncf'
     ierror = nf_close(idncm)
!     if (ierror.ne.0) print*,'Error closing Zmod.out.ncf'
  endif
  if (mype.eq.0) then
     ierror = nf_close(idnsp)
     if (ierror.ne.0) print*,'Error closing WVspc3D.out.ncf'
  endif
#endif
  
  if (mype.eq.0) then
     close (58)
     close (59)
     print*,'    '
     write(6,5000) time3, time3/60., time3/3600., time3/86400.
     print*,'    '
  endif

#if (MPI == 1)
  call mpi_finalize(ierror)
#endif

5000 format(1x,'CPU time required for main loop = ',F7.0,' s = ', F7.1,' m = ',F7.2,' h =',F7.3,' d.')
end PROGRAM MAIN


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! INITIALIZATION
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine initialize (zx,zy,zz,tt,uu,vv,ww,ge,g1,g2,zxr,zyr,zzr,ur,vr,wr,tr,delt,      & 
                 ampk,ampp,L,irest,kxa,kya,kza,ktrunc_x,ktrunc_y,ktrunc_z,        & 
                 aj,bj,cor,ts,L1,L2,L3,iwfr,kxwv,kywv,khwv,kzwv,iweng)

! Initialises stuff like wavenumbers, indices, spectra, phases, etc.

  implicit none
  include 'param.inc'

  integer :: ikx,iky,ikz,ikza
  integer :: L(iktx,ikty,iktzp),irest
  integer :: jw,ntime
  integer :: i,j,k,ja,jj,ilc,inormke,inormve
  complex, dimension(iktx,ikty,iktzp) :: zx, zy, zz, tt, uu, vv, ww, ge, g1, g2
  real,    dimension(n1d,n3d,n2dp)       :: zxr,zyr,zzr,tr, ur, vr, wr
  complex, dimension(iktx,ikty,iktzp) :: zxwv, zywv, zzwv, ttwv
  complex :: zi
  real ::  kx,ky,kz,wk,ek,ampk,ampp,aj,bj,ki,aaa
  real :: iwfr,kxwv,kywv,khwv,kzwv,iweng
  real ::  phase,twopi,kinen,poten,cor,vorten,waven
  real ::  delt,ranno,time1,ktrunc_x,ktrunc_y,ktrunc_z,noiseamp
  real ::  kxa(iktx),kya(ikty),kza(iktz),ts,kh,L1,L2,L3,khn,wkn,kzn
  real ::  r0,u0,x,y,z,r,dxtan,dytan,tanhx,tanhy
  real ::  bessj0,bessj1,x0,y0,pert,kzpert,xr,yr
  real, parameter :: mu1=3.83170597020751
  external :: ranno,proj,realit

  integer :: mype
  common/mpi/mype

  zi = cmplx(0.,1.)
  twopi = 4.*asin(1.)

! Initialize wavenumber arrays
  do  ikx = 1,iktx
     kxa(ikx) = float(ikx-1) * twopi/L1
  enddo
  do iky = 1,ikty
     jj = iky - 1
     if (iky.gt.kty)   jj = jj - 2*kty
     if (iky.eq.kty+1) jj = 0
     if (iky.gt.2*kty) jj = 0
     kya(iky) = float(jj) * twopi/L2
  enddo
  do ikz = 1,iktz
     jj = ikz - 1
     if (ikz.gt.ktz)   jj = jj - 2*ktz
     if (ikz.eq.ktz+1) jj = 0
     if (ikz.gt.2*ktz) jj = 0
     kza(ikz) = float(jj) * twopi/L3
  enddo
  

! L(ikx,iky,ikz) is unity for retained modes and zero beyond the truncation, 
! at k=0 and for half the modes on the plane kx=0.

  if (mype.eq.0) then
#if (SPHTRUNC == 1)
     print*,' '
     print*,'Using spherical truncation.'
#else
     print*,' '
     print*,'Using cylindrical truncation.'
#endif
  endif

  L  = 1
  zx = cmplx(0.,0.)
  zy = cmplx(0.,0.)
  zz = cmplx(0.,0.)
  tt = cmplx(0.,0.)
  uu = cmplx(0.,0.)
  vv = cmplx(0.,0.)
  ww = cmplx(0.,0.)
  ge = cmplx(0.,0.)
  g1 = cmplx(0.,0.)
  g2 = cmplx(0.,0.)

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx = 1,iktx
           kx = kxa(ikx)
           kh = sqrt(kx*kx + ky*ky)
           wk = kx*kx + ky*ky + kz*kz
           wk = sqrt( wk )
           
! Set L=0 where necessary:
           if (iky.eq.kty+1)                         L(ikx,iky,ikz) = 0
           if (ikza.eq.ktz+1)                        L(ikx,iky,ikz) = 0
           if (kx.lt.0)                              L(ikx,iky,ikz) = 0
           if (kx.eq.0 .and. ky.lt.0)                L(ikx,iky,ikz) = 0
           if (kx.eq.0 .and. ky.eq.0 .and. kz.lt.0)  L(ikx,iky,ikz) = 0
           if (iky.gt.2*kty)                         L(ikx,iky,ikz) = 0
           if (ikza.gt.2*ktz)                        L(ikx,iky,ikz) = 0

#if (SPHTRUNC == 1)
           if (wk.eq.0. .or. wk.gt.ifix(KTRUNC_X + 0.5) - 0.5)                L(ikx,iky,ikz) = 0
#else
           if (wk.eq.0. .or. (kh*L1/twopi).gt.ifix(float(n1)/3.+0.5)-0.5)     L(ikx,iky,ikz) = 0         
           if (abs(kz*L3/twopi).gt.ifix(float(n3)/3. + 0.5)-0.5)              L(ikx,iky,ikz) = 0          
#endif
        enddo
     enddo
  enddo


  !  Restart from output file then add the wave
#if (NETCDF == 1)
      call ncreadrst(zx,zy,zz,tt,ur,vr,irest,ts)
#else
      call binreadrst(zx,zy,zz,tt,irest,ts)
#endif
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
                     wk = sqrt(kx*kx + ky*ky + kz*kz)
                     if (abs(kz - kzwv).lt.(twopi/2/L3)) then
                        g1(ikx,iky,ikz) = kh*sqrt(iweng)
                        print*,'g1 at kx ky kz',g1(ikx,iky,ikz), kx, ky, kz
!                        g2(ikx,iky,ikz) = kh*sqrt(iweng)
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

      call atowb(ge,g1,g2,zxwv,zywv,zzwv,ttwv,uu,vv,ww,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

      zx = zx + zxwv
      zy = zy + zywv
      zz = zz + zzwv
      tt = tt + ttwv

      call proj(zx,zy,zz,L,kxa,kya,kza)      
      
end subroutine initialize


subroutine initksp(kxa,kya,kza,kxsp,kysp,kzsp,ikxsp,ikysp,ikzsp,L1,L2,L3,dumpsp,totnzsp)
  
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: ikx,iky,ikz,ikza,iikz,jj,istatus,ikz1,ikzsp1,nkzsp,stflag,totnzsp
  integer :: ikxsp(iktxsp),ikysp(iktysp-1),ikzsp(iktzsp-1),dumpsp(4),ismypein,ierr
  real    :: kxa(iktx),kya(ikty),kza(iktz),L1,L2,L3,twopi
  real    :: kxsp(iktxsp),kysp(iktysp-1),kzsp(iktzsp-1)

!! 3D spectrum of waves (and potentially other variables)
  integer :: idnsp,idkxsp,idkysp,idkzsp,idtsp
  integer, dimension(4) :: ncdimsp,ncstartsp,nccountsp
  integer, dimension(3) :: ncdimkk,ncstartkk,nccountkk
  integer :: idtimesp,idkkxx,idkkyy,idkkzz,idwvp,idwvm
  
  common/netcdfwvsp/idnsp,idkxsp,idkysp,idkzsp,idtsp,ncstartsp,nccountsp, &
       ncstartkk,nccountkk,idtimesp,idkkxx,idkkyy,idkkzz,idwvp,idwvm
  
  integer :: mype
  common/mpi/mype
  
  twopi = 4.*asin(1.)
  
  ! Initialize wavenumbers used for 3D spectra
  do  ikx = 1,iktxsp
     kxsp(ikx) = float(ikx-1) * twopi/L1
  enddo
  do iky = 1,iktysp - 1
     jj = iky - 1
     if (iky.gt.ktysp)   jj = jj - 2*ktysp + 1
     if (iky.gt.2*ktysp) print*, 'something screwed up in wv sp stuff'
     kysp(iky) = float(jj) * twopi/L2
  enddo
  do ikz = 1,iktzsp - 1
     jj = ikz - 1
     if (ikz.gt.ktzsp)   jj = jj - 2*ktzsp + 1
     if (ikz.gt.2*ktzsp) print*, 'something screwed up in wv sp stuff'
     kzsp(ikz) = float(jj) * twopi/L3
  enddo

  jj = 1
  do  ikx = 1,iktx
     if (kxsp(jj).eq.kxa(ikx)) then
        ikxsp(jj) = ikx
        jj = jj + 1
     endif
  enddo
  jj = 1
  do  iky = 1,ikty
     if (kysp(jj).eq.kya(iky)) then
        ikysp(jj) = iky
        jj = jj + 1
     endif
  enddo
  jj = 1
  do  ikz = 1,iktz
     if (kzsp(jj).eq.kza(ikz)) then
        ikzsp(jj) = ikz
        jj = jj + 1
     endif
  enddo


  ikzsp1 = -1 ! starting index of kzsp that communicates with the current process
  ikz1   = -1 ! starting index between 1:iktzp that communicates kzsp
  nkzsp  = 0  ! number of vert. modes in this thread/proc that are used for 3D spec
  stflag = 0  ! a flag to show if we started taking from this proc

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     if (ikza.le.iktz/2) then
        do iikz = 1,iktzsp/2
           if (kzsp(iikz).eq.kza(ikza)) then
              if (stflag.eq.0) then       
                 nkzsp = nkzsp + 1
                 ikz1 = ikz
                 ikzsp1 = iikz
                 stflag = 1
                 exit
              else
                 nkzsp = nkzsp + 1
                 exit
              endif
           endif
        enddo
     else
        do iikz = iktzsp/2+1,iktzsp - 1
           if (kzsp(iikz).eq.kza(ikza)) then
              if (stflag.eq.0) then       
                 nkzsp = nkzsp + 1
                 ikz1 = ikz
                 ikzsp1 = iikz
                 stflag = 1
                 exit
              else
                 nkzsp = nkzsp + 1
                 exit
              endif
           endif
        enddo
     endif
  enddo
  
  dumpsp(1) = mype   ! pass the number of proc 
  dumpsp(2) = ikzsp1 ! index in kzsp from which we start writing the info of proc mype
  dumpsp(3) = ikz1
  dumpsp(4) = nkzsp  ! total number of vertical modes in mype that is going to be dumped

  if (nkzsp.ne.0) then
     ismypein = 1
  else
     ismypein = 0
  endif
  
  call mpi_reduce(ismypein,totnzsp,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 
  
  nccountkk(1) = 1
  nccountkk(2) = 1
  nccountkk(3) = 1
  if (mype.eq.0) then
     print*,'totnzsp', totnzsp
     do  ikz = 1,iktzsp - 1 
        do  iky = 1,iktysp - 1
           do  ikx = 1,iktxsp 
              ncstartkk(1) = ikx
              ncstartkk(2) = iky
              ncstartkk(3) = ikz
              istatus = nf_put_vara_real(idnsp,idkkxx,ncstartkk,nccountkk,kxsp(ikx))
              istatus = nf_put_vara_real(idnsp,idkkyy,ncstartkk,nccountkk,kysp(iky))
              istatus = nf_put_vara_real(idnsp,idkkzz,ncstartkk,nccountkk,kzsp(ikz))
           enddo
        enddo
     enddo
  endif
  
end subroutine initksp
    

subroutine initkmodes(kxa,kya,kza,kxmode,kymode,kzmode,khmode,imyp,ikxmo,ikymo,ikzmo) 
 
  implicit none
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif
  real    :: khmode(nkhm)         ! the horizontal radious of each ring 
  real    :: kxmode(nmodeh,nkhm), kymode(nmodeh,nkhm), kzmode(nkhm,nkzm) !picked wavenumbers
  integer :: ikxmo(nmodeh,nkhm), ikymo(nmodeh,nkhm), ikzmo(nkhm,nkzm) !index of picked wavenumbers
  integer :: imyp(nkhm,nkzm)           !the number of processor for each vertical level 
  integer  :: ikx,iky,ikz,ikza
  real     :: kxa(iktx),kya(ikty),kza(iktz)
  integer  :: itet,irad,ized
  real     :: theta
  
  integer :: mype
  common/mpi/mype
  
  do irad = 1 , nkhm 
     do itet = 1 , nmodeh
        theta = real(itet-1)/nmodeh*((9.0/10.0)*3.14159265)-3.14159265/2+3.14159265/10
        kxmode(itet,irad) = NINT(khmode(irad)*cos(theta))
        kymode(itet,irad) = NINT(khmode(irad)*sin(theta))
        do ikx = 1,iktx
           if (abs(kxmode(itet,irad)-kxa(ikx)) .lt. 0.001) then
              ikxmo(itet,irad)=ikx
              exit
           endif
        enddo
        
        do iky = 1,ikty
           if (abs(kymode(itet,irad)-kya(iky)) .lt. 0.001) then
              ikymo(itet,irad)=iky
              exit
           endif
        enddo
     enddo
  enddo
  do ized = 1, nkzm
     do irad = 1, nkhm 
        do ikz = 1,iktzp
           ikza = mype*iktzp+ikz
           if (abs(kzmode(irad,ized)-kza(ikza)).lt.0.001) then
              ikzmo(irad,ized)=ikz
              imyp(irad,ized)=mype
              exit
           endif
        enddo
     enddo
  enddo
end subroutine initkmodes

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! FORCING
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine force(nzxk,nzyk,nzzk,nttk,AMPV,AMPW,gtau,L, & 
                 nfmax,alpha,kf,kxa,kya,kza,aj,bj,cor,time,L1,L2,L3)

! Use this after call convol.  Forcing is specified to be divergence-free

  implicit none
  include 'param.inc'

  integer :: nfmax,ikx,iky,ikz,ikza,ik
  integer ::  L(iktx,ikty,iktzp)
  complex, dimension(iktx,ikty,iktzp) :: nzxk,nzyk,nzzk,nttk
  complex :: gtau(nfmax,4), zi, fv,fw,ft
  real :: kxa(iktx),kya(ikty),kza(iktz)
  real :: rang,twopi,aj,bj,cor,L1,L2,L3
  real :: alpha,beta,KF,kx,ky,kz,wkh,time,cor2,bf
  real :: ampv,ampw,sqrt2,wk,wkh2,C,G,TF,theta,sk,bf2
  external :: proj,rang

  integer :: mype
  common/mpi/mype
  
  beta  = sqrt(1.-alpha**2)
  twopi = 4.*asin(1.)
  sqrt2 = sqrt(2.)
  zi    = cmplx(0.,1.)
  c     = sqrt(bj/(2.*aj))
  bf2   = aj*bj
  bf    = sqrt(bf2)
  cor2  = cor**2
  ik    = 0

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           if (L(ikx,iky,ikz).eq.1) then 
              kx = kxa(ikx)
              wkh2 = kx**2 + ky**2
              wkh  = sqrt(wkh2) 
              wk   = sqrt(kx**2 + ky**2 + kz**2)
              theta = wkh/wk
              tf    = 1./sqrt(2.)
              
              if (abs(wkh*L1/twopi-KF).le.1. .and. abs(kz*L3/twopi   ).lt.0.5) then
!             if (abs(wkh*L1/twopi-KF).le.1. .and. abs(kz*L3/twopi-1.).lt.0.5) then
                 ik=ik+1                   
                 gtau(ik,1) = alpha*gtau(ik,1) + beta*cmplx(rang(0),rang(0))
                 gtau(ik,2) = alpha*gtau(ik,2) + beta*cmplx(rang(0),rang(0))
                 gtau(ik,3) = alpha*gtau(ik,3) + beta*cmplx(rang(0),rang(0))
                 gtau(ik,4) = alpha*gtau(ik,4) + beta*cmplx(rang(0),rang(0))

                 G = (wkh*L1/twopi-(kf-1.))*(kf+1.-wkh*L1/twopi) 
!                G = 100.*(theta-(TF-0.1))*(TF+0.1-theta)
!                G = 1.
                       
                 Fv = ampv * G * gtau(ik,1) 
                 Fw = ampw * G * ( gtau(ik,2) + gtau(ik,3))
                 Ft = ampw * G * (-gtau(ik,2) + gtau(ik,3))

                 sk = sqrt(BF2*wkh2 + cor2*kz**2)
                 nzxk(ikx,iky,ikz) = nzxk(ikx,iky,ikz) - kx*kz*BF/sk*Fv & 
                                   - wk*ky/wkh/sqrt2*Fw + ZI*cor*kx*kz**2/sqrt2/sk/wkh*Ft
                 nzyk(ikx,iky,ikz) = nzyk(ikx,iky,ikz) - ky*kz*BF/sk*Fv & 
                                   + wk*kx/wkh/sqrt2*Fw + ZI*cor*ky*kz**2/sqrt2/sk/wkh*Ft
                 nzzk(ikx,iky,ikz) = nzzk(ikx,iky,ikz) + BF*wkh**2/sk*Fv - ZI*cor*wkh*kz/sqrt2/sk*Ft
                 nttk(ikx,iky,ikz) = nttk(ikx,iky,ikz) - C*sqrt2*ZI*cor*kz/sk*Fv + C*BF*wkh/sk*Ft
                 if (ik.eq.nfmax) then
                    print*,'Forcing too many modes!'
                    stop
                 endif
              endif
           endif
        enddo
     enddo
  enddo
        
  return
end subroutine force


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! DIAGNOSTICS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine spec (zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,ispec,iu,kh_g)


! Calculates spectra: 
! ispec = 1: k
! ispec = 2: kh
! ispec = 3: kz

  implicit none 
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: ikx,iky,ikz,ikza,ispec,fkh_g
  integer :: j,i,L(iktx,ikty,iktzp),iu,j0,nspz,nspz8,ierr,ilap
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,ge,g1,g2,ux,uy,uz
  real    :: kx,ky,kz,wk,wkh2,wkh,wkh2n,kzn,kz2
  real    :: aj,bj,vz,vt,time,visch,viscz,kvisc,kh_g
  real    :: vzx,vzy,vzz,kxa(iktx),kya(ikty),kza(iktz),L1,L2,L3,twopi
  real, dimension(0:kts,8) :: spz,spztot
  integer, dimension(0:kts) :: n,ntot
  external :: velo

  integer :: mype
  common/mpi/mype

  nspz  = kts+1
  nspz8 = (kts+1)*8

  if (ispec.eq.1) then
     j0 = 1
  elseif (ispec.eq.2) then
     j0 = 0
  elseif (ispec.eq.3) then
     j0 = 0
  else
     print*,"ispec error"
     stop
  endif

  twopi = 4.*asin(1.)
  N     = 0
  spz   = 0.

  call velo (zx,zy,zz,ux,uy,uz,L,kxa,kya,kza)

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
	   kvisc = visch*wkh2**ilap+viscz*kz2**ilap
           
           if (ispec.eq.1) then
              j = ifix(wk*L1/twopi+0.5)
           elseif (ispec.eq.2) then
              j = ifix(wkh*L1/twopi+0.5)
           elseif (ispec.eq.3) then
              j = ifix(abs(kz)*L3/twopi+0.5)
           endif

           if (L(ikx,iky,ikz).eq.1) then
              if (j.lt.j0 .or. j.gt.kts) print*,'SPEC: SCREW-UP.   k= ',j,'    ktx= ',ktx
              wk   = max(wk,  1.e-15)
              wkh2 = max(wkh2,1.e-15)

! Kinetic and potential energy.
              vzx      = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
              spz(j,1) = spz(j,1) + vzx/wk**2
              vzy      = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
              spz(j,1) = spz(j,1) + vzy/wk**2
              vzz      = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
              spz(j,1) = spz(j,1) + vzz/wk**2
              vt       = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )
              spz(j,2) = spz(j,2) + vt*aj/bj

! Geo, ageo decompostition.
! k \in R_k
              if (wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 vz       = real( ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)) )
                 spz(j,4) = spz(j,4) + vz/wkh2
                 spz(j,7) = spz(j,7) + kvisc*vz/wkh2
                 vz       = real( g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)) )
                 spz(j,5) = spz(j,5) + vz/wkh2
                 spz(j,8) = spz(j,8) + kvisc*vz/wkh2
                 vz       = real( g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)) )
                 spz(j,5) = spz(j,5) + vz/wkh2
                 spz(j,8) = spz(j,8) + kvisc*vz/wkh2
              endif

! Special cases: i) k_h=0, ii) k_z=0.
              vzx=real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              vzy=real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              vzz=real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )

! k \in V_k
              if(wkh2n.lt.1.e-10.and.abs(kzn).gt.1.e-10) then
                 spz(j,4) = spz(j,4) + vzz + vt*aj/bj
                 spz(j,5) = spz(j,5) + vzx + vzy
                 spz(j,7) = spz(j,7) + kvisc*(vzz + vt*aj/bj)
                 spz(j,8) = spz(j,8) + kvisc*(vzx + vzy)     
              endif

! k \in B_k
              if(abs(kzn).lt.1.e-10.and.wkh2n.gt.1.e-10) then
                 spz(j,4) = spz(j,4) + vzx + vzy
                 spz(j,5) = spz(j,5) + vzz + vt*aj/bj
                 spz(j,7) = spz(j,7) + kvisc*(vzx + vzy)     
                 spz(j,8) = spz(j,8) + kvisc*(vzz + vt*aj/bj)
              endif

! Buoyancy Flux
              if (aj.gt.0.) then 
                 spz(j,6) = spz(j,6) + aj*real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
              else
                 spz(j,6) = spz(j,6) +    real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
              endif

              n(j) = n(j) + 2 
           endif
        enddo
     enddo
  enddo




  spz(:,3) = spz(:,1) + spz(:,2)

#if ( MPI == 1) 
  call mpi_reduce(spz,spztot,nspz8,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(  n,  ntot, nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  spztot=spz
  ntot=n
#endif

  fkh_g=0
  if (mype.eq.0) then
     do j=j0,kts-1
        if (ntot(j).ne.0) then
           write(iu,5000) float(j),(spztot(j,i),i=1,8),ntot(j)
        endif
        if (ispec.eq.2) then
           if (j.lt.kts/3 .and. fkh_g.eq.0) then
              if ( (spztot(j+1,4)-spztot(j,4))*(spztot(j+2,4)-spztot(j+1,4)).lt.0.) then
!                 print*,'j =', j
!                 print*,'y(j) ', spztot(j,4)
!                 print*,'y(j+1)', spztot(j+1,4)
!                 print*,'y(j+2)', spztot(j+2,4)
                 kh_g=float(j)+0.5-(spztot(j+1,4)-spztot(j,4))/(spztot(j+2,4)-2*spztot(j+1,4)+spztot(j,4))
                 fkh_g=1
              endif
           endif
        endif
     enddo
  write(iu,*) '           '
!     write(iu,*) '           '
     call flush(iu)
  endif

  return
5000  format(1X,F4.0,4X,8(E14.6,1x),6X,I10)
5043  format(2X,' k',9X,'KE',11X,'PE',10X,' E',12X,'GE',11X,'AE',11X,'BF')
end subroutine spec

subroutine specTimeScales(dGdt,ge,g1,g2,dAb1dt,Ab1,dAb2dt,Ab2,dAu1dt,dAu2dt, &
     L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,ispec,iu)


! Calculates spectra: 
! ispec = 1: k
! ispec = 2: kh
! ispec = 3: kz

  implicit none 
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: ikx,iky,ikz,ikza,ispec
  integer :: j,i,L(iktx,ikty,iktzp),iu,j0,nspz,nspz8,ierr,ilap
  complex, dimension(iktx,ikty,iktzp) :: dGdt,ge,g1,g2,dAb1dt,Ab1,dAb2dt,Ab2,dAu1dt,dAu2dt
  complex :: TimeGe, TimeAb, TimeAu
  real    :: kx,ky,kz,wk,wkh2,wkh,wkh2n,kzn,kz2,r1
  real    :: aj,bj,vz,vt,time,visch,viscz,kvisc,kh_g
  real    :: vzx,vzy,vzz,kxa(iktx),kya(ikty),kza(iktz),L1,L2,L3,twopi
  real, dimension(0:kts,3) :: spz,spztot
  ! SPZ(:,1)= G/(dG/dt) Geostrophic Time Scale,
  ! SPZ(:,2)= Ab1/(dAb1/dt)+Ab2/(dAb2/dt) Balanced Ageostrophic Time Scale
  ! SPZ(:,3)= Au1/(dAu1/dt)+Au2/(dAu2/dt) Unbalanced Ageostrophic Time Scale
  integer, dimension(0:kts) :: n,ntot

  integer :: mype
  common/mpi/mype

  nspz  = kts+1
  nspz8 = (kts+1)*8

  if (ispec.eq.1) then
     j0 = 1
  elseif (ispec.eq.2) then
     j0 = 0
  elseif (ispec.eq.3) then
     j0 = 0
  else
     print*,"ispec =",ispec
     print*,"ispec error in spcTimeScales"
     stop
  endif

  twopi = 4.*asin(1.)
  N     = 0
  spz   = 0.

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
           
           if (ispec.eq.1) then
              j = ifix(wk*L1/twopi+0.5)
           elseif (ispec.eq.2) then
              j = ifix(wkh*L1/twopi+0.5)
           elseif (ispec.eq.3) then
              j = ifix(abs(kz)*L3/twopi+0.5)
           endif

           if (L(ikx,iky,ikz).eq.1) then
              if (j.lt.j0 .or. j.gt.kts) print*,'SPEC: SCREW-UP.   k= ',j,'    ktx= ',ktx
#if (SPHTRUNC == 1)
                    r1 = visch/delt*wk**2**ilap
#else
                    r1 = visch*wkh2**ilap + viscz*kz2**ilap
#endif 

! Derive the time scales
              TimeGe   = ge(ikx,iky,ikz)/(dGdt(ikx,iky,ikz)-r1*ge(ikx,iky,ikz)) 
              spz(j,1) = spz(j,1) + sqrt(real(TimeGe*conjg(TimeGe)))

              TimeAb   = Ab1(ikx,iky,ikz)/dAb1dt(ikx,iky,ikz) 
              spz(j,2) = spz(j,2) + sqrt(real(TimeAb*conjg(TimeAb)))
              TimeAb   = Ab2(ikx,iky,ikz)/dAb2dt(ikx,iky,ikz) 
              spz(j,2) = spz(j,2) + sqrt(real(TimeAb*conjg(TimeAb)))

              TimeAu   = (g1(ikx,iky,ikz)-Ab1(ikx,iky,ikz))/dAu1dt(ikx,iky,ikz) 
              spz(j,3) = spz(j,3) + sqrt(real(TimeAu*conjg(TimeAu)))
              TimeAu   = (g2(ikx,iky,ikz)-Ab2(ikx,iky,ikz))/dAu2dt(ikx,iky,ikz) 
              spz(j,3) = spz(j,3) + sqrt(real(TimeAu*conjg(TimeAu)))

              n(j) = n(j) + 2 
           endif
        enddo
     enddo
  enddo

#if ( MPI == 1) 
  call mpi_reduce(spz,spztot,nspz8,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(  n,  ntot, nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  spztot=spz
  ntot=n
#endif

  if (mype.eq.0) then
     do j=j0,kts-1
        if (ntot(j).ne.0) then
           write(iu,5023) float(j),(spztot(j,i),i=1,3),ntot(j)
        endif
     enddo
  write(iu,*) '           '
     call flush(iu)
  endif

  return
5023  format(1X,F4.0,4X,3(E14.6,1x),6X,I10)
end subroutine specTimeScales

subroutine spec2d (zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,L,time,kxa,kya,kza,aj,bj,visch,viscz,ilap,L1,L2,L3,iu)

! Calculates spectra: 
! ispec = 1: k
! ispec = 2: kh
! ispec = 3: kz

  implicit none 
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: ikx,iky,ikz,ikza,ispec
  integer :: jh,jv,i,L(iktx,ikty,iktzp),iu,j0,nspz,nspz8,ierr,ilap
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,ge,g1,g2,ux,uy,uz
  real    :: kx,ky,kz,wk,wkh2,wkh,wkh2n,kzn,kz2
  real    :: aj,bj,vz,vt,time,visch,viscz,kvisc
  real    :: vzx,vzy,vzz,kxa(iktx),kya(ikty),kza(iktz),L1,L2,L3,twopi
  real, dimension(0:kts,0:ktz,8) :: spz2d,spztot2d
  integer, dimension(0:kts,0:ktz) :: nn2d,ntot2d
  external :: velo

  integer :: mype
  common/mpi/mype

  nspz  = (kts+1)*(ktz+1)
  nspz8 = (kts+1)*(ktz+1)*8

  

  twopi = 4.*asin(1.)
  nn2d    = 0
  spz2d   = 0.

  call velo (zx,zy,zz,ux,uy,uz,L,kxa,kya,kza)


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
	   kvisc = visch*wkh2**ilap+viscz*kz2**ilap
           
           jh = ifix(wkh*L1/twopi+0.5)
           jv = ifix(abs(kz)*L3/twopi+0.5)
              
           if (L(ikx,iky,ikz).eq.1) then
              if (jv.lt.j0 .or. jv.gt.ktz) print*,'SPEC: SCREW-UP.   kz= ',jv,'    ktz= ',ktz
              if (jh.lt.j0 .or. jh.gt.ktx) print*,'SPEC: SCREW-UP.   kh= ',jh,'    ktx= ',ktx

              wk   = max(wk,  1.e-15)
              wkh2 = max(wkh2,1.e-15)

! Kinetic and potential energy.
              vzx      = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
              spz2d(jh,jv,1) = spz2d(jh,jv,1) + vzx/wk**2
              vzy      = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
              spz2d(jh,jv,1) = spz2d(jh,jv,1) + vzy/wk**2
              vzz      = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
              spz2d(jh,jv,1) = spz2d(jh,jv,1) + vzz/wk**2
              vt       = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )
              spz2d(jh,jv,2) = spz2d(jh,jv,2) + vt*aj/bj

! Geo, ageo decompostition.
! k \in R_k
              if (wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 vz       = real( ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)) )
                 spz2d(jh,jv,4) = spz2d(jh,jv,4) + vz/wkh2
                 spz2d(jh,jv,7) = spz2d(jh,jv,7) + kvisc*vz/wkh2
                 vz       = real( g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)) )
                 spz2d(jh,jv,5) = spz2d(jh,jv,5) + vz/wkh2
                 spz2d(jh,jv,8) = spz2d(jh,jv,8) + kvisc*vz/wkh2
                 vz       = real( g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)) )
                 spz2d(jh,jv,5) = spz2d(jh,jv,5) + vz/wkh2
                 spz2d(jh,jv,8) = spz2d(jh,jv,8) + kvisc*vz/wkh2
              endif

! Special cases: i) k_h=0, ii) k_z=0.
              vzx=real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              vzy=real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              vzz=real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )

! k \in V_k
              if(wkh2n.lt.1.e-10.and.abs(kzn).gt.1.e-10) then
                 spz2d(jh,jv,4) = spz2d(jh,jv,4) + vzz + vt*aj/bj
                 spz2d(jh,jv,5) = spz2d(jh,jv,5) + vzx + vzy
                 spz2d(jh,jv,7) = spz2d(jh,jv,7) + kvisc*(vzz + vt*aj/bj)
                 spz2d(jh,jv,8) = spz2d(jh,jv,8) + kvisc*(vzx + vzy)     
              endif

! k \in B_k
              if(abs(kzn).lt.1.e-10.and.wkh2n.gt.1.e-10) then
                 spz2d(jh,jv,4) = spz2d(jh,jv,4) + vzx + vzy
                 spz2d(jh,jv,5) = spz2d(jh,jv,5) + vzz + vt*aj/bj
                 spz2d(jh,jv,7) = spz2d(jh,jv,7) + kvisc*(vzx + vzy)     
                 spz2d(jh,jv,8) = spz2d(jh,jv,8) + kvisc*(vzz + vt*aj/bj)
              endif

! Buoyancy Flux
              if (aj.gt.0.) then 
                 spz2d(jh,jv,6) = spz2d(jh,jv,6) + aj*real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
              else
                 spz2d(jh,jv,6) = spz2d(jh,jv,6) +    real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
              endif

              nn2d(jh,jv) = nn2d(jh,jv) + 2 
           endif
        enddo
     enddo
  enddo
  
  spz2d(:,:,3) = spz2d(:,:,1) + spz2d(:,:,2)

#if ( MPI == 1) 
  call mpi_reduce(spz2d,spztot2d,nspz8,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(  nn2d,  ntot2d, nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  spztot2d=spz2d
  ntot2d=nn2d
#endif

  if (mype.eq.0) then
     do jh=0,kts-1
        do jv=0,ktz-1
           if (ntot2d(jh,jv).ne.0) then
              write(iu,5000) float(jh),float(jv),(spztot2d(jh,jv,i),i=1,8),ntot2d(jh,jv)
           endif
        enddo
     enddo
     write(iu,*) '           '
     call flush(iu)
  endif

  return

5000  format(1X,F4.0,2X,F4.0,4X,8(E14.6,1x),6X,I6)

end subroutine spec2d


!!$subroutine wn4spc3d(kxa,kya,kza,ktxsp,ktysp,ktzsp,iu)
!!$
!!$   implicit none 
!!$  include 'param.inc'
!!$#if (MPI == 1)
!!$  include 'mpif.h'
!!$#endif
!!$
!!$  integer :: ikx,iky,ikz,ikza,nk3d,iu,ierr
!!$  real    :: kx,ky,kz
!!$  real  :: ktxsp,ktysp,ktzsp
!!$  real    :: kxa(iktx),kya(ikty),kza(iktz)
!!$  real, dimension(iktx,ikty,iktz) :: kx3d,kx3dtot,ky3d,ky3dtot,kz3d,kz3dtot
!!$  integer, dimension(iktx,ikty,iktz) :: nmpty,nmptytot
!!$
!!$  integer :: mype
!!$  common/mpi/mype
!!$
!!$  nk3d = (iktx)*(ikty)*(iktz)
!!$  kx3d = 0.
!!$  ky3d = 0.
!!$  kz3d = 0.
!!$  nmpty = 0
!!$  
!!$  do ikz=1,iktzp
!!$     ikza = mype*iktzp+ikz
!!$     kz = kza(ikza)
!!$     if (abs(kz).le.abs(ktzsp)) then
!!$        do iky=1,ikty
!!$           ky = kya(iky)
!!$           if (abs(ky).le.abs(ktysp)) then
!!$              do ikx=1,ktxsp
!!$                 kx = kxa(ikx)
!!$                 if (abs(kx).le.abs(ktxsp)) then
!!$                    kx3d(ikx,iky,ikza) = kx
!!$                    ky3d(ikx,iky,ikza) = ky
!!$                    kz3d(ikx,iky,ikza) = kz
!!$                    nmpty(ikx,iky,ikza) = nmpty(ikx,iky,ikza) + 1
!!$                 endif
!!$              enddo
!!$           endif
!!$        enddo
!!$     endif
!!$  enddo
!!$
!!$#if ( MPI == 1) 
!!$  call mpi_reduce(kx3d,kx3dtot,nk3d,MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!$  call mpi_reduce(ky3d,ky3dtot,nk3d,MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!$  call mpi_reduce(kz3d,kz3dtot,nk3d,MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!$  call mpi_reduce(nmpty,nmptytot,nk3d,MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!$#else
!!$  kx3dtot = kx3d
!!$  ky3dtot = ky3d
!!$  kz3dtot = kz3d
!!$  nmptytot = nmpty
!!$#endif
!!$
!!$  if (mype.eq.0) then
!!$     do ikz=1,iktz
!!$        do iky=1,ikty
!!$           do ikx=1,ktxsp
!!$              if (nmpty(ikx,iky,ikz).ne.0) then
!!$                 write(iu,5013) kx3dtot(ikx,iky,ikz),ky3dtot(ikx,iky,ikz),kz3dtot(ikx,iky,ikz)
!!$              endif
!!$           enddo
!!$        enddo
!!$     enddo
!!$     write(iu,*) '           '
!!$     call flush(iu)
!!$  endif ! mype = 0
!!$  
!!$  return
!!$
!!$5013  format(1X,3(F8.2,1x),6X)
!!$ end subroutine wn4spc3d 

subroutine wvspc3D(g1,g2,time,nspdump,kxa,kya,kza,ikxsp,ikysp,ikzsp,dumpsp,totnzsp)
  
  implicit none 
  include 'param.inc'
  include 'netcdf.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: ikx,iky,ikz,ikza,istatus,ikz1,ikzsp1,nkzsp,nspdump
  integer :: nbufsend,nbufrec,nsends,totnzsp,iproc,sourcid
  integer :: ikxsp(iktxsp),ikysp(iktysp-1),ikzsp(iktzsp-1),dumpsp(4),dumpinf(4)
  real    :: kxa(iktx),kya(ikty),kza(iktz),kx,ky,kz,wkh2,time 
  complex, dimension(iktx,ikty,iktzp) :: g1,g2
  real,dimension(iktxsp,iktysp-1,iktzp) :: ew1redc,ew2redc

  !! 3D spectrum of waves (and potentially other variables)
  integer :: idnsp,idkxsp,idkysp,idkzsp,idtsp
  integer, dimension(4) :: ncdimsp,ncstartsp,nccountsp
  integer, dimension(3) :: ncdimkk,ncstartkk,nccountkk
  integer :: idtimesp,idkkxx,idkkyy,idkkzz,idwvp,idwvm
  
  common/netcdfwvsp/idnsp,idkxsp,idkysp,idkzsp,idtsp,ncstartsp,nccountsp, &
       ncstartkk,nccountkk,idtimesp,idkkxx,idkkyy,idkkzz,idwvp,idwvm

#if ( MPI == 1 )
  integer :: status(MPI_STATUS_SIZE)
#endif 
  integer :: mype
  common/mpi/mype

  ew1redc = 0.0
  ew2redc = 0.0
  
  nkzsp = dumpsp(4) ! number of vert modes in current proc that is going to be dumped
  ikz1  = dumpsp(3) ! index of kza in the current proc from which
                    ! we start sending variables to root to write in the output file


  if (nkzsp.ne.0) then
     
     nbufsend = iktxsp*(iktysp-1)*nkzsp
     nsends   = (totnzsp-1)
     ncstartsp(1) = 1
     ncstartsp(2) = 1
     ncstartsp(4) = nspdump
     nccountsp(1) = iktxsp
     nccountsp(2) = iktysp-1
     nccountsp(4) = 1

     do ikx = 1, iktxsp
        kx = kxa(ikxsp(ikx))
        do iky = 1,iktysp-1
           ky = kya(ikysp(iky))
           wkh2  = kx*kx+ky*ky
           if (wkh2.gt.1.e-10) then
              ew1redc(ikx,iky,1:nkzsp)=real(g1(ikxsp(ikx),ikysp(iky),ikz1:ikz1+nkzsp-1) &
                   *conjg(g1(ikxsp(ikx),ikysp(iky),ikz1:ikz1+nkzsp-1)))/wkh2
              ew2redc(ikx,iky,1:nkzsp)=real(g2(ikxsp(ikx),ikysp(iky),ikz1:ikz1+nkzsp-1) &
                   *conjg(g2(ikxsp(ikx),ikysp(iky),ikz1:ikz1+nkzsp-1)))/wkh2
           endif
        enddo
     enddo

     if (mype.eq.0) then
        istatus = nf_put_vara_real(idnsp,idtimesp,nspdump,1,time) 
     endif

     if (mype.gt.0) then
        call mpi_send(dumpsp,4,MPI_INTEGER,0,701,MPI_COMM_WORLD,istatus)
        call mpi_send(ew1redc,nbufsend,MPI_REAL,0,731,MPI_COMM_WORLD,istatus)
        call mpi_send(ew2redc,nbufsend,MPI_REAL,0,732,MPI_COMM_WORLD,istatus)
     endif

     if (mype.eq.0) then
        ncstartsp(3) = dumpsp(2)
        nccountsp(3) = nkzsp
        istatus = nf_put_vara_real(idnsp,idwvp,ncstartsp,nccountsp,ew1redc)
        istatus = nf_put_vara_real(idnsp,idwvm,ncstartsp,nccountsp,ew2redc)

        do iproc=1,nsends
           call mpi_recv(dumpinf,4,MPI_INTEGER,MPI_ANY_SOURCE,701,MPI_COMM_WORLD,status,istatus)
           sourcid = dumpinf(1)
           nbufrec = iktxsp*(iktysp-1)*dumpinf(4)
           ncstartsp(3) = dumpinf(2)
           nccountsp(3) = dumpinf(4)
           call mpi_recv(ew1redc,nbufrec,MPI_REAL,sourcid,731,MPI_COMM_WORLD,status,istatus)
           istatus = nf_put_vara_real(idnsp,idwvp,ncstartsp,nccountsp,ew1redc)
           call mpi_recv(ew2redc,nbufrec,MPI_REAL,sourcid,732,MPI_COMM_WORLD,status,istatus)
           istatus = nf_put_vara_real(idnsp,idwvm,ncstartsp,nccountsp,ew2redc)
        enddo
     endif ! mype = 0
           
  endif ! nkzsp > 0
  
 

end subroutine wvspc3D




subroutine transf(zx,zy,zz,tt,geok,gw1k,gw2k,nzx,nzy,nzz,ntt,ngeok,ngw1k,ngw2k, &
                  nk1,nk2,nk3,L,time,kxa,kya,kza,aj,bj,f,L1,L2,L3,ispec,iu)

! Calculates transfer spectra.
! ispec = 1: k
! ispec = 2: kh
! ispec = 3: kz

  implicit none 
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: ikx,iky,ikz,ikza,j,iu
  integer :: L(iktx,ikty,iktzp),ispec,j0,nspz,nspz6,ierr
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,geok,gw1k,gw2k,nk1,nk2,nk3
  complex, dimension(iktx,ikty,iktzp) :: nzx,nzy,nzz,ntt,ngeok,ngw1k,ngw2k
  complex :: u,v,w,c1,c2,c3,zi
  real :: kx,ky,kz,k,k2,kh2,wkh,kh2n,kzn
  real :: aj,bj,F,time
  real :: kxa(iktx),kya(ikty),kza(iktz)
  real :: vzx,vzy,vzz,vtt,L1,L2,L3,twopi
  real, dimension(0:kts,6)  :: spz,spztot
  integer, dimension(0:kts) :: n,ntot

  integer :: mype
  common/mpi/mype

  nspz  = kts+1
  nspz6 = (kts+1)*6

  if (ispec.eq.1) then
     j0 = 1
  elseif (ispec.eq.2) then
     j0 = 0
  elseif (ispec.eq.3) then
     j0 = 0
  else
     print*,"ispec error"
     stop
  endif
  
  zi = cmplx(0.,1.)
  twopi = 4.*asin(1.)

  call wtoab( zx, zy, zz, tt, geok, gw1k, gw2k,nk1,nk2,nk3,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)
  call wtoab(nzx,nzy,nzz,ntt,ngeok,ngw1k,ngw2k,nk1,nk2,nk3,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)
  call velo(nzx,nzy,nzz,nk1,nk2,nk3,L,kxa,kya,kza)

  N   = 0
  spz = 0.

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kzn = kz*L3/twopi
     do iky = 1,ikty
        ky  = kya(iky)
        do ikx = 1,iktx
           kx   = kxa(ikx)
           kh2  = kx*kx + ky*ky
           kh2n = kh2 * (L1/twopi)**2
           kh2  = max(1.e-15,kh2)
           wkh  = sqrt(kh2)
           k2   = kx*kx + ky*ky + kz*kz
           k    = sqrt(k2)
           k2   = max(1.e-15,k2)
           
           if (ispec.eq.1) then
              j = ifix(k*L1/twopi+0.5)
           elseif (ispec.eq.2) then
              j = ifix(wkh*L1/twopi+0.5)
           elseif (ispec.eq.3) then
              j = ifix(abs(kz)*L3/twopi+0.5)
           endif

           if (L(ikx,iky,ikz).eq.1)  then
              if (j.lt.j0 .or. j.gt.kts) then
                 print*,'transf: screw-up.   k= ',j,kx,ky,kz,L(ikx,iky,ikz)
              endif

              c1 = ky*zz(ikx,iky,ikz) - kz*zy(ikx,iky,ikz)
              c2 = kz*zx(ikx,iky,ikz) - kx*zz(ikx,iky,ikz)
              c3 = kx*zy(ikx,iky,ikz) - ky*zx(ikx,iky,ikz)
              u = zi * c1 / k2
              v = zi * c2 / k2
              w = zi * c3 / k2

! k \in R_k
              if (kh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 spz(j,1) = spz(j,1) + real(geok(ikx,iky,ikz)*conjg(ngeok(ikx,iky,ikz)) )/kh2
                 spz(j,2) = spz(j,2) + real(gw1k(ikx,iky,ikz)*conjg(ngw1k(ikx,iky,ikz)) )/kh2
                 spz(j,2) = spz(j,2) + real(gw2k(ikx,iky,ikz)*conjg(ngw2k(ikx,iky,ikz)) )/kh2
              endif

! Special cases: i) k_h=0, ii) k_z=0.
              vzx = real(u              *conjg( nk1(ikx,iky,ikz)) )
              vzy = real(v              *conjg( nk2(ikx,iky,ikz)) )
              vzz = real(w              *conjg( nk3(ikx,iky,ikz)) )
              vtt = real(tt(ikx,iky,ikz)*conjg(ntt(ikx,iky,ikz)) )

!  k \in V_k
              if (kh2n.lt.1.e-10.and.abs(kzn).gt.1.e-10) then
                 spz(j,1) = spz(j,1) + vtt
                 spz(j,2) = spz(j,2) + vzx + vzy            
              endif

!  k \in B_k
              if (kh2n.gt.1.e-10.and.abs(kzn).lt.1.e-10) then
                 spz(j,1) = spz(j,1) + vzx + vzy 
                 spz(j,2) = spz(j,2) + vzz + vtt*aj/bj
              endif

!  Now, compute KE and PE transfer
              spz(j,4) = spz(j,4) + vzx + vzy + vzz
              spz(j,5) = spz(j,5) + vtt*aj/bj

              n(j)   = n(j) + 2
           endif
        enddo
     enddo
  enddo

  spz(:,3) = spz(:,1) + spz(:,2)
  spz(:,6) = spz(:,4) + spz(:,5)

#if ( MPI == 1 )
  call mpi_reduce(spz,spztot,nspz6,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(  n,  ntot, nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  spztot=spz
  ntot=n
#endif

  if (mype.eq.0) then  
     do j = j0,kts-1	
        if (ntot(j).ne.0) then
           write(iu,5000) float(j),(spztot(j,iky),iky=1,6),ntot(j)
        endif
     enddo
     write(iu,*) '     '
     write(iu,*) '     '
  endif

  return
5000 format(1X,F4.0,4X,6(E12.4,1x),10X,I6)
5010 format(1X, 6(E12.4,1x))
5020 format(1x, 'Isotropic transfer totals at time t = ',f12.6, ":")
end subroutine transf


subroutine transf2unbal(zxunbalk,zyunbalk,zzunbalk,ttunbalk,geok,gw1k,gw2k,dAudt1,dAudt2, &
                  nk1,nk2,nk3,L,time,kxa,kya,kza,aj,bj,f,delt,v2h,k2h,v2z,k2z,ilap,ilap2,L1,L2,L3,ispec,iu)

! Calculates transfer spectra.
! ispec = 1: k
! ispec = 2: kh
! ispec = 3: kz

  implicit none 
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: ikx,iky,ikz,ikza,j,iu
  integer :: L(iktx,ikty,iktzp),ispec,j0,nspz,nspz6,ierr
  complex, dimension(iktx,ikty,iktzp) :: zxunbalk,zyunbalk,zzunbalk,ttunbalk,geok,gw1k,gw2k,nk1,nk2,nk3
  complex, dimension(iktx,ikty,iktzp) :: dAudt1,dAudt2
  complex :: u,v,w,c1,c2,c3,zi
  real :: kx,ky,kz,k,k2,kh2,wkh,kh2n,kzn,freq,wk2,khn,kh
  real :: aj,bj,f,time,delt,v2h,k2h,v2z,k2z,r1,r2
  real :: kxa(iktx),kya(ikty),kza(iktz)
  real :: vzx,vzy,vzz,vtt,L1,L2,L3,twopi
  real, dimension(0:kts)  :: spz,spztot
  integer, dimension(0:kts) :: n,ntot
  integer :: ilap,ilap2

  integer :: mype
  common/mpi/mype

  nspz  = kts+1
  nspz6 = (kts+1)*6

  if (ispec.eq.1) then
     j0 = 1
  elseif (ispec.eq.2) then
     j0 = 0
  elseif (ispec.eq.3) then
     j0 = 0
  else
     print*,"ispec error"
     stop
  endif
  
  zi = cmplx(0.,1.)
  twopi = 4.*asin(1.)

  call wtoab(zxunbalk,zyunbalk,zzunbalk,ttunbalk,geok,gw1k,gw2k,nk1,nk2,nk3,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)
 
  N   = 0
  spz = 0.

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kzn = kz*L3/twopi
     do iky = 1,ikty
        ky  = kya(iky)
        do ikx = 1,iktx
           kx   = kxa(ikx)
           kh2  = kx*kx + ky*ky
           kh  = sqrt( kx*kx + ky*ky )
           kh2n = kh2 * (L1/twopi)**2
           khn = kh * L1/twopi
           wk2  = kx*kx + ky*ky + kz*kz
           k2   = kx*kx + ky*ky + kz*kz
           wkh  = sqrt(kh2)
           k  = sqrt(wk2)
           k2   = max(1.e-15,k2)
           freq = sqrt(kz*kz*f**2 + kh*kh*AJ*BJ)/k   
#if (SPHTRUNC == 1)
                    r1 = v2h/delt*wk2**ilap
                    r2 = k2h/delt*wk2**ilap
#else
                    r1 = v2h/delt*kh**ilap2 + v2z/delt*kz**ilap2
                    r2 = k2h/delt*kh**ilap2 + k2z/delt*kz**ilap2
#endif  
           if (ispec.eq.1) then
              j = ifix(k*L1/twopi+0.5)
           elseif (ispec.eq.2) then
              j = ifix(wkh*L1/twopi+0.5)
           elseif (ispec.eq.3) then
              j = ifix(abs(kz)*L3/twopi+0.5)
           endif

           if (L(ikx,iky,ikz).eq.1)  then
              if (j.lt.j0 .or. j.gt.kts) then
                 print*,'transf: screw-up.   k= ',j,kx,ky,kz,L(ikx,iky,ikz)
              endif
! k \in R_k
              if (kh2n.gt.1.e-10) then
                 spz(j) = spz(j) + real(gw1k(ikx,iky,ikz)*(conjg(dAudt1(ikx,iky,ikz)-(+zi*freq-r1)*gw1k(ikx,iky,ikz))))/kh2
                 spz(j) = spz(j) + real(gw2k(ikx,iky,ikz)*(conjg(dAudt2(ikx,iky,ikz)-(-zi*freq-r1)*gw2k(ikx,iky,ikz))))/kh2
              endif
              n(j)   = n(j) + 2
           endif
        enddo
     enddo
  enddo


#if ( MPI == 1 )
  call mpi_reduce(spz,spztot,nspz6,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(  n,  ntot, nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  spztot=spz
  ntot=n
#endif

  if (mype.eq.0) then  
     do j = j0,kts-1	
        if (ntot(j).ne.0) then
           write(iu,5000) float(j),spztot(j),ntot(j)
        endif
     enddo
     write(iu,*) '     '
     write(iu,*) '     '
  endif

  return
5000 format(1X,F4.0,4X,1(E12.4,1x),10X,I6)
5010 format(1X, 1(E12.4,1x))
end subroutine transf2unbal


subroutine transtrms2g(NGGk,NGAbk,NGAuk,NAbAbk,NAbAuk,NAuAuk,zx,zy,zz,tt,geok,gw1k,gw2k, &
     uk,vk,wk,L,time,kxa,kya,kza,aj,bj,f,L1,L2,L3,ispec,iu)

! Calculates transfer spectra.
! ispec = 1: k
! ispec = 2: kh
! ispec = 3: kz

  implicit none 
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: ikx,iky,ikz,ikza,j,iu
  integer :: L(iktx,ikty,iktzp),ispec,j0,nspz,nspz6,ierr
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,geok,gw1k,gw2k,uk,vk,wk
  complex, dimension(iktx,ikty,iktzp) :: NGGk,NGAbk,NGAuk,NAbAbk,NAbAuk,NAuAuk
  real :: kx,ky,kz,k,k2,kh2,wkh,kh2n,kzn
  real :: aj,bj,F,time
  real :: kxa(iktx),kya(ikty),kza(iktz)
  real :: L1,L2,L3,twopi
  complex :: zi
  real, dimension(0:kts,6)  :: spz,spztot
  integer, dimension(0:kts) :: n,ntot

  integer :: mype
  common/mpi/mype

  nspz  = kts+1
  nspz6 = (kts+1)*6

  if (ispec.eq.1) then
     j0 = 1
  elseif (ispec.eq.2) then
     j0 = 0
  elseif (ispec.eq.3) then
     j0 = 0
  else
     print*,"ispec error"
     stop
  endif
  
  zi = cmplx(0.,1.)
  twopi = 4.*asin(1.)

  call wtoab( zx, zy, zz, tt, geok, gw1k, gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)
  
  N   = 0
  spz = 0.

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kzn = kz*L3/twopi
     do iky = 1,ikty
        ky  = kya(iky)
        do ikx = 1,iktx
           kx   = kxa(ikx)
           kh2  = kx*kx + ky*ky
           kh2n = kh2 * (L1/twopi)**2
           kh2  = max(1.e-15,kh2)
           wkh  = sqrt(kh2)
           k2   = kx*kx + ky*ky + kz*kz
           k    = sqrt(k2)
           k2   = max(1.e-15,k2)
           
           if (ispec.eq.1) then
              j = ifix(k*L1/twopi+0.5)
           elseif (ispec.eq.2) then
              j = ifix(wkh*L1/twopi+0.5)
           elseif (ispec.eq.3) then
              j = ifix(abs(kz)*L3/twopi+0.5)
           endif

           if (L(ikx,iky,ikz).eq.1)  then
              if (j.lt.j0 .or. j.gt.kts) then
                 print*,'transf: screw-up.   k= ',j,kx,ky,kz,L(ikx,iky,ikz)
              endif        
! k \in R_k
              if (kh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 spz(j,1) = spz(j,1) + real(geok(ikx,iky,ikz)*conjg(NGGk(ikx,iky,ikz)) )
                 spz(j,2) = spz(j,2) + real(geok(ikx,iky,ikz)*conjg(NGAbk(ikx,iky,ikz)) )
                 spz(j,3) = spz(j,3) + real(geok(ikx,iky,ikz)*conjg(NGAuk(ikx,iky,ikz)) )
                 spz(j,4) = spz(j,4) + real(geok(ikx,iky,ikz)*conjg(NAbAbk(ikx,iky,ikz)) )
                 spz(j,5) = spz(j,5) + real(geok(ikx,iky,ikz)*conjg(NAbAuk(ikx,iky,ikz)) )
                 spz(j,6) = spz(j,6) + real(geok(ikx,iky,ikz)*conjg(NAuAuk(ikx,iky,ikz)) )
              endif
              n(j)   = n(j) + 2
           endif
        enddo
     enddo
  enddo

#if ( MPI == 1 )
  call mpi_reduce(spz,spztot,nspz6,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(  n,  ntot, nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  spztot=spz
  ntot=n
#endif

  if (mype.eq.0) then  
     do j = j0,kts-1	
        if (ntot(j).ne.0) then
           write(iu,5202) float(j),(spztot(j,iky),iky=1,6),ntot(j)
        endif
     enddo
     write(iu,*) '     '
     write(iu,*) '     '
  endif

  return
5202 format(1X,F4.0,4X,6(E12.4,1x),10X,I6)
end subroutine transtrms2g


subroutine out(zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,nt,delt,ts,  &
               NSTOP,io,L,kxa,kya,kza,aj,bj,f,visch,viscz,ilap,L1,L2,L3,kh_g)

  implicit none
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif
  
  integer :: ikx,iky,ikz,ikza,ilap
  integer, dimension(iktx,ikty,iktzp) :: L
  integer :: nt,nstop,io,ierr
  complex, dimension(iktx,ikty,iktzp) ::  zx,zy,zz,ux,uy,uz,tt,ge,g1,g2
  complex :: der,zi
  real :: kxa(iktx),kya(ikty),kza(iktz)
  real :: aj,bj,ts,delt,e,kx,ky,kz,time,wk,wkh2,wkh2n,kzn
  real :: vh,vz,L1,L2,L3,twopi
  real :: visch,viscz,vzx,vzy,vzz,f,n,bv2,ffac,n2fac,sigma
  real :: vort2,div2,rossby,fr_z,fr_h,ke,pe,eg,ea,rossbyVel
  real :: pv,h,v,epsk,epsp,eps,tmp
  real :: zero_kz_geo,zero_kz_grv,zero_kh_grv,zero_kh_geo,kh_g
  external :: proj, velo

  integer :: mype
  common/mpi/mype

  twopi = 4.*asin(1.)
  io    = io + 1
  zi    = cmplx(0.,1.)
  time  = ts + nt*delt
  bv2   = aj*bj
  n     = sqrt(bv2)
  ffac  = max(f,1.e-15)
  n2fac = max(bv2,1.e-15)
  
  call proj(zx,zy,zz,L,kxa,kya,kza)
  call velo(zx,zy,zz,ux,uy,uz,L,kxa,kya,kza)

  vort2  = 0.
  div2   = 0.
  Rossby = 0.
  Fr_z   = 0.
  Fr_h   = 0.
  zero_kz_geo = 0.
  zero_kz_grv = 0.
  zero_kh_geo = 0.
  zero_kh_grv = 0.
  ke = 0.
  pe = 0.
  eg = 0.
  ea = 0.
  pv = 0.
  h  = 0.
  v  = 0.
  epsk = 0.
  epsp = 0.
  eps  = 0.

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kzn = kz * (L3/twopi)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx = 1,iktx
           kx = kxa(ikx)
           wkh2 = kx*kx + ky*ky
           wk = kx*kx + ky*ky + kz*kz
           if (L(ikx,iky,ikz).eq.1) then
              wkh2n = wkh2 * (L1/twopi)**2
              
              vzx   = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
              ke    = ke + vzx/wk
              
              epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vzx/wk
              v     = v  + vzx
              vzy   = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
              ke    = ke + vzy/wk

              epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vzy/wk
              v     = v + vzy
              vzz   = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
              ke    = ke + vzz/wk

              epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vzz/wk
              v     = v  + vzz
              vort2 = vort2 + vzz

              vh    = real( zx(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              h     = h + vh
              vh    = real( zy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              h     = h + vh
              vh    = real( zz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )
              h     = h + vh
              vh    = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )
              pe    = pe  + vh
              epsp  = epsp + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vh*aj/bj

              vzx   = real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              vzy   = real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              vzz   = real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )

              if(wkh2n.lt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 zero_kh_grv = zero_kh_grv + vzx + vzy
                 zero_kh_geo = zero_kh_geo + vh*aj/bj
                 pv = pv + F*F*kz*kz*vh*aj/bj
              endif
              
              if(wkh2n.gt.1.e-10 .and. abs(kzn).lt.1.e-10) then
                 zero_kz_geo = zero_kz_geo + vzx + vzy
                 zero_kz_grv = zero_kz_grv + vzz + vh*aj/bj
                 pv = pv + aj*bj*wkh2*(vzx + vzy)
              endif
              
              if(wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 vz    = real(ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)))
                 eg    = eg  + vz/wkh2
                 sigma = wkh2*bv2 + kz*kz*F*F
                 pv    = pv + sigma*vz/wkh2
                 vz    = real(g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)))
                 ea    = ea  + vz/wkh2
                 vz    = real(g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)))
                 ea    = ea  + vz/wkh2
              endif
              
              rossby = rossby + real(zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)))/(ffac*ffac)
              der    =          real(zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)))
              der    = (der +   real(zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz))))
              fr_z   = fr_z + der/N2FAC
              der    = zi*(kx*ux(ikx,iky,ikz) + ky*uy(ikx,iky,ikz))
              div2   = div2 + real(der*conjg(der))
           endif ! L
        enddo
     enddo
  enddo

  v      = 2.*v
  epsk   = 2.*epsk
  epsp   = 2.*epsp
  eg     = eg + zero_kz_geo + zero_kh_geo
  ea     = ea + zero_kz_grv + zero_kh_grv
  if (aj.ne.0. .and. bj.ne.0.) pe = aj*pe/bj

#if ( MPI == 1 )
  call mpi_reduce(ke,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ke=tmp
  call mpi_reduce(pe,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  pe=tmp
  call mpi_reduce(eg,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  eg=tmp
  call mpi_reduce(ea,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ea=tmp
  call mpi_reduce(pv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  pv=tmp
  call mpi_reduce(rossby,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rossby=tmp
  call mpi_reduce(fr_z,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  fr_z=tmp
  call mpi_reduce(vort2,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  vort2=tmp
  call mpi_reduce(v,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  v=tmp
  call mpi_reduce(zero_kz_grv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  zero_kz_grv=tmp
  call mpi_reduce(zero_kz_geo,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  zero_kz_geo=tmp
  call mpi_reduce(zero_kh_grv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  zero_kh_grv=tmp
  call mpi_reduce(zero_kh_geo,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  zero_kh_geo=tmp
  call mpi_reduce(epsk,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  epsk=tmp
  call mpi_reduce(epsp,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  epsp=tmp
#endif

  if (mype.eq.0) then ! prep for output

     if (f.gt.1.e-8) then
        rossby = sqrt(2.*rossby)
     else
        rossby = - 999.
     endif

     rossbyVel = sqrt(eg*2.)/(twopi/kh_g)/ffac

     if (N.gt.1.e-8) then
        fr_z = sqrt(fr_z)
        fr_h = sqrt(2.*vort2)/N
     else
        fr_z = - 999.
        fr_h = - 999.
     endif
     eps    = epsk + epsp
     if (aj.ne.0. .and. bj.ne.0.) then
        e  = (pe + ke)
     else
        e = 99999999.
     endif

     write(79,5046) time, epsk, epsp, eps
     write(46,5045) time,ke,pe,e,eg,ea,pv,rossby,fr_z,fr_h,zero_kh_grv,zero_kh_geo,v,rossbyVel,kh_g
     if (nt.eq.0) then
!        write(6,*) 'Rossby Velocity = ', rossbyVel
!        write(6,*) 'Geo spectrum peaks at  ', kh_g
        write(6,5043) 
        write(6,5042) '------------------------------------------------------------------------------------------'
     endif
     write( 6,5044) time,ke,pe,e,eg,ea,rossby,fr_z,fr_h,v,rossbyVel,kh_g
     call flush(46)
     call flush(79)
  endif
  return

5042  format(1x,a91)
5043  format(7x,' T',9x,'KE',8x,'PE',8x,'E',9x,'GE',8x,'AE',8x,'Ro',8x,'Fr_z',7x,'Fr_h',5x,'V',9x,'RoVel',8x,'kh_peak')
5044  format(1x,f12.2,2x,9(f8.3,2x),1x,f10.6,1x,f8.3,1x)
5045  format(1x,f12.2,2x,12(e12.4,1x),1x,f10.6,1x,f8.3,1x)
5046  format(1x,f12.2,2x,3(e22.14,1x))

end subroutine out

subroutine engeps(zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,nt,delt,ts,  &
                     L,kxa,kya,kza,aj,bj,f,visch,viscz,ilap,L1,L2,L3,iueng,iueps,kh_g)

  implicit none
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif
  
  integer :: ikx,iky,ikz,ikza,ilap
  integer, dimension(iktx,ikty,iktzp) :: L
  integer :: nt,ierr,iueng,iueps
  complex, dimension(iktx,ikty,iktzp) ::  zx,zy,zz,ux,uy,uz,tt,ge,g1,g2
  complex :: der,zi
  real :: kxa(iktx),kya(ikty),kza(iktz)
  real :: aj,bj,ts,delt,e,kx,ky,kz,time,wk,wkh2,wkh2n,kzn
  real :: vh,vz,L1,L2,L3,twopi
  real :: visch,viscz,vzx,vzy,vzz,f,n,bv2,ffac,n2fac,sigma
  real :: vort2,div2,rossby,fr_z,fr_h,ke,pe,eg,ea,rossbyVel
  real :: pv,h,v,epsk,epsp,eps,tmp
  real :: zero_kz_geo,zero_kz_grv,zero_kh_grv,zero_kh_geo,kh_g
  external :: proj, velo

  integer :: mype
  common/mpi/mype

  twopi = 4.*asin(1.)
  zi    = cmplx(0.,1.)
  time  = ts + nt*delt
  bv2   = aj*bj
  n     = sqrt(bv2)
  ffac  = max(f,1.e-15)    !Coriolis parameter makes it very small if it is zero 
  n2fac = max(bv2,1.e-15)  !Brunt Vaisala makes it very small if it is zero
  
  call proj(zx,zy,zz,L,kxa,kya,kza)
  call velo(zx,zy,zz,ux,uy,uz,L,kxa,kya,kza)

  vort2  = 0.
  div2   = 0.
  Rossby = 0.
  Fr_z   = 0.
  Fr_h   = 0.
  zero_kz_geo = 0.
  zero_kz_grv = 0.
  zero_kh_geo = 0.
  zero_kh_grv = 0.
  ke = 0.
  pe = 0.
  eg = 0.
  ea = 0.
  pv = 0.
  h  = 0.
  v  = 0.
  epsk = 0.
  epsp = 0.
  eps  = 0.

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kzn = kz * (L3/twopi)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx = 1,iktx
           kx = kxa(ikx)
           wkh2 = kx*kx + ky*ky
           wk = kx*kx + ky*ky + kz*kz
           if (L(ikx,iky,ikz).eq.1) then
              wkh2n = wkh2 * (L1/twopi)**2
              
              vzx   = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
              ke    = ke + vzx/wk
              
              epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vzx/wk
              v     = v  + vzx
              vzy   = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
              ke    = ke + vzy/wk

              epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vzy/wk
              v     = v + vzy
              vzz   = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
              ke    = ke + vzz/wk

              epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vzz/wk
              v     = v  + vzz
              vort2 = vort2 + vzz

              vh    = real( zx(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              h     = h + vh
              vh    = real( zy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              h     = h + vh
              vh    = real( zz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )
              h     = h + vh
              vh    = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )
              pe    = pe  + vh
              epsp  = epsp + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vh*aj/bj

              vzx   = real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              vzy   = real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              vzz   = real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )

              if(wkh2n.lt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 zero_kh_grv = zero_kh_grv + vzx + vzy
                 zero_kh_geo = zero_kh_geo + vh*aj/bj
                 pv = pv + F*F*kz*kz*vh*aj/bj
              endif
              
              if(wkh2n.gt.1.e-10 .and. abs(kzn).lt.1.e-10) then
                 zero_kz_geo = zero_kz_geo + vzx + vzy
                 zero_kz_grv = zero_kz_grv + vzz + vh*aj/bj
                 pv = pv + aj*bj*wkh2*(vzx + vzy)
              endif
              
              if(wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 vz    = real(ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)))
                 eg    = eg  + vz/wkh2
                 sigma = wkh2*bv2 + kz*kz*F*F
                 pv    = pv + sigma*vz/wkh2
                 vz    = real(g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)))
                 ea    = ea  + vz/wkh2
                 vz    = real(g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)))
                 ea    = ea  + vz/wkh2
              endif
              
              rossby = rossby + real(zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)))/(ffac*ffac)
              der    =          real(zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)))
              der    = (der +   real(zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz))))
              fr_z   = fr_z + der/N2FAC
              der    = zi*(kx*ux(ikx,iky,ikz) + ky*uy(ikx,iky,ikz))
              div2   = div2 + real(der*conjg(der))
           endif ! L
        enddo
     enddo
  enddo

  v      = 2.*v
  epsk   = 2.*epsk
  epsp   = 2.*epsp
  eg     = eg + zero_kz_geo + zero_kh_geo
  ea     = ea + zero_kz_grv + zero_kh_grv
  if (aj.ne.0. .and. bj.ne.0.) pe = aj*pe/bj

#if ( MPI == 1 )
  call mpi_reduce(ke,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ke=tmp
  call mpi_reduce(pe,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  pe=tmp
  call mpi_reduce(eg,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  eg=tmp
  call mpi_reduce(ea,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ea=tmp
  call mpi_reduce(pv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  pv=tmp
  call mpi_reduce(rossby,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rossby=tmp
  call mpi_reduce(fr_z,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  fr_z=tmp
  call mpi_reduce(vort2,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  vort2=tmp
  call mpi_reduce(v,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  v=tmp
  call mpi_reduce(zero_kz_grv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  zero_kz_grv=tmp
  call mpi_reduce(zero_kz_geo,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  zero_kz_geo=tmp
  call mpi_reduce(zero_kh_grv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  zero_kh_grv=tmp
  call mpi_reduce(zero_kh_geo,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  zero_kh_geo=tmp
  call mpi_reduce(epsk,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  epsk=tmp
  call mpi_reduce(epsp,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  epsp=tmp
#endif

  if (mype.eq.0) then ! prep for output

     if (f.gt.1.e-8) then
        rossby = sqrt(2.*rossby)
     else
        rossby = - 999.
     endif

     rossbyVel = sqrt(eg*2.)/(twopi/kh_g)/ffac

     if (N.gt.1.e-8) then
        fr_z = sqrt(fr_z)
        fr_h = sqrt(2.*vort2)/N
     else
        fr_z = - 999.
        fr_h = - 999.
     endif
     eps    = epsk + epsp
     if (aj.ne.0. .and. bj.ne.0.) then
        e  = (pe + ke)
     else
        e = 99999999.
     endif

     write(iueps,5046) time, epsk, epsp, eps
     write(iueng,5045) time,ke,pe,e,eg,ea,pv,rossby,fr_z,fr_h,zero_kh_grv,zero_kh_geo,v,rossbyVel,kh_g
     call flush(iueps)
     call flush(iueng)
  endif
  return
5045  format(1x,f12.2,2x,12(e14.5,1x),1x,f10.6,1x,f8.3,1x)
5046  format(1x,f12.2,2x,3(e22.14,1x))

end subroutine engeps


subroutine nltrmsgrms(NGG,NGAu,NGAb,NAbAb,NAbAu,NAuAu,nt,delt,ts,L,itrnrms)

! derives the rms of of nonlinear terms in dG/dt
! nltrmsgrms = Non-Linear TeRMS in the evolution of Geo modes RMS

  implicit none
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif
  
  integer :: ikx,iky,ikz,ikza
  integer, dimension(iktx,ikty,iktzp) :: L
  integer :: nt,ierr,itrnrms
  complex, dimension(iktx,ikty,iktzp) ::  NGG,NGAu,NGAb,NAbAb,NAbAu,NAuAu
  real    :: rmsNGG,rmsNGAu,rmsNGAb,rmsNAbAb,rmsNAbAu,rmsNAuAu
  real    :: time,ts,delt,tmp
  

  integer :: mype
  common/mpi/mype

  time  = ts + nt*delt
  rmsNGG    = 0
  rmsNGAu   = 0
  rmsNGAb   = 0
  rmsNAbAb  = 0
  rmsNAbAu  = 0
  rmsNAuAu  = 0
 

  do ikz = 1,iktzp
     do iky = 1,iktz
        do ikx = 1,iktx
           if (L(ikx,iky,ikz).eq.1) then             
              rmsNGG   = rmsNGG  + real( NGG(ikx,iky,ikz)*conjg(NGG(ikx,iky,ikz)) )
              rmsNGAb  = rmsNGAb + real( NGAb(ikx,iky,ikz)*conjg(NGAb(ikx,iky,ikz)) )
              rmsNGAu  = rmsNGAu + real( NGAu(ikx,iky,ikz)*conjg(NGAu(ikx,iky,ikz)) )
              rmsNAbAb = rmsNAbAb+ real( NAbAb(ikx,iky,ikz)*conjg(NAbAb(ikx,iky,ikz)) )
              rmsNAbAu = rmsNAbAu+ real( NAbAu(ikx,iky,ikz)*conjg(NAbAu(ikx,iky,ikz)) )
              rmsNAuAu = rmsNAuAu+ real( NAuAu(ikx,iky,ikz)*conjg(NAuAu(ikx,iky,ikz)) )
           endif ! L
        enddo
     enddo
  enddo
 

#if ( MPI == 1 )
  call mpi_reduce(rmsNGG,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rmsNGG=sqrt(tmp)
  call mpi_reduce(rmsNGAb,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rmsNGAb=sqrt(tmp)
  call mpi_reduce(rmsNGAu,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rmsNGAu=sqrt(tmp)
  call mpi_reduce(rmsNAbAb,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rmsNAbAb=sqrt(tmp)
  call mpi_reduce(rmsNAbAu,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rmsNAbAu=sqrt(tmp)
  call mpi_reduce(rmsNAuAu,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rmsNAuAu=sqrt(tmp)
#endif

  if (mype.eq.0) then ! prep for output
     write(itrnrms,5065) time,rmsNGG,rmsNGAu,rmsNGAb,rmsNAbAb,rmsNAbAu,rmsNAuAu
     call flush(itrnrms)
  endif

  return
5065  format(1x,f12.2,2x,6(e14.5,1x))

end subroutine nltrmsgrms

subroutine vartimeseries(zz,tt,ge,g1,L,time,kxmode,kymode,kzmode,khmode,imyp,ikxmo,ikymo,ikzmo,iu)


  implicit none 
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  real     :: kxmode(nmodeh,nkhm), kymode(nmodeh,nkhm), kzmode(nkhm,nkzm) !picked wavenumbers
  real     :: khmode(nkhm)         ! the horizontal radious of each ring 
  integer  :: ikxmo(nmodeh,nkhm), ikymo(nmodeh,nkhm), ikzmo(nkhm,nkzm) !index of picked wavenumbers
  integer  :: imyp(nkhm,nkzm)           !the number of processor for each vertical level
  integer  :: itet,irad,ized
  real     :: time
  integer  :: j,jj,kk,L(iktx,ikty,iktzp),iu
  complex, dimension(iktx,ikty,iktzp) :: zz,tt,ge,g1
  real, dimension(nkhm*nkzm*nmodeh,8) :: vars
  real,    dimension(nkhm*nkzm*nmodeh,4) :: wavnums  
  
  integer :: mype
  common/mpi/mype

  do ized = 1, nkzm
     jj=0
     do irad = 1, nkhm
        if (imyp(irad,ized).eq.mype) then
           do itet = 1 , nmodeh
              jj = jj + 1
              j = (ized-1)*nmodeh*nkhm+jj
              wavnums(j,1)=khmode(irad)
              wavnums(j,2)=kxmode(itet,irad)
              wavnums(j,3)=kymode(itet,irad) 
              wavnums(j,4)=kzmode(irad,ized)
              vars(j,1)=real(zz(ikxmo(itet,irad),ikymo(itet,irad),ikzmo(irad,ized)))
              vars(j,2)=AIMAG(zz(ikxmo(itet,irad),ikymo(itet,irad),ikzmo(irad,ized)))
              vars(j,3)=real(tt(ikxmo(itet,irad),ikymo(itet,irad),ikzmo(irad,ized)))
              vars(j,4)=AIMAG(tt(ikxmo(itet,irad),ikymo(itet,irad),ikzmo(irad,ized)))
              vars(j,5)=real(ge(ikxmo(itet,irad),ikymo(itet,irad),ikzmo(irad,ized)))
              vars(j,6)=AIMAG(ge(ikxmo(itet,irad),ikymo(itet,irad),ikzmo(irad,ized)))
              vars(j,7)=real(g1(ikxmo(itet,irad),ikymo(itet,irad),ikzmo(irad,ized)))
              vars(j,8)=AIMAG(g1(ikxmo(itet,irad),ikymo(itet,irad),ikzmo(irad,ized)))
           enddo
        endif
     enddo     
  enddo

  write(iu,*),time
  if (mype.eq.0) then
     do j= 1 , nkhm*nkzm*nmodeh
           write(iu,5000) (wavnums(j,jj),jj=1,4),(vars(j,kk),kk=1,8)
     enddo
     call flush(iu)
  endif
  return
5000  format(1X,4(F6.1,3X),4X,8(E14.6,1x))
end subroutine vartimeseries


subroutine energy_full(zx,zy,zz,tt,ux,uy,uz,ge,g1,g2, &
                        L,kxa,kya,kza,aj,bj,f,ke,pe,eg,ea,L1,L2,L3)

! Computes total energy (KE, PE, VE, WE) 

  implicit none 
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: ikx,iky,ikz,ikza
  integer :: L(iktx,ikty,iktzp),ierr
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,ux,uy,uz,tt,ge,g1,g2
  complex :: zi
  real :: aj,bj,ke,pe,eg,ea,vh,L1,L2,L3
  real :: kx,ky,kz,wk,wkh2,wkh2n,kzn
  real :: vzx,vzy,vzz,vz,f,n,bf2,ffac,n2fac,twopi
  real :: kxa(iktx),kya(ikty),kza(iktz)
  real :: zero_kz_geo,zero_kz_grv,zero_kh_grv,zero_kh_geo,tmp
  external proj, velo

  integer :: mype
  common/mpi/mype

  zi   = cmplx(0.,1.)
  bf2  = aj*bj
  n    = sqrt(bf2)
  ffac = max(F,1.e-15)
  n2fac= max(bf2,1.e-15)
  twopi = 4.*asin(1.)   

  call velo(zx,zy,zz,ux,uy,uz,L,kxa,kya,kza)
  call wtoab(zx,zy,zz,tt,ge,g1,g2,ux,uy,uz,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)

  zero_kz_geo = 0.
  zero_kz_grv = 0.
  zero_kh_geo = 0.
  zero_kh_grv = 0.
  ke = 0.
  pe = 0.
  eg = 0.
  ea  = 0.

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kzn = kz*L3/twopi
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           wkh2 = kx*kx + ky*ky
           wkh2n = wkh2*(L1/twopi)**2
           wk = kx*kx + ky*ky + kz*kz

           if (L(ikx,iky,ikz).eq.1) then
              wk   = max(wk,  1.e-15)
              wkh2 = max(wkh2,1.e-15)
              
              vzx = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
              ke = ke + vzx/wk
              vzy = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
              ke = ke + vzy/wk
              vzz = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
              ke = ke + vzz/wk
              vh = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )
              pe = pe  + vh
            
              vzx=real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              vzy=real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              vzz=real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )
            
              if(wkh2n.lt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 zero_kh_grv = zero_kh_grv+vzx+vzy
                 zero_kh_geo = zero_kh_geo+VH*aj/bj
              endif
            
              if(wkh2n.gt.1.e-10 .and. abs(kzn).lt.1.e-10) then
                 zero_kz_geo = zero_kz_geo+VZX+VZY
                 zero_kz_grv = zero_kz_grv+VZZ+VH*aj/bj
              endif
              
              if(wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 vz = real(ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)))
                 eg = eg  + vz/wkh2
                 vz = real(g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)))
                 ea = ea  + vz/wkh2
                 vz = real(g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)))
                 ea = ea  + vz/wkh2
              endif

           endif
        enddo
     enddo
  enddo

  if (aj.ne.0. .and. bj.ne.0.) pe = aj*pe/bj
  eg = eg + zero_kz_geo + zero_kh_geo
  ea = ea + zero_kz_grv + zero_kh_grv

#if ( MPI == 1 )
  call mpi_allreduce(ke,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  ke=tmp
  call mpi_allreduce(pe,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  pe=tmp
  call mpi_allreduce(eg,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  eg=tmp
  call mpi_allreduce(ea,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  ea=tmp
#endif

end subroutine energy_full

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Baer Tribbia Balancing 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine BalancePart(zxok,zyok,zzok,ttok,zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,  &
            ur,vr,wr,zxnr,zynr,zznr,ttnr,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3,ilap,visch,viscz)

! Inside this subroutine z*ok = fields BEFORE balancing and z*nk = fields AFTER balancing
! which is different from the main subroutine (Do not get confused!)

   implicit none 
   include 'param.inc'

   complex, dimension(iktx,ikty,iktzp) :: zxok,zyok,zzok,ttok,zxnk,zynk,zznk,ttnk
   real,    dimension(n1d,n3d,n2dp)    :: zxor,zyor,zzor,ttor,zxnr,zynr,zznr,ttnr
   complex, dimension(iktx,ikty,iktzp) :: nzxk,nzyk,nzzk,nttk,uk,  vk,  wk
   real,    dimension(n1d,n3d,n2dp)    :: nzxr,nzyr,nzzr,nttr,ur,  vr,  wr
   complex, dimension(iktx,ikty,iktzp) :: rhtt,geok,gw1k,gw2k,NNN1,NNN2 
   complex, dimension(iktx,ikty,iktzp) :: MGGk,NGG1k,NGG2k,NAA1k,NAA2k,NN1k,NN2k,NdGdG1,NdGdG2
   
   complex :: Ord1gw1k,Ord1gw2k,Ord2gw1k,Ord2gw2k
   complex :: N10gw1k,N10gw2k,N01gw1k,N01gw2k
   complex :: zi
   integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp)
   integer :: ilap
   real    :: aj,bj,cor,L1,L2,L3,r1,r2,visch,viscz
   real    :: wkh,wkh2,kx,ky,kz,kh,khn,wk2,wkhn,k,freq
   real    :: kxa(iktx),kya(ikty),kza(iktz)
   real    :: twopi
   external :: proj,velo,fftwrk,fftwkr

   integer :: mype
   common/mpi/mype

   zi = cmplx(0.,1.)
   twopi    = 4.*asin(1.)

   zxnk(:,:,:) = zxok(:,:,:)
   zynk(:,:,:) = zyok(:,:,:)
   zznk(:,:,:) = zzok(:,:,:)  
   ttnk(:,:,:) = ttok(:,:,:)

! Baer-Tribbia

! ---------------- finding NGG
   call killwaves(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
   call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L,  &
        uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
   call wtoab(zxnk,zynk,zznk,ttnk,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
   call wtoab(nzxk,nzyk,nzzk,nttk,MGGk,NGG1k,NGG2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
! ---------------- finding NAA

   do ikz = 1,iktzp
      ikza = mype*iktzp+ikz
      kz = kza(ikza)
      do iky = 1,IKTY
         ky = kya(iky)
         do ikx=1,IKTX
            if ( L(ikx,iky,ikz).eq.1 ) then
               kx  = kxa(ikx)
               kh  = sqrt( kx*kx + ky*ky )
               khn = kh * L1/twopi
               wk2 = kx*kx + ky*ky + kz*kz
               k  = sqrt(wk2)
               freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k

#if (SPHTRUNC == 1)
               r1 = visch*wk2**ilap
               r2 = visch*wk2**ilap
#else
               r1 = visch*kh**(ilap*2) + viscz*kz**(ilap*2)
               r2 = visch*kh**(ilap*2) + viscz*kz**(ilap*2)
#endif
              
               gw1k(ikx,iky,ikz)=NGG1k(ikx,iky,ikz)/(-zi*freq+r1)
               gw2k(ikx,iky,ikz)=NGG2k(ikx,iky,ikz)/(+zi*freq+r1)
            endif ! Ly
         enddo
      enddo
   enddo
  
   geok=cmplx(0.,0.)

   call atowb(geok,gw1k,gw2k,zxnk,zynk,zznk,ttnk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
   call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L,  &
        uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
   call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NAA1k,NAA2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

! ------------ finding N=NGG+NGA+NAA
  call wtoab(zxok,zyok,zzok,ttok,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,IKTY
        ky = kya(iky)
        do ikx=1,IKTX
           if ( L(ikx,iky,ikz).eq.1 ) then
              kx  = kxa(ikx)
              kh  = sqrt( kx*kx + ky*ky )
              khn = kh * L1/twopi
              wk2 = kx*kx + ky*ky + kz*kz
              k  = sqrt(wk2)
              freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k
#if (SPHTRUNC == 1)
               r1 = visch*wk2**ilap
               r2 = visch*wk2**ilap
#else
               r1 = visch*kh**(ilap*2) + viscz*kz**(ilap*2)
               r2 = visch*kh**(ilap*2) + viscz*kz**(ilap*2)
#endif             
              gw1k(ikx,iky,ikz)=NGG1k(ikx,iky,ikz)/(-zi*freq+r1)
              gw2k(ikx,iky,ikz)=NGG2k(ikx,iky,ikz)/(+zi*freq+r1)
           endif ! Ly
        enddo
     enddo
  enddo
  
  call atowb(geok,gw1k,gw2k,zxnk,zynk,zznk,ttnk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
  call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L,  &
       uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
  call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NN1k,NN2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

! ---------------- Calculating N01

! ---------------- STEP#1 calculating NdGdG

  call wtoab(zxok,zyok,zzok,ttok,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,IKTY
        ky = kya(iky)
        do ikx=1,IKTX
           if ( L(ikx,iky,ikz).eq.1 ) then
              geok(ikx,iky,ikz)=MGGk(ikx,iky,ikz)
           endif ! Ly
        enddo
     enddo
  enddo
  
  gw1k=cmplx(0.,0.)
  gw2k=cmplx(0.,0.)
  
  call atowb(geok,gw1k,gw2k,zxnk,zynk,zznk,ttnk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
  call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L,  &
       uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
  call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NdGdG1,NdGdG2,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

! ---------------- STEP#2 calculating NNN=N(G+dG)(G+dG)

  call wtoab(zxok,zyok,zzok,ttok,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
  
  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,IKTY
        ky = kya(iky)
        do ikx=1,IKTX
           if ( L(ikx,iky,ikz).eq.1 ) then
              geok(ikx,iky,ikz)=MGGk(ikx,iky,ikz)+geok(ikx,iky,ikz)
           endif ! Ly
        enddo
     enddo
  enddo
  
  gw1k=cmplx(0.,0.)
  gw2k=cmplx(0.,0.)
  
  call atowb(geok,gw1k,gw2k,zxnk,zynk,zznk,ttnk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
  call constr (zxnk,zynk,zznk,ttnk,nzxk,nzyk,nzzk,nttk,L,  &
       uk,vk,wk,ur,vr,wr,zxnr,zynr,zznr,ttnr,nzxr,nzyr,nzzr,nttr,kxa,kya,kza)
  call wtoab(nzxk,nzyk,nzzk,nttk,rhtt,NNN1,NNN2,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)
  
  call wtoab(zxok,zyok,zzok,ttok,geok,gw1k,gw2k,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

  

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,IKTY
        ky = kya(iky)
        do ikx=1,IKTX
           if ( L(ikx,iky,ikz).eq.1 ) then
              kx  = kxa(ikx)
              kh  = sqrt( kx*kx + ky*ky )
              khn = kh * L1/twopi
              wk2 = kx*kx + ky*ky + kz*kz
              k  = sqrt(wk2)
              freq = sqrt(kz*kz*COR**2 + kh*kh*AJ*BJ)/k
#if (SPHTRUNC == 1)
               r1 = visch*wk2**ilap
               r2 = visch*wk2**ilap
#else
               r1 = visch*kh**(ilap*2) + viscz*kz**(ilap*2)
               r2 = visch*kh**(ilap*2) + viscz*kz**(ilap*2)
#endif
    
              Ord1gw1k=NGG1k(ikx,iky,ikz)/(-zi*freq+r1)
              Ord1gw2k=NGG2k(ikx,iky,ikz)/(+zi*freq+r1)
                  
              N10gw1k=NN1k(ikx,iky,ikz)-NGG1k(ikx,iky,ikz)-NAA1k(ikx,iky,ikz)
              N10gw2k=NN2k(ikx,iky,ikz)-NGG2k(ikx,iky,ikz)-NAA2k(ikx,iky,ikz)

              N01gw1k=NNN1(ikx,iky,ikz)-NdGdG1(ikx,iky,ikz)-NGG1k(ikx,iky,ikz)
              N01gw2k=NNN2(ikx,iky,ikz)-NdGdG2(ikx,iky,ikz)-NGG2k(ikx,iky,ikz)
                  
              Ord2gw1k=N10gw1k/(-zi*freq+r1)-(N01gw1k-2*r1*NGG1k(ikx,iky,ikz))/(-zi*freq+r1)**2
              Ord2gw2k=N10gw2k/(+zi*freq+r1)-(N01gw2k-2*r1*NGG1k(ikx,iky,ikz))/(+zi*freq+r1)**2
          
              gw1k(ikx,iky,ikz)=Ord1gw1k+Ord2gw1k
              gw2k(ikx,iky,ikz)=Ord1gw2k+Ord2gw2k
              
           endif ! Ly
        enddo
     enddo
  enddo
  
  call atowb(geok,gw1k,gw2k,zxnk,zynk,zznk,ttnk,uk,vk,wk,L,kxa,kya,kza,aj,bj,cor,L1,L2,L3)

end subroutine BalancePart

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! NORMAL MODE CONVERSION
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine wtoab(zx,zy,zz,tt,geok,gw1k,gw2k,u,v,w,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)

! Converts from (vorticity, buoyancy) to (geok,grav.wave_1k,grav.wave_2k)
                  
   implicit none 
   include 'param.inc'
   
   integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp),ikz0
   complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,u,v,w,geok,gw1k,gw2k
   complex :: div,zk,tk,dk,bk,zi
   real :: aj,bj,f,omega,bv2,bv,norm,f2,sqr2,L1,L2,L3,twopi
   real :: wkh,wkh2,kx,ky,kz,wk,wk2,wkhn
   real :: kxa(iktx),kya(ikty),kza(iktz)
   
   integer :: mype
   common/mpi/mype
   
   geok = cmplx(0.,0.)
   gw1k = cmplx(0.,0.)
   gw2k = cmplx(0.,0.)
   zi   = cmplx(0.,1.)
   bv2  = aj*bj
   f2   = f*f
   bv   = sqrt(bv2)
   sqr2 = sqrt(2.)
   twopi = 4.*asin(1.)

   call velo(zx,zy,zz,u,v,w,L,kxa,kya,kza)

   if (mype.eq.0) then 
      ikz0 = 2
   else
      ikz0 = 1
   endif

   do iky=1,ikty
      ky = kya(iky)
      do ikx=1,iktx
         kx = kxa(ikx)
         wkh2 = kx*kx + ky*ky
         wkh  = sqrt(wkh2)
         wkhn = wkh*L1/twopi
         if (wkhn.lt.1.e-10) then
            geok(ikx,iky,ikz0:iktzp) = aj*tt(ikx,iky,ikz0:iktzp)
            gw1k(ikx,iky,ikz0:iktzp) = u(ikx,iky,ikz0:iktzp)-zi*v(ikx,iky,ikz0:iktzp)
            gw2k(ikx,iky,ikz0:iktzp) = u(ikx,iky,ikz0:iktzp)+zi*v(ikx,iky,ikz0:iktzp)
         else
            do ikz=ikz0,iktzp
               if (L(ikx,iky,ikz).eq.1) then
                  ikza  = mype*iktzp+ikz
                  kz    = kza(ikza)
                  wk2   = wkh2 + kz*kz
                  wk    = sqrt(wk2)                
                  wk    = max(wk,1.e-15)
                  omega = sqrt(f2*kz*kz + bv2*wkh2)/wk 
                  div   = zi*(kx*u(ikx,iky,ikz)+ky*v(ikx,iky,ikz))
                  bk    = aj*tt(ikx,iky,ikz)
                  
                  zk = zz(ikx,iky,ikz)
                  dk = (wk/kz) * div
                  tk = (wkh/bv) * bk
                  
                  norm  = omega*wk
                  geok(ikx,iky,ikz) = (bv*wkh*zk + zi*f*kz*tk)/norm
                  
                  norm  = sqr2*omega*wk
                  gw1k(ikx,iky,ikz) = (-zi*f*kz*zk + omega*wk*dk - bv*wkh*tk)/norm
                  
                  norm  = sqr2*omega*wk
                  gw2k(ikx,iky,ikz) = (+zi*f*kz*zk + omega*wk*dk + bv*wkh*tk)/norm
               endif
            enddo
         endif
      enddo
   enddo

   if (mype.eq.0) then  ! do kz=0 modes
      do iky=1,ikty
         ky = kya(iky)
         do ikx=1,iktx
            kx = kxa(ikx)
            if (L(ikx,iky,1).eq.1) then
               bk=aj*tt(ikx,iky,1)
               geok(ikx,iky,1)=zz(ikx,iky,1)
               gw1k(ikx,iky,1)=w(ikx,iky,1)-zi*bk/bv
               gw2k(ikx,iky,1)=w(ikx,iky,1)+zi*bk/bv
            endif
         enddo
      enddo
   endif
   
   return
 end subroutine wtoab
 

 subroutine atowb(geok,gw1k,gw2k,zx,zy,zz,tt,u,v,w,L,kxa,kya,kza,aj,bj,f,L1,L2,L3)

! Converts from (geok,grav.wave_1k,grav.wave_2k) to (zeta,d,t).

   implicit none 
   include 'param.inc'

   integer :: ikx,iky,ikz,L(iktx,ikty,iktzp),ikz0,ikza
   complex, dimension(iktx,ikty,iktzp) ::zx,zy,zz,tt,geok,gw1k,gw2k,u,v,w
   complex :: zi,zk,dk,tk,bk,gn,psi,div
   complex :: zgk,z1k,z2k,dgk,d1k,d2k,tgk,t1k,t2k
   real :: aj,bj,f,omega,bv2,bv,sqr2,f2,L1,L2,L3,twopi
   real :: wkh,wkh2,kx,ky,kz,wk,wk2,wkhn
   real :: kxa(iktx),kya(ikty),kza(iktz)

   integer :: mype
   common/mpi/mype

   u  = cmplx(0.,0.)
   v  = cmplx(0.,0.)
   w  = cmplx(0.,0.)
   tt = cmplx(0.,0.)
   zx = cmplx(0.,0.)
   zy = cmplx(0.,0.)
   zz = cmplx(0.,0.)
   zi = cmplx(0.,1.)
   bv2 = aj*bj
   bv = sqrt(bv2)
   sqr2 = sqrt(2.)
   f2 = f*f
   twopi = 4.*asin(1.)

   if (mype.eq.0) then 
      ikz0 = 2
   else
      ikz0 = 1
   endif
   
   do ikz=ikz0,iktzp
      ikza = mype*iktzp+ikz
      kz   = kza(ikza)
      do iky=1,ikty
         ky = kya(iky)
         do ikx=1,iktx
            kx = kxa(ikx)
            wkh2 = kx*kx + ky*ky
            wkh  = sqrt(wkh2)
            wkhn = wkh*L1/twopi

            if (L(ikx,iky,ikz).eq.1) then 
               wk2 = wkh2 + kz*kz
               wk  = sqrt(wk2)
               wk  = max(wk,1.e-15)
               omega = sqrt(f2*kz*kz + bv2*wkh2)/wk
               
               gn  = geok(ikx,iky,ikz) / (omega*wk)
               zgk = bv * wkh       * gn
               dgk = cmplx(0.,0.)
               tgk = - zi * f * kz * gn
               
               gn  = gw1k(ikx,iky,ikz) / (sqr2*omega*wk)
               z1k = + zi  * f * kz * gn
               d1k = omega * wk     * gn
               t1k = - bv   * wkh    * gn

               gn  = gw2k(ikx,iky,ikz) / (sqr2*omega*wk)
               z2k = - zi  * f * kz * gn
               d2k = omega * wk     * gn
               t2k = + bv   * wkh   * gn
               
               zk = zgk + z1k + z2k
               dk = dgk + d1k + d2k
               tk = tgk + t1k + t2k
               
               if (wkhn.gt.1.e-10) then
                  div             = dk*kz/wk
                  u(ikx,iky,ikz)  = +zi*(ky*zk-kx*div)/wkh2             
                  v(ikx,iky,ikz)  = -zi*(kx*zk+ky*div)/wkh2             
                  w(ikx,iky,ikz)  = zi*div/kz
                  tt(ikx,iky,ikz) = bv*tk/(aj*wkh)
               else
                  u(ikx,iky,ikz) =     0.5*(gw1k(ikx,iky,ikz)+gw2k(ikx,iky,ikz))
                  v(ikx,iky,ikz) = -zi*0.5*(gw1k(ikx,iky,ikz)-gw2k(ikx,iky,ikz))
                  w(ikx,iky,ikz) = cmplx(0.,0.)
                 tt(ikx,iky,ikz) = geok(ikx,iky,ikz)/aj
              endif

           endif
        enddo
     enddo
  enddo

  if (mype.eq.0) then  ! do kz=0 modes
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           if (L(ikx,iky,1).eq.1) then
              psi = -geok(ikx,iky,1)/(kx*kx+ky*ky)
              u(ikx,iky,1) = -zi*ky*psi
              v(ikx,iky,1) = +zi*kx*psi
              w(ikx,iky,1) = 0.5*(gw1k(ikx,iky,1)+gw2k(ikx,iky,1))
              bk = zi*bv*0.5*(gw1k(ikx,iky,1)-gw2k(ikx,iky,1))
              tt(ikx,iky,1) = bk/aj
           endif
        enddo
     enddo
  endif
     
  call vort(u,v,w,zx,zy,zz,L,kxa,kya,kza)
  
  return
end subroutine atowb


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! VELOCITY <--> VORTICITY
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine velo(zx,zy,zz,u,v,w,L,kxa,kya,kza)

 ! Calculates k-space velocity from k-space vorticity.
 ! curl (vorticity) = - laplacian (velocity) if velocity is solenoidal.

  implicit none
  include 'param.inc'

  integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp)
  real :: kx,ky,kz,k2,kxa(iktx),kya(ikty),kza(iktz)
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,u,v,w
  complex :: c1,c2,c3,i

  integer :: mype
  common/mpi/mype

  i = cmplx(0.,1.)
  u = cmplx(0.,0.)
  v = cmplx(0.,0.)
  w = cmplx(0.,0.)

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           if (L(ikx,iky,ikz).eq.1) then
              kx = kxa(IKX)
              k2 = max(kx*kx+ky*ky+kz*kz,1.e-15)
              c1 = ky*zz(ikx,iky,ikz) - kz*zy(ikx,iky,ikz)               
              c2 = kz*zx(ikx,iky,ikz) - kx*zz(ikx,iky,ikz)
              c3 = kx*zy(ikx,iky,ikz) - ky*zx(ikx,iky,ikz)               
              u(ikx,iky,ikz) = i * c1 / k2
              v(ikx,iky,ikz) = i * c2 / k2
              w(ikx,iky,ikz) = i * c3 / k2
           endif
        enddo
     enddo
  enddo

  return
end subroutine velo


subroutine vort(u,v,w,zx,zy,zz,L,kxa,kya,kza)

! Calculates k-space vortcity from k-space velocity.

  implicit none
  include 'param.inc'

  integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp)
  real ::  kx,ky,kz,kxa(iktx),kya(ikty),kza(iktz)
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,u,v,w
  complex :: c1,c2,c3,i

  integer :: mype
  common/mpi/mype

  i = cmplx(0.,1.)

  zx = cmplx(0.,0.)
  zy = cmplx(0.,0.)
  zz = cmplx(0.,0.)

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           if (L(ikx,iky,ikz).eq.1) then
              c1 = ky*w(ikx,iky,ikz) - kz*v(ikx,iky,ikz)               
              c2 = kz*u(ikx,iky,ikz) - kx*w(ikx,iky,ikz)
              c3 = kx*v(ikx,iky,ikz) - ky*u(ikx,iky,ikz)               
              zx(ikx,iky,ikz) = i * c1
              zy(ikx,iky,ikz) = i * c2
              zz(ikx,iky,ikz) = i * c3
           endif
        enddo
     enddo
  enddo
  
  return
end subroutine vort


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c NONLINEAR TERMS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine convol(zxk,zyk,zzk,ttk,nxk,nyk,nzk,ntk,L,uk,vk,wk,ur,vr,wr, &
                  wxk,wyk,wzk,wtk,zxr,zyr,zzr,ttr,nxr,nyr,nzr,ntr,kxa,kya,kza)

! Calculates convolution sums, calls ffts, etc.

  implicit none 
  include 'param.inc'

  integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp)
  real    :: kxa(iktx),kya(ikty),kza(iktz),kx,ky,kz
  real,    dimension(n1d,n3d,n2dp)       :: ur,vr,wr,nxr,nyr,nzr,ntr,zxr,zyr,zzr,ttr
  complex, dimension(iktx,ikty,iktzp) :: nxk,nyk,nzk,ntk,zxk,zyk,zzk,ttk
  complex, dimension(iktx,ikty,iktzp) :: wxk,wyk,wzk,wtk,uk,vk,wk
  complex :: c1,c2,c3,zi
  external :: fftwrk,fftwkr,velo

  integer :: mype
  common/mpi/mype

  zi  = cmplx(0.,1.)
  ntk = cmplx(0.,0.)
  wxk = zxk
  wyk = zyk
  wzk = zzk
  wtk = ttk

! Nonlinear term in the temperature equation.

  call velo(zxk,zyk,zzk,uk,vk,wk,L,kxa,kya,kza)
  call fftwkr(ur,uk)
  call fftwkr(vr,vk)
  call fftwkr(wr,wk)
  call fftwkr(ttr,ttk)

  nxr = ur*ttr
  nyr = vr*ttr
  nzr = wr*ttr

  call fftwrk(nxr,nxk)
  call fftwrk(nyr,nyk)
  call fftwrk(nzr,nzk)
  
  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           ntk(ikx,iky,ikz) = - zi * ( kx*nxk(ikx,iky,ikz) &
                                     + ky*nyk(ikx,iky,ikz) &
                                     + kz*nzk(ikx,iky,ikz) )*L(ikx,iky,ikz)
        enddo
     enddo
  enddo

! Nonlinear term in the vorticity equation.

  call fftwkr(zxr,zxk)
  call fftwkr(zyr,zyk)
  call fftwkr(zzr,zzk)

  nxr = wr*zyr - vr*zzr
  nyr = ur*zzr - wr*zxr
  nzr = vr*zxr - ur*zyr

  call fftwrk(nxr,nxk)
  call fftwrk(nyr,nyk)
  call fftwrk(nzr,nzk)

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           c1 = ky*nzk(ikx,iky,ikz) - kz*nyk(ikx,iky,ikz)               
           c2 = kz*nxk(ikx,iky,ikz) - kx*nzk(ikx,iky,ikz)
           c3 = kx*nyk(ikx,iky,ikz) - ky*nxk(ikx,iky,ikz)               
           nxk(ikx,iky,ikz) = - zi*c1*L(ikx,iky,ikz)
           nyk(ikx,iky,ikz) = - zi*c2*L(ikx,iky,ikz)
           nzk(ikx,iky,ikz) = - zi*c3*L(ikx,iky,ikz)
        enddo
     enddo
  enddo

  zxk = wxk
  zyk = wyk
  zzk = wzk
  ttk = wtk
  
  return
end subroutine convol


subroutine constr(zxk,zyk,zzk,ttk,nxk,nyk,nzk,ntk,L,uk,vk,wk,ur,vr,wr, &
                  zxr,zyr,zzr,ttr,nxr,nyr,nzr,ntr,kxa,kya,kza)

! Calculates convolution sums, calls ffts, etc. without scratch arrays

  implicit none 
  include 'param.inc'

  integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp)

  real :: kxa(iktx),kya(ikty),kza(iktz),kx,ky,kz
  real,    dimension(n1d,n3d,n2dp)       :: ur,vr,wr,zxr,zyr,zzr,ttr,nxr,nyr,nzr,ntr
  complex, dimension(iktx,ikty,iktzp) :: uk,vk,wk,zxk,zyk,zzk,ttk,nxk,nyk,nzk,ntk
  complex :: c1,c2,c3,zi

  integer :: mype
  common/mpi/mype

  external :: fftwrk,fftwkr,velo
  
  zi  = cmplx(0.,1.)
  ntr = cmplx(0.,0.)

! Nonlinear term in temperature equation.

  call velo(zxk,zyk,zzk,uk,vk,wk,L,kxa,kya,kza)
  call fftwkr(ur,uk)
  call fftwkr(vr,vk)
  call fftwkr(wr,wk)
  call fftwkr(ttr,ttk)

  ur = ur*ttr
  vr = vr*ttr
  wr = wr*ttr

! these calls not in convol:
  call fftwrk(ur,uk)
  call fftwrk(vr,vk)
  call fftwrk(wr,wk)
  call fftwrk(ttr,ttk)

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           ntk(ikx,iky,ikz) = - ZI * ( kx*uk(ikx,iky,ikz)  &
                                     + ky*vk(ikx,iky,ikz)  &
                                     + kz*wk(ikx,iky,ikz)  ) * L(ikx,iky,ikz)
        enddo
     enddo
  enddo
  
 
! Nonlinear term in vorticity equation.

  call velo(zxk,zyk,zzk,uk,vk,wk,L,kxa,kya,kza)
  call fftwkr(ur,uk)
  call fftwkr(vr,vk)
  call fftwkr(wr,wk)
  call fftwkr(zxr,zxk)
  call fftwkr(zyr,zyk)
  call fftwkr(zzr,zzk)

  nxr = wr*zyr - vr*zzr
  nyr = ur*zzr - wr*zxr
  nzr = vr*zxr - ur*zyr
  
  call fftwrk(nxr,nxk)
  call fftwrk(nyr,nyk)
  call fftwrk(nzr,nzk)
  call fftwrk(zxr,zxk)
  call fftwrk(zyr,zyk)
  call fftwrk(zzr,zzk)

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           c1 = ky*nzk(ikx,iky,ikz) - kz*nyk(ikx,iky,ikz)               
           c2 = kz*nxk(ikx,iky,ikz) - kx*nzk(ikx,iky,ikz)
           c3 = kx*nyk(ikx,iky,ikz) - ky*nxk(ikx,iky,ikz)               
           nxk(ikx,iky,ikz) = - zi * c1 * L(ikx,iky,ikz)
           nyk(ikx,iky,ikz) = - zi * c2 * L(ikx,iky,ikz)
           nzk(ikx,iky,ikz) = - zi * c3 * L(ikx,iky,ikz)
        enddo
     enddo
  enddo

  return
end subroutine constr


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c FFTS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine fftwrk(zr,zk)

  implicit none
  include 'param.inc'

  complex, dimension(iktx,ikty,iktzp) :: zk
  real,    dimension(n1d,n3d,n2dp)    :: zr
  real :: norm
  integer :: i,j,k,ikstart
  integer, parameter :: n1n3=n1d*n3d
  integer, parameter :: iktxz=iktx*iktz,iktxy=iktx*ikty,iktyp=ikty/npe
  integer, parameter :: iktxyzp=iktx*ikty*iktzp

  complex, dimension(iktx,iktz,iktyp) :: zkt
  complex, dimension(iktx,ikty,iktzp) :: zk1,zk2 
  common/scratch/zkt,zk1,zk2

  integer*8 :: plan2_rk,plan2_kr,plan1_rk,plan1_kr
  common/FTW2D/plan2_rk,plan2_kr
  common/FTW1D/plan1_rk,plan1_kr

  integer :: mype
  common/mpi/mype

  norm = float(n1*n2*n3)

  ! Do 2D (x,z) transforms at n2 levels all at once 
  ! Note, output data has dimensions iktx,iktz,iktyp 
  ! zk1 is scratch
  call rfftwnd_f77_real_to_complex(plan2_rk,n2p,zr,1,n1n3,zk1,1,iktxz)

  ! 2D transforms are in-place so output is in zk
  ! but has dimensions iktx,iktz,iktyp
  ! copy to zkt in prep for transpose
!! zkt=zk ! this doesn't work on sharcnet!  Do this instead:
   do i=1,iktxyzp
     zkt(i,1,1)=zk(i,1,1)
   enddo

  ! Transpose zkt(iktx,iktz,iktyp) -> zk(iktx,ikty,iktzp)
  ! Transpose output from previous step to have dimensions iktx,ikty,iktzp
#if (MPI == 1 ) 
  call mpitranspose(zkt,iktx,iktz,iktyp,zk,ikty,iktzp,npe,zk1,zk2)
#else
  call serialtranspose(zkt,zk,iktx,iktz,ikty)
#endif

  ! Do remaining 1D (y) transforms at iktx*iktz rows
  do k=1,iktzp
     ikstart=1+(k-1)*iktxy
     call fftw_f77(plan1_rk,iktx,zk(ikstart,1,1),iktx,1,zk1,iktx,1)
  enddo
  
  ! Normalize
  zk=zk/norm
end subroutine fftwrk


subroutine fftwkr(zr,zk)

  implicit none
  include 'param.inc'

  complex, dimension(iktx,ikty,iktzp) :: zk
  real,    dimension(n1d,n3d,n2dp)    :: zr
  integer :: i,k,ikstart
  integer, parameter :: n1n3=n1d*n3d
  integer, parameter :: iktxz=iktx*iktz,iktxy=iktx*ikty,iktyp=ikty/npe
  integer, parameter :: iktxyzp=iktx*ikty*iktzp

  complex, dimension(iktx,iktz,iktyp) :: zkt
  complex, dimension(iktx,ikty,iktzp) :: zk1,zk2 
  common/scratch/zkt,zk1,zk2

  integer*8 :: plan2_rk,plan2_kr,plan1_rk,plan1_kr
  common/FTW2D/plan2_rk,plan2_kr
  common/FTW1D/plan1_rk,plan1_kr

  integer :: mype
  common/mpi/mype

  call realit(zk)

  ! Do 1D (y) transforms at iktx*iktz rows
  do k=1,iktzp
     ikstart=1+(k-1)*iktxy
     call fftw_f77(plan1_kr,iktx,zk(ikstart,1,1),iktx,1,zk1,iktx,1)
  enddo

  ! Transpose zk(iktx,ikty,iktzp) -> zkt(iktx,iktz,iktyp)
  ! Note, output of transpose has dimensions iktx,iktz,iktyp but store it in zk
#if (MPI == 1)
  call mpitranspose(zk,iktx,ikty,iktzp,zkt,iktz,iktyp,npe,zk1,zk2)
#else
  call serialtranspose(zk,zkt,iktx,ikty,iktz)
#endif

  ! Copy transposed array back into zk in preparation for in-place 2D transforms
!! zk=zkt ! this doesnt work on sharcnet
  do i=1,iktxyzp
     zk(i,1,1)=zkt(i,1,1)
  enddo

  ! Finally do 2D (x,z) transforms at n2 levels all at once
  ! Note, input data in zk has dimensions iktx,iktz,iktyp
  ! zk1 is scrtatch
  call rfftwnd_f77_complex_to_real(plan2_kr,n2p,zk,1,iktxz,zk1,1,n1n3)

end subroutine fftwkr


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c RANDOM NUMBERS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function rang(i)

! If i ne 0 then initializes RANNO with seed = i
! If i eq 0 then draws a random GAUSSIAN number with 
! mean and std = 1

  implicit none
  integer :: i 
  real :: v1,v2,R,FAC,twopi,ranno,rang
  external :: ranno

  twopi = 4.*asin(1.)
  
  if (i.ne.0) then
     v1 = ranno(i)
  else
200  v1 = 2.*(ranno(0)+twopi/2.)/twopi -1.
     v2 = 2.*(ranno(0)+twopi/2.)/twopi -1.
     r = v1**2. + v2**2.
     if (r.gt.1.) goto 200
     fac = sqrt( -2.*log(r)/r)
     rang = v1*fac
  endif
  return
end function RANG


function ranno (i)

! Controls random number generator.
!-----------------------------------------
! - If argument i.ne.0 it performs initialization with i=seed no.
! - If argument i.eq.0 it draws a random no.
!-----------------------------------------
  implicit none
  integer :: i,junk,ihold
  real :: ranno,twopi,ran1
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
  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ntab=32,ndiv=1+(im-1)/ntab
  real :: ran1
  real, parameter :: am=1./im,eps=1.2e-7,rnmx=1.-eps
  integer :: idum,j,k,iv(32),iy,junk
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



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c PHYSICAL SPACE DIAGNOSTICS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dumprealHCVC(zzk,zzr,nRS,noutput,time,iuH,iuV,nskipH,nskipV,nHmax,nVmax,justHor)


! dump real space fields into HC.out VC.out
! if  justHor = 1 just keep the horizontal slice otherwise keep both horizontal and vertical slices 
  implicit none
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  real, dimension(n1d,n3d,n2dp) :: zzr,junk
  complex, dimension(iktx,ikty,iktzp) :: zzk
  integer :: iproc,nbuf,nsends,istatus,istart,noutput
  integer :: ikx,iky,ikz,ikza,i,j,k,nRS,justHor
  real    :: time
  integer :: horcutindex,vercutindex,jj,iuH,iuV,nskipH,nskipV,nHmax,nVmax

  real, dimension(n1,n2d) :: horvars
  external :: fftwrk,fftwkr

#if ( MPI == 1 )
  integer :: status(MPI_STATUS_SIZE)
#endif
  integer :: mype
  common/mpi/mype

  ! Move vorticity and velocity to physical space.
  call fftwkr(zzr,zzk)
#if ( MPI == 1 )
  nbuf = n1d*n3d*n2dp
  nsends = npe-1
#endif
 
  nRS=nRS+1
  if (justHor.ne.1) vercutindex = 1
  horcutindex = 1
  
  
  if (mype.eq.0) then
     if (nRS.eq.1) then
        if (justHor.ne.1) then
           WRITE(iuV,'(I8)') noutput
           WRITE(iuV,'(I8,2X,I8)') nHmax, nVmax
           WRITE(iuV,'(I8,2X,I8)') nskipH, nskipV
        endif
        WRITE(iuH,'(I8)') noutput
        WRITE(iuH,'(I8,2X,I8)') nHmax, nHmax
        WRITE(iuH,'(I8,2X,I8)') nskipH, nskipH
     endif

     if (justHor.ne.1) then
        WRITE(iuV,'(f12.2)') TIME
        do j=1,nVmax,nskipV
           WRITE(iuV,'(<n3>(E14.6,1x))') (ZZR(I,J,vercutindex),I=1,nHmax,nskipH)
        enddo
     endif
  endif

#if ( MPI == 1 )
  if (mype.gt.0) then
     call mpi_send(zzr,nbuf,MPI_REAL,0,137,MPI_COMM_WORLD,istatus)
  endif
#endif


  if (mype.eq.0) then
     jj=jj+1
     do k=1,n2dp
        horvars(1:n1,k) = zzr(1:n1,horcutindex,k)
     enddo

    
#if ( MPI == 1 )
     do iproc=1,nsends
        call mpi_recv(junk,nbuf,MPI_REAL,MPI_ANY_SOURCE,137,MPI_COMM_WORLD,status,istatus)
        do k=1,n2dp
           istart=status(MPI_SOURCE)*n2p
           horvars(1:n1,istart+k)=junk(1:n1,horcutindex,k)
        enddo  
     enddo
#endif

     WRITE(iuH,'(f12.2)') TIME
     do j=1,nHmax,nskipH
        WRITE(iuH,'(<n3>(E14.6,1x))') (horvars(I,J),I=1,nHmax,nskipH)
     enddo

  endif

! Put everything back to fourier for future timestepping
  call fftwrk(zzr,zzk)
  return 

end subroutine dumprealHCVC



subroutine bindumprst(zx,zy,zz,tt,junk,ntdump,time)

! dump restart into
  
  implicit none
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  real :: time
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  complex, dimension(iktx,ikty,iktzp) :: junk
  integer :: iproc,nbuf,ntdump
  integer ikx, iky, ikz, istatus


#if ( MPI == 1 )
  integer :: status(MPI_STATUS_SIZE)
#endif
  integer :: mype
  common/mpi/mype

#if ( MPI == 1 )
  nbuf = iktx*ikty*iktzp
#endif




#if ( MPI == 1 )  
  if (mype.gt.0) then
     ! send to mype=0; use netcdf var ids as tags
     call mpi_send(zx,nbuf,MPI_COMPLEX,0,1201,MPI_COMM_WORLD,istatus)
     call mpi_send(zy,nbuf,MPI_COMPLEX,0,1202,MPI_COMM_WORLD,istatus)
     call mpi_send(zz,nbuf,MPI_COMPLEX,0,1203,MPI_COMM_WORLD,istatus)
     call mpi_send(tt,nbuf,MPI_COMPLEX,0,1204,MPI_COMM_WORLD,istatus)
  endif
#endif

  if (mype.eq.0) then

     open (90,file='ZkX.out',   form='unformatted')
     open (91,file='ZkY.out',   form='unformatted')
     open (92,file='ZkZ.out',   form='unformatted')
     open (93,file='ZkT.out',   form='unformatted')

     write(90) iktx
     write(90) ikty
     write(90) iktz
     write(90) time

     write(90) (((zx(ikx,iky,ikz),ikx=1,iktx),iky=1,ikty),ikz=1,iktzp)             
     write(91) (((zy(ikx,iky,ikz),ikx=1,iktx),iky=1,ikty),ikz=1,iktzp)
     write(92) (((zz(ikx,iky,ikz),ikx=1,iktx),iky=1,ikty),ikz=1,iktzp)
     write(93) (((tt(ikx,iky,ikz),ikx=1,iktx),iky=1,ikty),ikz=1,iktzp)
     
        
#if ( MPI == 1 )
     do iproc=1,npe-1
        call mpi_recv(junk,nbuf,MPI_COMPLEX,iproc,1201,MPI_COMM_WORLD,status,istatus)
        write(90) (((junk(ikx,iky,ikz),ikx=1,iktx),iky=1,ikty),ikz=1,iktzp)
        call mpi_recv(junk,nbuf,MPI_COMPLEX,iproc,1202,MPI_COMM_WORLD,status,istatus)
        write(91) (((junk(ikx,iky,ikz),ikx=1,iktx),iky=1,ikty),ikz=1,iktzp)
        call mpi_recv(junk,nbuf,MPI_COMPLEX,iproc,1203,MPI_COMM_WORLD,status,istatus)
        write(92) (((junk(ikx,iky,ikz),ikx=1,iktx),iky=1,ikty),ikz=1,iktzp)
        call mpi_recv(junk,nbuf,MPI_COMPLEX,iproc,1204,MPI_COMM_WORLD,status,istatus)
        write(93) (((junk(ikx,iky,ikz),ikx=1,iktx),iky=1,ikty),ikz=1,iktzp)
     enddo
#endif
  endif

 
  if (mype.eq.0) then
     close (90)
     close (91)
     close (92)
     close (93)
     print*,'T = ',time
     print*,'Just wrote on restart file.'
  endif


end subroutine bindumprst

subroutine binreadrst(zx,zy,zz,tt,irest,time)

  implicit none
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  real :: time, time1
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,junkk
  integer irest,jw,kx1,ky1,kz1,ntime,nbuf,istatus,iproc
  integer ikx, iky, ikz

#if ( MPI == 1 )
  integer :: status(MPI_STATUS_SIZE)
#endif
  integer :: mype
  common/mpi/mype



  nbuf = iktx*ikty*iktzp
  
  if (mype.eq.0) then
     print*,'  '
     print*, 'Restarting from output data Zk.in '
     open(81,file='ZkX.in',   form='unformatted')
     open(82,file='ZkY.in',   form='unformatted')
     open(83,file='ZkZ.in',   form='unformatted')
     open(84,file='ZkT.in',   form='unformatted')
!     jw = 1
!     kx1 = iktx
!     ky1 = ikty
!     kz1 = iktz
!     time1 = 5.0
!    read(81) jw
     read(81) kx1
     read(81) ky1
     read(81) kz1
     read(81) time1
     jw = 1

     print*, '  IKTX,Y,Z: ',kx1,ky1,kz1
     
     if (jw.lt.IREST) then
        print*, 'Files not long enough.'
        stop
     endif
     
     if (kx1.ne.iktx .or. ky1.ne.ikty .or. kz1.ne.iktz) then
        print*,'The resolution does not match!'
        stop  
     else
        print*,'Resolution unchanged. '
     endif


#if ( MPI == 1 ) 
     read(81) (((zx(ikx,iky,ikz),ikx=1,kx1),iky=1,ky1),ikz=1,iktzp)
     read(82) (((zy(ikx,iky,ikz),ikx=1,kx1),iky=1,ky1),ikz=1,iktzp)
     read(83) (((zz(ikx,iky,ikz),ikx=1,kx1),iky=1,ky1),ikz=1,iktzp)
     read(84) (((tt(ikx,iky,ikz),ikx=1,kx1),iky=1,ky1),ikz=1,iktzp)
     do iproc = 1, npe-1
        read(81) (((junkk(ikx,iky,ikz),ikx=1,kx1),iky=1,ky1),ikz=1,iktzp)
        call mpi_send(junkk,nbuf,MPI_COMPLEX,iproc,120,MPI_COMM_WORLD,istatus)
        read(82) (((junkk(ikx,iky,ikz),ikx=1,kx1),iky=1,ky1),ikz=1,iktzp)
        call mpi_send(junkk,nbuf,MPI_COMPLEX,iproc,220,MPI_COMM_WORLD,istatus)
        read(83) (((junkk(ikx,iky,ikz),ikx=1,kx1),iky=1,ky1),ikz=1,iktzp)
        call mpi_send(junkk,nbuf,MPI_COMPLEX,iproc,320,MPI_COMM_WORLD,istatus)
        read(84) (((junkk(ikx,iky,ikz),ikx=1,kx1),iky=1,ky1),ikz=1,iktzp)
        call mpi_send(junkk,nbuf,MPI_COMPLEX,iproc,420,MPI_COMM_WORLD,istatus)
        if (istatus.ne.0) print*,'Error sending Zk.in'
     enddo
#endif
  endif !mype=0

#if ( MPI == 1 )     
     if (mype.gt.0) then
        ! send to mype=0; use netcdf var ids as tags
        call mpi_recv(zx,nbuf,MPI_COMPLEX,0,120,MPI_COMM_WORLD,status,istatus)
        call mpi_recv(zy,nbuf,MPI_COMPLEX,0,220,MPI_COMM_WORLD,status,istatus)
        call mpi_recv(zz,nbuf,MPI_COMPLEX,0,320,MPI_COMM_WORLD,status,istatus)
        call mpi_recv(tt,nbuf,MPI_COMPLEX,0,420,MPI_COMM_WORLD,status,istatus)
     endif
#endif

     if (mype.eq.0) then  
        print*, 'Just read time= ', time1
        print*, 'Integration begins at t= ', time1
        time = time1
        close (81)
        close (82)
        close (83)
        close (84)
     endif

end subroutine binreadrst



subroutine dumpreal(zxk,zyk,zzk,ttk,zxr,zyr,zzr,ttr,uk,vk,wk,ur,vr,wr,junk, &
                    L,kxa,kya,kza,ntdump,time)

  implicit none
  include 'param.inc'

  integer :: ikx,iky,ikz,ikza,i,j,k,L(iktx,ikty,iktzp),ntdump
  real :: kxa(iktx),kya(ikty),kza(iktz),kx,ky,kz,time
  real,    dimension(n1d,n3d,n2dp)       :: ur,vr,wr,zxr,zyr,zzr,ttr,junk
  complex, dimension(iktx,ikty,iktzp) :: zxk,zyk,zzk,ttk,uk,vk,wk
  external :: fftwrk,fftwkr,velo

  integer :: mype
  common/mpi/mype

  call velo(zxk,zyk,zzk,uk,vk,wk,L,kxa,kya,kza)

! Move vorticity and velocity to physical space.
  call fftwkr(zxr,zxk)
  call fftwkr(zyr,zyk)
  call fftwkr(zzr,zzk)
  call fftwkr(ttr,ttk)
  call fftwkr( ur, uk)
  call fftwkr( vr, vk)
  call fftwkr( wr, wk)

#if (NETCDF == 1)
  call ncdumprsp(ur,vr,wr,zxr,zyr,zzr,ttr,junk,ntdump,time)
#endif

! Put everything back to fourier for future timestepping
  call fftwrk(zxr,zxk)
  call fftwrk(zyr,zyk)
  call fftwrk(zzr,zzk)
  call fftwrk(ttr,ttk)

  uk = cmplx(0.,0.)
  vk = cmplx(0.,0.)
  wk = cmplx(0.,0.)

  return
end subroutine dumpreal


#if (NETCDF == 1)

subroutine ncprep(nrsp,nrst,nspwv)
  
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: istatus

!! real space output vars
  integer :: nrsp,idnc,rs1slic
  integer :: iputu,iputv,iputw,iputth,iputzx,iputzy,iputzz
  integer :: idtimes,idu,idv,idw,idth,idzx,idzy,idzz
  integer :: idx,idy,idz,idt
  integer, dimension(4) :: ncstart,nccount,ncdims
  common/netcdfrsp/idnc,ncstart,nccount, &
         iputu,iputv,iputw,iputth,iputzx,iputzy,iputzz,&
         idtimes,idu,idv,idw,idth,idzx,idzy,idzz


!! restart ourput vars
  integer :: nrst
  integer :: idzxk,idzyk,idzzk,idttk
  integer, dimension(5) :: ncstartk,nccountk,ncdimstk,ncdimsxk,ncdimsyk,ncdimszk
  integer :: idncxk
  integer :: idxkx,idxky,idxkz,idxkri,idxkt,idxtimesk
  integer :: idncyk
  integer :: idykx,idyky,idykz,idykri,idykt,idytimesk
  integer :: idnczk
  integer :: idzkx,idzky,idzkz,idzkri,idzkt,idztimesk
  integer :: idnctk
  integer :: idtkx,idtky,idtkz,idtkri,idtkt,idttimesk
  common/netcdfrst/idncxk,idncyk,idnczk,idnctk,ncstartk,nccountk, &
       idxtimesk,idytimesk,idztimesk,idttimesk,idzxk,idzyk,idzzk,idttk

!! 3D spectrum of waves (and potentially other variables)
  integer :: nspwv
  integer :: idnsp,idkxsp,idkysp,idkzsp,idtsp
  integer, dimension(4) :: ncdimsp,ncstartsp,nccountsp
  integer, dimension(3) :: ncdimkk,ncstartkk,nccountkk
  integer :: idtimesp,idkkxx,idkkyy,idkkzz,idwvp,idwvm
  
  common/netcdfwvsp/idnsp,idkxsp,idkysp,idkzsp,idtsp,ncstartsp,nccountsp, &
       ncstartkk,nccountkk,idtimesp,idkkxx,idkkyy,idkkzz,idwvp,idwvm

  integer :: mype
  common/mpi/mype
  !!! first, prep physical space output file out.ncf

  if (nrsp.gt.0) then

  if (mype.eq.0) then

     print*,'Creating netcdf file for output'
     istatus = nf_create("out.ncf",IOR(NF_NOCLOBBER,NF_64BIT_OFFSET),idnc)
     if (istatus.ne.0) print*,'error in nf_create'
     if (n2rs.gt.n2dp) then
        if (mod(n2rs,n2dp).ne.0) then
           print*, 'n2rs should be either less than n2dp or integer*n2dp'
           stop
        endif
     endif

     if (rs1slic.eq.1) then
        istatus = nf_def_dim(idnc,"X",n1rs,idx)
        if (istatus.ne.0) print*,'error in nf_def_dim x'
        istatus = nf_def_dim(idnc,"Y",n2rs,idy)
        if (istatus.ne.0) print*,'error in nf_def_dim y'
        istatus = nf_def_dim(idnc,"Z",n3rs,idz)
        if (istatus.ne.0) print*,'error in nf_def_dim z'
     else
        istatus = nf_def_dim(idnc,"X",n1,idx)
        if (istatus.ne.0) print*,'error in nf_def_dim x'
        istatus = nf_def_dim(idnc,"Y",n2,idy)
        if (istatus.ne.0) print*,'error in nf_def_dim y'
        istatus = nf_def_dim(idnc,"Z",n3,idz)
        if (istatus.ne.0) print*,'error in nf_def_dim z'
     endif     
     istatus = nf_def_dim(idnc,"T",nrsp,idt)
     if (istatus.ne.0) print*,'error in nf_def_dim t'

     ncdims(1) = idx
     ncdims(2) = idz
     ncdims(3) = idy
     ncdims(4) = idt

     istatus = nf_def_var(idnc,"TIMES",NF_FLOAT,1,idt,idtimes)
     if (istatus.ne.0) print*,'error nf_def_var TIMES'

     if (iputu.eq.1) then
        istatus = nf_def_var(idnc,"U",NF_FLOAT,4,ncdims,idu)
        if (istatus.ne.0) print*,'error nf_def_var U'
     endif
     if (iputv.eq.1) then
        istatus = nf_def_var(idnc,"V",NF_FLOAT,4,ncdims,idv)
        if (istatus.ne.0) print*,'error nf_def_var V'
     endif
     if (iputw.eq.1) then
        istatus = nf_def_var(idnc,"W",NF_FLOAT,4,ncdims,idw)
        if (istatus.ne.0) print*,'error nf_def_var W'
     endif
     if (iputth.eq.1) then
        istatus = nf_def_var(idnc,"TH",NF_FLOAT,4,ncdims,idth)
        if (istatus.ne.0) print*,'error nf_def_var TH'
     endif
     if (iputzx.eq.1) then
        istatus = nf_def_var(idnc,"ZX",NF_FLOAT,4,ncdims,idzx)
        if (istatus.ne.0) print*,'error nf_def_var ZX'
     endif
     if (iputzy.eq.1) then
        istatus = nf_def_var(idnc,"ZY",NF_FLOAT,4,ncdims,idzy)
        if (istatus.ne.0) print*,'error nf_def_var ZY'
     endif
     if (iputzz.eq.1) then
        istatus = nf_def_var(idnc,"ZZ",NF_FLOAT,4,ncdims,idzz)
        if (istatus.ne.0) print*,'error nf_def_var ZZ'
     endif
     istatus = nf_enddef(idnc)
     if (istatus.ne.0) print*,'error enddef'

     ncstart(1) = 1
     ncstart(2) = 1

     if (rs1slic.eq.1) then
        nccount(1) = n1rs
        nccount(2) = n3rs
        if (n2rs.le.n2dp) then
           nccount(3) = n2rs
        else
           nccount(3) = n2dp
        endif
     else
        nccount(1) = n1
        nccount(3) = n2p
        nccount(2) = n3
     end if
     
     
     nccount(4) = 1

  endif ! mype

#if ( MPI == 1) 
  ! broadcast variable IDs to all procs
  call mpi_bcast(idu, 1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus) 
  call mpi_bcast(idv, 1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idw, 1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idth,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idzx,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idzy,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idzz,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
#endif

  endif ! nrsp


  !!! next, prep restart file

  if (mype.eq.0) then

     print*,'Creating netcdf restart file'
     istatus = nf_create("ZkX.out.ncf",ior(nf_clobber,nf_64bit_offset),idncxk)
     if (istatus.ne.0) print*,'error nf_create ZkX.out.ncf'
     istatus = nf_def_dim(idncxk,"KX",iktx,idxkx)
     if (istatus.ne.0) print*,'error nf_def_dim kx'
     istatus = nf_def_dim(idncxk,"KY",ikty,idxky)
     if (istatus.ne.0) print*,'error nf_def_dim ky'
     istatus = nf_def_dim(idncxk,"KZ",iktz,idxkz)
     if (istatus.ne.0) print*,'error nf_def_dim kz'
     istatus = nf_def_dim(idncxk,"RI",2,idxkri)
     if (istatus.ne.0) print*,'error nf_def_dim ri'
     istatus = nf_def_dim(idncxk,"T", nrst,idxkt)
     if (istatus.ne.0) print*,'error nf_def_dim t'

     ncdimsxk(1) = idxkx
     ncdimsxk(2) = idxky
     ncdimsxk(3) = idxkz
     ncdimsxk(4) = idxkri  ! dim 4 is real/imag part
     ncdimsxk(5) = idxkt

     istatus = nf_def_var(idncxk,"TIMES",NF_FLOAT,1,idxkt,idxtimesk)
     if (istatus.ne.0) print*,'error nf_def_var TIMES'

     istatus = nf_def_var(idncxk,"ZXK",NF_FLOAT,5,ncdimsxk,idzxk)
     if (istatus.ne.0) print*,'error nf_def_var ZX'
    
     istatus = nf_enddef(idncxk)
     if (istatus.ne.0) print*,'enddef error'

    

     istatus = nf_create("ZkY.out.ncf",ior(nf_clobber,nf_64bit_offset),idncyk)
     if (istatus.ne.0) print*,'error nf_create ZkY.out.ncf' 
     istatus = nf_def_dim(idncyk,"KX",iktx,idykx)
     if (istatus.ne.0) print*,'error nf_def_dim kx'
     istatus = nf_def_dim(idncyk,"KY",ikty,idyky)
     if (istatus.ne.0) print*,'error nf_def_dim ky'
     istatus = nf_def_dim(idncyk,"KZ",iktz,idykz)
     if (istatus.ne.0) print*,'error nf_def_dim kz'
     istatus = nf_def_dim(idncyk,"RI",2,idykri)
     if (istatus.ne.0) print*,'error nf_def_dim ri'
     istatus = nf_def_dim(idncyk,"T", nrst,idykt)
     if (istatus.ne.0) print*,'error nf_def_dim t'
     
     ncdimsyk(1) = idykx
     ncdimsyk(2) = idyky
     ncdimsyk(3) = idykz
     ncdimsyk(4) = idykri  ! dim 4 is real/imag part
     ncdimsyk(5) = idykt

     istatus = nf_def_var(idncyk,"TIMES",NF_FLOAT,1,idykt,idytimesk)
     if (istatus.ne.0) print*,'error nf_def_var TIMES'

   
     istatus = nf_def_var(idncyk,"ZYK",NF_FLOAT,5,ncdimsyk,idzyk)
     if (istatus.ne.0) print*, 'error nf_def_var ZY'
     
     istatus = nf_enddef(idncyk)
     if (istatus.ne.0) print*,'enddef error'

  
     istatus = nf_create("ZkZ.out.ncf",ior(nf_clobber,nf_64bit_offset),idnczk)
     if (istatus.ne.0) print*,'error nf_create ZkZ.out.ncf'
     istatus = nf_def_dim(idnczk,"KX",iktx,idzkx)
     if (istatus.ne.0) print*,'error nf_def_dim kx'
     istatus = nf_def_dim(idnczk,"KY",ikty,idzky)
     if (istatus.ne.0) print*,'error nf_def_dim ky'
     istatus = nf_def_dim(idnczk,"KZ",iktz,idzkz)
     if (istatus.ne.0) print*,'error nf_def_dim kz'
     istatus = nf_def_dim(idnczk,"RI",2,idzkri)
     if (istatus.ne.0) print*,'error nf_def_dim ri'
     istatus = nf_def_dim(idnczk,"T", nrst,idzkt)
     if (istatus.ne.0) print*,'error nf_def_dim t'

     ncdimszk(1) = idzkx
     ncdimszk(2) = idzky
     ncdimszk(3) = idzkz
     ncdimszk(4) = idzkri  ! dim 4 is real/imag part
     ncdimszk(5) = idzkt

     istatus = nf_def_var(idnczk,"TIMES",NF_FLOAT,1,idzkt,idztimesk)
     if (istatus.ne.0) print*,'error nf_def_var TIMES'

    
     istatus = nf_def_var(idnczk,"ZZK",NF_FLOAT,5,ncdimszk,idzzk)
     if (istatus.ne.0) print*,'error nf_def_var ZZ'
  
     istatus = nf_enddef(idnczk)
     if (istatus.ne.0) print*,'enddef error'


     istatus = nf_create("ZkT.out.ncf",ior(nf_clobber,nf_64bit_offset),idnctk)
     if (istatus.ne.0) print*,'error nf_create'
     istatus = nf_def_dim(idnctk,"KX",iktx,idtkx)
     if (istatus.ne.0) print*,'error nf_def_dim kx'
     istatus = nf_def_dim(idnctk,"KY",ikty,idtky)
     if (istatus.ne.0) print*,'error nf_def_dim ky'
     istatus = nf_def_dim(idnctk,"KZ",iktz,idtkz)
     if (istatus.ne.0) print*,'error nf_def_dim kz'
     istatus = nf_def_dim(idnctk,"RI",2,idtkri)
     if (istatus.ne.0) print*,'error nf_def_dim ri'
     istatus = nf_def_dim(idnctk,"T", nrst,idtkt)
     if (istatus.ne.0) print*,'error nf_def_dim t'

     ncdimstk(1) = idtkx
     ncdimstk(2) = idtky
     ncdimstk(3) = idtkz
     ncdimstk(4) = idtkri  ! dim 4 is real/imag part
     ncdimstk(5) = idtkt

     istatus = nf_def_var(idnctk,"TIMES",NF_FLOAT,1,idtkt,idttimesk)
     if (istatus.ne.0) print*,'error nf_def_var TIMES'

     istatus = nf_def_var(idnctk,"TTK",NF_FLOAT,5,ncdimstk,idttk)
     if (istatus.ne.0) print*,'error nf_def_var TT'

     istatus = nf_enddef(idnctk)
     if (istatus.ne.0) print*,'enddef error'

     ncstartk(1) = 1
     ncstartk(2) = 1

     nccountk(1) = iktx
     nccountk(2) = ikty
     nccountk(3) = iktzp
     nccountk(4) = 1
     nccountk(5) = 1

  endif ! mype



#if ( MPI == 1 )
  ! broadcast variable IDs to all procs
  call mpi_bcast(idzxk,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idzyk,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idzzk,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idttk,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
#endif
  

  !! Creating netcdf file for 3D spectrum of waves
  if (mype.eq.0) then
     
     print*,'Creating netcdf 3D spectrum file'
     istatus = nf_create("WVspc3D.out.ncf",ior(nf_clobber,nf_64bit_offset),idnsp)
     if (istatus.ne.0) print*,'error nf_create WVspc3D.out.ncf'
     istatus = nf_def_dim(idnsp,"KX",iktxsp,idkxsp)
     if (istatus.ne.0) print*,'error nf_def_dim kx in waves spectrum'
     istatus = nf_def_dim(idnsp,"KY",iktysp-1,idkysp)
     if (istatus.ne.0) print*,'error nf_def_dim ky in waves spectrum'
     istatus = nf_def_dim(idnsp,"KZ",iktzsp-1,idkzsp)
     if (istatus.ne.0) print*,'error nf_def_dim kz in waves spectrum'
     istatus = nf_def_dim(idnsp,"T", nspwv,idtsp)
     if (istatus.ne.0) print*,'error nf_def_dim t in waves spectrum'

     ncdimsp(1) = idkxsp
     ncdimsp(2) = idkysp
     ncdimsp(3) = idkzsp
     ncdimsp(4) = idtsp
     ncdimkk(1) = idkxsp
     ncdimkk(2) = idkysp
     ncdimkk(3) = idkzsp

     istatus = nf_def_var(idnsp,"TIMES",NF_FLOAT,1,idtsp,idtimesp)
     if (istatus.ne.0) print*,'error nf_def_var TIMES in waves spectrum'

     istatus = nf_def_var(idnsp,"KKXX",NF_FLOAT,3,ncdimkk,idkkxx)
     if (istatus.ne.0) print*,'error nf_def_var KKXX in waves spectrum'

     istatus = nf_def_var(idnsp,"KKYY",NF_FLOAT,3,ncdimkk,idkkyy)
     if (istatus.ne.0) print*,'error nf_def_var KKYY in waves spectrum'

     istatus = nf_def_var(idnsp,"KKZZ",NF_FLOAT,3,ncdimkk,idkkzz)
     if (istatus.ne.0) print*,'error nf_def_var KKZZ in waves spectrum'

     istatus = nf_def_var(idnsp,"EWVp",NF_FLOAT,4,ncdimsp,idwvp)
     if (istatus.ne.0) print*,'error nf_def_var WVp'

     istatus = nf_def_var(idnsp,"EWVm",NF_FLOAT,4,ncdimsp,idwvm)
     if (istatus.ne.0) print*,'error nf_def_var WVm'
    
     istatus = nf_enddef(idnsp)
     if (istatus.ne.0) print*,'enddef error'

!!$     ncstartsp(1) = 1
!!$     ncstartsp(2) = 1
!!$     ncstartkk(1) = 1
!!$     ncstartkk(2) = 1
!!$
!!$     nccountsp(1) = iktxsp
!!$     nccountsp(2) = iktysp
!!$     nccountk(3) = iktzp
!!$     nccountk(4) = 
!!$     nccountk(5) = 1

  endif ! mype

#if ( MPI == 1 )
  ! broadcast variable IDs to all procs
  call mpi_bcast(idkkxx,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idkkyy,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idkkzz,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idwvp,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
  call mpi_bcast(idwvm,1,MPI_INTEGER,0,MPI_COMM_WORLD,istatus)
#endif
  
  
end subroutine ncprep

subroutine ncdumprsp(ur,vr,wr,zxr,zyr,zzr,ttr,junk,ntdump,time,rs1slic)

! dump real space fields into out.ncf
  
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  real :: time
  real, dimension(n1d,n3d,n2dp) :: ur,vr,wr,zxr,zyr,zzr,ttr,junk
  integer :: iproc,nbuf,nsends,ntdump,rs1slic,n2end
  integer :: idnc,istatus,ip,idvar
  integer :: iputu,iputv,iputw,iputth,iputzx,iputzy,iputzz
  integer :: idtimes,idu,idv,idw,idth,idzx,idzy,idzz
  integer, dimension(4) :: ncstart,nccount
  common/netcdfrsp/idnc,ncstart,nccount, &
       iputu,iputv,iputw,iputth,iputzx,iputzy,iputzz,  &
       idtimes,idu,idv,idw,idth,idzx,idzy,idzz

#if ( MPI == 1 )
  integer :: status(MPI_STATUS_SIZE)
#endif
  integer :: mype
  common/mpi/mype

#if ( MPI == 1 )
  nbuf       = n1d*n3d*n2dp
  nsends     = (npe-1)*(iputu+iputv+iputw+iputth+iputzx+iputzy+iputzz)
#endif

  ncstart(4) = ntdump

  if (mype.eq.0) then
!     print*,'Writing to out.ncf at time ',time
     istatus    = nf_put_vara_real(idnc,idtimes,ntdump,1,time)
  endif

  
  if ((rs1slic.eq.0).or.(n2rs.gt.n2dp)) then
#if ( MPI == 1 )
     if (mype.gt.0) then
        ! send to mype=0; use netcdf var ids as tags
        if (iputu.eq.1)  call mpi_send(ur, nbuf,MPI_REAL,0,idu, MPI_COMM_WORLD,istatus)
        if (iputv.eq.1)  call mpi_send(vr, nbuf,MPI_REAL,0,idv, MPI_COMM_WORLD,istatus)
        if (iputw.eq.1)  call mpi_send(wr, nbuf,MPI_REAL,0,idw, MPI_COMM_WORLD,istatus)
        if (iputth.eq.1) call mpi_send(ttr,nbuf,MPI_REAL,0,idth,MPI_COMM_WORLD,istatus)
        if (iputzx.eq.1) call mpi_send(zxr,nbuf,MPI_REAL,0,idzx,MPI_COMM_WORLD,istatus)
        if (iputzy.eq.1) call mpi_send(zyr,nbuf,MPI_REAL,0,idzy,MPI_COMM_WORLD,istatus)
        if (iputzz.eq.1) call mpi_send(zzr,nbuf,MPI_REAL,0,idzz,MPI_COMM_WORLD,istatus)
     endif
#endif
  endif


  if (mype.eq.0) then
     if (rs1slic.eq.1) then 
        ncstart(3) = 1
        if (n2rs.le.n2dp) then
           n2end = n2rs
        else
           n2end = n2dp
        endif

        if (iputu.eq.1)  istatus = nf_put_vara_real(idnc,idu, ncstart,nccount, ur(1:n1rs,1:n3rs,1:n2end))
        if (iputv.eq.1)  istatus = nf_put_vara_real(idnc,idv, ncstart,nccount, vr(1:n1rs,1:n3rs,1:n2end))
        if (iputw.eq.1)  istatus = nf_put_vara_real(idnc,idw, ncstart,nccount, wr(1:n1rs,1:n3rs,1:n2end))
        if (iputth.eq.1) istatus = nf_put_vara_real(idnc,idth,ncstart,nccount,ttr(1:n1rs,1:n3rs,1:n2end))
        if (iputzx.eq.1) istatus = nf_put_vara_real(idnc,idzx,ncstart,nccount,zxr(1:n1rs,1:n3rs,1:n2end))
        if (iputzy.eq.1) istatus = nf_put_vara_real(idnc,idzy,ncstart,nccount,zyr(1:n1rs,1:n3rs,1:n2end))
        if (iputzz.eq.1) istatus = nf_put_vara_real(idnc,idzz,ncstart,nccount,zzr(1:n1rs,1:n3rs,1:n2end))
     else
        ncstart(3) = 1
        if (iputu.eq.1)  istatus = nf_put_vara_real(idnc,idu, ncstart,nccount, ur(1:n1,:,:))
        if (iputv.eq.1)  istatus = nf_put_vara_real(idnc,idv, ncstart,nccount, vr(1:n1,:,:))
        if (iputw.eq.1)  istatus = nf_put_vara_real(idnc,idw, ncstart,nccount, wr(1:n1,:,:))
        if (iputth.eq.1) istatus = nf_put_vara_real(idnc,idth,ncstart,nccount,ttr(1:n1,:,:))
        if (iputzx.eq.1) istatus = nf_put_vara_real(idnc,idzx,ncstart,nccount,zxr(1:n1,:,:))
        if (iputzy.eq.1) istatus = nf_put_vara_real(idnc,idzy,ncstart,nccount,zyr(1:n1,:,:))
        if (iputzz.eq.1) istatus = nf_put_vara_real(idnc,idzz,ncstart,nccount,zzr(1:n1,:,:))
     endif

#if ( MPI == 1 )
     if ((rs1slic.eq.0).or.(n2rs.gt.n2dp)) then
        do iproc=1,nsends
           call mpi_recv(junk,nbuf,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
           ncstart(3) = status(MPI_SOURCE)*n2p+1
           idvar = status(MPI_TAG) 
           if ((rs1slic.eq.1).and.(n2rs.gt.n2dp)) then
              if (ncstart(3).lt.n2rs) then
                 istatus = nf_put_vara_real(idnc,idvar,ncstart,nccount,junk(1:n1rs,1:n3rs,:))
              endif
           else
              istatus = nf_put_vara_real(idnc,idvar,ncstart,nccount,junk(1:n1,:,:))
           endif
        enddo
     endif
#endif 
  endif

  
!  if (mype.eq.0) then
!     istatus = nf_close(idnc)
!     if (istatus.ne.0) print*,'Error closing out.ncf'
!  endif
  

end subroutine ncdumprsp

subroutine ncdumprst(zx,zy,zz,tt,junk,ntdump,time)

! dump restart into out.ncf
  
  implicit none
  include 'param.inc'
  include 'netcdf.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  real :: time
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  complex, dimension(iktx,ikty,iktzp) :: junk
  integer :: iproc,nbuf,nsends,ntdump
  integer :: idncxk,idncyk,idnczk,idnctk,istatus,ip,idvar
  integer :: idzxk,idzyk,idzzk,idttk
  integer :: idxtimesk,idytimesk,idztimesk,idttimesk

  integer, dimension(5) :: ncstartk,nccountk
  common/netcdfrst/idncxk,idncyk,idnczk,idnctk,ncstartk,nccountk, &
       idxtimesk,idytimesk,idztimesk,idttimesk,idzxk,idzyk,idzzk,idttk



#if ( MPI == 1 )
  integer :: status(MPI_STATUS_SIZE)
#endif
  integer :: mype
  common/mpi/mype

#if ( MPI == 1 )
  nbuf = iktx*ikty*iktzp
  nsends = (npe-1)*4
#endif

  ncstartk(5) = ntdump

  if (mype.eq.0) then
     print*,'Writing to ZkX.out at time ',time
     istatus    = nf_put_vara_real(idncxk,idxtimesk,ntdump,1,time)
     print*,'Writing to ZkY.out at time ',time
     istatus    = nf_put_vara_real(idncyk,idytimesk,ntdump,1,time)
     print*,'Writing to ZkZ.out at time ',time
     istatus    = nf_put_vara_real(idnczk,idztimesk,ntdump,1,time)
     print*,'Writing to ZkT.out at time ',time
     istatus    = nf_put_vara_real(idnctk,idttimesk,ntdump,1,time)
  endif

#if ( MPI == 1 )  
  if (mype.gt.0) then
     ! send to mype=0; use netcdf var ids as tags
     call mpi_send(zx,nbuf,MPI_COMPLEX,0,1201,MPI_COMM_WORLD,istatus)
     call mpi_send(zy,nbuf,MPI_COMPLEX,0,1202,MPI_COMM_WORLD,istatus)
     call mpi_send(zz,nbuf,MPI_COMPLEX,0,1203,MPI_COMM_WORLD,istatus)
     call mpi_send(tt,nbuf,MPI_COMPLEX,0,1204,MPI_COMM_WORLD,istatus)
  endif
#endif

  if (mype.eq.0) then
     ncstartk(3) = 1
     ncstartk(5) = 1
     ncstartk(4) = 1 ! (real part)
     istatus = nf_put_vara_real(idncxk,idzxk,ncstartk,nccountk,real(zx))
     istatus = nf_put_vara_real(idncyk,idzyk,ncstartk,nccountk,real(zy))
     istatus = nf_put_vara_real(idnczk,idzzk,ncstartk,nccountk,real(zz))
     istatus = nf_put_vara_real(idnctk,idttk,ncstartk,nccountk,real(tt))
     ncstartk(4) = 2 ! (imag part)
     istatus = nf_put_vara_real(idncxk,idzxk,ncstartk,nccountk,aimag(zx))
     istatus = nf_put_vara_real(idncyk,idzyk,ncstartk,nccountk,aimag(zy))
     istatus = nf_put_vara_real(idnczk,idzzk,ncstartk,nccountk,aimag(zz))
     istatus = nf_put_vara_real(idnctk,idttk,ncstartk,nccountk,aimag(tt))

        
#if ( MPI == 1 )
     do iproc=1,nsends
        call mpi_recv(junk,nbuf,MPI_COMPLEX,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
        idvar = status(MPI_TAG)
        if (idvar.eq.1201) then
           ncstartk(3) = status(MPI_SOURCE)*iktzp+1
           ncstartk(4) = 1 ! (real part)
           istatus = nf_put_vara_real(idncxk,idzxk,ncstartk,nccountk,real(junk))
           if (istatus.ne.0) print*, 'Error dumpingZxReal'
           ncstartk(4) = 2 ! (imag part)
           istatus = nf_put_vara_real(idncxk,idzxk,ncstartk,nccountk,aimag(junk))
        endif
        if (idvar.eq.1202) then
           ncstartk(3) = status(MPI_SOURCE)*iktzp+1
           ncstartk(4) = 1 ! (real part)
           istatus = nf_put_vara_real(idncyk,idzyk,ncstartk,nccountk,real(junk))
           ncstartk(4) = 2 ! (imag part)
           istatus = nf_put_vara_real(idncyk,idzyk,ncstartk,nccountk,aimag(junk))
        endif
        if (idvar.eq.1203) then
           ncstartk(3) = status(MPI_SOURCE)*iktzp+1
           ncstartk(4) = 1 ! (real part)
           istatus = nf_put_vara_real(idnczk,idzzk,ncstartk,nccountk,real(junk))
           ncstartk(4) = 2 ! (imag part)
           istatus = nf_put_vara_real(idnczk,idzzk,ncstartk,nccountk,aimag(junk))
        endif
        if (idvar.eq.1204) then
           ncstartk(3) = status(MPI_SOURCE)*iktzp+1
           ncstartk(4) = 1 ! (real part)
           istatus = nf_put_vara_real(idnctk,idttk,ncstartk,nccountk,real(junk))
           ncstartk(4) = 2 ! (imag part)
           istatus = nf_put_vara_real(idnctk,idttk,ncstartk,nccountk,aimag(junk))
        endif
     enddo
#endif
  endif

  if (mype.eq.0) then
     istatus = nf_close(idncxk)
     if (istatus.ne.0) print*,'Error closing ZkX.out.ncf'
     istatus = nf_close(idncyk)
     if (istatus.ne.0) print*,'Error closing ZkY.out.ncf'
     istatus = nf_close(idnczk)
     if (istatus.ne.0) print*,'Error closing ZkZ.out.ncf'
     istatus = nf_close(idnctk)
     if (istatus.ne.0) print*,'Error closing ZkT.out.ncf'
  endif


end subroutine ncdumprst

subroutine ncreadrst(zx,zy,zz,tt,wr,wi,irest,time)

  implicit none
  include 'param.inc'
  include 'netcdf.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  real :: time
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  real, dimension(iktx,ikty,iktzp) :: wr,wi
  integer iktx1,ikty1,iktz1,nwrites,irest,iread,nbuf
  complex zi

  integer istatus,idkx,idky,idkz,idkri,idkt
  integer idncx,idncy,idncz,idnct
  integer idxtimesk,idytimesk,idztimesk,idttimesk,idzxk,idzyk,idzzk,idttk
  integer, dimension(5) :: ncstartrk,ncstartik,nccountk

#if ( MPI == 1 )
  integer :: status(MPI_STATUS_SIZE)
#endif
  integer :: mype
  common/mpi/mype

  zi  = cmplx(0.,1.)
  nbuf = iktx*ikty*iktzp

  if (mype.eq.0) then

     print*, 'Restarting from output data ZkX.in '
     istatus = nf_open('ZkX.in.ncf',NF_NOWRITE,idncx)
     if (istatus.ne.0) print*,'ZkX.in bad idnc'
     print*, 'Restarting from output data ZkY.in '
     istatus = nf_open('ZkY.in.ncf',NF_NOWRITE,idncy)
     if (istatus.ne.0) print*,'ZkY.in bad idnc'
     print*, 'Restarting from output data ZkZ.in '
     istatus = nf_open('ZkZ.in.ncf',NF_NOWRITE,idncz)
     if (istatus.ne.0) print*,'ZkZ.in bad idnc'
     print*, 'Restarting from output data ZkT.in '
     istatus = nf_open('ZkT.in.ncf',NF_NOWRITE,idnct)
     if (istatus.ne.0) print*,'ZkT.in bad idnc'
     

     ! get dimension IDs
     istatus = nf_inq_dimid(idncy,'KX',idkx)
     if (istatus.ne.0) print*,'Error getting kx ID'
     istatus = nf_inq_dimid(idncy,'KY',idky)
     if (istatus.ne.0) print*,'Error getting ky ID'
     istatus = nf_inq_dimid(idncy,'KZ',idkz)
     if (istatus.ne.0) print*,'Error getting kz ID'
     istatus = nf_inq_dimid(idncy,'RI',idkri)
     if (istatus.ne.0) print*,'Error getting real/imag ID'
     istatus = nf_inq_dimid(idncy,'T', idkt)
     if (istatus.ne.0) print*,'Error getting time ID'

     ! get dimension sizes
     istatus = nf_inq_dimlen(idncy,idkx,iktx1)  
     if (istatus.ne.0) print*,'Error getting idkx1'
     istatus = nf_inq_dimlen(idncy,idky,ikty1)  
     if (istatus.ne.0) print*,'Error getting idky1'
     istatus = nf_inq_dimlen(idncy,idkz,iktz1)  
     if (istatus.ne.0) print*,'Error getting idkz1'
     istatus = nf_inq_dimlen(idncy,idkt,nwrites)
     if (istatus.ne.0) print*,'Error getting nwrites'
     if (nwrites.lt.irest) then
        print*,'Not enough outputs in Zk.in'
        stop
     endif

     if (iktx1.ne.iktx .or. ikty1.ne.ikty .or. iktz1.ne.iktz) then
        print*,'Sorry, do not know how to change resolution.  Use restartres.F90'
        stop
     endif

     ! get variable IDs
     istatus = nf_inq_varid(idncy,'TIMES',idytimesk)
     if (istatus.ne.0) print*,'Error getting times ID'
     istatus = nf_inq_varid(idncx,'ZXK',idzxk)
     if (istatus.ne.0) print*,'Error getting zx ID'
     istatus = nf_inq_varid(idncy,'ZYK',idzyk)
     if (istatus.ne.0) print*,'Error getting zy ID'
     istatus = nf_inq_varid(idncz,'ZZK',idzzk)
     if (istatus.ne.0) print*,'Error getting zz ID'
     istatus = nf_inq_varid(idnct,'TTK',idttk)
     if (istatus.ne.0) print*,'Error getting tt ID'
     
     ! prep netcdf read
     ncstartrk(1) = 1
     ncstartik(1) = 1
     ncstartrk(2) = 1
     ncstartik(2) = 1
     ncstartrk(4) = 1 ! for real part
     ncstartik(4) = 2 ! for imaginary part
     ncstartrk(5) = irest
     ncstartik(5) = irest

     nccountk(1) = iktx
     nccountk(2) = ikty
     nccountk(3) = iktzp
     nccountk(4) = 1
     nccountk(5) = 1

     ! read in time
     istatus = nf_get_vara_real(idncy,idytimesk,irest,1,time)
     if (istatus.ne.0) print*,'Error reading time'

     do iread=npe-1,0,-1  ! read one block at a time; read first block last

        ncstartrk(3) = iread*iktzp+1 
        ncstartik(3) = iread*iktzp+1 

        istatus = nf_get_vara_real(idncx,idzxk,ncstartrk,nccountk,wr)
        if (istatus.ne.0) print*,'Error reading real(zx)'
        istatus = nf_get_vara_real(idncx,idzxk,ncstartik,nccountk,wi)
        if (istatus.ne.0) print*,'Error reading imag(zx)'
        zx = wr + zi*wi

        istatus = nf_get_vara_real(idncy,idzyk,ncstartrk,nccountk,wr)
        if (istatus.ne.0) print*,'Error reading real(zy)'
        istatus = nf_get_vara_real(idncy,idzyk,ncstartik,nccountk,wi)
        if (istatus.ne.0) print*,'Error reading imag(zy)'
        zy = wr + zi*wi
        
        istatus = nf_get_vara_real(idncz,idzzk,ncstartrk,nccountk,wr)
        if (istatus.ne.0) print*,'Error reading real(zz)'
        istatus = nf_get_vara_real(idncz,idzzk,ncstartik,nccountk,wi)
        if (istatus.ne.0) print*,'Error reading imag(zz)'
        zz = wr + zi*wi
        
        istatus = nf_get_vara_real(idnct,idttk,ncstartrk,nccountk,wr)
        if (istatus.ne.0) print*,'Error reading real(tt)'
        istatus = nf_get_vara_real(idnct,idttk,ncstartik,nccountk,wi)
        if (istatus.ne.0) print*,'Error reading imag(tt)'
        tt = wr + zi*wi

#if ( MPI == 1 )        
        if (iread.gt.0) then  ! send to processor iread
           call mpi_send(zx,nbuf,MPI_COMPLEX,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending zx to proc', iread
           call mpi_send(zy,nbuf,MPI_COMPLEX,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending zy to proc', iread
           call mpi_send(zz,nbuf,MPI_COMPLEX,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending zz to proc', iread
           call mpi_send(tt,nbuf,MPI_COMPLEX,iread,0,MPI_COMM_WORLD,istatus)
           if (istatus.ne.0) print*,'Error sending tt to proc', iread
        endif
#endif
     enddo
  endif

#if ( MPI == 1 )
  if (mype.gt.0) then ! get block from proc 0
     call mpi_recv(zx,nbuf,MPI_COMPLEX,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving zx'
     call mpi_recv(zy,nbuf,MPI_COMPLEX,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving zy'
     call mpi_recv(zz,nbuf,MPI_COMPLEX,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving zz'
     call mpi_recv(tt,nbuf,MPI_COMPLEX,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,istatus)
     if (istatus.ne.0) print*,mype,'Error receiving tt'
  endif
#endif

  if (mype.eq.0) then
     istatus = nf_close(idncx)
     if (istatus.ne.0) print*,'Error closing ZkX.in'
     istatus = nf_close(idncy)
     if (istatus.ne.0) print*,'Error closing ZkY.in'
     istatus = nf_close(idncz)
     if (istatus.ne.0) print*,'Error closing ZkZ.in'
     istatus = nf_close(idnct)
     if (istatus.ne.0) print*,'Error closing ZkT.in'
  endif

end subroutine ncreadrst

#endif



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c KILLING SUBROUTINES
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine killwaves(zx,zy,zz,tt,geok,gw1k,gw2k,u,v,w,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)

! Killing the gravity waves

   implicit none 
   include 'param.inc'
   
   integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp),ikz0
   complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,u,v,w,geok,gw1k,gw2k
   complex :: div,zk,tk,dk,bk,zi
   real :: aj,bj,f,omega,bv2,bv,norm,f2,sqr2,L1,L2,L3,twopi
   real :: wkh,wkh2,kx,ky,kz,wk,wk2,wkhn
   real :: kxa(iktx),kya(ikty),kza(iktz)
   
   integer :: mype
   common/mpi/mype

   call wtoab(zx,zy,zz,tt,geok,gw1k,gw2k,u,v,w,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)
 

   do ikx = 1,IKTX
      do iky = 1,IKTY
         do ikz = 1,IKTZP
            gw1k(ikx,iky,ikz) = cmplx(0.,0.)
            gw2k(ikx,iky,ikz) = cmplx(0.,0.)
         enddo
      enddo
   enddo

   call atowb(geok,gw1k,gw2k,zx,zy,zz,tt,u,v,w,L,kxa,kya,kza,aj,bj,f,L1,L2,L3)
 
     return
end subroutine  killwaves

subroutine killgeo(zx,zy,zz,tt,geok,gw1k,gw2k,u,v,w,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)

! Killing the gravity waves

   implicit none 
   include 'param.inc'
   
   integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp),ikz0
   complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,u,v,w,geok,gw1k,gw2k
   complex :: div,zk,tk,dk,bk,zi
   real :: aj,bj,f,omega,bv2,bv,norm,f2,sqr2,L1,L2,L3,twopi
   real :: wkh,wkh2,kx,ky,kz,wk,wk2,wkhn
   real :: kxa(iktx),kya(ikty),kza(iktz)
   
   integer :: mype
   common/mpi/mype

   call wtoab(zx,zy,zz,tt,geok,gw1k,gw2k,u,v,w,L,kxa,kya,kza,aj,bj,F,L1,L2,L3)
 

   do ikx = 1,IKTX
      do iky = 1,IKTY
         do ikz = 1,IKTZP
            geok(ikx,iky,ikz) = cmplx(0.,0.)
         enddo
      enddo
   enddo

   call atowb(geok,gw1k,gw2k,zx,zy,zz,tt,u,v,w,L,kxa,kya,kza,aj,bj,f,L1,L2,L3)
 
     return
end subroutine  killgeo


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c MISC SUBROUTINES
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine proj(zx,zy,zz,L,kxa,kya,kza)

! Fourier-space determination of the solenoidal part of a vector zx,y,z.

  implicit none
  include 'param.inc'

  integer :: ikx,iky,ikz,ikza,L(iktx,ikty,iktzp)
  real :: kx,ky,kz,k2,kxa(iktx),kya(ikty),kza(iktz)
  complex, dimension(iktx,ikty,iktzp) :: zx,zy,zz
  complex :: c1,c2,c3

  integer :: mype
  common/mpi/mype

  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx = 1,iktx
	           kx = kxa(ikx)     
           if (L(ikx,iky,ikz).eq.1) then
              k2 = max(kx*kx + ky*ky + kz*kz, 1.e-15)
              c1 =  (k2-kx*kx)*zx(ikx,iky,ikz) - kx*ky*zy(ikx,iky,ikz) - kx*kz*zz(ikx,iky,ikz)
              c2 = -ky*kx*zx(ikx,iky,ikz) + (k2-ky*ky)*zy(ikx,iky,ikz) - ky*kz*zz(ikx,iky,ikz)
              c3 = -kz*kx*zx(ikx,iky,ikz) - kz*ky*zy(ikx,iky,ikz) + (k2-kz*kz)*zz(ikx,iky,ikz)
              zx(ikx,iky,ikz) = c1 / k2
              zy(ikx,iky,ikz) = c2 / k2
              zz(ikx,iky,ikz) = c3 / k2
           endif
        enddo
     enddo
  enddo
  return
end subroutine proj


subroutine realit(zk)

! Enforces the reality condition on the plane kx=0
! by writing on modes with L(ikx,iky,ikz) = 0.

  implicit none
  include 'param.inc'
#if (MPI == 1)
  include 'mpif.h'
#endif

  integer :: n2h,n3h
  integer :: ikx,iky,ikz,kz,inkz,ky,inky 
  complex :: zk(iktx,ikty,iktzp)

  integer :: ikza
  real :: kyy,kzz
  integer, parameter :: iktyh = ikty/2
  integer :: mype,ierr
  common/mpi/mype

#if ( MPI == 1 )
  integer :: status(MPI_STATUS_SIZE)
  complex :: buf1(iktyh,iktzp-1),buf2(iktyh)
  integer :: nto,nfrom,nbuf1,nbuf2,nph
#endif

#if ( MPI == 1 )
  nph   = npe/2
  nbuf1 = iktyh*(iktzp-1)
  nbuf2 = iktyh
#endif

  n2h = n2/2
  n3h = n3/2

! First, negative ky axis; no communication required
! Set zk(0,-ky,0) = conjg(zk(0,ky,0)) 
  if (mype.eq.0) then
     ikx = 1
     ikz = 1
     do ky=1,n2h-1
        iky  = n2h+1-ky
        inky = n2h+1+ky
        zk(ikx,inky,ikz) = conjg( zk(ikx,iky,ikz) )           
     enddo
  endif


#if ( MPI == 0 )
  ikx = 1
  iky = 1
  do kz=1,n3h-1
     ikz  = n3h+1-kz
     inkz = n3h+1+kz
     zk(ikx,iky,inkz) = conjg( zk(ikx,iky,ikz) )           
  enddo

  ikx = 1
  do kz=1,n3h-1
     ikz  = n3h+1-kz
     inkz = n3h+1+kz
     do ky=1,n2h-1
        iky  = n2h+1-ky 
        inky = n2h+1+ky 
        zk(ikx,inky,inkz) = conjg( zk(ikx,iky,ikz) )
        zk(ikx,inky,ikz)  = conjg( zk(ikx,iky,inkz) )
     enddo
  enddo
#else /* MPI */
! Next, send kz>0,ky>=0 to kz<0,ky<=0
! Write zk(0,-ky,-kz) = conjg(zk(0,ky, kz)) 
! Have to do it in 2 parts

  if (mype.le.nph-1) then  
     nto  = npe-mype-1
     buf1 = zk(1,1:ikty/2,2:iktzp)
     call mpi_send(buf1,nbuf1,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,ierr)
  else
     nfrom = npe-mype-1
     call mpi_recv(buf1,nbuf1,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,ierr)
     ikx = 1
     do kz=2,iktzp
        ikz  = kz-1
        inkz = iktzp+2-kz
        iky = 1
        zk(ikx,iky,inkz) = conjg( buf1(iky,ikz) )  ! do ky=0
        do ky=1,n2h-1
           iky  = n2h+1-ky 
           inky = n2h+1+ky 
           zk(ikx,inky,inkz) = conjg( buf1(iky,ikz) )
        enddo
     enddo
  endif

  if ((mype.gt.0).and.(mype.le.nph-1)) then  
     nto  = npe-mype
     buf2 = zk(1,1:ikty/2,1)
     call mpi_send(buf2,nbuf2,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,ierr)
  endif
  if (mype.gt.nph) then
     nfrom = npe-mype
     call mpi_recv(buf2,nbuf2,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,ierr)
     ikx  = 1
     inkz = 1
     iky  = 1
     zk(ikx,iky,inkz) = conjg( buf2(iky) )  ! do ky=0
     do ky=1,n2h-1
        iky  = n2h+1-ky 
        inky = n2h+1+ky 
        zk(ikx,inky,inkz) = conjg( buf2(iky) )
     enddo
  endif

! Finally, send kz<0,ky>=0 to kz>0,ky<=0
! Write zk(0,-ky,kz) = conjg(zk(0,ky,-kz)) 
! Have to do it in 2 parts

  if (mype.ge.nph) then  
     nto  = npe-mype-1
     buf1 = zk(1,1:ikty/2,2:iktzp)
     call mpi_send(buf1,nbuf1,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,ierr)
  else
     nfrom = npe-mype-1
     call mpi_recv(buf1,nbuf1,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,ierr)
     ikx = 1
     do kz=2,iktzp
        ikz  = kz-1
        inkz = iktzp+2-kz
        iky = 1
        do ky=1,n2h-1
           iky  = n2h+1-ky 
           inky = n2h+1+ky 
           zk(ikx,inky,inkz) = conjg( buf1(iky,ikz) )
        enddo
     enddo
  endif

  if (mype.gt.nph) then
     nto  = npe-mype
     buf2 = zk(1,1:ikty/2,1)
     call mpi_send(buf2,nbuf2,MPI_COMPLEX,nto,1,MPI_COMM_WORLD,ierr)
  endif
  if ((mype.gt.0).and.(mype.le.nph-1)) then
     nfrom = npe-mype
     call mpi_recv(buf2,nbuf2,MPI_COMPLEX,nfrom,1,MPI_COMM_WORLD,status,ierr)
     ikx  = 1
     inkz = 1
     iky  = 1
     do ky=1,n2h-1
        iky  = n2h+1-ky 
        inky = n2h+1+ky 
        zk(ikx,inky,inkz) = conjg( buf2(iky) )
     enddo
  endif
#endif

  return
end subroutine realit



#if ( MPI == 1 )
subroutine mpitranspose(x,n1,n2,n3p,xt,n3,n2p,npe,xsend,xrecv)

  implicit none
  include 'mpif.h'

  integer :: status(MPI_STATUS_SIZE),npe
  integer :: n1,n2,n3,n2p,n3p
  integer :: i,j,k,ip,ito,ifrom,jblock,kblock,koff,n123p
  complex, dimension(n1,n2,n3p) :: x
  complex, dimension(n1,n3,n2p) :: xt
  complex, dimension(n1,n2p,n3p,npe) :: xsend,xrecv

  integer :: mype,ierror
  common/mpi/mype

  n123p=n1*n2p*n3p

  do ip=1,npe
     do k=1,n3p
        do j=1,n2p
           jblock=j+(ip-1)*n2p
           do i=1,n1
              xsend(i,j,k,ip) = x(i,jblock,k)
           enddo
        enddo
     enddo
  enddo
  
  call mpi_alltoall(xsend,n123p,MPI_COMPLEX,xrecv,n123p,MPI_COMPLEX,MPI_COMM_WORLD,ierror)
 
  do ip=1,npe
     do j=1,n2p
        do k=1,n3p
           do i=1,n1
              xt(i,k+(ip-1)*n3p,j) = xrecv(i,j,k,ip)
           enddo
        enddo
     enddo
  enddo
  
end subroutine mpitranspose

#else /* MPI == 0 */

subroutine serialtranspose(x,xt,n1,n2,n3)

  implicit none

  integer :: n1,n2,n3
  integer :: i,j,k
  complex, dimension(n1,n2,n3) :: x
  complex, dimension(n1,n3,n2) :: xt

  do j=1,n2
     do k=1,n3
        do i=1,n1
           xt(i,k,j) = x(i,j,k)
        enddo
     enddo
  enddo

end subroutine serialtranspose
#endif 

