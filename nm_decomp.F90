!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module contains Normal Mode Decomposition as described in Bartello, 1995, J. Atmos. Sci.  !!
!! The subroutines in this module: wtoab, atowb                                                   !!
!! Modules used : param.F90 and velvorproj.F90                                                    !!
!! Required for: ncf2Dspc.F90 boussinesq.F90 (main file)                                          !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module nm_decomp
  use param
  use velvorproj
  implicit none

CONTAINS

subroutine wtoab(zx,zy,zz,tt,geok,gw1k,gw2k,u,v,w)
! Converts from (vorticity, buoyancy) to (geok,grav.wave_1k,grav.wave_2k)
  implicit none 
  complex, intent(in),  dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: u,v,w
  complex, intent(out), dimension(iktx,ikty,iktzp) :: geok,gw1k,gw2k
  integer :: ikx,iky,ikz,ikza,ikz0
  real :: omega,norm
  real :: wkh,wkh2,kx,ky,kz,wk,wk2,wkhn
  complex :: div,zk,tk,dk,bk
   
  geok = cmplx(0.,0.)
  gw1k = cmplx(0.,0.)
  gw2k = cmplx(0.,0.)
  call velo(zx,zy,zz,u,v,w)
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
            omega = sqrt(cor2*kz*kz + bf2*wkh2)/wk 
            div   = zi*(kx*u(ikx,iky,ikz)+ky*v(ikx,iky,ikz))
            bk    = aj*tt(ikx,iky,ikz)
            
            zk = zz(ikx,iky,ikz)
            dk = (wk/kz) * div
            tk = (wkh/bf) * bk
            
            norm  = omega*wk
            geok(ikx,iky,ikz) = (bf*wkh*zk + zi*cor*kz*tk)/norm
            
            norm  = sqrt2*omega*wk
            gw1k(ikx,iky,ikz) = (-zi*cor*kz*zk + omega*wk*dk - bf*wkh*tk)/norm
            
            norm  = sqrt2*omega*wk
            gw2k(ikx,iky,ikz) = (+zi*cor*kz*zk + omega*wk*dk + bf*wkh*tk)/norm
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
          gw1k(ikx,iky,1)=w(ikx,iky,1)-zi*bk/bf
          gw2k(ikx,iky,1)=w(ikx,iky,1)+zi*bk/bf
        endif
      enddo
    enddo
  endif
   
  return
end subroutine wtoab

subroutine atowb(geok,gw1k,gw2k,zx,zy,zz,tt,u,v,w)
! Converts from (geok,grav.wave_1k,grav.wave_2k) to (zeta,d,t).
  implicit none 
  complex, intent(in), dimension(iktx,ikty,iktzp) :: geok,gw1k,gw2k
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: u,v,w
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt

  integer :: ikx,iky,ikz,ikz0,ikza
  real :: omega
  real :: wkh,wkh2,kx,ky,kz,wk,wk2,wkhn
  complex :: zk,dk,tk,bk,gn,psi,div
  complex :: zgk,z1k,z2k,dgk,d1k,d2k,tgk,t1k,t2k

  u  = cmplx(0.,0.)
  v  = cmplx(0.,0.)
  w  = cmplx(0.,0.)
  tt = cmplx(0.,0.)
  zx = cmplx(0.,0.)
  zy = cmplx(0.,0.)
  zz = cmplx(0.,0.)

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
          omega = sqrt(cor2*kz*kz + bf2*wkh2)/wk
          
          gn  = geok(ikx,iky,ikz) / (omega*wk)
          zgk = bf * wkh       * gn
          dgk = cmplx(0.,0.)
          tgk = - zi * cor * kz * gn
          
          gn  = gw1k(ikx,iky,ikz) / (sqrt2*omega*wk)
          z1k = + zi  * cor * kz * gn
          d1k = omega * wk     * gn
          t1k = - bf   * wkh    * gn
          gn  = gw2k(ikx,iky,ikz) / (sqrt2*omega*wk)
          z2k = - zi  * cor * kz * gn
          d2k = omega * wk     * gn
          t2k = + bf   * wkh   * gn
          
          zk = zgk + z1k + z2k
          dk = dgk + d1k + d2k
          tk = tgk + t1k + t2k
               
          if (wkhn.gt.1.e-10) then
            div             = dk*kz/wk
            u(ikx,iky,ikz)  = +zi*(ky*zk-kx*div)/wkh2             
            v(ikx,iky,ikz)  = -zi*(kx*zk+ky*div)/wkh2             
            w(ikx,iky,ikz)  = zi*div/kz
            tt(ikx,iky,ikz) = bf*tk/(aj*wkh)
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
          bk = zi*bf*0.5*(gw1k(ikx,iky,1)-gw2k(ikx,iky,1))
          tt(ikx,iky,1) = bk/aj
        endif
      enddo
    enddo
  endif     
  call vort(u,v,w,zx,zy,zz)
  
  return
end subroutine atowb
 
end module nm_decomp
