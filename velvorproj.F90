!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module contains the conversion between vorticity and velocity (VELOCITY <--> VORTICITY)   !!
!! And also determines the solenoidal part of vorticity field                                     !!
!! -->  Note that vorticity has to solenoindal since it is the curl of velocity
!! The subroutines in this module: velo, vort, proj                                               !!
!! One of the essential parts of the code that cannot be eliminated even in the lightest version  !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module velvorproj
  use param
  implicit none

CONTAINS
  subroutine velo(zx,zy,zz,u,v,w)
    ! Calculates k-space velocity from k-space vorticity.
    ! curl (vorticity) = - laplacian (velocity) if velocity is solenoidal.
    implicit none
    include 'mpif.h'
    complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz
    complex, intent(out), dimension(iktx,ikty,iktzp) :: u,v,w
    integer :: ikx,iky,ikz,ikza
    real :: kx,ky,kz,k2
    complex :: c1,c2,c3

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
                u(ikx,iky,ikz) = zi * c1 / k2
                v(ikx,iky,ikz) = zi * c2 / k2
                w(ikx,iky,ikz) = zi * c3 / k2
             endif
          enddo
       enddo
    enddo
    
    return
  end subroutine velo

  subroutine vort(u,v,w,zx,zy,zz)
    ! Calculates k-space vortcity from k-space velocity.
    implicit none
    include 'mpif.h'
    complex, intent(in), dimension(iktx,ikty,iktzp) :: u,v,w
    complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz
    integer :: ikx,iky,ikz,ikza
    real ::  kx,ky,kz
    complex :: c1,c2,c3

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
                zx(ikx,iky,ikz) = zi * c1
                zy(ikx,iky,ikz) = zi * c2
                zz(ikx,iky,ikz) = zi * c3
             endif
          enddo
       enddo
    enddo
    
    return
  end subroutine vort

  subroutine proj(zx,zy,zz)
    ! Fourier-space determination of the solenoidal part of a vector zx,y,z.
    implicit none
    include 'mpif.h'
    complex, intent(inout), dimension(iktx,ikty,iktzp) :: zx,zy,zz
    
    integer :: ikx,iky,ikz,ikza
    real :: kx,ky,kz,k2
    complex :: c1,c2,c3
    
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
 
end module velvorproj
