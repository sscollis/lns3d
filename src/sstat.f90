!=============================================================================!
        module stats
        
        real, allocatable :: prms(:)
        
        end module stats
!=============================================================================!
        subroutine sstat(vl,vml,xl,yl) 
!  
! Output for the wave problem
!  
!=============================================================================!
        use global
        use stats
        implicit none

        real :: vl(ny,nx,ndof), vml(ny,nx,ndof), xl(ny,nx), yl(ny,nx)
        integer :: i, j
        
        real :: rhom(nx), tm(nx), cm(nx)
        real :: rho(nx), u1(nx), u2(nx), u3(nx), t(nx), p(nx)
        real :: c3(nx)
!=============================================================================!

        if (istep .eq. 1) then
          open(unit=53,file='wave.dat',form='unformatted',status='unknown')
        end if

        if ( mod(istep,itout) .eq. 0 ) then
          j = 1

          rhom = vml(j,:,1)
          tm   = vml(j,:,5)
          cm   = one / Ma * sqrt( tm )

          rho = vl(j,:,1)
          u1  = vl(j,:,2)
          u2  = vl(j,:,3)
          u3  = vl(j,:,4)
          t   = vl(j,:,5)
          p   = one / (gamma * Ma**2) * ( rhom * t + tm * rho )

          do i = 1, nx
            c3 =  rhom * cm * u1 + p
          end do

          write(53) istep, time, c3
        end if

        return
        end
!=============================================================================!
        subroutine sstat_1(vl,vml,xl,yl) 
!  
! Used for the cylinder scattering calculation
!  
!=============================================================================!
        use global
        use stats
        implicit none

        real :: vl(ny,nx,ndof), vml(ny,nx,ndof), xl(ny,nx), yl(ny,nx)
        integer :: i, j, ix
        
        real :: rhom(ny), tm(ny), cm(ny)
        real :: rho(ny), u1(ny), u2(ny), u3(ny), t(ny), p(ny)
        real :: c1(ny), c2(ny), c3(ny), c4(ny)
        
        real :: rr, rho0, c0, lambda, d, kk, a
!=============================================================================!
        ix = nx / 2

        if (istep .eq. 1) then
          allocate( prms(ny) )
          prms = zero
        end if
        
        rhom = vml(:,ix,1)
        tm   = vml(:,ix,5)
        cm   = one / Ma * sqrt( tm )

        rho = vl(:,ix,1)
        u1  = vl(:,ix,2)
        u2  = vl(:,ix,3)
        u3  = vl(:,ix,4)
        t   = vl(:,ix,5)
        p   = one / (gamma * Ma**2) * ( rhom * t + tm * rho )

        rr     = 1.0
        rho0   = 1.0
        c0     = 1.0
!       x0     = -11.0
        lambda = 2.5
        d      = (one/Re) * pt5 / rho0 * ( onept33 + gamma1 / Pr )
        kk     = two * pi / lambda
        a      = kk**2 * d / c0

        c1 = zero
        c2 = zero
        c3 = cos( kk * ( xl(:,ix) - cm * time) ) * &
             exp( -a * (xl(:,ix) - x0) )
        c4 = zero

        rho = rho - ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
        u1  = u1 - ( c3 - c4 ) * pt5 / ( rhom * cm )
        u2  = u2 -  one / ( rhom * cm ) * c2
        p   = p - ( c3 + c4 ) * pt5
        t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

        p   = p * exp( a * (abs(xl(:,ix)) - one) )      ! viscosity correction
        prms = prms + p**2
        
        if (istep .eq. nstep) then
          prms = sqrt( prms/float(nstep) )
          do j = 1, ny
            write(80,"(2(1pe13.6,1x))") xl(j,ix), prms(j)
          end do
        end if
                
        return
        end
