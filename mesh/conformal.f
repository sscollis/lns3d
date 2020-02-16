c=============================================================================c
      subroutine conformal( npts, s, xb, yb, cur )
c=============================================================================c
      parameter (zero=0.0d0,    pt5=0.5d0,    pt25=0.25d0,    one=1.0d0, 
     &           onept5=1.5d0,  two=2.0d0,    twopt5=2.5d0, 
     &           three=3.0d0,   pi=3.1415926535897932385d+0,  four=4.0,
     &           tol=1.0e-7)

      dimension s(npts), xb(npts), yb(npts), cur(npts)

      parameter ( mpts=2048 )
      dimension x(mpts), y(mpts), th(mpts), psi(mpts)
      dimension phi(mpts), eps(mpts), tmp(mpts)

      parameter ( nint = 800 )
      dimension aa(nint), bb(nint)
      complex   cc(nint), z, zeta, rk(nint), rh(nint), ra(nint)

      common /tmp1/ n, korder, xknot(mpts), bs(mpts), phip

      external integrand, RKQCR, BSITG
c=============================================================================c
      write(*,"(/,'Conformal is not tested, therefore skipping...')")
      return

      korder = 5
      n = npts
      if (n .gt. mpts) then
         write(*,"('ERROR:  Increase mpts in conformal ')")
         call exit(1)
      end if
c
c.... Map the airfoil to an approximate circle
c
      open(20,file='map.out')
      rn = one / cur(1)
c     a  = (one - pt5 * rn) * pt25
c     x0 = pt5 * rn + two * a
      a  = (one - rn) * pt25
      x0 = rn + two * a
      y0 = zero
      do i = 1, npts
         x(i) = -(xb(i) - x0)
         y(i) = yb(i) - y0
         p = one - (x(i)/(two*a))**2 - (y(i)/(two*a))**2
         th(i) = asin( sqrt( pt5*(p + sqrt(p**2 + (y(i)/a)**2)) ) )
         if (x(i) .lt. zero) th(i) = pi - th(i)
c.error  psi(i) = asinh( sqrt( pt5*(-p + sqrt(p**2 + (y(i)/a)**2)) ) )
         write(20,100) s(i), x(i), y(i), th(i), psi(i), 
     &                 a*exp(psi(i))*cos(th(i)), 
     &                 a*exp(psi(i))*sin(th(i))
      end do
      close(20)
c
c.... correct for finite precision
c
      th(1)  = zero
      th(n)  = pi
      psi(n) = zero
c
c.... initialize
c
      do i = 1, n
         eps(i) = zero
      end do

      do k = 1, nint-1
         aa(k) = log( abs( sin(pi*(two*k+one)/(four*nint)) / 
     &                     sin(pi*(two*k-one)/(four*nint)) ) )
      end do
c
c.... compute the mapping to a circle
c
      niter = 1
      write(*,"('Enter niter ==> ',$)")
      read(*,*) niter
      do iter = 1, niter
         write(*,"('Interation = ',i4)") iter
         do i = 1, n
            phi(i) = th(i) + eps(i)
         end do
         call BSNAK(n, phi, korder, xknot)
         call BSINT(n, phi, psi, korder, xknot, bs )

         do i = 1, n
            phip   = phi(i)
            eps(i) = pi/nint * BSDER( 1, phip, korder, xknot, n, bs )
            do k = 1, nint-1
               
               phia = phip + k * pi / nint
               if (phia .ge. two*pi) phia = phia - two*pi
               
               phib = phip - k * pi / nint
               if (phib .le. zero) phib = phib + two*pi
               
               if (phia .ge. zero .and. phia .le. pi) then
                  psia = BSDER( 0, phia, korder, xknot, n, bs )
               else
                  phia = two * pi - phia
                  psia = BSDER( 0, phia, korder, xknot, n, bs )
               end if
               
               if (phib .ge. zero .and. phib .le. pi) then
                  psib = BSDER( 0, phib, korder, xknot, n, bs )
               else
                  phib = two * pi - phib
                  psib = BSDER( 0, phib, korder, xknot, n, bs )
               end if
               
               eps(i) = eps(i) + aa(k) * ( psia - psib )
            end do
            eps(i) = -one/pi * eps(i)
         end do
c
c.... enforce symmetry
c
         eps(1) = zero
         eps(n) = zero

      end do

      do i = 1, n
         phi(i) = th(i) + eps(i)
      end do
      call BSNAK(n, phi, korder, xknot)
      call BSINT(n, phi, psi, korder, xknot, bs )
      
      open(20,file='conf.out')
      do i = 1, n
         write(20,100) phi(i), th(i), psi(i), eps(i),
     &   BSDER( 0, phi(i), korder, xknot, n, bs ),
     &   BSDER( 1, phi(i), korder, xknot, n, bs )
      end do
      close(20)
c
c.... compute psi0
c
      psi0 = one/pi * BSITG( zero, pi, korder, xknot, n, bs )
      write(*,"('Psi0 = ',1pe14.7)") psi0
c
c.... check to make sure that they are congugate
c
      call BSNAK(n, phi, korder, xknot)
      call BSINT(n, phi, eps, korder, xknot, bs )

      do i = 1, n
         phip   = phi(i)
         psi(i) = pi/nint * BSDER( 1, phip, korder, xknot, n, bs )
         do k = 1, nint-1
            
            phia = phip + k * pi / nint
            if (phia .ge. two*pi) phia = phia - two*pi
            
            phib = phip - k * pi / nint
            if (phib .le. zero) phib = phib + two*pi
            
            if (phia .ge. zero .and. phia .le. pi) then
               epsa = BSDER( 0, phia, korder, xknot, n, bs )
            else
               phia = two * pi - phia
               epsa = -BSDER( 0, phia, korder, xknot, n, bs )
            end if
            
            if (phib .ge. zero .and. phib .le. pi) then
               epsb = BSDER( 0, phib, korder, xknot, n, bs )
            else
               phib = two * pi - phib
               epsb = -BSDER( 0, phib, korder, xknot, n, bs )
            end if
            
            psi(i) = psi(i) + aa(k) * ( epsa - epsb )
         end do
         psi(i) = one/pi * psi(i)
      end do

      open(20,file='check.out')
      do i = 1, n
         write(20,100) phi(i), th(i), psi(i), eps(i),
     &                 BSDER( 1, phi(i), korder, xknot, n, bs )
      end do
      close(20)
c
c.... compute the coefficients of the mapping 
c.... (The b are zero due to symmetry)
c
      R = a * exp( psi0 )
      write(*,"('a = ',1pe13.7)") a
      write(*,"('R = ',1pe13.7)") R
      do k = 1, 10
         do i = 1, n
            tmp(i) =  psi(i) * cos( real(k) * phi(i) )
         end do
         call BSINT(n, phi, tmp, korder, xknot, bs )
         aa(k) = R**k * two / pi * BSITG( zero, pi, korder, xknot, n, bs )
         write(*,"('A(',i2,') = ',1pe14.7)") k, aa(k)
      end do

      do k = 1, 10
         bb(k) = zero
         cc(k) = cmplx(aa(k),bb(k))
      end do

      rk(1) = cc(1)
      rk(2) = 1.0 / 2.0 * ( rk(1) * cc(1) + 2.0 * cc(2) )
      rk(3) = 1.0 / 3.0 * ( rk(2) * cc(1) + 2.0 * rk(1) * cc(2) + 
     &                      3.0 * cc(3) )
      rk(4) = 1.0 / 4.0 * ( rk(3) * cc(1) + 2.0 * rk(2) * cc(2) + 
     &                      3.0 * rk(1) * cc(3) + 4.0 * rk(4) )
      rk(5) = 1.0 / 5.0 * ( rk(4) * cc(1) + 2.0 * rk(3) * cc(2) + 
     &                      3.0 * rk(2) * cc(3) + 4.0 * rk(1) * cc(4) +
     &                      5.0 * cc(5) )
      rk(6) = 1.0 / 6.0 * ( rk(5) * cc(1) + 2.0 * rk(4) * cc(2) + 
     &                      3.0 * rk(3) * cc(3) + 4.0 * rk(2) * cc(4) +
     &                      5.0 * rk(1) * cc(5) + 6.0 * cc(6) )

      rh(1) = -cc(1)
      rh(2) = 1.0 / 2.0 * ( rk(1) * (-cc(1)) + 2.0 * (-cc(2)) )
      rh(3) = 1.0 / 3.0 * ( rk(2) * (-cc(1)) + 2.0 * rk(1) * (-cc(2)) + 
     &                      3.0 * (-cc(3)) )
      rh(4) = 1.0 / 4.0 * ( rk(3) * (-cc(1)) + 2.0 * rk(2) * (-cc(2)) + 
     &                      3.0 * rk(1) * (-cc(3)) + 4.0 * rk(4) )
      rh(5) = 1.0 / 5.0 * ( rk(4) * (-cc(1)) + 2.0 * rk(3) * (-cc(2)) + 
     &                      3.0 * rk(2) * (-cc(3)) + 4.0 * rk(1) * (-cc(4)) +
     &                      5.0 * (-cc(5)) )

      ra(1) = rk(2) + a**2
      ra(2) = rk(3) + a**2 * rh(1)
      ra(3) = rk(4) + a**2 * rh(2)
      ra(4) = rk(5) + a**2 * rh(3)
      ra(5) = rk(6) + a**2 * rh(4)

      open(20,file='body.new')
      do i = 1, n
         z = a * exp(psi0) * exp( cmplx(zero,one) * phi(i) )

         zeta = z + cc(1) + ra(1) / z + ra(2) / z**2 + ra(3) / z**3  
     &        + ra(4) / z**4 
c    &        + ra(5) / z**5

         write(20,100) real(zeta), aimag(zeta), real(z), aimag(z)
      end do
      close(20)

      return
 100  format(8(1pe13.6,1x))
      end

c=============================================================================c
      subroutine integrand( nvar, epsl, phil, depsdphi )
c=============================================================================c
      parameter (zero=0.0d0,    pt5=0.5d0,    pt25=0.25d0,    one=1.0d0, 
     &           onept5=1.5d0,  two=2.0d0,    twopt5=2.5d0, 
     &           three=3.0d0,   big=1.0e10,   pi=3.1415926535897932385d+0,
     &           tol=1.0e-7)

      parameter ( mpts=2048 )
      common /tmp1/ n, korder, xknot(mpts), bs(mpts), phip
c=============================================================================c

      dpsidphi = BSDER( 1, phil, korder, xknot, n, bs )

      depsdphi = one/pi * ( dpsidphi * log(abs(sin(pt5*(phil - phip)))) ) 

      return
      end
