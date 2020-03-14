!=============================================================================c
        program ogrid
!=============================================================================c
!
!  Purpose:  Given a closed body, this program generates an O-grid 
!            suitable for Euler calculations
!
!  Author:   S. Scott Collis
!
!  Revised:  6-17-97
!
!=============================================================================c
        implicit none

        real, parameter :: zero=0.0d0, pt5=0.5d0, one=1.0d0, onept5=1.5d0, &
                           two=2.0d0, twopt5=2.5d0, three=3.0d0,           &
                           pi=3.1415926535897932385d+0

        integer, parameter :: iin=10, iout=11, ibs=12

        integer :: i, j, im, jm

        real, allocatable ::  x(:,:), y(:,:)
        real :: dx, theta, dtheta, r = 10.0, d, dd

        integer :: nmax
        real :: eps

        real, external :: sharp_foil
        real, external :: rounded_foil

        character(80) :: name

        integer :: narg, iarg
        character(128) arg
        logical :: airfoil=.false.
        logical :: rounded=.true.
!=============================================================================c
 
!.... parse the argument list
 
        narg = iargc()
        do iarg = 1, narg
          call getarg(iarg,arg)
          select case (arg(1:2))
          case ('-a')
            airfoil = .true.
            rounded = .true.
          case ('-s')
            airfoil = .true.
            rounded = .false.
          case ('-c')
            airfoil = .false.
          case ('-h')
            write(*,"('-----------------------------------------------')")
            write(*,"('Usage:  ogrid [options]                        ')")
            write(*,"('-----------------------------------------------')")
            write(*,"('   -h:  this help                              ')")
            write(*,"('-----------------------------------------------')")
            write(*,"('   -a:  use rounded TE airfoil geometry        ')")
            write(*,"('   -s:  use sharp TE airfoil geometry        ')")
            write(*,"('   -c:  use circle geometry (default)          ')")
            write(*,"('-----------------------------------------------')")
            call exit(0)
          case default
            write(*,"('Argument ',i2,' ignored.')") iarg
          end select
        end do
 
!.... get the maximum indices in ij directions
 
        write(*,"('Enter im, jm ==> ',$)")
        read(*,*) im, jm

        allocate( x(im,jm), y(im,jm) )

!.... make the mesh along the body

        dx = one / real(im - 1)
        dtheta = two * pi / real(im-1)
        if (airfoil) then
          name = "foil.xyz"
          write(*,*) "Building mesh for an airfoil"
          do i = 1, im
#if 0
            theta = pi - real(i-1) * dtheta
            x(i,1) = real(i-1) * dx
#else
            theta = -real(i-1) * dtheta
            x(i,1) = ( cos( theta ) + one) * pt5
#endif
            if (rounded) then
              y(i,1) = rounded_foil( x(i,1), 0.1 )
            else
              y(i,1) = sharp_foil( x(i,1), 0.1 )
            endif
            if (theta .gt. -pi) y(i,1) = -y(i,1)
            write (10,*) x(i,1), y(i,1), theta
          end do
        else
          name = "circ.xyz"
          write(*,*) "Building mesh for a circular cylinder"
          do i = 1, im
            theta = - real(i-1) * dtheta
            x(i,1) = ( cos( theta ) + one) * pt5
            y(i,1) = ( sin( theta ) ) * pt5
            write (10,*) x(i,1), y(i,1), theta
          end do
        end if

!.... make the mesh along the top circle

        dtheta = two * pi / real(im-1)
        do i = 1, im
!          theta = pi - real(i-1) * dtheta
           theta = - real(i-1) * dtheta
           x(i,jm) = r * cos(theta) + 0.5
           y(i,jm) = r * sin(theta)
        end do

!.... now connect the boundary points by straight lines

        do i = 1, im
           dd  = sqrt( (x(i,jm)-x(i,1))**2 + (y(i,jm)-y(i,1))**2 ) / real(jm-1)
           theta = atan2( y(i,jm)-y(i,1), x(i,jm)-x(i,1) )
           do j = 2, jm-1
              d = real(j-1) * dd
              x(i,j) = x(i,1) + d * cos(theta)
              y(i,j) = y(i,1) + d * sin(theta)
           end do
        end do

!.... elliptic grid generator

        write(*,"(/,'Using elliptic smoother...')")
        write(*,"('Enter eps, nmax ==> ',$)")
        read(*,*) eps, nmax
        call elliptic( im, jm, x, y, eps, nmax)

!.... output grid file

        open(unit=iout, file=name, form='unformatted')
        write(iout) im, jm, 1
        write(iout) ((x(i,j), i=1,im), j=1,jm), &
                    ((y(i,j), i=1,im), j=1,jm), &
                    ((  zero, i=1,im), j=1,jm)
        close(iout)

        call exit(0)
        end 

!=============================================================================c
        function rounded_foil(x,t)
!=============================================================================c
        implicit none
        real foil, x, t
        real r, x0, y0, y, term, xp
        real :: rounded_foil
        real, external :: sharp_foil
!=============================================================================c
        r  = 0.00255       ! SSC: hand tuned
        x0 = 1.0 - r
        y0 = 0.0
        xp = x/(1.0+0.02)  ! SSC:  hand tuned
        foil = sharp_foil(xp,t)
        if (x.ge.x0) then
          term = r**2 - (x-x0)**2
          if (term.le.0) then
            y = 0.0
          else
            y = sqrt(term) + y0
          end if
          write (*,*) x, y, foil 
          if (y.lt.foil) foil = y
        end if 
        rounded_foil = foil
        return
        end

!=============================================================================c
        function sharp_foil(x,t)
!=============================================================================c
        implicit none
        real foil, x, t
        real sharp_foil
!=============================================================================c
        foil = t/0.2 * ( 0.2969 * sqrt(x) - 0.1281 * x - 0.3516 * x**2 + &
               0.2843 * x**3 - 0.1015 * x**4 )
        sharp_foil = foil
        return
        end

!=============================================================================c
      subroutine elliptic( nx, ny, x, y, eps, nmax )
!
!     This routine performs an elliptic smoothing of the input grid
!     
!     Author:  Scott Collis
!
!     Date:    6-17-97
!
!=============================================================================c
      implicit none

      integer :: nx, ny, nmax
      real :: x(nx,ny), y(nx,ny), eps

      real, parameter :: zero=0.0, pt5=0.5, one=1.0, two=2.0

      integer :: n, i, j, ibc=1
      real :: dxi, deta, dxdeta, dydeta, dxdxi, dydxi
      real :: a, b, c, xnew, ynew, err, jac

      real :: w, rx(nx), ry(nx), ad(nx), bd(nx), cd(nx), xl(nx), yl(nx)
      real :: xi(nx), eta(ny), p(nx), q(nx), dn(nx)
      real :: px, py
!=============================================================================c

      dxi = one / float(nx-1)
      deta = one / float(ny-1)

      do j = 1, ny
        eta(j) = real(j-1) * deta
      end do

      do i = 1, nx
        xi(i) = real(i-1) * dxi
        dn(i) = sqrt( (x(i,2)-x(i,1))**2 + (y(i,2)-y(i,1))**2 )
      end do

!.... Line Gauss-Siedel with SOR (because of the nonlinearity
!.... SOR is unstable for w > 1?)

      w = one
!     write(*,"(/,'Enter the SOR parmameter (0 < w < 2, w = 1.78) ==> ',$)")
!     read(*,*) w
      n = 0
 20   continue
         n = n + 1
         err = zero

         if (.false. .and. n .ge. 1) then
           call ortho(nx, ny, dxi, deta, x, y, dn, p, q)
         else
           p = zero
           q = zero
         end if

         do j = 2, ny-1

            do i = 2, nx-1
               dxdeta = (x(i,j+1)-x(i,j-1)) / (two * deta)
               dydeta = (y(i,j+1)-y(i,j-1)) / (two * deta)
               dxdxi  = (x(i+1,j)-x(i-1,j)) / (two * dxi)
               dydxi  = (y(i+1,j)-y(i-1,j)) / (two * dxi)

               a = (dxdeta**2 + dydeta**2) / dxi**2
               b = (dxdxi*dxdeta + dydxi*dydeta) / ( two * dxi * deta )
               c = (dxdxi**2 + dydxi**2) / deta**2
               
               rx(i) = -two*(one-w)*(a + c) * x(i,j) + w*b*(x(i+1,j+1) - &
                       x(i+1,j-1) + x(i-1,j-1) - x(i-1,j+1)) -           &
                       w*c*(x(i,j+1) + x(i,j-1))

               ry(i) = -two*(one-w)*(a + c) * y(i,j) + w*b*(y(i+1,j+1) - &
                       y(i+1,j-1) + y(i-1,j-1) - y(i-1,j+1)) -           &
                       w*c*(y(i,j+1) + y(i,j-1))

               jac = one / ( dxdxi * dydeta - dydxi * dxdeta )
               if (p(i) .gt. zero) then
                 dxdxi = (x(i,j)-x(i-1,j)) / dxi
                 dydxi = (y(i,j)-y(i-1,j)) / dxi
               else
                 dxdxi = (x(i+1,j)-x(i,j)) / dxi
                 dydxi = (y(i+1,j)-y(i,j)) / dxi
               end if
               if (q(i) .gt. zero) then
                 dxdeta = (x(i,j)-x(i,j-1)) / deta
                 dydeta = (y(i,j)-y(i,j-1)) / deta
               else
                 dxdeta = (x(i,j+1)-x(i,j)) / deta
                 dydeta = (y(i,j+1)-y(i,j)) / deta
               end if
               px = -one/jac**2 * ( p(i) * exp( -(eta(j)-eta(1)) ) * dxdxi + &
                                    q(i) * exp( -(eta(j)-eta(1)) ) * dxdeta )
               py = -one/jac**2 * ( p(i) * exp( -(eta(j)-eta(1)) ) * dydxi + &
                                    q(i) * exp( -(eta(j)-eta(1)) ) * dydeta )

               rx(i) = rx(i) + w * px
               ry(i) = ry(i) + w * py

               ad(i) = w * a
               bd(i) = -two * (a + c)
               cd(i) = w * a
            end do

!.... apply the boundary conditions

            if (ibc .eq. 0) then             ! fixed boundaries
               rx(1)  = x(1,j)
               ry(1)  = y(1,j)
               rx(nx) = x(nx,j)
               ry(nx) = y(nx,j)
            
               ad(1) = zero
               bd(1) = one
               cd(1) = zero
            
               ad(nx) = zero
               bd(nx) = one
               cd(nx) = zero

            else if (ibc .eq. 1) then        ! symmetry boundaries
               i = 1
               dxdeta = (x(i,j+1)-x(i,j-1)) / (two * deta)
               dydeta = (y(i,j+1)-y(i,j-1)) / (two * deta)
               dxdxi  = (x(i+1,j)-x(i+1,j)) / (two * dxi)
               dydxi  = (y(i+1,j)+y(i+1,j)) / (two * dxi)

               a = (dxdeta**2 + dydeta**2) / dxi**2
               b = (dxdxi*dxdeta + dydxi*dydeta) / ( two * dxi * deta )
               c = (dxdxi**2 + dydxi**2) / deta**2
               
               rx(i) = -two*(one-w)*(a + c) * x(i,j) + w*b*(x(i+1,j+1) - &
                       x(i+1,j-1) + x(i+1,j-1) - x(i+1,j+1)) -           &
                       w*c*(x(i,j+1) + x(i,j-1))

               ry(i) = -two*(one-w)*(a + c) * y(i,j) + w*b*(y(i+1,j+1) - &
                       y(i+1,j-1) - y(i+1,j-1) + y(i+1,j+1)) -           &
                       w*c*(y(i,j+1) + y(i,j-1))

               jac = one / ( dxdxi * dydeta - dydxi * dxdeta )
               if (p(i) .gt. zero) then
                 dxdxi = (x(i,j)-x(i+1,j)) / dxi
                 dydxi = (y(i,j)-y(i+1,j)) / dxi
               else
                 dxdxi = (x(i+1,j)-x(i,j)) / dxi
                 dydxi = (y(i+1,j)-y(i,j)) / dxi
               end if
               if (q(i) .gt. zero) then
                 dxdeta = (x(i,j)-x(i,j-1)) / deta
                 dydeta = (y(i,j)-y(i,j-1)) / deta
               else
                 dxdeta = (x(i,j+1)-x(i,j)) / deta
                 dydeta = (y(i,j+1)-y(i,j)) / deta
               end if
               px = -one/jac**2 * ( p(i) * exp( -(eta(j)-eta(1)) ) * dxdxi + &
                                    q(i) * exp( -(eta(j)-eta(1)) ) * dxdeta )
               py = -one/jac**2 * ( p(i) * exp( -(eta(j)-eta(1)) ) * dydxi + &
                                    q(i) * exp( -(eta(j)-eta(1)) ) * dydeta )

               rx(i) = rx(i) + w * px
               ry(i) = ry(i) + w * py

               ad(i) = zero
               bd(i) = -two * (a + c)
               cd(i) = two * w * a

               i = nx
               dxdeta = ( x(i,j+1)-x(i,j-1)) / (two * deta)
               dydeta = ( y(i,j+1)-y(i,j-1)) / (two * deta)
               dxdxi  = ( x(i-1,j)-x(i-1,j)) / (two * dxi)
               dydxi  = (-y(i-1,j)-y(i-1,j)) / (two * dxi)

               a = (dxdeta**2 + dydeta**2) / dxi**2
               b = (dxdxi*dxdeta + dydxi*dydeta) / ( two * dxi * deta )
               c = (dxdxi**2 + dydxi**2) / deta**2
               
               rx(i) = -two*(one-w)*(a + c) * x(i,j) + w*b*(x(i-1,j+1) - &
                       x(i-1,j-1) + x(i-1,j-1) - x(i-1,j+1)) -           &
                       w*c*(x(i,j+1) + x(i,j-1))

               ry(i) = -two*(one-w)*(a + c) * y(i,j) + w*b*(-y(i-1,j+1) + &
                       y(i-1,j-1) + y(i-1,j-1) - y(i-1,j+1)) -            &
                       w*c*(y(i,j+1) + y(i,j-1))

               jac = one / ( dxdxi * dydeta - dydxi * dxdeta )
               if (p(i) .gt. zero) then
                 dxdxi = (x(i,j)-x(i-1,j)) / dxi
                 dydxi = (y(i,j)-y(i-1,j)) / dxi
               else
                 dxdxi = (x(i-1,j)-x(i,j)) / dxi
                 dydxi = (-y(i-1,j)-y(i,j)) / dxi
               end if
               if (q(i) .gt. zero) then
                 dxdeta = (x(i,j)-x(i,j-1)) / deta
                 dydeta = (y(i,j)-y(i,j-1)) / deta
               else
                 dxdeta = (x(i,j+1)-x(i,j)) / deta
                 dydeta = (y(i,j+1)-y(i,j)) / deta
               end if

               px = -one/jac**2 * ( p(i) * exp( -(eta(j)-eta(1)) ) * dxdxi + &
                                    q(i) * exp( -(eta(j)-eta(1)) ) * dxdeta )
               py = -one/jac**2 * ( p(i) * exp( -(eta(j)-eta(1)) ) * dydxi + &
                                    q(i) * exp( -(eta(j)-eta(1)) ) * dydeta )

               rx(i) = rx(i) + w * px
               ry(i) = ry(i) + w * py

               ad(i) = two * w * a
               bd(i) = -two * (a + c)
               cd(i) = zero
            end if
            
            call tridag( ad, bd, cd, rx, xl, nx )

            if (ibc .eq. 1) then
               cd(1) = zero
               ad(nx) = zero
            end if

            call tridag( ad, bd, cd, ry, yl, nx )

            do i = 1, nx
               err = err + abs(x(i,j)-xl(i)) + abs(y(i,j)-yl(i))
               x(i,j) = xl(i)
               y(i,j) = yl(i)
            end do

         end do

         err = err / real( (nx-2)*(ny-2) )
         write(*,"(i5,1x,1pe13.6)") n, err
      if ( err .gt. eps .and. n .lt. nmax) goto 20

      return
      end

!=============================================================================c
      subroutine ortho (nx, ny, dxi, deta, x, y, dn, p, q)
!=============================================================================c
      implicit none

      integer, intent(in) :: nx, ny
      real, intent(in)    :: dxi, deta, x(nx,ny), y(nx,ny), dn(nx)
      real, intent(out)   :: p(nx), q(nx)

      real, parameter :: zero=0.0d0, pt5=0.5d0, one=1.0d0, onept5=1.5d0, &
                         two=2.0d0, twopt5=2.5d0, three=3.0d0,           &
                         pi=3.1415926535897932385d+0

      integer :: i
      real :: dxdxi(nx), dydxi(nx), d2xdxi2(nx), d2ydxi2(nx)
      real :: dxdeta(nx), dydeta(nx), d2xdeta2(nx), d2ydeta2(nx)
      real :: d2xdxideta(nx), d2ydxideta(nx)

      real :: jac, r1, r2, a, b, c
!=============================================================================c
      i = 1
      dxdxi(i) = ( x(i+1,1) - x(i+1,1) ) / (two * dxi)
      dydxi(i) = ( y(i+1,1) + y(i+1,1) ) / (two * dxi)
      d2xdxi2(i) = ( x(i+1,1) - two * x(i,1) + x(i+1,1) ) / dxi**2      
      d2ydxi2(i) = ( y(i+1,1) - two * y(i,1) - y(i+1,1) ) / dxi**2
      do i = 2, nx-1
        dxdxi(i) = ( x(i+1,1) - x(i-1,1) ) / (two * dxi)
        dydxi(i) = ( y(i+1,1) - y(i-1,1) ) / (two * dxi)
        d2xdxi2(i) = ( x(i+1,1) - two * x(i,1) + x(i-1,1) ) / dxi**2    
        d2ydxi2(i) = ( y(i+1,1) - two * y(i,1) + y(i-1,1) ) / dxi**2
      end do
      i = nx
      dxdxi(i) = ( x(i-1,1) - x(i-1,1) ) / (two * dxi)
      dydxi(i) = (-y(i-1,1) - y(i-1,1) ) / (two * dxi)
      d2xdxi2(i) = ( x(i-1,1) - two * x(i,1) + x(i-1,1) ) / dxi**2      
      d2ydxi2(i) = (-y(i-1,1) - two * y(i,1) + y(i-1,1) ) / dxi**2

      do i = 1, nx
        dxdeta(i) = dn(i)/deta * (-dydxi(i)) / sqrt(dxdxi(i)**2 + dydxi(i)**2)
        dydeta(i) = dn(i)/deta * ( dxdxi(i)) / sqrt(dxdxi(i)**2 + dydxi(i)**2)
      end do

      do i = 1, nx
        d2xdeta2(i) = (-7.0 * x(i,1) + 8.0 * x(i,2) - x(i,3)) / &
                      (2.0 * deta**2) - 3.0 * dxdeta(i) / deta
        d2ydeta2(i) = (-7.0 * y(i,1) + 8.0 * y(i,2) - y(i,3)) / &
                      (2.0 * deta**2) - 3.0 * dydeta(i) / deta
      end do

      i = 1
      d2xdxideta(i) = ( dxdeta(i+1) - dxdeta(i+1) ) / (two * dxi)
      d2ydxideta(i) = ( dydeta(i+1) + dydeta(i+1) ) / (two * dxi)
      do i = 2, nx-1
        d2xdxideta(i) = ( dxdeta(i+1) - dxdeta(i-1) ) / (two * dxi)
        d2ydxideta(i) = ( dydeta(i+1) - dydeta(i-1) ) / (two * dxi)
      end do
      i = nx
      d2xdxideta(i) = ( dxdeta(i-1) - dxdeta(i-1) ) / (two * dxi)
      d2ydxideta(i) = (-dydeta(i-1) - dydeta(i-1) ) / (two * dxi)

      do i = 1, nx
        a = dxdeta(i)**2 + dydeta(i)**2
        b = dxdxi(i)*dxdeta(i) + dydxi(i)*dydeta(i)
        c = dxdxi(i)**2 + dydxi(i)**2

        jac = one / ( dxdxi(i) * dydeta(i) - dydxi(i) * dxdeta(i) )

        r1 = -(jac**2) * ( a * d2xdxi2(i) - two * b * d2xdxideta(i) + &
                           c * d2xdeta2(i) )

        r2 = -(jac**2) * ( a * d2ydxi2(i) - two * b * d2ydxideta(i) + &
                           c * d2ydeta2(i) )

        p(i) = jac * ( dydeta(i) * r1 - dxdeta(i) * r2 )
        q(i) = jac * ( -dydxi(i) * r1 + dxdxi(i) * r2 )
      end do

      return
      end

!=============================================================================c
      SUBROUTINE TRIDAG(A,B,C,R,U,N)
!=============================================================================c
      DIMENSION GAM(N),A(N),B(N),C(N),R(N),U(N)
      IF(B(1).EQ.0.)THEN
        WRITE(*,*) 'ERROR IN TRIDAG:  B(1) = 0'
        CALL EXIT(1)
      END IF
      BET=B(1)
      U(1)=R(1)/BET
      DO 11 J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF(BET.EQ.0.)THEN
          WRITE(*,*) 'ERROR IN TRIDAG: B = 0', J
          CALL EXIT(1)
        END IF
        U(J)=(R(J)-A(J)*U(J-1))/BET
11    CONTINUE
      DO 12 J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
12    CONTINUE
      RETURN
      END
