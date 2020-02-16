c=============================================================================c
      subroutine elliptic( nx, ny, x, ldx, y, ldy, eps, nmax )
c
c     This routine performs an elliptic smoothing of the input grid
c     
c     Author:  Scott Collis
c
c     Date:    6-17-97
c
c=============================================================================c
      implicit none

      integer nx, ny, ldx, ldy, nmax
      real x(ldx,ny), y(ldy,ny), eps

      real zero, pt5, one, two
      parameter (zero=0.0, pt5=0.5, one=1.0, two=2.0)

      integer n, i, j
      real dxi, deta, xeta, yeta, xxi, yxi, a, b, c, xnew, ynew, err

      integer mx
      parameter (mx=1000)
      real w, rx(mx), ry(mx), ad(mx), bd(mx), cd(mx), xl(mx), yl(mx)
c=============================================================================c
      if (nx.gt.mx) then
         write(*,*) 'Increase mx in elliptic'
         stop
      end if

      dxi = one / float(nx-1)
      deta = one / float(ny-1)
c
c.... Gauss-Siedel
c
      if (.false.) then

      n = 0
 10   continue
         n = n + 1
         err = zero
         do i = 2, nx-1
            do j = 2, ny-1
               xeta = (x(i,j+1)-x(i,j-1)) / (two * deta)
               yeta = (y(i,j+1)-y(i,j-1)) / (two * deta)
               xxi  = (x(i+1,j)-x(i-1,j)) / (two * dxi)
               yxi  = (y(i+1,j)-y(i-1,j)) / (two * dxi)

               a = (xeta**2 + yeta**2) / dxi**2
               b = (xxi*xeta + yxi*yeta) / ( two * dxi * deta )
               c = (xxi**2 + yxi**2) / deta**2
               
               xnew = (a * (x(i+1,j) + x(i-1,j)) + c * (x(i,j+1) + x(i,j-1))
     &                - b * (x(i+1,j+1) - x(i+1,j-1) + x(i-1,j-1) -
     &                x(i-1,j+1)) ) / ( two * (a + c) )
               ynew = (a * (y(i+1,j) + y(i-1,j)) + c * (y(i,j+1) + y(i,j-1))
     &                - b * (y(i+1,j+1) - y(i+1,j-1) + y(i-1,j-1) -
     &                y(i-1,j+1)) ) / ( two * (a + c) )

               err = err + abs(x(i,j)-xnew) + abs(y(i,j)-ynew)

               x(i,j) = xnew
               y(i,j) = ynew
            end do
         end do
         err = err / real( (nx-2)*(ny-2) )
         write(*,"(i5,1x,1pe13.6)") n, err
      if ( err .gt. eps .and. n .lt. nmax) goto 10
c
c.... Line Gauss-Siedel with SOR (because of the nonlinearity, SOR is unstable
c.... for w > 1.
c
      else

      write(*,"(/,'Enter the SOR parmameter (0 < w < 2, w = 1.78) ==>',$)")
      read(*,*) w
      n = 0
 20   continue
         n = n + 1
         err = zero
         do j = 2, ny-1
            do i = 2, nx-1
               xeta = (x(i,j+1)-x(i,j-1)) / (two * deta)
               yeta = (y(i,j+1)-y(i,j-1)) / (two * deta)
               xxi  = (x(i+1,j)-x(i-1,j)) / (two * dxi)
               yxi  = (y(i+1,j)-y(i-1,j)) / (two * dxi)

               a = (xeta**2 + yeta**2) / dxi**2
               b = (xxi*xeta + yxi*yeta) / ( two * dxi * deta )
               c = (xxi**2 + yxi**2) / deta**2
               
               rx(i) = -two*(one-w)*(a + c) * x(i,j) + w*b*(x(i+1,j+1) -
     &                 x(i+1,j-1) + x(i-1,j-1) - x(i-1,j+1)) -
     &                 w*c*(x(i,j+1) + x(i,j-1))

               ry(i) = -two*(one-w)*(a + c) * y(i,j) + w*b*(y(i+1,j+1) -
     &                 y(i+1,j-1) + y(i-1,j-1) - y(i-1,j+1)) -
     &                 w*c*(y(i,j+1) + y(i,j-1))

               ad(i) = w * a
               bd(i) = -two * (a + c)
               cd(i) = w * a
            end do
c
c.... apply the boundary conditions
c
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
            
            call tridag( ad, bd, cd, rx, xl, nx )
            call tridag( ad, bd, cd, ry, yl, nx )
            
            do i = 2, nx-1
               err = err + abs(x(i,j)-xl(i)) + abs(y(i,j)-yl(i))
               x(i,j) = xl(i)
               y(i,j) = yl(i)
            end do
         end do

         err = err / real( (nx-2)*(ny-2) )
         write(*,"(i5,1x,1pe13.6)") n, err
      if ( err .gt. eps .and. n .lt. nmax) goto 20

      end if

      return
      end

      SUBROUTINE TRIDAG(A,B,C,R,U,N)
      DIMENSION A(N),B(N),C(N),R(N),U(N)
      INTEGER NMAX
      PARAMETER (NMAX=1000)
      DIMENSION GAM(NMAX)
      IF(N.GT.NMAX)PAUSE
      IF(B(1).EQ.0.)PAUSE
      BET=B(1)
      U(1)=R(1)/BET
      DO 11 J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF(BET.EQ.0.)PAUSE
        U(J)=(R(J)-A(J)*U(J-1))/BET
11    CONTINUE
      DO 12 J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
12    CONTINUE
      RETURN
      END
