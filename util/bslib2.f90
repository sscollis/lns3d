!
!.... Driver routines that mimic the IMSL Bspline routine
!
!.... Author:  S. Scott Collis
!
!.... Copyright:  (c)2020 S. Scott Collis
!
!.... Be sure to compile with -r8 otherwise you'll get roundoff errors for
!.... large nx
!
      subroutine BSNAK( nx, x, korder, knot )
!
!.... Make a Not-a-Knot sequence.  Modeled on IMSL routine BSNAK
!
!.... S. Collis
!
      integer nx, korder
      real  x(nx), knot(nx+korder)

      do i = 1, korder
         knot(i) = x(1)
         knot(nx+i) = x(nx)
      end do
      if (mod(korder,2).eq.0) then
         do i = korder+1, nx
            knot(i) = x(i-korder/2)
         end do
      else
         do i = korder+1, nx
            knot(i) = 0.5 * ( x(i-(korder+1)+korder/2 + 1) + &
                              x(i-(korder+1)+korder/2 + 2) )
         end do
      end if

      return
      end

      subroutine BSINT( nx, x, u, korder, knot, bs )
      implicit none
      integer :: nx, korder
      real :: x(nx), u(nx), knot(nx+korder), bs(nx) 
      real :: work(nx*(2*korder-1)), work2(2*korder)
      integer :: iflag
#if 1 
      call BINTK( x, u, knot, nx, korder, bs, work, work2 )
#else
      call splint( x, u, knot, nx, korder, work, bs, iflag )
      if (iflag.ne.1) then
        write(*,*) 'Error in BSINT'
        stop
      end if
#endif
      return
      end

      subroutine BS2IN( nx, x, ny, y, u, ldu, kxord, kyord, xknot, yknot, bs )
      implicit none
      integer :: nx, ny, kxord, kyord, ldu
      real :: x(nx), y(ny), u(nx,ny), xknot(nx+kxord), yknot(ny+kyord), &
              bs(nx,ny), ts(nx,ny)
      real :: xwork(nx*(2*kxord-1)), xwork2(2*kxord)
      real :: ywork(ny*(2*kyord-1)), ywork2(2*kyord)
      integer :: iflag, i, j 
#if 0
      do j = 1, ny
        call BINTK( x, u(:,j), xknot, nx, kxord, ts(:,j), xwork, xwork2 )
      end do
      do i = 1, nx
        call BINTK( y, ts(i,:), yknot, ny, kyord, bs(i,:), ywork, ywork2 )
      end do
#else
      do i = 1, nx
        call BINTK( y, u(i,:), yknot, ny, kyord, ts(i,:), ywork, ywork2 )
      end do
      do j = 1, ny
        call BINTK( x, ts(:,j), xknot, nx, kxord, bs(:,j), xwork, xwork2 )
      end do
#endif
      return
      end

      subroutine BS2GD( ider, jder, nx, x, ny, y, kxord, kyord, xknot, &
                        yknot, nx1, ny1, bs, v, ldv )
      implicit none
      integer :: ider, jder, kxord, kyord, nx, ny, nx1, ny1, ldv, i, j
      real :: x(nx), y(ny), v(nx,ny)
      real :: xknot(nx1+kxord), yknot(ny1+kyord), bs(nx1,ny1)
      real :: BS2DR
      do j=1,ny
        do i=1,nx
          v(i,j) = BS2DR(ider, jder, x(i), y(j), kxord, kyord, xknot, &
                         yknot, nx1, ny1, bs)
         end do
      end do
      return
      end subroutine BS2GD

      real function BS2DR( ider, jder, x, y, kxord, kyord, xknot, yknot, &
                           nx, ny, bs )
      implicit none
      integer :: ider, jder, kxord, kyord, nx, ny
      real :: x, y, xknot(nx+kxord), yknot(ny+kyord), bs(nx,ny)
      real, external :: BVALU
      real :: xwork(3*kxord), ywork(3*kyord)
      integer :: inbv=1, i, j
#if 0
      real :: tmp(nx)
      do i = 1, nx
        tmp(i) = BVALU( yknot, bs(i,:), ny, kyord, jder, y, inbv, ywork )
      end do
      BS2DR = BVALU( xknot, tmp, nx, kxord, ider, x, inbv, xwork )
#else
      real :: tmp(ny)
      do j = 1, ny
        tmp(j) = BVALU( xknot, bs(:,j), nx, kxord, ider, x, inbv, xwork )
      end do
      BS2DR = BVALU( yknot, tmp, ny, kyord, jder, y, inbv, ywork )
#endif
      return
      end

      real function BSDER( ider, x, korder, knot, nx, bs )
      implicit none
      integer :: ider, korder, nx
      real :: x, knot(nx+korder), bs(nx)
      real, external :: bvalue, BVALU
      real :: work(3*korder)
      integer :: inbv=1

      BSDER = BVALU( knot, bs, nx, korder, ider, x, inbv, work)

      return
      end

      real function BSVAL( x, korder, knot, nx, bs )
      implicit none
      integer :: ider=1, korder, nx
      real :: x, knot(nx+korder), bs(nx)
      real, external :: bvalue, BVALU
      real :: work(3*korder)
      integer :: inbv=1

      BSVAL = BVALU( knot, bs, nx, korder, ider, x, inbv, work)
!     BSVAL = BVALUE( knot, bs, nx, korder, x, ider )

      return
      end

      real function BSITG( x1, x2, korder, knot, nx, bs )
      implicit none
      integer :: korder, nx
      real :: x1, x2, knot(nx+korder), bs(nx)
      real :: work(3*korder)
      integer :: inbv=1

      call BSQAD( knot, bs, nx, korder, x1, x2, BSITG, work)

      return
      end

      subroutine QDAG( fun, x1, x2, err1, err2, idum, val, errest)
      implicit none
      external fun
      real :: x1, x2, err1, err2, errest, val
      integer :: idum, ierr

      errest = min(err1,err2)
      call GAUS8( fun, x1, x2, errest, val, ierr)

      return
      end
