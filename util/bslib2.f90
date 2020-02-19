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
      
      call BINTK( x, u, knot, nx, korder, bs, work, work2 )

!      call splint( x, u, knot, nx, korder, work, bs, iflag )
!      if (iflag.ne.1) then
!              write(*,*) 'Error in BSINT'
!              stop
!      end if

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

!     BSDER = bvalue( knot, bs, nx, korder, x, ider )

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

!     BSVAL = bvalue( knot, bs, nx, korder, x, ider )

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
