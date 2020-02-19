!=============================================================================!
module fmax

      private :: fun
      integer, private :: n, kord = 5
      real, allocatable, private :: knot(:), bs(:)
      
      real, external :: BSDER
        
contains

!=============================================================================!
      function getval( nl, x, f, xmax)
!=============================================================================!
            implicit none
        
            integer :: nl, i, u, imax
            real :: getval, x(nl), f(nl), xmax
!=============================================================================!
            n = nl
            
            allocate( knot(n+kord), bs(n) )
            call BSNAK( n, x, kord, knot )
            call BSINT( n, x, f, kord, knot, bs )
            
            getval = BSDER( 0, xmax, kord, knot, n, bs )
            
            deallocate( knot, bs )
            
            return
      end function getval

!=============================================================================!
      subroutine fun( x, g, d )
!=============================================================================!
            implicit none
            
            real :: x, f, g, d
!=============================================================================!
            
!.... maximum value of function f
            
            f = BSDER( 0, x, kord, knot, n, bs )
            g = BSDER( 1, x, kord, knot, n, bs )
            d = BSDER( 2, x, kord, knot, n, bs )
            
            return
      end subroutine fun

!=============================================================================!
      function findmax( nl, x, f, xmax)
!=============================================================================!
            implicit none
            
            integer :: nl, i, u, imax
            real :: findmax, x(nl), f(nl), xmax
            
            real, external :: RTSAFE
!=============================================================================!
            n = nl
            
            allocate( knot(n+kord), bs(n) )
            call BSNAK( n, x, kord, knot )
            call BSINT( n, x, f, kord, knot, bs )
            
            do i = 1, n-1
              if ( f(i+1) .lt. f(i) ) goto 10
            end do
            write(*,*) 'Error in findmax'
            findmax = 1.0
            goto 20
            
10          continue
            imax = i
            xmax = RTSAFE( fun, x(i-2), x(i+2), 1.0e-14 )
            findmax = BSDER( 0, xmax, kord, knot, n, bs )
            
20          deallocate( knot, bs )
            
            return
      end function findmax

end module fmax

