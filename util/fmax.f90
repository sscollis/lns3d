!=============================================================================!
module fmax

      private :: fun
      integer, private :: n, kord = 5
      real, allocatable, private :: knot(:), bs(:)
      real, external :: bsder
        
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
            
        getval = bsder( 0, xmax, kord, knot, n, bs )
          
        deallocate( knot, bs )
            
        return
      end function getval
!=============================================================================!
      function f_un(x)
        implicit none
        real :: x, f_un, f, d
        call fun(x,f,d)
        f_un = f
        return
      end function
!=============================================================================!
      subroutine fun( x, g, d )
!=============================================================================!
        implicit none
            
        real :: x, f, g, d
!=============================================================================!
            
!.... maximum value of function f
            
        f = bsder( 0, x, kord, knot, n, bs )
        g = bsder( 1, x, kord, knot, n, bs )
        d = bsder( 2, x, kord, knot, n, bs )
            
        return
      end subroutine fun

!=============================================================================!
      function findmax( nl, x, f, xmax)
!=============================================================================!
        implicit none
        
        integer :: nl, i, u, imax
        real :: findmax, x(nl), f(nl), xmax
#ifdef USE_NR    
        real, external :: rtsafe 
#else
        real, external :: zeroin 
#endif
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
            
10      continue
        imax = i
#ifdef USE_NR
        xmax = rtsafe( fun, x(i-2), x(i+2), 1.0e-14 )
#else
        xmax = zeroin( x(i-2), x(i+2), f_un, 1.0e-14 )
#endif
        findmax = bsder( 0, xmax, kord, knot, n, bs )
            
20      deallocate( knot, bs )
            
        return
     end function findmax

end module fmax
