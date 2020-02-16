!=============================================================================!
        subroutine cgradbc(g1v, g2v, g11v, g12v, g22v) 
!  
!       Enforce the gradient boundary conditions 
!
!=============================================================================!
        use global
        implicit none

        complex :: g1v(ndof,nx,ny),  g2v(ndof,nx,ny)
        complex :: g11v(ndof,nx,ny), g12v(ndof,nx,ny), g22v(ndof,nx,ny)
        
!.... enforce adiabatic wall BC on the mean field

        if ( wallt.eq.1 .or. wallt.eq.2 ) g2v(ndof,:,1)  = zero

        return
        end
