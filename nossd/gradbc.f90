!=============================================================================!
	subroutine gradbc(g1v, g2v, g11v, g12v, g22v) 
!  
!       Enforce the gradient boundary conditions 
!
!=============================================================================!
	use global
	implicit none

	real :: g1v(ny,nx,ndof),  g2v(ny,nx,ndof)
	real :: g11v(ny,nx,ndof), g12v(ny,nx,ndof), g22v(ny,nx,ndof)
	
!.... enforce adiabatic wall BC on the mean field

	if ( wallt.eq.1 .or. wallt.eq.2 ) g2v(1,:,ndof)  = zero

	return
	end
