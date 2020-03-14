!=============================================================================!
	subroutine lhsbc2f3d( mat, vl, vml)
!  
!  Correct the LHS for boundary conditions in the \eta direction
!
!  This version supports the 4th order LHS
!
!  Revised: 10-17-95
!=============================================================================!
	use global
	use stencil
	implicit none
	
	complex :: mat(ny,nx,ndof,ndof,5), vl(ny,nx,ndof)
	real    :: vml(ny,nx,ndof)

	complex :: matl(5,2)
	
	integer :: i, j, ij
!=============================================================================!
!	L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
	if (.not. yper) then

	if (Navier) then
	
!.... wall boundary condition

	if (wall.eq.1) then
	  mat(1,:,1,:,:) = zero
	  mat(1,:,1,1,3) = -gc1 
	  mat(1,:,1,1,4) = -gc2
	  mat(1,:,1,1,5) = -gc3
	  mat(1,:,1,1,1) = -gc4
	  mat(1,:,1,1,2) = -gc5
	end if

!.... no-slip

	mat(1,:,2:4,:,:) = zero
	mat(1,:,2,2,3)   = one
	mat(1,:,3,3,3)   = one
	mat(1,:,4,4,3)   = one

!.... isothermal

	if (wallt.eq.0) then
	  mat(1,:,5,:,:) = zero
	  mat(1,:,5,5,3) = one
	end if

!.... adiabatic

	if (wallt.eq.1) then
	  mat(1,:,5,:,:) = zero
	  mat(1,:,5,5,3) = -gc1
	  mat(1,:,5,5,4) = -gc2
	  mat(1,:,5,5,5) = -gc3
	  mat(1,:,5,5,1) = -gc4
	  mat(1,:,5,5,2) = -gc5
	end if

	else		! inviscid:  rotate to body normal coordinates
	
	  do i = 1, nx
	    matl(1,1) = mat(1,i,1,2,3) * bnb(i,2) - mat(1,i,1,3,3) * bnb(i,1)
	    matl(1,2) = mat(1,i,1,2,3) * bnb(i,1) + mat(1,i,1,3,3) * bnb(i,2)
	    matl(2,1) = mat(1,i,2,2,3) * bnb(i,2) - mat(1,i,2,3,3) * bnb(i,1)
	    matl(2,2) = mat(1,i,2,2,3) * bnb(i,1) + mat(1,i,2,3,3) * bnb(i,2)
	    matl(3,1) = mat(1,i,3,2,3) * bnb(i,2) - mat(1,i,3,3,3) * bnb(i,1)
	    matl(3,2) = mat(1,i,3,2,3) * bnb(i,1) + mat(1,i,3,3,3) * bnb(i,2)
	    matl(4,1) = mat(1,i,4,2,3) * bnb(i,2) - mat(1,i,4,3,3) * bnb(i,1)
	    matl(4,2) = mat(1,i,4,2,3) * bnb(i,1) + mat(1,i,4,3,3) * bnb(i,2)
	    matl(5,1) = mat(1,i,5,2,3) * bnb(i,2) - mat(1,i,5,3,3) * bnb(i,1)
	    matl(5,2) = mat(1,i,5,2,3) * bnb(i,1) + mat(1,i,5,3,3) * bnb(i,2)
	    mat(1,i,:,2:3,3) = matl
	    
	    matl(1,1) = mat(2,i,1,2,2) * bnb(i,2) - mat(2,i,1,3,2) * bnb(i,1)
	    matl(1,2) = mat(2,i,1,2,2) * bnb(i,1) + mat(2,i,1,3,2) * bnb(i,2)
	    matl(2,1) = mat(2,i,2,2,2) * bnb(i,2) - mat(2,i,2,3,2) * bnb(i,1)
	    matl(2,2) = mat(2,i,2,2,2) * bnb(i,1) + mat(2,i,2,3,2) * bnb(i,2)
	    matl(3,1) = mat(2,i,3,2,2) * bnb(i,2) - mat(2,i,3,3,2) * bnb(i,1)
	    matl(3,2) = mat(2,i,3,2,2) * bnb(i,1) + mat(2,i,3,3,2) * bnb(i,2)
	    matl(4,1) = mat(2,i,4,2,2) * bnb(i,2) - mat(2,i,4,3,2) * bnb(i,1)
	    matl(4,2) = mat(2,i,4,2,2) * bnb(i,1) + mat(2,i,4,3,2) * bnb(i,2)
	    matl(5,1) = mat(2,i,5,2,2) * bnb(i,2) - mat(2,i,5,3,2) * bnb(i,1)
	    matl(5,2) = mat(2,i,5,2,2) * bnb(i,1) + mat(2,i,5,3,2) * bnb(i,2)
	    mat(2,i,:,2:3,2) = matl
	    
	    matl(1,1) = mat(3,i,1,2,1) * bnb(i,2) - mat(3,i,1,3,1) * bnb(i,1)
	    matl(1,2) = mat(3,i,1,2,1) * bnb(i,1) + mat(3,i,1,3,1) * bnb(i,2)
	    matl(2,1) = mat(3,i,2,2,1) * bnb(i,2) - mat(3,i,2,3,1) * bnb(i,1)
	    matl(2,2) = mat(3,i,2,2,1) * bnb(i,1) + mat(3,i,2,3,1) * bnb(i,2)
	    matl(3,1) = mat(3,i,3,2,1) * bnb(i,2) - mat(3,i,3,3,1) * bnb(i,1)
	    matl(3,2) = mat(3,i,3,2,1) * bnb(i,1) + mat(3,i,3,3,1) * bnb(i,2)
	    matl(4,1) = mat(3,i,4,2,1) * bnb(i,2) - mat(3,i,4,3,1) * bnb(i,1)
	    matl(4,2) = mat(3,i,4,2,1) * bnb(i,1) + mat(3,i,4,3,1) * bnb(i,2)
	    matl(5,1) = mat(3,i,5,2,1) * bnb(i,2) - mat(3,i,5,3,1) * bnb(i,1)
	    matl(5,2) = mat(3,i,5,2,1) * bnb(i,1) + mat(3,i,5,3,1) * bnb(i,2)
	    mat(3,i,:,2:3,1) = matl

	    mat(1,i,2,:,:) = bnb(i,2) * mat(1,i,2,:,:) - &
			     bnb(i,1) * mat(1,i,3,:,:)
	    mat(1,i,3,:,:) = bnb(i,1) * mat(1,i,2,:,:) + &
			     bnb(i,2) * mat(1,i,3,:,:)

!.... constrain the normal velocity to be zero

	    mat(1,i,3,:,:) = zero
	    mat(1,i,3,3,3) = one

	  end do

	end if 		! Navier

!.... freestream boundary conditions

	if (top.eq.0 .or. top.eq.1) then
	  mat(ny,:,:,:,:) = zero
	  mat(ny,:,1,1,3) = one
	  mat(ny,:,2,2,3) = one
	  mat(ny,:,3,3,3) = one
	  mat(ny,:,4,4,3) = one
	  mat(ny,:,5,5,3) = one
	end if
	
	end if 		! yper

!=============================================================================!

	if (.not. xper) then

!.... Left boundary

	  if (left.eq.0 .or. left.eq.4 .or. left.eq.5) then

	    mat(:,1,:,:,:) = zero
	    mat(:,1,1,1,3) = one
	    mat(:,1,2,2,3) = one
	    mat(:,1,3,3,3) = one
	    mat(:,1,4,4,3) = one
	    mat(:,1,5,5,3) = one

          end if

!.... Right boundary

          if (right.eq.0 .or. right.eq.4 .or. right.eq.5) then

	    mat(:,nx,:,:,:) = zero
	    mat(:,nx,1,1,3) = one
	    mat(:,nx,2,2,3) = one
	    mat(:,nx,3,3,3) = one
	    mat(:,nx,4,4,3) = one
	    mat(:,nx,5,5,3) = one
	
	  end if

	end if		! xper

	return
	end
