!=============================================================================!
	subroutine rhsBC(rl, vl, vml) 
!  
!  Satisfy the boundary conditions 
!  
!=============================================================================!
	use global
	use stencil
	use pot
	implicit none

	real :: rl(ny, nx, ndof), vl(ny, nx, ndof), vml(ny, nx, ndof)
	
        real, allocatable :: p(:), pnorm(:)

	integer :: i, j
!=============================================================================!
!	L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
	if (linear.eq.1) then

	if (yper) then
	
	  rl(ny,:,:) = zero
	
	else			! yper

!=============================================================================!
!	W a l l
!=============================================================================!
	if (Navier) then

!.... density boundary condition

	  if (wall.eq.1) then
	    rl(1,:,1) = gc1 * vl(1,:,1) + gc2 * vl(2,:,1) + &
	                gc3 * vl(3,:,1) + gc4 * vl(4,:,1) + &
		        gc5 * vl(5,:,1)
	  end if
	  
!.... linear extrapolation of rho at the wall

	  if (wall.eq.2) then
	    rl(1,:,1) = vl(1,:,1) - two * vl(2,:,1) + vl(3,:,1)
	  end if

!.... no-slip boundary condition

	  rl(1,:,2:ndof-1)  = zero

!.... isothermal wall

	  if (wallt.eq.0) rl(1,:,ndof) = zero

!.... adiabatic boundary condition

	  if (wallt.eq.1) then
	    rl(1,:,ndof) = gc1 * vl(1,:,ndof) + gc2 * vl(2,:,ndof) + &
	                   gc3 * vl(3,:,ndof) + gc4 * vl(4,:,ndof) + &
		           gc5 * vl(5,:,ndof)
	  end if
	
	else		! inviscid wall
	
	  allocate( p(nx), pnorm(nx) )
	  call lwallbc(vl, vml, p, pnorm)
	  do i = 1, nx
	    rl(1,i,1) = ( gc1 * vl(1,i,1) + gc2 * vl(2,i,1) + &
			  gc3 * vl(3,i,1) + gc4 * vl(4,i,1) + &
			  gc5 * vl(5,i,1) ) / deta - Pnorm(i)
	  end do
	  deallocate( p, pnorm )

	  do i = 1, nx
	    if ( abs(bnb(i,2)) .ge. abs(bnb(i,1)) ) then
	      rl(1,i,2) = bnb(i,2) * rl(1,i,2) - bnb(i,1) * rl(1,i,3)
	      rl(1,i,3) = -( bnb(i,1) * vl(1,i,2) + bnb(i,2) * vl(1,i,3) )
	    else
	      rl(1,i,3) = bnb(i,2) * rl(1,i,2) - bnb(i,1) * rl(1,i,3)
	      rl(1,i,2) = -( bnb(i,1) * vl(1,i,2) + bnb(i,2) * vl(1,i,3) )
	    end if
	  end do

	end if		! Navier
!=============================================================================!
!	T o p
!=============================================================================!

!.... freestream zero disturbance boundary conditions

	  if (top.eq.0) then
	    rl(ny,:,:) = zero
	  end if
	
	end if			! yper

	if (xper) then
	  rl(:,nx,:) = zero
	else			! xper
	
!=============================================================================!
!	L e f t   B o u n d a r y
!=============================================================================!
        if (left.eq.0) then

          rl(:,1,:)  = zero

        else if (left.eq.1) then

	  rl(nbl+1:ny,1,1) = vl(nbl+1:ny,1,1) - &
	          (gamma * Ma * vml(nbl+1:ny,1,1) * &
	          sqrt(vml(nbl+1:ny,1,5)) * vl(nbl+1:ny,1,2) - &
	          vml(nbl+1:ny,1,1) * vl(nbl+1:ny,1,5)) / vml(nbl+1:ny,1,5)

!	  rl(1:nbl,1,:) = vl(1:nbl,1,:) - vl(1:nbl,2,:)

	else if (left.eq.2) then

	  do j = 1, ny
!	    rl(j,1,1) = gc1 * vl(j,1,1) + gc2 * vl(j,2,1)  + &
!			gc3 * vl(j,3,1) + gc4 * vl(j,4,1)  + &
!			gc5 * vl(j,5,1)

	    rl(j,1,2) = gc1 * vl(j,1,2) + gc2 * vl(j,2,2)  + &
			gc3 * vl(j,3,2) + gc4 * vl(j,4,2)  + &
			gc5 * vl(j,5,2)
			    
	    rl(j,1,3) = zero

	    rl(j,1,4) = gc1 * vl(j,1,4) + gc2 * vl(j,2,4)  + &
			gc3 * vl(j,3,4) + gc4 * vl(j,4,4)  + &
			gc5 * vl(j,5,4)

	    rl(j,1,5) = gc1 * vl(j,1,5) + gc2 * vl(j,2,5)  + &
			gc3 * vl(j,3,5) + gc4 * vl(j,4,5)  + &
			gc5 * vl(j,5,5)
	  end do

	else if (left.eq.7) then           ! symmetry boundary
	  rl(:,1,3) = zero
	end if

!=============================================================================!
!	R i g h t   B o u n d a r y
!=============================================================================!
        if (right.eq.0) then

          rl(:,nx,:) = zero

        else if (right.eq.1) then

	  rl(nbl+1:ny,nx,1) = vl(nbl+1:ny,nx,1) - &
	          (gamma * Ma * vml(nbl+1:ny,nx,1) * &
	          sqrt(vml(nbl+1:ny,nx,5)) * vl(nbl+1:ny,nx,2) - &
	          vml(nbl+1:ny,nx,1) * vl(nbl+1:ny,nx,5)) /vml(nbl+1:ny,nx,5)

!	  rl(1:nbl,nx,:) = vl(1:nbl,nx,:) - vl(1:nbl,nx-1,:)

	end if

	end if			! xper
			
!=============================================================================!
!	N O N L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
	else			! linear = 0	

	if (yper) then
	
	  rl(ny,:,:) = zero
	
	else			! yper

!=============================================================================!
!	W a l l
!=============================================================================!
	if (Navier) then

!.... density boundary condition

	  if (wall.eq.1) then
	    rl(1,:,1) = gc1 * vl(1,:,1) + gc2 * vl(2,:,1) + &
	                gc3 * vl(3,:,1) + gc4 * vl(4,:,1) + &
		        gc5 * vl(5,:,1) 
	  end if

!.... no-slip boundary condition

	  rl(1,:,2:ndof-1) = zero

!.... isothermal wall

	  if (wallt.eq.0) then
	    rl(1,:,ndof) = vl(1,:,ndof) - ( one + pt5 * gamma1 * Ma**2 * &
                           (one + tan(alpha)**2) * sqrt(Pr) )
	  end if

!.... adiabatic boundary condition

	  if (wallt.eq.1) then
	    rl(1,:,ndof) = gc1 * vl(1,:,ndof) + gc2 * vl(2,:,ndof) + &
			   gc3 * vl(3,:,ndof) + gc4 * vl(4,:,ndof) + &
			   gc5 * vl(5,:,ndof)
	  end if

	  if (wall.eq.2) then
	    allocate( p(nx), pnorm(nx) )
	    call wallbc(vl, p, pnorm)
	    do i = 1, nx
	      rl(1,i,1) = ( gc1 * vl(1,i,1) + gc2 * vl(2,i,1) + &
			    gc3 * vl(3,i,1) + gc4 * vl(4,i,1) + &
			    gc5 * vl(5,i,1) ) / deta -          &
			    gamma * Ma**2 * Pnorm(i) / vl(1,i,5)
	    end do
	    deallocate( p, pnorm )
	  end if
	  
	else    ! inviscid
	
!....     rl(1,:,2) = wall tangent momentum
!....     rl(1,:,3) = wall normal velocity
	
	  rl(1,:,2) =  ( bnb(:,2) * rl(1,:,2) - bnb(:,1) * rl(1,:,3) )	  
	  rl(1,:,3) = -( bnb(:,1) * vl(1,:,2) + bnb(:,2) * vl(1,:,3) )

	end if	! Navier
!=============================================================================!
!	T o p
!=============================================================================!
	if (Ma.gt.one) then

	  rl(ny,:,1) = zero		! set all variables
	  rl(ny,:,2) = zero
	  rl(ny,:,3) = zero
	  rl(ny,:,4) = zero
	  rl(ny,:,5) = zero

	else	! Ma

	  if (top.eq.0) then
	    call ReimannRHS( nx, vl(ny,:,:), vl(ny-1,:,:), vl(ny-2,:,:), &
                             bnt, x(ny,:), y(ny,:), rl(ny,:,:), &
			     rhobt, ubt, vbt, wbt, tbt, pbt, cbt )
	  end if
	
	end if 			! Ma
	
	end if			! yper
	
!=============================================================================!

	if (xper) then
	
	  rl(:,nx,:) = zero
	
	else			! xper
	
!=============================================================================!
!	L e f t  B o u n d a r y
!=============================================================================!

	if (Ma.lt.one) then

!.... first-order Riemann

	  if (left.eq.0) then
	    call ReimannRHS( ny, vl(:,1,:), vl(:,2,:), vl(:,3,:), &
                             bnl, x(:,1), y(:,1), rl(:,1,:), &
			     rhobl, ubl, vbl, wbl, tbl, pbl, cbl )
	  end if

!.... Symmetry plane

	  if (left.eq.2) then
	    do j = 1, ny
	      rl(j,1,1) = gc1 * vl(j,1,1) + gc2 * vl(j,2,1)  + &
			  gc3 * vl(j,3,1) + gc4 * vl(j,4,1)  + &
			  gc5 * vl(j,5,1)
  
	      rl(j,1,2) = gc1 * vl(j,1,2) + gc2 * vl(j,2,2)  + &
			  gc3 * vl(j,3,2) + gc4 * vl(j,4,2)  + &
			  gc5 * vl(j,5,2)
			      
	      rl(j,1,3) = zero
	      rl(j,1,4) = gc1 * vl(j,1,4) + gc2 * vl(j,2,4)  + &
			  gc3 * vl(j,3,4) + gc4 * vl(j,4,4)  + &
			  gc5 * vl(j,5,4)
	      rl(j,1,5) = gc1 * vl(j,1,5) + gc2 * vl(j,2,5)  + &
			  gc3 * vl(j,3,5) + gc4 * vl(j,4,5)  + &
			  gc5 * vl(j,5,5)
	    end do
	  end if

!.... Symmetry boundary

	  if (left.eq.7) then
	    rl(:,1,3) = zero
	  end if
	  
	end if 	! Ma

	if (.false.) then
	  
!.... Zero'th order extrapolation in the viscous layers

	  if (extrap.eq.0) then
	    rl(1:nbl,1,1) = vl(1:nbl,1,1) - vl(1:nbl,2,1)
	    rl(1:nbl,1,2) = vl(1:nbl,1,2) - vl(1:nbl,2,2)
	    rl(1:nbl,1,3) = vl(1:nbl,1,3) - vl(1:nbl,2,3)
	    rl(1:nbl,1,4) = vl(1:nbl,1,4) - vl(1:nbl,2,4)
	    rl(1:nbl,1,5) = vl(1:nbl,1,5) - vl(1:nbl,2,5)
	  
!.... First-order extrapolation in the viscous layers

	  else if (extrap.eq.1) then
	    rl(1:nbl,1,:) = vl(1:nbl,1,:) - &
			    ( two * vl(1:nbl,2,:) - vl(1:nbl,3,:) )
	  
!.... Second-order extrapolation in the viscous layers (no good)

	  else if (extrap.eq.2) then
	    rl(1:nbl,1,:) = vl(1:nbl,1,:) - &
			    ( three * vl(1:nbl,2,:) - three * vl(1:nbl,3,:) + &
			      vl(1:nbl,4,:) )
	  end if
	  
	end if
!=============================================================================!
!	R i g h t  B o u n d a r y
!=============================================================================!

	if (Ma.lt.one) then

!.... First-order Riemann

	  if (right.eq.0) then
	    call ReimannRHS( ny, vl(:,nx,:), vl(:,nx-1,:), vl(:,nx-2,:), &
                             bnr, x(:,nx), y(:,nx), rl(:,nx,:), &
			     rhobr, ubr, vbr, wbr, tbr, pbr, cbr )
	  end if

!.... Symmetry boundary

          if (right.eq.7) then
            rl(:,nx,3) = zero
          end if

	end if 		! Ma

!.... Zero'th order extrapolation in the viscous layers

	if (extrap.eq.0) then
	  rl(1:nbl,nx,:) = vl(1:nbl,nx,:) - vl(1:nbl,nx-1,:)
	  
!.... First-order extrapolation in the viscous layers

	else if (extrap.eq.1) then
	  rl(1:nbl,nx,:) = vl(1:nbl,nx,:) - &
			  ( two * vl(1:nbl,nx-1,:) - vl(1:nbl,nx-2,:) )
	end if
	
	if (right.eq.8) then		! hold initial condition
	  rl(:,nx,:) = zero
	end if

	end if			! xper
	
	end if			! linear

	return
	end
