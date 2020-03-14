!=============================================================================!
	subroutine rhsBC3D(rl, vl, vml) 
!  
!  Satisfy the boundary conditions for the 3D residual
!  
!=============================================================================!
	use global
	use stencil
	use bump
	implicit none

	complex :: rl(ny, nx, ndof), vl(ny, nx, ndof)
	real    :: vml(ny, nx, ndof)
	
!.... local variables

	real :: tmp
	
	real, allocatable :: rhom(:), cm(:), tm(:), um(:)
	complex, allocatable :: c3(:)
	
	complex :: ub(nx)
	
	integer :: i, j

!.... forcing parameters

	real :: a, d, kk
!=============================================================================!
!	L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
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

	  if (wall.eq.3) then
	    rl(1,:,2) = -(vl(1,:,2) - u1w)
	    rl(1,:,3) = -(vl(1,:,3) - u2w)
	    rl(1,:,4) = -(vl(1,:,4) - u3w)
	  else
	    rl(1,:,2) = -(vl(1,:,2) - zero)
	    rl(1,:,3) = -(vl(1,:,3) - zero)
	    rl(1,:,4) = -(vl(1,:,4) - zero)
	  end if

!.... isothermal wall

	  if (wallt.eq.0) then
	    if (wall.eq.3) then
	      rl(1,:,ndof) = -(vl(1,:,ndof) - tw)
	    else
	      rl(1,:,ndof) = zero
	    end if
	  end if

!.... adiabatic boundary condition

	  if (wallt.eq.1) then
	    if (wall.eq.3) then
	      rl(1,:,ndof) = gc1 * vl(1,:,ndof) + gc2 * vl(2,:,ndof) + &
			     gc3 * vl(3,:,ndof) + gc4 * vl(4,:,ndof) + &
			     gc5 * vl(5,:,ndof) - twp
	    else
	      rl(1,:,ndof) = gc1 * vl(1,:,ndof) + gc2 * vl(2,:,ndof) + &
			     gc3 * vl(3,:,ndof) + gc4 * vl(4,:,ndof) + &
			     gc5 * vl(5,:,ndof)
	    end if
	  end if

	else	 ! inviscid
	
!....     rl(1,:,2) = wall tangent momentum
!....     rl(1,:,3) = wall normal velocity
	
	  rl(1,:,2) =  ( bnb(:,2) * rl(1,:,2) - bnb(:,1) * rl(1,:,3) )	  
	  rl(1,:,3) = -( bnb(:,1) * vl(1,:,2) + bnb(:,2) * vl(1,:,3) )

	end if 	 ! Navier
!=============================================================================!
!	T o p
!=============================================================================!

!.... freestream zero disturbance boundary conditions

	if (top.eq.0) then

	  rl(ny,:,:) = zero

	else if (top.eq.1) then

!.... compute the characteristic amplitudes on the top boundary
	  
	  allocate( rhom(nx), tm(nx), cm(nx), um(nx), c3(nx) )

	  rhom = vml(ny,:,1)
	  tm   = vml(ny,:,5)
	  cm   = sqrt( tm ) / Ma
	  um   = vml(ny,:,2)

	  d  = pt5 * ( onept33 + gamma1 / Pr ) / Re

	  do i = 1, nx
	    kk = omega / (cm(i)+um(i))
	    a  = omega**2 * d / (cm(i)+um(i))**3
	    c3(i) = exp( -a * (x(ny,i) - x0) ) * exp( im * kk * x(ny,i) )
!	    c3(i) = wamp(i) * exp( im * kk * x(ny,i) )
	  end do

	  rl(ny,:,1) = -(vl(ny,:,1) - pt5 * c3 / cm**2)
	  rl(ny,:,2) = -(vl(ny,:,2) - c3 * pt5 / ( rhom * cm ))
	  rl(ny,:,3) = -(vl(ny,:,3) - zero)
	  rl(ny,:,4) = -(vl(ny,:,4) - zero)
	  rl(ny,:,5) = -(vl(ny,:,5) - (gamma*Ma**2 * c3 * pt5 - &
                                       tm * pt5 * c3 / cm**2) / rhom)

	  deallocate( rhom, tm, cm, um, c3 )

	end if
	
	end if			! yper

	if (xper) then
	  rl(:,nx,:) = zero
	else			! xper
	
!=============================================================================!
!	L e f t   B o u n d a r y
!=============================================================================!
	if (left.eq.0) then		! zero disturbance

          rl(:,1,:) = zero

	else if (left.eq.4) then	! eigenfunction inflow
	
	  rl(:,1,:) = zero
	
	else if (left.eq.5) then	! acoustic wave inflow
	  
	  allocate( rhom(ny), tm(ny), cm(ny), um(ny), c3(ny) )

	  rhom = vml(:,1,1)
	  tm   = vml(:,1,5)
	  cm   = sqrt( tm ) / Ma
	  um   = vml(:,1,2)

	  do j = 1, ny
	    kk = omega / ( cm(j)+um(j) )
	    c3(j) = exp( im * kk * x(j,1) )
	  end do

	  rl(:,1,1) = -(vl(:,1,1) - pt5 * c3 / cm**2)
	  rl(:,1,2) = -(vl(:,1,2) - c3 * pt5 / ( rhom * cm ))
	  rl(:,1,3) = -(vl(:,1,3) - zero)
	  rl(:,1,4) = -(vl(:,1,4) - zero)
	  rl(:,1,5) = -(vl(:,1,5) - (gamma*Ma**2 * c3 * pt5 - &
                                     tm * pt5 * c3 / cm**2) / rhom)

	  deallocate( rhom, tm, cm, um, c3 )

	else if (left.eq.7) then           ! symmetry boundary

	  rl(:,1,3) = zero
	    
	end if

!=============================================================================!
!	R i g h t   B o u n d a r y
!=============================================================================!
        if (right.eq.0) then              ! zero disturbance

          rl(:,nx,:) = zero

	else if (right.eq.4) then
	
	  rl(:,nx,:) = zero

	else if (right.eq.5) then

!.... compute the characteristic amplitudes on the right boundary
	  
	  allocate( rhom(ny), tm(ny), cm(ny), um(ny), c3(ny) )

	  rhom = vml(:,nx,1)
	  tm   = vml(:,nx,5)
	  cm   = sqrt( tm ) / Ma
	  um   = vml(:,nx,2)

	  do j = 1, ny
	    kk = omega / ( cm(j)+um(j) )
	    c3(j) = exp( im * kk * x(j,nx) )
	  end do

	  rl(:,nx,1) = -(vl(:,nx,1) - pt5 * c3 / cm**2)
	  rl(:,nx,2) = -(vl(:,nx,2) - c3 * pt5 / ( rhom * cm ))
	  rl(:,nx,3) = -(vl(:,nx,3) - zero)
	  rl(:,nx,4) = -(vl(:,nx,4) - zero)
	  rl(:,nx,5) = -(vl(:,nx,5) - (gamma*Ma**2 * c3 * pt5 - &
                                     tm * pt5 * c3 / cm**2) / rhom)

	  deallocate( rhom, tm, cm, um, c3 )

	else if (right.eq.8) then        ! hold initial condition

          rl(:,nx,:) = zero

	else if (right .eq. 9) then      ! extrapolation 

	  rl(:,nx,:) = -( vl(:,nx,:) - exp( two * log(vl(:,nx-1,:)) - &
                                                  log(vl(:,nx-2,:)) ) )

	end if

	end if			          ! xper
			
	return
	end
