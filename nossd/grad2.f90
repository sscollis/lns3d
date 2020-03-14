!=============================================================================!
	subroutine grad2 (ndof, nx, ny, v, g1v, g11v, g12v, g22v, dx, dy, &
	                  optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)
!
!  Take the second derivative of a 2-D field.
!  Updated to sixth-order accurate differencing on the interior with
!  the option for optimized fourth-order differencing.
!
!  Added g1v as input to improve efficiency
!
!  Revised: 6-28-95
!=============================================================================!
	use stencil
	implicit none
	
	integer :: ndof, nx, ny, optx, opty
	logical :: xper, yper
	logical :: lsym, rsym, bsym, tsym, carp
	real    :: v(ny,nx,ndof), g1v(ny,nx,ndof)
	real    :: g11v(ny,nx,ndof), g12v(ny,nx,ndof), g22v(ny,nx,ndof)
	real    :: dx, dy

	real, parameter :: zero = 0.0, one = 1.0, pt5 = 0.5, two = 2.0
	real :: dxinv, dyinv, dxsinv, dysinv
	real :: a, b, c, w
	real :: gx1, gx2, gx3, gx4, gx5, gx6
	real :: gy1, gy2, gy3, gy4, gy5, gy6
	real :: dx1, dx2, dx3, dx4, dx5, dx6, dx7
	real :: dy1, dy2, dy3, dy4, dy5, dy6, dy7
	
	integer :: i, j, idof
	real :: eps = 1.0e-12, isign
!=============================================================================!

	dxinv  = one / dx
	dyinv  = one / dy
	dxsinv = one / dx**2
	dysinv = one / dy**2

!.... seven point stencil in x

	if (optx.eq.0) then
	  c = 1.0 / 60.0
	else if (optx.eq.-1) then
	  c = 0.0
	else
	  w = 2.0 * 3.1415926535897932385e+0 / 12.0
	  c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
	      (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
	end if
	
	a = (2.0 + 15.0 * c) / 3.0
	b = -(1.0 + 48.0 * c) / 12.0

	gx1 =  -c * dxinv
	gx2 =  -b * dxinv
	gx3 =  -a * dxinv
	gx4 =   a * dxinv
	gx5 =   b * dxinv
	gx6 =   c * dxinv

	if (optx.eq.0) then
	  c = 1.0 / 90.0
	else if (optx.eq.-1) then
	  c = 0.0
	else
	  w = 2.0 * 3.1415926535897932385e+0 / 12.0
	  c = -(3.0*w**2 - 16.0*Sin(w/2.0)**2 + Sin(w)**2) / &
               (12.0*(-15.0*Sin(w/2.0)**2 + 6.0*Sin(w)**2 -  &
               Sin(3.0*w/2.0)**2))
	end if
	
	a = (4.0 + 45.0 * c) / 3.0
	b = -(1.0 + 72.0 * c) / 12.0

	dx1 =  c * dxsinv
	dx2 =  b * dxsinv
	dx3 =  a * dxsinv
	dx4 = -2.0 * ( a + b + c ) * dxsinv
	dx5 =  a * dxsinv
	dx6 =  b * dxsinv
	dx7 =  c * dxsinv

!.... seven point stencil in y

	if (opty.eq.0) then
	  c = 1.0 / 60.0
	else if (opty.eq.-1) then
	  c = 0.0
	else
	  w = 2.0 * 3.1415926535897932385e+0 / 12.0
	  c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
	      (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
	end if
	
	a = (2.0 + 15.0 * c) / 3.0
	b = -(1.0 + 48.0 * c) / 12.0

	gy1 =  -c * dyinv
	gy2 =  -b * dyinv
	gy3 =  -a * dyinv
	gy4 =   a * dyinv
	gy5 =   b * dyinv
	gy6 =   c * dyinv

	if (opty.eq.0) then
	  c = 1.0 / 90.0
	else if (opty.eq.-1) then
	  c = 0.0
	else
	  w = 2.0 * 3.1415926535897932385e+0 / 12.0
	  c = -(3.0*w**2 - 16.0*Sin(w/2.0)**2 + Sin(w)**2) / &
               (12.0*(-15.0*Sin(w/2.0)**2 + 6.0*Sin(w)**2 -  &
               Sin(3.0*w/2.0)**2))
	end if
	
	a = (4.0 + 45.0 * c) / 3.0
	b = -(1.0 + 72.0 * c) / 12.0

	dy1 =  c * dysinv
	dy2 =  b * dysinv
	dy3 =  a * dysinv
	dy4 = -2.0 * ( a + b + c ) * dysinv
	dy5 =  a * dysinv
	dy6 =  b * dysinv
	dy7 =  c * dysinv

!=============================================================================!
!.... compute the second derivative in x
!=============================================================================!

	if (xper) then

	  g11v(:,1,:)       = ( dx1 * v(:,nx-3,:)   + &
				dx2 * v(:,nx-2,:)   + &
				dx3 * v(:,nx-1,:)   + &
				dx4 * v(:,1,:)      + &
				dx5 * v(:,2,:)	    + &
				dx6 * v(:,3,:)	    + &
				dx7 * v(:,4,:)      ) 
  
	  g11v(:,2,:)       = ( dx1 * v(:,nx-2,:)   + &
				dx2 * v(:,nx-1,:)   + &
				dx3 * v(:,1,:)      + &
				dx4 * v(:,2,:)      + &
				dx5 * v(:,3,:)      + &
				dx6 * v(:,4,:)      + &
				dx7 * v(:,5,:)      ) 
  
	  g11v(:,3,:)       = ( dx1 * v(:,nx-1,:)   + &
				dx2 * v(:,1,:)      + &
				dx3 * v(:,2,:)      + &
				dx4 * v(:,3,:)      + &
				dx5 * v(:,4,:)      + &
				dx6 * v(:,5,:)      + &
				dx7 * v(:,6,:)      ) 
  
	  g11v(:,nx-2,:)    = ( dx1 * v(:,nx-5,:)   + &
				dx2 * v(:,nx-4,:)   + &
				dx3 * v(:,nx-3,:)   + &
				dx4 * v(:,nx-2,:)   + &
				dx5 * v(:,nx-1,:)   + &
				dx6 * v(:,1,:)      + &
				dx7 * v(:,2,:)      ) 
  
	  g11v(:,nx-1,:)    = ( dx1 * v(:,nx-4,:)   + &
				dx2 * v(:,nx-3,:)   + &
				dx3 * v(:,nx-2,:)   + &
				dx4 * v(:,nx-1,:)   + &
				dx5 * v(:,1,:)      + &
				dx6 * v(:,2,:)      + &
				dx7 * v(:,3,:)   ) 
  
	  g11v(:,nx,:) = g11v(:,1,:)
	  
	else
	
	  g11v(:,1,:)       = ( dd1 * v(:,1,:) + &
				dd2 * v(:,2,:) + &
				dd3 * v(:,3,:) + &
				dd4 * v(:,4,:) + &
				dd5 * v(:,5,:) ) * dxsinv
  
	  g11v(:,2,:)       = ( db1 * v(:,1,:) + &
				db2 * v(:,2,:) + &
				db3 * v(:,3,:) + &
				db4 * v(:,4,:) + &
				db5 * v(:,5,:) ) * dxsinv
  
	  g11v(:,3,:)       = ( da1 * v(:,1,:) + &
				da2 * v(:,2,:) + &
				da3 * v(:,3,:) + &
				da4 * v(:,4,:) + &
				da5 * v(:,5,:) ) * dxsinv
  
	  g11v(:,nx-2,:)    = ( da1 * v(:,nx-4,:) + &
				da2 * v(:,nx-3,:) + &
				da3 * v(:,nx-2,:) + &
				da4 * v(:,nx-1,:) + &
				da5 * v(:,nx,:)   ) * dxsinv
  
	  g11v(:,nx-1,:)    =   ( db1 * v(:,nx,:)   + &
				  db2 * v(:,nx-1,:) + &
				  db3 * v(:,nx-2,:) + &
				  db4 * v(:,nx-3,:) + &
				  db5 * v(:,nx-4,:) ) * dxsinv
  
	  g11v(:,nx,:)      =   ( dd1 * v(:,nx,:)   + &
				  dd2 * v(:,nx-1,:) + &
				  dd3 * v(:,nx-2,:) + &
				  dd4 * v(:,nx-3,:) + &
				  dd5 * v(:,nx-4,:) ) * dxsinv

	end if
	
!.... implement symmetry conditions

	if (lsym) then
	  do idof = 1, ndof
	    if (idof .eq. 3) then
	      isign = -one
	    else
	      isign = one
	    end if
	    
	    g11v(:,1,idof)    = ( isign * dx1 * v(:,4,idof)   	+ &
				  isign * dx2 * v(:,3,idof)   	+ &
				  isign * dx3 * v(:,2,idof)   	+ &
				  dx4 * v(:,1,idof)   		+ &
				  dx5 * v(:,2,idof)		+ &
				  dx6 * v(:,3,idof)		+ &
				  dx7 * v(:,4,idof)   		) 
    
	    g11v(:,2,idof)    = ( isign * dx1 * v(:,3,idof)   	+ &
				  isign * dx2 * v(:,2,idof)   	+ &
				  dx3 * v(:,1,idof)   		+ &
				  dx4 * v(:,2,idof)   		+ &
				  dx5 * v(:,3,idof)   		+ &
				  dx6 * v(:,4,idof)   		+ &
				  dx7 * v(:,5,idof)   		) 
    
	    g11v(:,3,idof)    = ( isign * dx1 * v(:,2,idof)   	+ &
				  dx2 * v(:,1,idof)   		+ &
				  dx3 * v(:,2,idof)   		+ &
				  dx4 * v(:,3,idof)   		+ &
				  dx5 * v(:,4,idof)   		+ &
				  dx6 * v(:,5,idof)   		+ &
				  dx7 * v(:,6,idof)   		) 
	  end do
	end if

	if (rsym) then
	  do idof = 1, ndof
	    if (idof .eq. 3) then
	      isign = -one
	    else
	      isign = one
	    end if
  
	    g11v(:,nx-2,idof) = ( dx1 * v(:,nx-5,idof)	        + &
				  dx2 * v(:,nx-4,idof)	        + &
				  dx3 * v(:,nx-3,idof)	        + &
				  dx4 * v(:,nx-2,idof)	        + &
				  dx5 * v(:,nx-1,idof)		+ &
				  dx6 * v(:,nx,idof)		+ &
				  isign * dx7 * v(:,nx-1,idof)	) 
    
	    g11v(:,nx-1,idof) = ( dx1 * v(:,nx-4,idof)	        + &
				  dx2 * v(:,nx-3,idof)	        + &
				  dx3 * v(:,nx-2,idof)		+ &
				  dx4 * v(:,nx-1,idof)		+ &
				  dx5 * v(:,nx,idof)		+ &
				  isign * dx6 * v(:,nx-1,idof)	+ &
				  isign * dx7 * v(:,nx-2,idof)	) 
    
	    g11v(:,nx,idof)   = ( dx1 * v(:,nx-3,idof)	        + &
				  dx2 * v(:,nx-2,idof)		+ &
				  dx3 * v(:,nx-1,idof)		+ &
				  dx4 * v(:,nx,idof)		+ &
				  isign * dx5 * v(:,nx-1,idof)	+ &
				  isign * dx6 * v(:,nx-2,idof)	+ &
				  isign * dx7 * v(:,nx-3,idof)	) 
	  end do
	end if
	
!.... interior

	g11v(:,4:nx-3,:) = ( dx1 * v(:,1:nx-6,:)   + &
	                     dx2 * v(:,2:nx-5,:)   + &
	                     dx3 * v(:,3:nx-4,:)   + &
	                     dx4 * v(:,4:nx-3,:)   + &
			     dx5 * v(:,5:nx-2,:)   + &
			     dx6 * v(:,6:nx-1,:)   + &
			     dx7 * v(:,7:nx  ,:)   ) 

!=============================================================================!
!.... compute the second derivative in y
!=============================================================================!

	if (yper) then
	
	  g22v(1,:,:)       = ( dy1 * v(ny-3,:,:)   + &
				dy2 * v(ny-2,:,:)   + &
				dy3 * v(ny-1,:,:)   + &
				dy4 * v(1,:,:)      + &
				dy5 * v(2,:,:)	    + &
				dy6 * v(3,:,:)	    + &
				dy7 * v(4,:,:)      ) 
  
	  g22v(2,:,:)       = ( dy1 * v(ny-2,:,:)   + &
				dy2 * v(ny-1,:,:)   + &
				dy3 * v(1,:,:)      + &
				dy4 * v(2,:,:)      + &
				dy5 * v(3,:,:)      + &
				dy6 * v(4,:,:)      + &
				dy7 * v(5,:,:)      ) 
  
	  g22v(3,:,:)       = ( dy1 * v(ny-1,:,:)   + &
				dy2 * v(1,:,:)      + &
				dy3 * v(2,:,:)      + &
				dy4 * v(3,:,:)      + &
				dy5 * v(4,:,:)      + &
				dy6 * v(5,:,:)      + &
				dy7 * v(6,:,:)      ) 
  
	  g22v(ny-2,:,:)    = ( dy1 * v(ny-5,:,:)   + &
				dy2 * v(ny-4,:,:)   + &
				dy3 * v(ny-3,:,:)   + &
				dy4 * v(ny-2,:,:)   + &
				dy5 * v(ny-1,:,:)   + &
				dy6 * v(1,:,:)      + &
				dy7 * v(2,:,:)      ) 
  
	  g22v(ny-1,:,:)    = ( dy1 * v(ny-4,:,:)   + &
				dy2 * v(ny-3,:,:)   + &
				dy3 * v(ny-2,:,:)   + &
				dy4 * v(ny-1,:,:)   + &
				dy5 * v(1,:,:)      + &
				dy6 * v(2,:,:)      + &
				dy7 * v(3,:,:)   ) 
  
	  g22v(ny,:,:) = g22v(1,:,:)

	else
	
	  g22v(1,:,:)      =  ( dd1 * v(1,:,:) + &
				dd2 * v(2,:,:) + &
				dd3 * v(3,:,:) + &
				dd4 * v(4,:,:) + &
				dd5 * v(5,:,:) ) * dysinv

!.... adiabatic wall  (doesn't work!)

!	  g22v(1,:,ndof)   =  ( dc1 * v(1,:,ndof) + &
!				dc2 * v(2,:,ndof) + &
!				dc3 * v(3,:,ndof) + &
!				dc4 * v(4,:,ndof) + &
!				dc5 * v(5,:,ndof) ) * dysinv

	  g22v(2,:,:)      =  ( db1 * v(1,:,:) + &
				db2 * v(2,:,:) + &
				db3 * v(3,:,:) + &
				db4 * v(4,:,:) + &
				db5 * v(5,:,:) ) * dysinv
  
	  g22v(3,:,:)       = ( da1 * v(1,:,:) + &
				da2 * v(2,:,:) + &
				da3 * v(3,:,:) + &
				da4 * v(4,:,:) + &
				da5 * v(5,:,:) ) * dysinv
  
	  g22v(ny-2,:,:)    = ( da1 * v(ny-4,:,:) + &
				da2 * v(ny-3,:,:) + &
				da3 * v(ny-2,:,:) + &
				da4 * v(ny-1,:,:) + &
				da5 * v(ny,:,:)   ) * dysinv
  
	  g22v(ny-1,:,:)   =  ( db1 * v(ny,:,:)   + &
				db2 * v(ny-1,:,:) + &
				db3 * v(ny-2,:,:) + &
				db4 * v(ny-3,:,:) + &
				db5 * v(ny-4,:,:) ) * dysinv
  
	  g22v(ny,:,:)     =  ( dd1 * v(ny,:,:)   + &
				dd2 * v(ny-1,:,:) + &
				dd3 * v(ny-2,:,:) + &
				dd4 * v(ny-3,:,:) + &
				dd5 * v(ny-4,:,:) ) * dysinv
			       
	endif

!.... interior

	g22v(4:ny-3,:,:) = ( dy1 * v(1:ny-6,:,:)   + &
	                     dy2 * v(2:ny-5,:,:)   + &
	                     dy3 * v(3:ny-4,:,:)   + &
	                     dy4 * v(4:ny-3,:,:)   + &
		             dy5 * v(5:ny-2,:,:)   + &
			     dy6 * v(6:ny-1,:,:)   + &
			     dy7 * v(7:ny  ,:,:)   ) 

!=============================================================================!
!.... compute the cross derivative
!=============================================================================!

	if (yper) then
	
	  g12v(1,:,:)       = ( gy1 * g1v(ny-3,:,:)	+ &
				gy2 * g1v(ny-2,:,:)	+ &
				gy3 * g1v(ny-1,:,:)	+ &
				gy4 * g1v(2,:,:)	+ &
				gy5 * g1v(3,:,:)	+ &
				gy6 * g1v(4,:,:)  	) 
  
	  g12v(2,:,:)       = ( gy1 * g1v(ny-2,:,:)	+ &
				gy2 * g1v(ny-1,:,:)	+ &
				gy3 * g1v(1,:,:)	+ &
				gy4 * g1v(3,:,:)	+ &
				gy5 * g1v(4,:,:)	+ &
				gy6 * g1v(5,:,:)	) 
  
	  g12v(3,:,:)       = ( gy1 * g1v(ny-1,:,:)	+ &
				gy2 * g1v(1,:,:)	+ &
				gy3 * g1v(2,:,:)	+ &
				gy4 * g1v(4,:,:)	+ &
				gy5 * g1v(5,:,:)	+ &
				gy6 * g1v(6,:,:)	) 
  
	  g12v(ny-2,:,:)    = ( gy1 * g1v(ny-5,:,:)	+ &
				gy2 * g1v(ny-4,:,:)	+ &
				gy3 * g1v(ny-3,:,:)	+ &
				gy4 * g1v(ny-1,:,:)	+ &
				gy5 * g1v(1,:,:)	+ &
				gy6 * g1v(2,:,:)	) 
  
	  g12v(ny-1,:,:)    = ( gy1 * g1v(ny-4,:,:)	+ &
				gy2 * g1v(ny-3,:,:)	+ &
				gy3 * g1v(ny-2,:,:)	+ &
				gy4 * g1v(1,:,:)	+ &
				gy5 * g1v(2,:,:)	+ &
				gy6 * g1v(3,:,:)	) 
  
	  g12v(ny,:,:) = g12v(1,:,:)

	else

	  g12v(1,:,:)     = ( gc1 * g1v(1,:,:)  + &
			      gc2 * g1v(2,:,:)  + &
			      gc3 * g1v(3,:,:)  + &
			      gc4 * g1v(4,:,:)  + &
			      gc5 * g1v(5,:,:)  ) * dyinv
  
	  g12v(2,:,:)     = ( gb1 * g1v(1,:,:)  + &
			      gb2 * g1v(2,:,:)  + &
			      gb3 * g1v(3,:,:)  + &
			      gb4 * g1v(4,:,:)  + &
			      gb5 * g1v(5,:,:)  ) * dyinv
  
	  g12v(3,:,:)     = ( ga1 * g1v(1,:,:)  + &
			      ga2 * g1v(2,:,:)  + &
			      ga3 * g1v(4,:,:)  + &
			      ga4 * g1v(5,:,:)  ) * dyinv
  
	  g12v(ny-2,:,:)  = ( ga1 * g1v(ny-4,:,:)  + &
			      ga2 * g1v(ny-3,:,:)  + &
			      ga3 * g1v(ny-1,:,:)  + &
			      ga4 * g1v(ny  ,:,:)  ) * dyinv
  
	  g12v(ny-1,:,:) = -( gb1 * g1v(ny  ,:,:)  + &
			      gb2 * g1v(ny-1,:,:)  + &
			      gb3 * g1v(ny-2,:,:)  + &
			      gb4 * g1v(ny-3,:,:)  + &
			      gb5 * g1v(ny-4,:,:)  ) * dyinv
  
	  g12v(ny,:,:)   = -( gc1 * g1v(ny  ,:,:)  + &
			      gc2 * g1v(ny-1,:,:)  + &
			      gc3 * g1v(ny-2,:,:)  + &
			      gc4 * g1v(ny-3,:,:)  + &
			      gc5 * g1v(ny-4,:,:)  ) * dyinv
			    
	end if
	
!.... interior

	g12v(4:ny-3,:,:) = ( gy1 * g1v(1:ny-6,:,:)	+ &
	                     gy2 * g1v(2:ny-5,:,:)	+ &
	                     gy3 * g1v(3:ny-4,:,:)	+ &
			     gy4 * g1v(5:ny-2,:,:)	+ &
			     gy5 * g1v(6:ny-1,:,:)	+ &
			     gy6 * g1v(7:ny  ,:,:)	) 

!.... Implement Carpenter's boundary stencil on the wall

	if (carp .and. (.not. yper) ) then
	  g12v(1,:,:)    = ( gg1 * g1v(1,:,:)	+ &
			     gg2 * g1v(2,:,:)  	+ &
			     gg3 * g1v(3,:,:)  	+ &
			     gg4 * g1v(4,:,:)  	+ &
			     gg5 * g1v(5,:,:)  	+ &
			     gg6 * g1v(6,:,:)  ) * dyinv
  
	  g12v(2,:,:)    = ( gh1 * g1v(1,:,:)  	+ &
			     gh2 * g1v(2,:,:)  	+ &
			     gh3 * g1v(3,:,:)  	+ &
			     gh4 * g1v(4,:,:)  	+ &
			     gh5 * g1v(5,:,:)  	+ &
			     gh6 * g1v(6,:,:)  ) * dyinv

	  g12v(3,:,:)    = ( gi1 * g1v(1,:,:)  	+ &
			     gi2 * g1v(2,:,:)  	+ &
			     gi3 * g1v(3,:,:)  	+ &
			     gi4 * g1v(4,:,:)  	+ &
			     gi5 * g1v(5,:,:)  	+ &
			     gi6 * g1v(6,:,:)  ) * dyinv

	  g12v(4,:,:)    = ( gj1 * g1v(1,:,:)  	+ &
			     gj2 * g1v(2,:,:)  	+ &
			     gj3 * g1v(3,:,:)  	+ &
			     gj4 * g1v(4,:,:)  	+ &
			     gj5 * g1v(5,:,:)  	+ &
			     gj6 * g1v(6,:,:)  ) * dyinv

	  g12v(ny-3,:,:) =-( gj1 * g1v(ny,:,:)    + &
			     gj2 * g1v(ny-1,:,:)  + &
			     gj3 * g1v(ny-2,:,:)  + &
			     gj4 * g1v(ny-3,:,:)  + &
			     gj5 * g1v(ny-4,:,:)  + &
			     gj6 * g1v(ny-5,:,:)  ) * dyinv

	  g12v(ny-2,:,:) =-( gi1 * g1v(ny,:,:)    + &
			     gi2 * g1v(ny-1,:,:)  + &
			     gi3 * g1v(ny-2,:,:)  + &
			     gi4 * g1v(ny-3,:,:)  + &
			     gi5 * g1v(ny-4,:,:)  + &
			     gi6 * g1v(ny-5,:,:)  ) * dyinv

	  g12v(ny-1,:,:) =-( gh1 * g1v(ny,:,:)    + &
			     gh2 * g1v(ny-1,:,:)  + &
			     gh3 * g1v(ny-2,:,:)  + &
			     gh4 * g1v(ny-3,:,:)  + &
			     gh5 * g1v(ny-4,:,:)  + &
			     gh6 * g1v(ny-5,:,:)  ) * dyinv

	  g12v(ny,:,:)   =-( gg1 * g1v(ny,:,:)    + &
			     gg2 * g1v(ny-1,:,:)  + &
			     gg3 * g1v(ny-2,:,:)  + &
			     gg4 * g1v(ny-3,:,:)  + &
			     gg5 * g1v(ny-4,:,:)  + &
			     gg6 * g1v(ny-5,:,:)  ) * dyinv
	end if

!.... implement a filter of roundoff noise

!	do idof = 1, ndof
!	  do i = 1, nx
!	    do j = 1, ny
!	      if ( abs(g11v(j,i,idof)) .lt. eps ) g11v(j,i,idof) = zero
!	      if ( abs(g12v(j,i,idof)) .lt. eps ) g12v(j,i,idof) = zero
!	      if ( abs(g22v(j,i,idof)) .lt. eps ) g22v(j,i,idof) = zero
!	    end do
!	  end do
!	end do

	return
	end
