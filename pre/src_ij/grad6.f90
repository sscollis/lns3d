!=============================================================================!
	module diff
!
!  Finite Difference coefficients of various orders
!
!=============================================================================!
!.... First derivatives
!=============================================================================!

!.... fourth order central difference ( 1 2 x 4 5 )

	real, parameter :: ga1 =  8.333333333333333333333E-02
	real, parameter :: ga2 = -6.666666666666666666667E-01
	real, parameter :: ga3 =  6.666666666666666666667E-01
	real, parameter :: ga4 = -8.333333333333333333333E-02

!.... fourth order one-sided ( x 2 3 4 5 )

	real, parameter :: gc1 = -2.083333333333333333333E+00
	real, parameter :: gc2 =  4.000000000000000000000E+00
	real, parameter :: gc3 = -3.000000000000000000000E+00
	real, parameter :: gc4 =  1.333333333333333333333E+00
	real, parameter :: gc5 = -2.500000000000000000000E-01

!.... fourth order biased difference ( 1 x 2 3 4 5 )

	real, parameter :: gb1 = -2.500000000000000000000E-01
	real, parameter :: gb2 = -8.333333333333333333333E-01
	real, parameter :: gb3 =  1.500000000000000000000E+00
	real, parameter :: gb4 = -5.000000000000000000000E-01
	real, parameter :: gb5 =  8.333333333333333333333E-02

!.... fifth order one-sided ( x 2 3 4 5 6 )

	real, parameter :: gd1 = -2.283333333333323E+00
	real, parameter :: gd2 =  4.999999999999969E+00
	real, parameter :: gd3 = -4.999999999999964E+00
	real, parameter :: gd4 =  3.333333333333314E+00
	real, parameter :: gd5 = -1.249999999999996E+00
	real, parameter :: gd6 =  2.000000000000000E-01

!.... sixth order one-sided ( x 2 3 4 5 6 7 )

	real, parameter :: ge1 = -2.449999999999916E+00
	real, parameter :: ge2 =  5.999999999999535E+00
	real, parameter :: ge3 = -7.499999999998948E+00
	real, parameter :: ge4 =  6.666666666665405E+00
	real, parameter :: ge5 = -3.749999999999149E+00
	real, parameter :: ge6 =  1.199999999999693E+00
	real, parameter :: ge7 = -1.666666666666202E-01

!.... sixth order one-pt biased ( 1 x 3 4 5 6 7 ) 

	real, parameter :: gf1 = -1.666666666666662E-01
	real, parameter :: gf2 = -1.283333333333328E+00
	real, parameter :: gf3 =  2.499999999999976E+00 
	real, parameter :: gf4 = -1.666666666666632E+00
	real, parameter :: gf5 =  8.333333333333081E-01
	real, parameter :: gf6 = -2.499999999999907E-01
	real, parameter :: gf7 =  3.333333333333190E-02

!.... sixth order two-pt biased ( 1 2 x 4 5 6 7 )

	real, parameter :: gg1 =  3.333333333333335E-02
	real, parameter :: gg2 = -3.999999999999997E-01
	real, parameter :: gg3 = -5.833333333333346E-01
	real, parameter :: gg4 =  1.333333333333335E+00
	real, parameter :: gg5 = -5.000000000000018E-01
	real, parameter :: gg6 =  1.333333333333340E-01
	real, parameter :: gg7 = -1.666666666666675E-02

!=============================================================================!
!.... Second derivatives
!=============================================================================!

!.... third order one-sided ( x 2 3 4 5 )

	real, parameter :: dd1 =  2.916666666666666666667E+00
	real, parameter :: dd2 = -8.666666666666666666667E+00
	real, parameter :: dd3 =  9.500000000000000000000E+00
	real, parameter :: dd4 = -4.666666666666666666667E+00
	real, parameter :: dd5 =  9.166666666666666666667E-01
	
!.... third order biased difference  ( 1 x 3 4 5 )

	real, parameter :: db1 =  9.166666666666666666667E-01
	real, parameter :: db2 = -1.666666666666666666667E+00
	real, parameter :: db3 =  5.000000000000000000000E-01
	real, parameter :: db4 =  3.333333333333333333333E-01
	real, parameter :: db5 = -8.333333333333333333333E-02
	
!.... fourth order central difference ( 1 2 x 4 5 )

	real, parameter :: da1 = -8.333333333333333333333E-02
	real, parameter :: da2 =  1.333333333333333333333E+00
	real, parameter :: da3 = -2.500000000000000000000E+00
	real, parameter :: da4 =  1.333333333333333333333E+00
	real, parameter :: da5 = -8.333333333333333333333E-02
	
!.... fourth order one-sided (assumes f'=0) [not smooth?]

	real, parameter :: dc1 = -5.763888888888888888889E+00
	real, parameter :: dc2 =  8.000000000000000000000E+00
	real, parameter :: dc3 = -3.000000000000000000000E+00
	real, parameter :: dc4 =  8.888888888888888888889E-01
	real, parameter :: dc5 = -1.250000000000000000000E-01
	
!.... fifth order one-sided ( x 2 3 4 5 6 7 )

	real, parameter :: de1 =  4.511111111110855E+00
	real, parameter :: de2 = -1.739999999999854E+01
	real, parameter :: de3 =  2.924999999999667E+01
	real, parameter :: de4 = -2.822222222221822E+01
	real, parameter :: de5 =  1.649999999999730E+01
	real, parameter :: de6 = -5.399999999999025E+00 
	real, parameter :: de7 =  7.611111111109641E-01

!.... fifth order one-pnt biased ( 1 x 3 4 5 6 7 )

	real, parameter :: df1 =  7.611111111111125E-01
	real, parameter :: df2 = -8.166666666666740E-01
	real, parameter :: df3 = -1.416666666666648E+00
	real, parameter :: df4 =  2.611111111111084E+00
	real, parameter :: df5 = -1.583333333333309E+00
	real, parameter :: df6 =  5.166666666666554E-01
	real, parameter :: df7 = -7.222222222222015E-02

!.... fifth order two-pnt biased ( 1 2 x 4 5 6 7 )

	real, parameter :: dg1 = -7.222222222222148E-02
	real, parameter :: dg2 =  1.266666666666662E+00 
	real, parameter :: dg3 = -2.333333333333324E+00
	real, parameter :: dg4 =  1.111111111111102E+00
	real, parameter :: dg5 =  8.333333333333835E-02
	real, parameter :: dg6 = -6.666666666666764E-02
	real, parameter :: dg7 =  1.111111111111120E-02

!.... sixth order one-sided ( x 2 3 4 5 6 7 8 )

	real, parameter :: dh1 =  5.211111111107453E+00
	real, parameter :: dh2 = -2.229999999997824E+01
	real, parameter :: dh3 =  4.394999999994423E+01
	real, parameter :: dh4 = -5.272222222214224E+01
	real, parameter :: dh5 =  4.099999999993068E+01
	real, parameter :: dh6 = -2.009999999996367E+01
	real, parameter :: dh7 =  5.661111111100436E+00
	real, parameter :: dh8 = -6.999999999986426E-01

!.... sixth order one-pnt biased ( 1 x 3 4 5 6 7 8 )

	real, parameter :: di1 =  6.999999999999785E-01
	real, parameter :: di2 = -3.888888888887693E-01
	real, parameter :: di3 = -2.700000000000279E+00
	real, parameter :: di4 =  4.750000000000371E+00
	real, parameter :: di5 = -3.722222222222532E+00
	real, parameter :: di6 =  1.800000000000160E+00
	real, parameter :: di7 = -5.000000000000471E-01
	real, parameter :: di8 =  6.111111111111715E-02

!.... sixth order two-pnt biased ( 1 2 x 4 5 6 7 8 )

	real, parameter :: dj1 = -6.111111111110915E-02
	real, parameter :: dj2 =  1.188888888888875E+00
	real, parameter :: dj3 = -2.099999999999959E+00
	real, parameter :: dj4 =  7.222222222221563E-01
	real, parameter :: dj5 =  4.722222222222845E-01
	real, parameter :: dj6 = -3.000000000000333E-01
	real, parameter :: dj7 =  8.888888888889816E-02
	real, parameter :: dj8 = -1.111111111111222E-02

	end module diff

!=============================================================================!
	subroutine grad2 (ndof, nx, ny, v, g11v, g12v, g22v, dx, dy, &
	                  optx, opty, xper, yper)
!
!  Take the second derivative of a 2-D field.
!  Updated to sixth order accurate differencing on the interior with
!  the option for optimized fourth order differencing.
!
!=============================================================================!
	use diff
	implicit none
	
	integer ndof, nx, ny, optx, opty
	logical :: xper, yper
	real v(ny,ndof,nx), g11v(ny,ndof,nx), g12v(ny,ndof,nx)
	real g22v(ny,ndof,nx), g1v(ny,ndof,nx)
	real dx, dy

	real, parameter :: one = 1.0
	real :: dxinv, dyinv, dxsinv, dysinv
	real :: a, b, c, w
	real :: gx1, gx2, gx3, gx4, gx5, gx6
	real :: gy1, gy2, gy3, gy4, gy5, gy6
	real :: dx1, dx2, dx3, dx4, dx5, dx6, dx7
	real :: dy1, dy2, dy3, dy4, dy5, dy6, dy7
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

	  g11v(:,:,1)       = ( dx1 * v(:,:,nx-3)   + &
				dx2 * v(:,:,nx-2)   + &
				dx3 * v(:,:,nx-1)   + &
				dx4 * v(:,:,1)      + &
				dx5 * v(:,:,2)	 + &
				dx6 * v(:,:,3)	 + &
				dx7 * v(:,:,4)      ) 
  
	  g11v(:,:,2)       = ( dx1 * v(:,:,nx-2)   + &
				dx2 * v(:,:,nx-1)   + &
				dx3 * v(:,:,1)      + &
				dx4 * v(:,:,2)      + &
				dx5 * v(:,:,3)      + &
				dx6 * v(:,:,4)      + &
				dx7 * v(:,:,5)      ) 
  
	  g11v(:,:,3)       = ( dx1 * v(:,:,nx-1)   + &
				dx2 * v(:,:,1)      + &
				dx3 * v(:,:,2)      + &
				dx4 * v(:,:,3)      + &
				dx5 * v(:,:,4)      + &
				dx6 * v(:,:,5)      + &
				dx7 * v(:,:,6)      ) 
  
	  g11v(:,:,nx-2)    = ( dx1 * v(:,:,nx-5)   + &
				dx2 * v(:,:,nx-4)   + &
				dx3 * v(:,:,nx-3)   + &
				dx4 * v(:,:,nx-2)   + &
				dx5 * v(:,:,nx-1)   + &
				dx6 * v(:,:,1)      + &
				dx7 * v(:,:,2)      ) 
  
	  g11v(:,:,nx-1)    = ( dx1 * v(:,:,nx-4)   + &
				dx2 * v(:,:,nx-3)   + &
				dx3 * v(:,:,nx-2)   + &
				dx4 * v(:,:,nx-1)   + &
				dx5 * v(:,:,1)      + &
				dx6 * v(:,:,2)      + &
				dx7 * v(:,:,3)   ) 
  
	  g11v(:,:,nx) = g11v(:,:,1)
	  
	else
	
	  g11v(:,:,1)       = ( dh1 * v(:,:,1) + &
				dh2 * v(:,:,2) + &
				dh3 * v(:,:,3) + &
				dh4 * v(:,:,4) + &
				dh5 * v(:,:,5) + &
				dh6 * v(:,:,6) + &
				dh7 * v(:,:,7) + &
				dh8 * v(:,:,8) ) * dxsinv
  
	  g11v(:,:,2)       = ( di1 * v(:,:,1) + &
				di2 * v(:,:,2) + &
				di3 * v(:,:,3) + &
				di4 * v(:,:,4) + &
				di5 * v(:,:,5) + &
				di6 * v(:,:,6) + &
				di7 * v(:,:,7) + &
				di8 * v(:,:,8) ) * dxsinv
  
	  g11v(:,:,3)       = ( dj1 * v(:,:,1) + &
				dj2 * v(:,:,2) + &
				dj3 * v(:,:,3) + &
				dj4 * v(:,:,4) + &
				dj5 * v(:,:,5) + &
				dj6 * v(:,:,6) + &
				dj7 * v(:,:,7) + &
				dj8 * v(:,:,8) ) * dxsinv
  
	  g11v(:,:,nx-2)   =  ( dj1 * v(:,:,nx) + &
				dj2 * v(:,:,nx-1) + &
				dj3 * v(:,:,nx-2) + &
				dj4 * v(:,:,nx-3) + &
				dj5 * v(:,:,nx-4) + &
				dj6 * v(:,:,nx-5) + &
				dj7 * v(:,:,nx-6) + &
				dj8 * v(:,:,nx-7)   ) * dxsinv
  
	  g11v(:,:,nx-1)    = ( di1 * v(:,:,nx)   + &
				di2 * v(:,:,nx-1) + &
				di3 * v(:,:,nx-2) + &
				di4 * v(:,:,nx-3) + &
				di5 * v(:,:,nx-4) + &
				di6 * v(:,:,nx-5) + &
				di7 * v(:,:,nx-6) + &
				di8 * v(:,:,nx-7) ) * dxsinv
  
	  g11v(:,:,nx)      = ( dh1 * v(:,:,nx)   + &
				dh2 * v(:,:,nx-1) + &
				dh3 * v(:,:,nx-2) + &
				dh4 * v(:,:,nx-3) + &
				dh5 * v(:,:,nx-4) + &
				dh6 * v(:,:,nx-5) + &
				dh7 * v(:,:,nx-6) + &
				dh8 * v(:,:,nx-7) ) * dxsinv
	end if
	
!.... interior

	g11v(:,:,4:nx-3) = ( dx1 * v(:,:,1:nx-6)   + &
	                     dx2 * v(:,:,2:nx-5)   + &
	                     dx3 * v(:,:,3:nx-4)   + &
	                     dx4 * v(:,:,4:nx-3)   + &
			     dx5 * v(:,:,5:nx-2)   + &
			     dx6 * v(:,:,6:nx-1)   + &
			     dx7 * v(:,:,7:nx  )   ) 

!=============================================================================!
!.... compute the second derivative in y
!=============================================================================!

	if (yper) then
	
	  g22v(1,:,:)       = ( dy1 * v(ny-3,:,:)   + &
				dy2 * v(ny-2,:,:)   + &
				dy3 * v(ny-1,:,:)   + &
				dy4 * v(1,:,:)      + &
				dy5 * v(2,:,:)	 + &
				dy6 * v(3,:,:)	 + &
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
	
	  g22v(1,:,:)      =  ( dh1 * v(1,:,:) + &
				dh2 * v(2,:,:) + &
				dh3 * v(3,:,:) + &
				dh4 * v(4,:,:) + &
				dh5 * v(5,:,:) + &
				dh6 * v(6,:,:) + &
				dh7 * v(7,:,:) + &
				dh8 * v(8,:,:) ) * dysinv
  
	  g22v(2,:,:)      =  ( di1 * v(1,:,:) + &
				di2 * v(2,:,:) + &
				di3 * v(3,:,:) + &
				di4 * v(4,:,:) + &
				di5 * v(5,:,:) + &
				di6 * v(6,:,:) + &
				di7 * v(7,:,:) + &
				di8 * v(8,:,:) ) * dysinv
  
	  g22v(3,:,:)       = ( dj1 * v(1,:,:) + &
				dj2 * v(2,:,:) + &
				dj3 * v(3,:,:) + &
				dj4 * v(4,:,:) + &
				dj5 * v(5,:,:) + &
				dj6 * v(6,:,:) + &
				dj7 * v(7,:,:) + &
				dj8 * v(8,:,:) ) * dysinv
  
	  g22v(ny-2,:,:)    = ( dj1 * v(ny  ,:,:) + &
				dj2 * v(ny-1,:,:) + &
				dj3 * v(ny-2,:,:) + &
				dj4 * v(ny-3,:,:) + &
				dj5 * v(ny-4,:,:) + &
				dj6 * v(ny-5,:,:) + &
				dj7 * v(ny-6,:,:) + &
				dj8 * v(ny-7,:,:)   ) * dysinv
  
	  g22v(ny-1,:,:)   =  ( di1 * v(ny,:,:)   + &
				di2 * v(ny-1,:,:) + &
				di3 * v(ny-2,:,:) + &
				di4 * v(ny-3,:,:) + &
				di5 * v(ny-4,:,:) + &
				di6 * v(ny-5,:,:) + &
				di7 * v(ny-6,:,:) + &
				di8 * v(ny-7,:,:) ) * dysinv
  
	  g22v(ny,:,:)     =  ( dh1 * v(ny,:,:)   + &
				dh2 * v(ny-1,:,:) + &
				dh3 * v(ny-2,:,:) + &
				dh4 * v(ny-3,:,:) + &
				dh5 * v(ny-4,:,:) + &
				dh6 * v(ny-5,:,:) + &
				dh7 * v(ny-6,:,:) + &
				dh8 * v(ny-7,:,:) ) * dysinv
			       
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

	if (xper) then

	  g1v(:,:,1)        = ( gx1 * v(:,:,nx-3)	+ &
				gx2 * v(:,:,nx-2)	+ &
				gx3 * v(:,:,nx-1)	+ &
				gx4 * v(:,:,2)	+ &
				gx5 * v(:,:,3)	+ &
				gx6 * v(:,:,4)	) 
  
	  g1v(:,:,2)        = ( gx1 * v(:,:,nx-2)	+ &
				gx2 * v(:,:,nx-1)	+ &
				gx3 * v(:,:,1)	+ &
				gx4 * v(:,:,3)	+ &
				gx5 * v(:,:,4)	+ &
				gx6 * v(:,:,5)	) 
  
	  g1v(:,:,3)        = ( gx1 * v(:,:,nx-1)	+ &
				gx2 * v(:,:,1)	+ &
				gx3 * v(:,:,2)	+ &
				gx4 * v(:,:,4)	+ &
				gx5 * v(:,:,5)	+ &
				gx6 * v(:,:,6)	) 
  
	  g1v(:,:,nx-2)     = ( gx1 * v(:,:,nx-5)	+ &
				gx2 * v(:,:,nx-4)	+ &
				gx3 * v(:,:,nx-3)	+ &
				gx4 * v(:,:,nx-1)	+ &
				gx5 * v(:,:,1)	+ &
				gx6 * v(:,:,2)	) 
  
	  g1v(:,:,nx-1)     = ( gx1 * v(:,:,nx-4)	+ &
				gx2 * v(:,:,nx-3)	+ &
				gx3 * v(:,:,nx-2)	+ &
				gx4 * v(:,:,1)	+ &
				gx5 * v(:,:,2)	+ &
				gx6 * v(:,:,3)	) 
  
	  g1v(:,:,nx) = g1v(:,:,1)
	  
	else
	
	  g1v(:,:,1)      = ( ge1 * v(:,:,1)  + &
			      ge2 * v(:,:,2)  + &
			      ge3 * v(:,:,3)  + &
			      ge4 * v(:,:,4)  + &
			      ge5 * v(:,:,5)  + &
			      ge6 * v(:,:,6)  + &
			      ge7 * v(:,:,7)  ) * dxinv
  
	  g1v(:,:,2)      = ( gf1 * v(:,:,1)  + &
			      gf2 * v(:,:,2)  + &
			      gf3 * v(:,:,3)  + &
			      gf4 * v(:,:,4)  + &
			      gf5 * v(:,:,5)  + &
			      gf6 * v(:,:,6)  + &
			      gf7 * v(:,:,7)  ) * dxinv
  
	  g1v(:,:,3)      = ( gg1 * v(:,:,1)  + &
			      gg2 * v(:,:,2)  + &
			      gg3 * v(:,:,3)  + &
			      gg4 * v(:,:,4)  + &
			      gg5 * v(:,:,5)  + &
			      gg6 * v(:,:,6)  + &
			      gg7 * v(:,:,7)  ) * dxinv
	  

	  g1v(:,:,nx-2)  = -( gg1 * v(:,:,nx)    + &
			      gg2 * v(:,:,nx-1)  + &
			      gg3 * v(:,:,nx-2)  + &
			      gg4 * v(:,:,nx-3)  + &
			      gg5 * v(:,:,nx-4)  + &
			      gg6 * v(:,:,nx-5)  + &
			      gg7 * v(:,:,nx-6)  ) * dxinv
  
	  g1v(:,:,nx-1)  = -( gf1 * v(:,:,nx  )  + &
			      gf2 * v(:,:,nx-1)  + &
			      gf3 * v(:,:,nx-2)  + &
			      gf4 * v(:,:,nx-3)  + &
			      gf5 * v(:,:,nx-4)  + &
			      gf6 * v(:,:,nx-5)  + &
			      gf7 * v(:,:,nx-6)  ) * dxinv
  
	  g1v(:,:,nx)    = -( ge1 * v(:,:,nx  )  + &
			      ge2 * v(:,:,nx-1)  + &
			      ge3 * v(:,:,nx-2)  + &
			      ge4 * v(:,:,nx-3)  + &
			      ge5 * v(:,:,nx-4)  + &
			      ge6 * v(:,:,nx-5)  + &
			      ge7 * v(:,:,nx-6)  ) * dxinv

	end if

!.... interior

	g1v(:,:,4:nx-3)  = ( gx1 * v(:,:,1:nx-6)	+ &
	                     gx2 * v(:,:,2:nx-5)	+ &
	                     gx3 * v(:,:,3:nx-4)	+ &
			     gx4 * v(:,:,5:nx-2)	+ &
			     gx5 * v(:,:,6:nx-1)	+ &
			     gx6 * v(:,:,7:nx  )	) 

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

	  g12v(1,:,:)    = ( ge1 * g1v(1,:,:)  + &
			     ge2 * g1v(2,:,:)  + &
			     ge3 * g1v(3,:,:)  + &
			     ge4 * g1v(4,:,:)  + &
			     ge5 * g1v(5,:,:)  + &
			     ge6 * g1v(6,:,:)  + &
			     ge7 * g1v(7,:,:)  ) * dyinv
  
	  g12v(2,:,:)    = ( gf1 * g1v(1,:,:)  + &
			     gf2 * g1v(2,:,:)  + &
			     gf3 * g1v(3,:,:)  + &
			     gf4 * g1v(4,:,:)  + &
			     gf5 * g1v(5,:,:)  + &
			     gf6 * g1v(6,:,:)  + &
			     gf7 * g1v(7,:,:)  ) * dyinv
  
	  g12v(3,:,:)    = ( gg1 * g1v(1,:,:)  + &
			     gg2 * g1v(2,:,:)  + &
			     gg3 * g1v(3,:,:)  + &
			     gg4 * g1v(4,:,:)  + &
			     gg5 * g1v(5,:,:)  + &
			     gg6 * g1v(6,:,:)  + &
			     gg7 * g1v(7,:,:)  ) * dyinv
  
	  g12v(ny-2,:,:)= -( gg1 * g1v(ny  ,:,:)  + &
			     gg2 * g1v(ny-1,:,:)  + &
			     gg3 * g1v(ny-2,:,:)  + &
			     gg4 * g1v(ny-3,:,:)  + &
			     gg5 * g1v(ny-4,:,:)  + &
			     gg6 * g1v(ny-5,:,:)  + &
			     gg7 * g1v(ny-6,:,:)  ) * dyinv
  
	  g12v(ny-1,:,:)= -( gf1 * g1v(ny  ,:,:)  + &
			     gf2 * g1v(ny-1,:,:)  + &
			     gf3 * g1v(ny-2,:,:)  + &
			     gf4 * g1v(ny-3,:,:)  + &
			     gf5 * g1v(ny-4,:,:)  + &
			     gf6 * g1v(ny-5,:,:)  + &
			     gf7 * g1v(ny-6,:,:)  ) * dyinv
  
	  g12v(ny,:,:)  = -( ge1 * g1v(ny  ,:,:)  + &
			     ge2 * g1v(ny-1,:,:)  + &
			     ge3 * g1v(ny-2,:,:)  + &
			     ge4 * g1v(ny-3,:,:)  + &
			     ge5 * g1v(ny-4,:,:)  + &
			     ge6 * g1v(ny-5,:,:)  + &
			     ge7 * g1v(ny-6,:,:)  ) * dyinv
			    
	end if
	
!.... interior

	g12v(4:ny-3,:,:) = ( gy1 * g1v(1:ny-6,:,:)	+ &
	                     gy2 * g1v(2:ny-5,:,:)	+ &
	                     gy3 * g1v(3:ny-4,:,:)	+ &
			     gy4 * g1v(5:ny-2,:,:)	+ &
			     gy5 * g1v(6:ny-1,:,:)	+ &
			     gy6 * g1v(7:ny  ,:,:)	) 

	return
	end

!=============================================================================!
	subroutine grad( ndof, nx, ny, v, g1v, g2v, dx, dy, optx, opty, &
	                 xper, yper)
!
!  Take the gradient of a 2-D field.
!  updated to fourth order accurate differencing
!
!=============================================================================!
	use diff
	implicit none
	
	integer ndof, nx, ny, optx, opty
	logical :: xper, yper
	real v(ny,ndof,nx), g1v(ny,ndof,nx), g2v(ny,ndof,nx)
	real dx, dy
	
	real,parameter :: one = 1.0
	real dxinv, dyinv
	real a, b, c, w
	real gx1, gx2, gx3, gx4, gx5, gx6
	real gy1, gy2, gy3, gy4, gy5, gy6
!=============================================================================!

	dxinv  = one / dx
	dyinv  = one / dy

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

!=============================================================================!
!.... compute the gradient in x
!=============================================================================!

	if (xper) then

	  g1v(:,:,1)        = ( gx1 * v(:,:,nx-3)	+ &
				gx2 * v(:,:,nx-2)	+ &
				gx3 * v(:,:,nx-1)	+ &
				gx4 * v(:,:,2)	+ &
				gx5 * v(:,:,3)	+ &
				gx6 * v(:,:,4)	) 
  
	  g1v(:,:,2)        = ( gx1 * v(:,:,nx-2)	+ &
				gx2 * v(:,:,nx-1)	+ &
				gx3 * v(:,:,1)	+ &
				gx4 * v(:,:,3)	+ &
				gx5 * v(:,:,4)	+ &
				gx6 * v(:,:,5)	) 
  
	  g1v(:,:,3)        = ( gx1 * v(:,:,nx-1)	+ &
				gx2 * v(:,:,1)	+ &
				gx3 * v(:,:,2)	+ &
				gx4 * v(:,:,4)	+ &
				gx5 * v(:,:,5)	+ &
				gx6 * v(:,:,6)	) 
  
	  g1v(:,:,nx-2)     = ( gx1 * v(:,:,nx-5)	+ &
				gx2 * v(:,:,nx-4)	+ &
				gx3 * v(:,:,nx-3)	+ &
				gx4 * v(:,:,nx-1)	+ &
				gx5 * v(:,:,1)	+ &
				gx6 * v(:,:,2)	) 
  
	  g1v(:,:,nx-1)     = ( gx1 * v(:,:,nx-4)	+ &
				gx2 * v(:,:,nx-3)	+ &
				gx3 * v(:,:,nx-2)	+ &
				gx4 * v(:,:,1)	+ &
				gx5 * v(:,:,2)	+ &
				gx6 * v(:,:,3)	) 
  
	  g1v(:,:,nx) = g1v(:,:,1)
	  
	else
	
	  g1v(:,:,1)      = ( gc1 * v(:,:,1)  + &
			      gc2 * v(:,:,2)  + &
			      gc3 * v(:,:,3)  + &
			      gc4 * v(:,:,4)  + &
			      gc5 * v(:,:,5)  ) * dxinv
  
	  g1v(:,:,2)      = ( gb1 * v(:,:,1)  + &
			      gb2 * v(:,:,2)  + &
			      gb3 * v(:,:,3)  + &
			      gb4 * v(:,:,4)  + &
			      gb5 * v(:,:,5)  ) * dxinv
  
	  g1v(:,:,3)      = ( ga1 * v(:,:,1)  + &
			      ga2 * v(:,:,2)  + &
			      ga3 * v(:,:,4)  + &
			      ga4 * v(:,:,5)  ) * dxinv
	  
	  g1v(:,:,nx-2)   = ( ga1 * v(:,:,nx-4)  + &
			      ga2 * v(:,:,nx-3)  + &
			      ga3 * v(:,:,nx-1)  + &
			      ga4 * v(:,:,nx  )  ) * dxinv
  
	  g1v(:,:,nx-1)  = -( gb1 * v(:,:,nx  )  + &
			      gb2 * v(:,:,nx-1)  + &
			      gb3 * v(:,:,nx-2)  + &
			      gb4 * v(:,:,nx-3)  + &
			      gb5 * v(:,:,nx-4)  ) * dxinv
  
	  g1v(:,:,nx)    = -( gc1 * v(:,:,nx  )  + &
			      gc2 * v(:,:,nx-1)  + &
			      gc3 * v(:,:,nx-2)  + &
			      gc4 * v(:,:,nx-3)  + &
			      gc5 * v(:,:,nx-4)  ) * dxinv

	end if

!.... interior

	g1v(:,:,4:nx-3)  = ( gx1 * v(:,:,1:nx-6)	+ &
	                     gx2 * v(:,:,2:nx-5)	+ &
	                     gx3 * v(:,:,3:nx-4)	+ &
			     gx4 * v(:,:,5:nx-2)	+ &
			     gx5 * v(:,:,6:nx-1)	+ &
			     gx6 * v(:,:,7:nx  )	) 

!=============================================================================!
!.... compute the gradient in y
!=============================================================================!

	if (yper) then
	
	  g2v(1,:,:)        = ( gy1 * v(ny-3,:,:)	+ &
				gy2 * v(ny-2,:,:)	+ &
				gy3 * v(ny-1,:,:)	+ &
				gy4 * v(2,:,:)	+ &
				gy5 * v(3,:,:)	+ &
				gy6 * v(4,:,:)  	) 
  
	  g2v(2,:,:)        = ( gy1 * v(ny-2,:,:)	+ &
				gy2 * v(ny-1,:,:)	+ &
				gy3 * v(1,:,:)	+ &
				gy4 * v(3,:,:)	+ &
				gy5 * v(4,:,:)	+ &
				gy6 * v(5,:,:)	) 
  
	  g2v(3,:,:)        = ( gy1 * v(ny-1,:,:)	+ &
				gy2 * v(1,:,:)	+ &
				gy3 * v(2,:,:)	+ &
				gy4 * v(4,:,:)	+ &
				gy5 * v(5,:,:)	+ &
				gy6 * v(6,:,:)	) 
  
	  g2v(ny-2,:,:)     = ( gy1 * v(ny-5,:,:)	+ &
				gy2 * v(ny-4,:,:)	+ &
				gy3 * v(ny-3,:,:)	+ &
				gy4 * v(ny-1,:,:)	+ &
				gy5 * v(1,:,:)	+ &
				gy6 * v(2,:,:)	) 
  
	  g2v(ny-1,:,:)     = ( gy1 * v(ny-4,:,:)	+ &
				gy2 * v(ny-3,:,:)	+ &
				gy3 * v(ny-2,:,:)	+ &
				gy4 * v(1,:,:)	+ &
				gy5 * v(2,:,:)	+ &
				gy6 * v(3,:,:)	) 
  
	  g2v(ny,:,:)       = g2v(1,:,:)

	else

	  g2v(1,:,:)     = ( gc1 * v(1,:,:)	+ &
			      gc2 * v(2,:,:)  	+ &
			      gc3 * v(3,:,:)  	+ &
			      gc4 * v(4,:,:)  	+ &
			      gc5 * v(5,:,:)  ) * dyinv
  
	  g2v(2,:,:)     = ( gb1 * v(1,:,:)  	+ &
			      gb2 * v(2,:,:)  	+ &
			      gb3 * v(3,:,:)  	+ &
			      gb4 * v(4,:,:)  	+ &
			      gb5 * v(5,:,:)  ) * dyinv
  
	  g2v(3,:,:)     = ( ga1 * v(1,:,:)  	+ &
			      ga2 * v(2,:,:)  	+ &
			      ga3 * v(4,:,:)  	+ &
			      ga4 * v(5,:,:)  ) * dyinv
  
	  g2v(ny-2,:,:)  = ( ga1 * v(ny-4,:,:)  	+ &
			      ga2 * v(ny-3,:,:)  	+ &
			      ga3 * v(ny-1,:,:)  	+ &
			      ga4 * v(ny  ,:,:)  ) 	* dyinv
  
	  g2v(ny-1,:,:) = -( gb1 * v(ny  ,:,:)  	+ &
			      gb2 * v(ny-1,:,:)  	+ &
			      gb3 * v(ny-2,:,:) 	 + &
			      gb4 * v(ny-3,:,:)  	+ &
			      gb5 * v(ny-4,:,:)  ) 	* dyinv
  
	  g2v(ny,:,:)   = -( gc1 * v(ny  ,:,:)  	+ &
			      gc2 * v(ny-1,:,:)  	+ &
			      gc3 * v(ny-2,:,:)  	+ &
			      gc4 * v(ny-3,:,:)  	+ &
			      gc5 * v(ny-4,:,:)  ) 	* dyinv
	end if
	
!.... interior

	g2v(4:ny-3,:,:) = ( gy1 * v(1:ny-6,:,:)	+ &
	                     gy2 * v(2:ny-5,:,:)	+ &
	                     gy3 * v(3:ny-4,:,:)	+ &
			     gy4 * v(5:ny-2,:,:)	+ &
			     gy5 * v(6:ny-1,:,:)	+ &
			     gy6 * v(7:ny  ,:,:)	) 

	return
	end

!=============================================================================!
	subroutine g1(v, g1v, nx, dx, optx, xper)
!
!  Take the gradient in the streamwise direction
!  updated to fourth order accurate differencing
!
!=============================================================================!
	use diff
	implicit none
	
	integer nx, optx
	logical :: xper
	real v(nx), g1v(nx), dx

	real,parameter :: one = 1.0
	real dxinv
	real a, b, c, w
	real gx1, gx2, gx3, gx4, gx5, gx6
!=============================================================================!

	dxinv  = one / dx

!.... seven point stencil in x

	if (optx.eq.0) then
	  c = 1.0 / 60.0
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

!=============================================================================!
!.... compute the gradient in x
!=============================================================================!

	if (xper) then

	  g1v(1)        = ( gx1 * v(nx-3)	+ &
				gx2 * v(nx-2)	+ &
				gx3 * v(nx-1)	+ &
				gx4 * v(2)	+ &
				gx5 * v(3)	+ &
				gx6 * v(4)	) 
  
	  g1v(2)        = ( gx1 * v(nx-2)	+ &
				gx2 * v(nx-1)	+ &
				gx3 * v(1)	+ &
				gx4 * v(3)	+ &
				gx5 * v(4)	+ &
				gx6 * v(5)	) 
  
	  g1v(3)        = ( gx1 * v(nx-1)	+ &
				gx2 * v(1)	+ &
				gx3 * v(2)	+ &
				gx4 * v(4)	+ &
				gx5 * v(5)	+ &
				gx6 * v(6)	) 
  
	  g1v(nx-2)     = ( gx1 * v(nx-5)	+ &
				gx2 * v(nx-4)	+ &
				gx3 * v(nx-3)	+ &
				gx4 * v(nx-1)	+ &
				gx5 * v(1)	+ &
				gx6 * v(2)	) 
  
	  g1v(nx-1)     = ( gx1 * v(nx-4)	+ &
				gx2 * v(nx-3)	+ &
				gx3 * v(nx-2)	+ &
				gx4 * v(1)	+ &
				gx5 * v(2)	+ &
				gx6 * v(3)	) 
  
	  g1v(nx) = g1v(1)
	  
	else
	
	  g1v(1)      = ( gc1 * v(1)  + &
			      gc2 * v(2)  + &
			      gc3 * v(3)  + &
			      gc4 * v(4)  + &
			      gc5 * v(5)  ) * dxinv
  
	  g1v(2)      = ( gb1 * v(1)  + &
			      gb2 * v(2)  + &
			      gb3 * v(3)  + &
			      gb4 * v(4)  + &
			      gb5 * v(5)  ) * dxinv
  
	  g1v(3)      = ( ga1 * v(1)  + &
			      ga2 * v(2)  + &
			      ga3 * v(4)  + &
			      ga4 * v(5)  ) * dxinv
	  
	  g1v(nx-2)   = ( ga1 * v(nx-4)  + &
			      ga2 * v(nx-3)  + &
			      ga3 * v(nx-1)  + &
			      ga4 * v(nx  )  ) * dxinv
  
	  g1v(nx-1)  = -( gb1 * v(nx  )  + &
			      gb2 * v(nx-1)  + &
			      gb3 * v(nx-2)  + &
			      gb4 * v(nx-3)  + &
			      gb5 * v(nx-4)  ) * dxinv
  
	  g1v(nx)    = -( gc1 * v(nx  )  + &
			      gc2 * v(nx-1)  + &
			      gc3 * v(nx-2)  + &
			      gc4 * v(nx-3)  + &
			      gc5 * v(nx-4)  ) * dxinv

	end if

!.... interior

	g1v(4:nx-3)  = ( gx1 * v(1:nx-6)	+ &
	                     gx2 * v(2:nx-5)	+ &
	                     gx3 * v(3:nx-4)	+ &
			     gx4 * v(5:nx-2)	+ &
			     gx5 * v(6:nx-1)	+ &
			     gx6 * v(7:nx  )	) 

	return
	end

!=============================================================================!
	subroutine g2(v, g2v, ny, dy, opty, yper)
!
!  Take the gradient in the wall normal direction
!  updated to fourth order accurate differencing
!
!=============================================================================!
	use diff
	implicit none

	integer ny, opty
	logical :: yper
	real v(ny), g2v(ny), dy

	real,parameter :: one = 1.0
	real dyinv
	real a, b, c, w
	real gy1, gy2, gy3, gy4, gy5, gy6
!=============================================================================!

	dyinv  = one / dy

!.... seven point stencil in y

	if (opty.eq.0) then
	  c = 1.0 / 60.0
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

!=============================================================================!
!.... compute the gradient in y
!=============================================================================!

	if (yper) then
	
	  g2v(1)       = ( gy1 * v(ny-3)	+ &
				gy2 * v(ny-2)	+ &
				gy3 * v(ny-1)	+ &
				gy4 * v(2)	+ &
				gy5 * v(3)	+ &
				gy6 * v(4)  	) 
  
	  g2v(2)       = ( gy1 * v(ny-2)	+ &
				gy2 * v(ny-1)	+ &
				gy3 * v(1)	+ &
				gy4 * v(3)	+ &
				gy5 * v(4)	+ &
				gy6 * v(5)	) 
  
	  g2v(3)       = ( gy1 * v(ny-1)	+ &
				gy2 * v(1)	+ &
				gy3 * v(2)	+ &
				gy4 * v(4)	+ &
				gy5 * v(5)	+ &
				gy6 * v(6)	) 
  
	  g2v(ny-2)    = ( gy1 * v(ny-5)	+ &
				gy2 * v(ny-4)	+ &
				gy3 * v(ny-3)	+ &
				gy4 * v(ny-1)	+ &
				gy5 * v(1)	+ &
				gy6 * v(2)	) 
  
	  g2v(ny-1)    = ( gy1 * v(ny-4)	+ &
				gy2 * v(ny-3)	+ &
				gy3 * v(ny-2)	+ &
				gy4 * v(1)	+ &
				gy5 * v(2)	+ &
				gy6 * v(3)	) 
  
	  g2v(ny) = g2v(1)

	else

	  g2v(1)     = ( gc1 * v(1)  + &
			      gc2 * v(2)  + &
			      gc3 * v(3)  + &
			      gc4 * v(4)  + &
			      gc5 * v(5)  ) * dyinv
  
	  g2v(2)     = ( gb1 * v(1)  + &
			      gb2 * v(2)  + &
			      gb3 * v(3)  + &
			      gb4 * v(4)  + &
			      gb5 * v(5)  ) * dyinv
  
	  g2v(3)     = ( ga1 * v(1)  + &
			      ga2 * v(2)  + &
			      ga3 * v(4)  + &
			      ga4 * v(5)  ) * dyinv
  
	  g2v(ny-2)  = ( ga1 * v(ny-4)  + &
			      ga2 * v(ny-3)  + &
			      ga3 * v(ny-1)  + &
			      ga4 * v(ny  )  ) * dyinv
  
	  g2v(ny-1) = -( gb1 * v(ny  )  + &
			      gb2 * v(ny-1)  + &
			      gb3 * v(ny-2)  + &
			      gb4 * v(ny-3)  + &
			      gb5 * v(ny-4)  ) * dyinv
  
	  g2v(ny)   = -( gc1 * v(ny  )  + &
			      gc2 * v(ny-1)  + &
			      gc3 * v(ny-2)  + &
			      gc4 * v(ny-3)  + &
			      gc5 * v(ny-4)  ) * dyinv
			    
	end if
	
!.... interior

	g2v(4:ny-3) = ( gy1 * v(1:ny-6)	+ &
	                     gy2 * v(2:ny-5)	+ &
	                     gy3 * v(3:ny-4)	+ &
			     gy4 * v(5:ny-2)	+ &
			     gy5 * v(6:ny-1)	+ &
			     gy6 * v(7:ny  )	) 

	return
	end
